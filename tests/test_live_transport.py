import json
import numpy as np
import pytest
from embryo import Config, Simulation
from embryo.transport import contact_transport, conservative_transport, transport_graph
from embryo.signaling import integrate, gm_jacobian, mode_growth
from embryo.live_transport_validation import run


def test_unequal_volume_flux_and_full_jacobian():
    graph = transport_graph(conservative_transport([[0, 2, 0], [2, 0, 1], [0, 1, 0]], [.2, .7, 1.3]))
    np.testing.assert_allclose(graph.masses @ graph.delta, 0, atol=1e-14)
    np.testing.assert_allclose(graph.delta @ np.ones(3), 0, atol=1e-14)
    j = gm_jacobian()
    full = np.kron(j, np.eye(3)) + np.kron(np.diag([.02, .4]), graph.delta)
    predicted = np.concatenate([np.linalg.eigvals(j - value * np.diag([.02, .4])) for value in graph.eigenvalues])
    np.testing.assert_allclose(np.sort_complex(np.linalg.eigvals(full)), np.sort_complex(predicted), atol=1e-12)


def test_live_step_matches_transport_then_amount_preserving_mechanics():
    sim = Simulation(Config(grid=16, interface_width=.16, max_cells=2))
    sim.divide(0, [1, 0, 0])
    daughters, neck, overlap, fraction = sim._cleavage_fields(0, sim.divisions[0])
    sim._complete_division(0, daughters, fraction, neck, overlap)
    sim.activator[:] = [.8, 1.2]
    sim.inhibitor[:] = [1.1, .9]
    volumes = sim.volumes()
    graph = sim.signaling_graph()
    expected = integrate(sim.activator, sim.inhibitor, graph, sim.config.dt, da=.02, dh=.4)
    sim.step()
    for name, concentration in zip(('activator', 'inhibitor'), expected):
        np.testing.assert_allclose(sim.volumes() * getattr(sim, name), volumes * concentration, rtol=1e-13)
    assert sim.graph_snapshot()['operator'] == 'conservative'


def test_division_conserves_measured_amount_even_if_daughter_volume_changes():
    sim = Simulation(Config(grid=16, interface_width=.16))
    sim.activator[:] = 1.7
    sim.inhibitor[:] = .8
    sim.divide(0, [1, 0, 0])
    daughters, neck, overlap, fraction = sim._cleavage_fields(0, sim.divisions[0])
    daughters[0] *= .98  # Force a volume discrepancy to test conservative remapping.
    amounts = [sim.volumes() @ sim.activator, sim.volumes() @ sim.inhibitor]
    sim._complete_division(0, daughters, fraction, neck, overlap)
    np.testing.assert_allclose([sim.volumes() @ sim.activator, sim.volumes() @ sim.inhibitor], amounts, rtol=1e-14)


def test_checkpoint_does_not_silently_migrate_legacy_transport(tmp_path):
    sim = Simulation(Config(grid=16, interface_width=.16))
    path = tmp_path / 'current.npz'
    sim.checkpoint(path)
    assert Simulation.restore(path).config.signal_transport == 'conservative'
    with np.load(path, allow_pickle=False) as data:
        payload = {key: data[key] for key in data.files}
    meta = json.loads(str(payload['metadata']))
    del meta['config']['signal_transport']
    meta['config'].update(signal_da=1., signal_dh=20.)
    payload['metadata'] = json.dumps(meta)
    np.savez(tmp_path / 'legacy.npz', **payload)
    restored = Simulation.restore(tmp_path / 'legacy.npz')
    assert restored.config.signal_transport == 'random_walk'
    assert restored.config.signal_da == 1.


def test_contact_geometry_rejects_ambiguous_coincident_centers():
    with pytest.raises(ValueError, match='distinct centers'):
        contact_transport([[0, 1], [1, 0]], [1, 1], np.zeros((2, 3)), .1)


def test_calibrated_contact_refinement_and_measured_mode_growth():
    assert run()['passed']
