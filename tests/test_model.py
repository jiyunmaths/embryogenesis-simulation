import numpy as np
import pytest

from embryo.model import Config, Simulation, occupancy


def small(**kwargs):
    return Simulation(Config(grid=20, interface_width=0.14, steps=20, **kwargs))


def finish_first_division(sim):
    for _ in range(400):
        sim.step()
        if len(sim.phi) > 1:
            return
    pytest.fail("cytokinesis did not resolve within the test interval")


def test_division_conserves_occupancy_targets_and_regulator_amount():
    sim = small(competence_cells=16, max_cells=2)
    before = occupancy(sim.phi).sum(axis=0)
    volume = sim.volumes().sum()
    sim.fate[:] = 0.3
    sim.divide(0, [1, 2, 3])
    assert len(sim.phi) == 1
    assert sim.lineage[0]["division"] is None
    np.testing.assert_allclose(occupancy(sim.phi).sum(axis=0), before, atol=3e-7)
    finish_first_division(sim)
    np.testing.assert_allclose(sim.volumes().sum(), volume, rtol=1e-6)
    np.testing.assert_allclose(sim.target.sum(), sim.initial_volume)
    np.testing.assert_allclose(np.dot(sim.fate, sim.target), 0.3 * sim.initial_volume)
    assert list(sim.parents) == [0, 0]
    assert sim.lineage[0]["division"] >= sim.config.cytokinesis_duration


def test_pair_contact_is_symmetric_and_self_contact_is_zero():
    sim = small()
    sim.divide(0)
    finish_first_division(sim)
    contact, exposure = sim.contacts()
    np.testing.assert_allclose(contact, contact.T)
    assert np.all(np.diag(contact) == 0)
    assert contact[0, 1] > 0
    assert np.all((exposure >= 0) & (exposure <= 1))


def test_single_cell_remains_centered_and_isotropic():
    sim = small(max_cells=1)
    for _ in range(40):
        sim.step()
    np.testing.assert_allclose(sim.centers(), 0, atol=1e-6)
    assert abs(sim.metrics()["axis_ratio"] - 1) < 1e-6
    assert abs(sim.metrics()["relative_volume_error"]) < .05


def test_checkpoint_continuation_is_exact(tmp_path):
    sim = small(division_interval=.15, max_cells=4, competence_cells=2)
    for _ in range(20):
        sim.step()
    sim.checkpoint(tmp_path / "state.npz")
    restored = Simulation.restore(tmp_path / "state.npz")
    assert sim.divisions  # Checkpoint in the middle of a live furrow.
    assert restored.metrics() == sim.metrics()
    for _ in range(180):
        sim.step()
        restored.step()
    np.testing.assert_array_equal(sim.phi, restored.phi)
    np.testing.assert_array_equal(sim.fate, restored.fate)
    np.testing.assert_array_equal(sim.activator, restored.activator)
    np.testing.assert_array_equal(sim.inhibitor, restored.inhibitor)
    np.testing.assert_array_equal(sim.polarity, restored.polarity)
    assert sim.graph_events == restored.graph_events
    assert sim.lineage == restored.lineage
    assert sim.divisions == restored.divisions
    assert len(sim.phi) > 1


def test_mechanics_control_has_no_differentiation():
    sim = small(differentiation=False, feedback=False, division_interval=.12, max_cells=4)
    for _ in range(300):
        sim.step()
    assert len(sim.phi) == 4
    np.testing.assert_array_equal(sim.fate, 0)
    assert np.all(np.isfinite(sim.phi))


def test_time_refinement_for_single_cell():
    first = small(max_cells=1, dt=.015)
    second = small(max_cells=1, dt=.0075)
    for _ in range(20):
        first.step()
        second.step()
        second.step()
    np.testing.assert_allclose(first.volumes(), second.volumes(), rtol=.002)


def test_control_preserves_division_sampling():
    enabled = small(division_interval=.12, max_cells=4, competence_cells=2, feedback=False)
    disabled = small(division_interval=.12, max_cells=4, differentiation=False, feedback=False)
    for _ in range(300):
        enabled.step()
        disabled.step()
    assert enabled.lineage == disabled.lineage
    np.testing.assert_array_equal(enabled.due, disabled.due)
    assert enabled.rng.bit_generator.state == disabled.rng.bit_generator.state


def test_hertwig_axis_follows_rotated_elongated_cell():
    sim = small()
    axis = np.array([1.0, 2.0, -1.0])
    axis /= np.linalg.norm(axis)
    along = np.einsum("d,dijk->ijk", axis, sim.xyz)
    radius2 = np.sum(sim.xyz**2, axis=0) - along**2
    level = np.sqrt(along**2 / 1.0**2 + radius2 / .4**2)
    sim.phi[0] = .5 * (1 - np.tanh((level - 1) * .4 / sim.config.interface_width))
    for _ in range(10):
        assert abs(sim.division_direction(0) @ axis) > .999


def test_isotropic_shape_does_not_select_grid_axis():
    sim = small()
    directions = np.array([sim.division_direction(0) for _ in range(50)])
    assert np.all(np.var(directions, axis=0) > .15)
    assert np.any(np.max(np.abs(directions), axis=1) < .8)


def test_cytokinesis_has_no_startup_contact_or_volume_jump():
    sim = small(max_cells=2)
    phi = sim.phi.copy()
    contacts, exposure = sim.contacts()
    sim.divide(0, [1, 0, 0])
    np.testing.assert_array_equal(sim.phi, phi)
    np.testing.assert_array_equal(sim.contacts()[0], contacts)
    np.testing.assert_array_equal(sim.contacts()[1], exposure)
    initial_volume = sim.volumes().sum()
    narrowed = False
    for _ in range(400):
        sim.step()
        # This includes the completion step, before daughter relaxation.
        assert abs(sim.volumes().sum() / initial_volume - 1) < 2e-6
        if len(sim.phi) == 2:
            break
        _, neck, _, _ = sim._cleavage_fields(0, sim.divisions[0])
        narrowed |= .2 < neck < .9
    assert len(sim.phi) == 2
    assert narrowed
    assert sim.lineage[0]["neck_at_abscission"] <= sim.config.neck_threshold
    assert sim.lineage[0]["overlap_at_abscission"] <= sim.config.division_overlap_tolerance
    from scipy.ndimage import label
    assert all(label(field >= .5)[1] == 1 for field in sim.phi)
    for _ in range(20):
        sim.step()
        assert np.max(np.abs(sim.volumes() / sim.target - 1)) < .025


def test_active_divisions_reserve_capacity_and_cannot_restart():
    sim = small(max_cells=3)
    assert sim.divide(0)
    assert not sim.divide(0)
    finish_first_division(sim)
    assert sim.divide(0)
    assert not sim.divide(1)
    for _ in range(250):
        sim.step()
    assert len(sim.phi) == 3


def test_timer_does_not_force_an_unresolved_neck_to_split():
    sim = small(max_cells=2)
    sim.divide(0)
    sim.time = 10 * sim.config.cytokinesis_duration
    sim._finish_divisions()
    assert len(sim.phi) == 1
    assert sim.metrics()["overdue_divisions"] == 1


def test_old_checkpoint_cannot_silently_resume_different_dynamics(tmp_path):
    import json
    sim = small()
    sim.checkpoint(tmp_path / "state.npz")
    with np.load(tmp_path / "state.npz", allow_pickle=False) as data:
        arrays = dict(data)
    meta = json.loads(str(arrays["metadata"]))
    del meta["schema_version"]
    arrays["metadata"] = json.dumps(meta)
    np.savez_compressed(tmp_path / "legacy.npz", **arrays)
    with pytest.raises(ValueError, match="previous model"):
        Simulation.restore(tmp_path / "legacy.npz")


@pytest.mark.parametrize("sign", [-1, 1])
def test_isolated_fate_switch_has_two_attractors(sign):
    sim = small(competence_cells=1, fate_noise=0, neighbor_inhibition=0, exposure_bias=0)
    sim.fate[:] = sign * .1
    for _ in range(1000):
        sim.update_fate(np.zeros((1, 1)), np.array([.5]))
    np.testing.assert_allclose(sim.fate, sign, atol=.001)


@pytest.mark.parametrize("kwargs", [{"dt": -1}, {"max_cells": 0}, {"grid": 8}, {"fate_tension": 1}, {"dt": 100},
                                   {"division_orientation": "unknown"}, {"cytokinesis_duration": .01},
                                   {"axis_degeneracy": 0}, {"neck_threshold": 1}, {"ring_strength": 100}])
def test_invalid_config_is_rejected(kwargs):
    with pytest.raises(ValueError):
        Config(**kwargs).validate()
