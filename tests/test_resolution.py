import json
import numpy as np
import pytest

from embryo import Config, Simulation
from embryo.resolution import manufactured, prepare, run, project_occupancy


def test_manufactured_cells_have_same_targets_signals_and_no_resampled_fields():
    coarse = manufactured(Config(grid=32, extent=2.24))
    fine = manufactured(Config(grid=40, extent=2.24))
    assert len(coarse.phi) == len(fine.phi) == 16
    np.testing.assert_array_equal(coarse.target, fine.target)
    np.testing.assert_array_equal(coarse.activator, fine.activator)
    np.testing.assert_array_equal(coarse.inhibitor, fine.inhibitor)
    assert not coarse.divisions and not fine.divisions
    assert np.isinf(coarse.due).all()
    assert coarse.phi.shape == (16, 32, 32, 32)
    assert fine.phi.shape == (16, 40, 40, 40)
    # Different quadrature is measured; it is not removed by projecting volumes.
    assert not np.array_equal(coarse.volumes(), fine.volumes())
    from scipy.special import expit
    center = np.array([-.87, -.29, -.29])
    expected = expit(np.sqrt(2) * (.3 - np.linalg.norm(fine.xyz - center[:, None, None, None], axis=0)) / .085)
    np.testing.assert_array_equal(fine.phi[0], expected.astype(np.float32))


def test_projection_is_diagnostic_and_preserves_same_grid_field():
    from embryo.model import occupancy
    sim = manufactured(Config(grid=32, extent=2.24))
    before = sim.phi.copy()
    projected = project_occupancy(sim, 32)
    np.testing.assert_allclose(projected, occupancy(sim.phi), atol=1e-7)
    np.testing.assert_array_equal(sim.phi, before)


def test_preparation_separates_space_and_time_and_refuses_overwrite(tmp_path):
    output = tmp_path/'study'
    protocol = prepare(output, grids=(32,36,40), duration=.015)
    spatial = [protocol['cases'][name] for name in protocol['space_cases']]
    temporal = [protocol['cases'][name] for name in protocol['time_cases']]
    assert {c['dt'] for c in spatial} == {.00375}
    assert {c['grid'] for c in temporal} == {36}
    assert {c['extent'] for c in protocol['cases'].values()} == {2.24}
    assert len(protocol['cases']) == 5  # Finest time case reused, not recomputed.
    with pytest.raises(FileExistsError):
        prepare(output)


def test_completed_smoke_screen_exports_all_cases_and_matches_physical_times(tmp_path):
    output = tmp_path/'study'
    protocol = prepare(output, grids=(32,36,40), duration=.015)
    run(output)
    status = json.loads((output/'status.json').read_text())
    assert status['state'] == 'completed'
    assert len(status['completed_cases']) == 5
    report = json.loads((output/'comparison.json').read_text())
    assert set(report['pairs']) == {'space','time'}
    assert (output/'RESULTS.md').exists()
    for name in protocol['cases']:
        sim = Simulation.restore(output/name/'final_state.npz')
        assert sim.time == pytest.approx(.015)
        assert len(sim.phi) == 16
        assert (output/name/'viewer.html').exists()
    with pytest.raises(ValueError, match='already'):
        run(output)


def test_upstream_failed_screen_blocks_new_simulations(tmp_path):
    output, dependency = tmp_path/'study', tmp_path/'domain'
    prepare(output, grids=(32,36,40), duration=.015)
    dependency.mkdir()
    (dependency/'status.json').write_text(json.dumps({'state':'completed'}))
    (dependency/'comparison.json').write_text(json.dumps({
        'largest_domain_boundary_screen_pass':False, 'largest_pair_shape_screen_pass':True}))
    run(output, after=dependency, poll_interval=.001)
    status = json.loads((output/'status.json').read_text())
    assert status['state'] == 'blocked'
    assert status['completed_cases'] == []
    assert not (output/'space-32').exists()
