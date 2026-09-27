from dataclasses import asdict
import json
import numpy as np
import pytest

from embryo import Config, Simulation
from embryo.domain import enlarge, prepare, run_branch, compare


def test_enlargement_preserves_state_voxels_coordinates_and_amounts():
    sim = Simulation(Config(grid=20, interface_width=.14))
    sim.activator[:] = 1.7
    sim.polarity[:] = [.2, 0, .1]
    sim.divide(0, [1, 0, 0])
    original_phi = sim.phi.copy()
    bigger = enlarge(sim, 28)
    assert bigger.config.extent == pytest.approx(2.24)
    assert bigger.dx == sim.dx
    np.testing.assert_array_equal(bigger.phi[:, 4:-4, 4:-4, 4:-4], sim.phi)
    np.testing.assert_array_equal(bigger.xyz[:, 4:-4, 4:-4, 4:-4], sim.xyz)
    np.testing.assert_array_equal(bigger.phi[:, :4], 0)
    np.testing.assert_allclose(bigger.volumes(), sim.volumes(), rtol=1e-14)
    np.testing.assert_allclose(bigger.volumes() @ bigger.activator, sim.volumes() @ sim.activator, rtol=1e-14)
    for name in ('fate', 'target', 'activator', 'inhibitor', 'polarity', 'ids', 'parents', 'due'):
        np.testing.assert_array_equal(getattr(bigger, name), getattr(sim, name))
    for name in ('rng', 'fate_rng', 'signal_rng'):
        assert getattr(bigger, name).bit_generator.state == getattr(sim, name).bit_generator.state
    assert bigger.divisions == sim.divisions
    assert bigger.lineage == sim.lineage
    assert bigger.step_number == sim.step_number
    bigger.phi[:] = 0
    np.testing.assert_array_equal(sim.phi, original_phi)


@pytest.mark.parametrize('grid', [18, 21, 22.5, True])
def test_enlargement_rejects_cropping_or_misaligned_voxels(grid):
    with pytest.raises(ValueError):
        enlarge(Simulation(Config(grid=20, interface_width=.14)), grid)


def test_fixed_domain_protocol_exports_comparison_and_cannot_overwrite(tmp_path):
    sim = Simulation(Config(grid=12, interface_width=.18, max_cells=1))
    checkpoint = tmp_path / 'source.npz'
    sim.checkpoint(checkpoint)
    output = tmp_path / 'study'
    protocol = prepare(output, grids=[12, 16], until=.03, sample_interval=.015, checkpoint=checkpoint)
    assert protocol['dx'] == sim.dx
    for grid in protocol['grids']:
        report = run_branch(output, grid)
        assert report['final']['time'] == pytest.approx(.03)
        restored = Simulation.restore(output / f'grid-{grid}' / 'final_state.npz')
        assert restored.config.signal_transport == 'conservative'
        assert restored.dx == pytest.approx(sim.dx)
    report = compare(output)
    assert report['largest_domain_boundary_screen_pass']
    assert report['pairs'][0]['cell_ids_align']
    assert (output / 'comparison.png').is_file()
    assert json.loads((output / 'comparison.json').read_text())['protocol']['source_sha256']
    with pytest.raises(FileExistsError):
        prepare(output)
    with pytest.raises(FileExistsError):
        run_branch(output, 12)


def test_large_domain_preset_preserves_default_dx():
    from pathlib import Path
    large = Config(**json.loads(Path('configs/large_domain.json').read_text()))
    small = Config()
    large.validate()
    assert 2 * large.extent / large.grid == pytest.approx(2 * small.extent / small.grid)
    assert large.extent > small.extent


def test_extension_preserves_each_branch_and_joins_history_without_duplicates(tmp_path):
    from embryo.domain import prepare_extension
    sim = Simulation(Config(grid=12, interface_width=.18, max_cells=1))
    checkpoint = tmp_path / 'source.npz'
    sim.checkpoint(checkpoint)
    source, output = tmp_path / 'early', tmp_path / 'later'
    prepare(source, grids=[12, 16], until=.03, sample_interval=.015, checkpoint=checkpoint)
    expected = {}
    for grid in (12, 16):
        run_branch(source, grid)
        reference = Simulation.restore(source / f'grid-{grid}' / 'final_state.npz')
        for _ in range(2):
            reference.step()
        expected[grid] = reference
    prepare_extension(source, output, until=.06)
    for grid in (12, 16):
        run_branch(output, grid)
        actual = Simulation.restore(output / f'grid-{grid}' / 'final_state.npz')
        for name in ('phi', 'activator', 'inhibitor', 'fate', 'polarity', 'target'):
            np.testing.assert_array_equal(getattr(actual, name), getattr(expected[grid], name))
        history = json.loads((output / f'grid-{grid}' / 'history.json').read_text())
        assert [row['time'] for row in history] == pytest.approx([0, .015, .03, .045, .06])
        frames = json.loads((output / f'grid-{grid}' / 'trajectory.json').read_text())['frames']
        assert len(frames) == len(history)
        assert frames[0]['metrics'] == history[0]
        assert (output / f'grid-{grid}' / 'progress_state.npz').exists()
    assert compare(output)['largest_domain_boundary_screen_pass']
    with pytest.raises(FileExistsError):
        prepare_extension(source, output, until=.09)


def test_monitor_publishes_final_report_for_completed_study(tmp_path):
    from embryo.domain_monitor import monitor
    sim = Simulation(Config(grid=12, interface_width=.18, max_cells=1))
    checkpoint = tmp_path / 'source.npz'
    sim.checkpoint(checkpoint)
    output = tmp_path / 'study'
    prepare(output, grids=[12, 16], until=.015, sample_interval=.015, checkpoint=checkpoint)
    for grid in (12, 16):
        run_branch(output, grid)
    report = monitor(output)
    assert report['largest_domain_boundary_screen_pass']
    assert json.loads((output / 'status.json').read_text())['state'] == 'completed'
    assert 'PASS' in (output / 'RESULTS.md').read_text()
