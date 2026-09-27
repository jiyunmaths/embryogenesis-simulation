"""Controlled domain enlargement at fixed voxel spacing and identical initial state."""
import argparse
from copy import deepcopy
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .model import Config, Simulation
from .shape import observe


def enlarge(simulation, grid):
    """Center-pad with empty space, without resampling or changing cell/RNG state.

    An even grid increment is required so old voxel centers remain aligned.
    This does not recover tails outside the old domain or undo prior wall effects.
    Prefer a common zygote, whose boundary occupancy is negligible.
    """
    old = simulation.config.grid
    if isinstance(grid, bool) or not isinstance(grid, int) or grid < old or (grid - old) % 2:
        raise ValueError('grid must be an integer >= old grid with an even increment')
    result = deepcopy(simulation)
    if grid == old:
        return result
    pad = (grid - old) // 2
    result.config = replace(simulation.config, grid=grid, extent=simulation.dx * grid / 2)
    result.config.validate()
    result.phi = np.pad(simulation.phi, ((0, 0), (pad, pad), (pad, pad), (pad, pad)))
    # Preserve the original coordinate centers bit-for-bit on the shared grid.
    axis = simulation.xyz[0, :, 0, 0]
    expanded = np.concatenate([axis[0] - simulation.dx * np.arange(pad, 0, -1),
                               axis, axis[-1] + simulation.dx * np.arange(1, pad + 1)])
    result.xyz = np.stack(np.meshgrid(expanded, expanded, expanded, indexing='ij'))
    return result


def observe_domain(sim):
    row = observe(sim)
    mask = np.any(sim.phi >= .5, axis=0)
    if np.any(mask):
        indices = np.argwhere(mask)
        gap_voxels = min(int(indices.min()), sim.config.grid - 1 - int(indices.max()))
        row['cell_surface_clearance'] = float((gap_voxels + .5) * sim.dx)
    else:
        row['cell_surface_clearance'] = None
    row['domain_half_width'] = sim.config.extent
    row['grid_spacing'] = sim.dx
    return row


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def prepare(output, grids=(40, 48, 56), until=30., sample_interval=.6, checkpoint=None):
    output = Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    sim = (Simulation.restore(checkpoint) if checkpoint else
           Simulation(Config(grid=40, extent=1.6, signal_transport='conservative')))
    if not np.isfinite([until, sample_interval]).all() or until <= sim.time or sample_interval <= 0:
        raise ValueError('until must exceed source time and sample_interval must be positive')
    if not np.isclose(until / sim.config.dt, round(until / sim.config.dt)) or not np.isclose(
            sample_interval / sim.config.dt, round(sample_interval / sim.config.dt)):
        raise ValueError('end time and sampling interval must be multiples of dt')
    if len(grids) < 2 or len(set(grids)) != len(grids) or sorted(grids) != list(grids):
        raise ValueError('provide at least two distinct increasing grid sizes')
    # Validate before writing anything; avoid constructing all large arrays at once.
    for grid in grids:
        if isinstance(grid, bool) or not isinstance(grid, int) or grid < sim.config.grid or (grid - sim.config.grid) % 2:
            raise ValueError('grids must be >= source grid with even increments')
    output.mkdir(parents=True)
    sim.checkpoint(output / 'source.npz')
    source_metrics = observe_domain(sim)
    protocol = {'schema': 1, 'source': str(Path(checkpoint).resolve()) if checkpoint else 'common_zygote',
                'source_sha256': hashlib.sha256((output / 'source.npz').read_bytes()).hexdigest(),
                'source_config': asdict(sim.config), 'source_time': sim.time,
                'source_boundary_occupancy': source_metrics['boundary_occupancy'],
                'grids': list(grids), 'half_widths': [grid * sim.dx / 2 for grid in grids],
                'dx': sim.dx, 'until': until, 'sample_interval': sample_interval,
                'criteria': {'boundary_occupancy_max': .01,
                             'largest_pair_axis_ratio_relative_difference_max': .01},
                'scope': 'Fixed voxel spacing, parameters, source state, and random streams; new space has phi=0. No interpolation. Boundary screening and domain sensitivity are distinct from spatial convergence.',
                'source_warning': 'Padding cannot remove prior boundary influence or recover omitted tails.'}
    write_json(output / 'protocol.json', protocol)
    return protocol


def prepare_extension(source, output, until=180.):
    """Continue every completed domain from its own checkpoint, retaining ancestry."""
    source, output = Path(source).resolve(), Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    previous = json.loads((source / 'protocol.json').read_text())
    if not np.isfinite(until) or until <= previous['until']:
        raise ValueError('until must exceed the completed source end time')
    checkpoints = {}
    for grid in previous['grids']:
        branch = source / f'grid-{grid}'
        checkpoint = branch / 'final_state.npz'
        sim = Simulation.restore(checkpoint)
        history = json.loads((branch / 'history.json').read_text())
        if (sim.config.grid != grid or not np.isclose(sim.dx, previous['dx'])
                or not np.isclose(sim.time, previous['until'])
                or not history or not np.isclose(history[-1]['time'], sim.time)):
            raise ValueError('source checkpoint and history do not match the completed protocol')
        if not np.isclose(until / sim.config.dt, round(until / sim.config.dt)):
            raise ValueError('end time must be a multiple of dt')
        checkpoints[str(grid)] = {'path': str(checkpoint),
                                  'sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest()}
    protocol = deepcopy(previous)
    protocol.update(until=until, continuation={'source': str(source),
                    'from_time': previous['until'], 'checkpoints': checkpoints,
                    'protocol_sha256': hashlib.sha256((source / 'protocol.json').read_bytes()).hexdigest()})
    output.mkdir(parents=True)
    write_json(output / 'protocol.json', protocol)
    return protocol


def run_branch(output, grid):
    output = Path(output)
    protocol = json.loads((output / 'protocol.json').read_text())
    if grid not in protocol['grids']:
        raise ValueError('grid is not in the prepared study')
    path = output / f'grid-{grid}'
    path.mkdir(exist_ok=False)
    continuation = protocol.get('continuation')
    old_history, old_frames = [], []
    if continuation:
        source = Path(continuation['source']) / f'grid-{grid}'
        entry = continuation['checkpoints'][str(grid)]
        checkpoint = Path(entry['path'])
        if hashlib.sha256(checkpoint.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError('source checkpoint changed after preparing continuation')
        sim = Simulation.restore(checkpoint)
        old_history = json.loads((source / 'history.json').read_text())[:-1]
        old_frames = json.loads((source / 'trajectory.json').read_text())['frames'][:-1]
    else:
        sim = enlarge(Simulation.restore(output / 'source.npz'), grid)
    sim.config = replace(sim.config, steps=round(protocol['until'] / sim.config.dt),
                         save_every=round(protocol['sample_interval'] / sim.config.dt))
    write_json(path / 'config.json', asdict(sim.config))
    write_json(path / 'preflight.json', sim.graph_snapshot())
    history, frames = old_history, old_frames
    started = time.monotonic()
    start_step = sim.step_number
    while True:
        if ((sim.step_number - start_step) % sim.config.save_every == 0
                or sim.step_number == sim.config.steps):
            row = observe_domain(sim)
            history.append(row)
            frames.append({'metrics': row, 'cells': sim.surfaces(max_points=450),
                           'graph': sim.graph_snapshot()})
            write_json(path / 'history.json', history)
            if len(history) % 10 == 1 or sim.step_number == sim.config.steps:
                sim.checkpoint(path / 'progress_state.npz')
                print(f'grid={grid}: t={sim.time:.2f}, N={len(sim.phi)}, boundary={row["boundary_occupancy"]:.3g}, R={row["axis_ratio"]:.5f}, elapsed={time.monotonic()-started:.1f}s', flush=True)
        if sim.step_number == sim.config.steps:
            break
        sim.step()
    sim.checkpoint(path / 'final_state.npz')
    payload = json.dumps({'config': asdict(sim.config), 'frames': frames}, separators=(',', ':'), allow_nan=False)
    (path / 'trajectory.json').write_text(payload)
    from .surface import viewer_template
    (path / 'viewer.html').write_text(viewer_template().replace('__SIMULATION_DATA__', payload))
    report = {'grid': grid, 'extent': sim.config.extent, 'dx': sim.dx,
              'elapsed_seconds': time.monotonic() - started, 'final': history[-1],
              'continuation': continuation,
              'first_boundary_threshold_crossing': next((row['time'] for row in history if row['boundary_occupancy'] >= protocol['criteria']['boundary_occupancy_max']), None),
              'max_boundary_occupancy': max(row['boundary_occupancy'] for row in history),
              'minimum_surface_clearance': min(row['cell_surface_clearance'] for row in history
                                                if row['cell_surface_clearance'] is not None)}
    write_json(path / 'analysis.json', report)
    return report


def compare(output):
    output = Path(output)
    protocol = json.loads((output / 'protocol.json').read_text())
    reports = [json.loads((output / f'grid-{grid}' / 'analysis.json').read_text()) for grid in protocol['grids']]
    histories = [json.loads((output / f'grid-{grid}' / 'history.json').read_text()) for grid in protocol['grids']]
    for history in histories[1:]:
        if [row['time'] for row in history] != [row['time'] for row in histories[0]]:
            raise ValueError('comparison requires matched observation times')
    pairs = []
    for index in range(len(reports) - 1):
        a, b = reports[index], reports[index + 1]
        small = Simulation.restore(output / f'grid-{a["grid"]}' / 'final_state.npz')
        large = Simulation.restore(output / f'grid-{b["grid"]}' / 'final_state.npz')
        aligned = np.array_equal(small.ids, large.ids) and np.array_equal(small.parents, large.parents)
        field_error = None
        if aligned:
            pad = (b['grid'] - a['grid']) // 2
            core = large.phi[:, pad:-pad, pad:-pad, pad:-pad]
            field_error = float(np.linalg.norm((small.phi - core).ravel()) / np.linalg.norm(core.ravel()))
        ratio_error = max(abs(x['axis_ratio'] / y['axis_ratio'] - 1)
                          for x, y in zip(histories[index], histories[index + 1]))
        pairs.append({'grids': [a['grid'], b['grid']], 'max_axis_ratio_relative_difference': ratio_error,
                      'final_common_box_phase_field_relative_l2': field_error,
                      'cell_ids_align': bool(aligned),
                      'final_activator_max_absolute_difference': float(np.max(abs(small.activator - large.activator))) if aligned else None})
    largest_clear = reports[-1]['max_boundary_occupancy'] < protocol['criteria']['boundary_occupancy_max']
    stable_ratio = pairs[-1]['max_axis_ratio_relative_difference'] < protocol['criteria']['largest_pair_axis_ratio_relative_difference_max']
    report = {'protocol': protocol, 'runs': reports, 'pairs': pairs,
              'largest_domain_boundary_screen_pass': largest_clear,
              'largest_pair_shape_screen_pass': stable_ratio,
              'interpretation': 'Finite-window, single-seed domain screen. Boundary occupancy below threshold does not prove zero boundary influence. Shape agreement does not establish signaling or full-field convergence. Cell resolution is unchanged.'}
    write_json(output / 'comparison.json', report)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for record, history in zip(reports, histories):
        label = f'{record["grid"]}³, half-width={record["extent"]:.2f}'
        times = [row['time'] for row in history]
        for ax, key in zip(axes.flat, ['boundary_occupancy', 'axis_ratio', 'activator_std', 'cell_surface_clearance']):
            ax.plot(times, [row[key] for row in history], label=label)
            ax.set(xlabel='Model time', ylabel=key.replace('_', ' '))
    axes[0, 0].axhline(.01, color='black', ls=':', label='Boundary screen')
    axes[0, 0].set_yscale('symlog', linthresh=1e-12)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f'Fixed-spacing domain study: dx={protocol["dx"]:.3g}, {protocol["source_config"]["signal_transport"]} transport')
    fig.savefig(output / 'comparison.png', dpi=150)
    plt.close(fig)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    create = sub.add_parser('prepare')
    create.add_argument('--output', type=Path, required=True)
    create.add_argument('--grids', type=int, nargs='+', default=[40, 48, 56])
    create.add_argument('--until', type=float, default=30.)
    create.add_argument('--sample-interval', type=float, default=.6)
    create.add_argument('--checkpoint', type=Path)
    extension = sub.add_parser('extend')
    extension.add_argument('--source', type=Path, required=True)
    extension.add_argument('--output', type=Path, required=True)
    extension.add_argument('--until', type=float, default=180.)
    branch = sub.add_parser('run')
    branch.add_argument('--output', type=Path, required=True)
    branch.add_argument('--grid', type=int, required=True)
    comparison = sub.add_parser('compare')
    comparison.add_argument('--output', type=Path, required=True)
    args = vars(parser.parse_args())
    command = args.pop('command')
    if command == 'prepare': prepare(**args)
    elif command == 'extend': prepare_extension(**args)
    elif command == 'run': run_branch(**args)
    else: compare(**args)


if __name__ == '__main__':
    main()
