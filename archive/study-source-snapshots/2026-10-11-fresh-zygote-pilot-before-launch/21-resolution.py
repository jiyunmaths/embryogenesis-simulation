"""Separate voxel/time refinement of a manufactured 16-cell conservative embryo.

The same analytic phase fields are sampled afresh on every grid. This tests
local coupled numerics, not spontaneous development or cleavage convergence.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from scipy.special import expit

from .model import Config, Simulation, occupancy
from .domain import observe_domain


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def manufactured(config):
    """An explicitly imposed aggregate and weak chemical perturbation, no cleavage."""
    if config.max_cells != 16 or config.signal_transport != 'conservative':
        raise ValueError('the controlled study requires 16 conservative compartments')
    sim = Simulation(config)
    centers = np.array([[x, y, z] for x in (-.87, -.29, .29, .87)
                        for y in (-.29, .29) for z in (-.29, .29)])
    radius, width = .3, config.interface_width
    fields = []
    for center in centers:
        distance = np.sqrt(np.sum((sim.xyz - center[:, None, None, None])**2, axis=0))
        fields.append(expit(np.sqrt(2) * (radius - distance) / width).astype(np.float32))
    sim.phi = np.stack(fields)
    # Infinite-space diffuse-sphere volume, with an exponentially small tail
    # beyond radius+30*width. This common target is independent of voxel count.
    volume = quad(lambda r: 4 * np.pi * r*r * occupancy(expit(np.sqrt(2)*(radius-r)/width)),
                  0, radius + 30 * width, epsabs=1e-12, epsrel=1e-12)[0]
    sim.target = np.full(16, volume)
    sim.initial_volume = 16 * volume
    sim.ids = np.arange(16)
    sim.parents = np.full(16, -1, dtype=int)
    sim.next_id = 16
    sim.due = np.full(16, np.inf)
    sim.lineage = [{'id': i, 'parent': -1, 'birth': 0., 'division': None} for i in range(16)]
    sim.divisions = {}
    sim.fate = np.zeros(16)
    sim.polarity = np.zeros((16, 3))
    sim.activator = 1 + .001 * (np.sin(np.pi * centers[:, 0]) + .5 * np.cos(np.pi * centers[:, 2]))
    sim.inhibitor = np.ones(16)
    sim.last_contact = np.zeros((16, 16))
    sim.last_exposure = np.ones(16)
    return sim


def _steps(duration, dt):
    if not np.isfinite([duration, dt]).all() or min(duration, dt) <= 0:
        raise ValueError('duration and dt must be positive and finite')
    steps = round(duration / dt)
    if steps < 1 or not np.isclose(steps * dt, duration, atol=1e-12, rtol=0):
        raise ValueError('duration must be an integer multiple of every dt')
    return steps


def prepare(output, grids=(56, 72, 88), dts=(.015, .0075, .00375), duration=.6, extent=2.24):
    output = Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    if len(grids) != 3 or list(grids) != sorted(set(grids)):
        raise ValueError('provide three distinct increasing grids')
    if len(dts) != 3 or list(dts) != sorted(set(dts), reverse=True):
        raise ValueError('provide three distinct decreasing time steps')
    cases = []
    for grid in grids:
        cases.append((f'space-{grid}', grid, dts[-1], 'space'))
    for dt in dts[:-1]:
        cases.append((f'time-{dt:g}', grids[1], dt, 'time'))
    configs = {}
    for name, grid, dt, _ in cases:
        config = Config(grid=grid, extent=extent, dt=dt, max_cells=16,
                        steps=_steps(duration, dt), save_every=1,
                        signal_transport='conservative', signal_partition_noise=0.)
        config.validate()
        configs[name] = config
    output.mkdir(parents=True)
    protocol = {'schema': 1, 'geometry': '16 analytic diffuse spheres, radius 0.30, 4x2x2 centers spaced 0.58',
                'duration': duration, 'extent': extent, 'interface_width': .085,
                'grids': list(grids), 'dts': list(dts),
                'space_cases': [name for name, _, _, kind in cases if kind == 'space'],
                'time_cases': [f'time-{dt:g}' for dt in dts[:-1]] + [f'space-{grids[1]}'],
                'cases': {name: asdict(config) for name, config in configs.items()},
                'criteria': {'finest_pair_field_relative_l2_max': .01,
                             'finest_pair_signal_relative_l2_max': .01,
                             'finest_pair_axis_ratio_relative_difference_max': .01,
                             'boundary_occupancy_max': .01,
                             'cell_volume_error_max': .05,
                             'minimum_radius_grid_cells': 4.},
                'scope': 'Manufactured initial aggregate and chemical perturbation; no cell division. Fixed physical geometry, interface width, diffusivities, and domain. Space and time vary separately. No claim of developmental or sharp-interface convergence.'}
    write_json(output / 'protocol.json', protocol)
    rows = []
    for name in protocol['space_cases']:
        sim = manufactured(configs[name])
        graph = sim.graph_snapshot()
        rows.append({'case': name, 'grid': sim.config.grid, 'dx': sim.dx,
                     'relative_volume_quadrature_error': float(np.max(abs(sim.volumes()/sim.target-1))),
                     'metrics': observe_domain(sim), 'graph': graph})
    changes = []
    for a, b in zip(rows[:-1], rows[1:]):
        ga, gb = np.array(a['graph']['weights']), np.array(b['graph']['weights'])
        ea, eb = np.array(a['graph']['eigenvalues'])[1:], np.array(b['graph']['eigenvalues'])[1:]
        changes.append({'grids': [a['grid'], b['grid']],
                        'conductance_relative_frobenius_difference': float(np.linalg.norm(ga-gb)/np.linalg.norm(gb)),
                        'max_positive_eigenvalue_relative_difference': float(np.max(abs(ea-eb)/np.maximum(eb, 1e-12))),
                        'unstable_mode_counts': [len(a['graph']['unstable_modes']), len(b['graph']['unstable_modes'])]})
    preflight = {'states': rows, 'changes': changes,
                 'interpretation': 'Frozen manufactured geometry, voxel quadrature only. This does not establish accuracy of the overlap-to-area physical closure or convergence of subsequent mechanics.'}
    write_json(output / 'frozen_geometry.json', preflight)
    write_json(output / 'status.json', {'state': 'prepared', 'completed_cases': []})
    return protocol


def project_occupancy(sim, grid):
    """Trilinear diagnostic projection, never used to initialize or evolve a case."""
    axis = (np.arange(grid)+.5) * (2*sim.config.extent/grid) - sim.config.extent
    coordinates = np.stack(np.meshgrid(*([axis]*3), indexing='ij'))
    index = (coordinates + sim.config.extent) / sim.dx - .5
    return np.stack([map_coordinates(occupancy(field), index, order=1, mode='nearest', prefilter=False)
                     for field in sim.phi])


def run_case(output, name):
    output = Path(output)
    protocol = json.loads((output/'protocol.json').read_text())
    config = Config(**protocol['cases'][name])
    path = output/name
    path.mkdir(exist_ok=False)
    sim = manufactured(config)
    initial = {'metrics': observe_domain(sim), 'cells': sim.surfaces(max_points=300), 'graph': sim.graph_snapshot()}
    history = [initial['metrics']]
    started = time.monotonic()
    for _ in range(config.steps):
        sim.step()
        # Every-step metrics catch clipping/volume excursions between endpoints.
        history.append(observe_domain(sim))
    sim.checkpoint(path/'final_state.npz')
    frame = {'metrics': history[-1], 'cells': sim.surfaces(max_points=300), 'graph': sim.graph_snapshot()}
    payload = json.dumps({'config': asdict(config), 'frames': [initial, frame]}, separators=(',', ':'), allow_nan=False)
    (path/'trajectory.json').write_text(payload)
    from .surface import viewer_template
    template = viewer_template().replace('One cell. Two possible identities.', 'Controlled cell-resolution test.')
    template = template.replace('A freely evolving 3D aggregate of deformable cells.', 'Prescribed initial geometry. Two endpoint snapshots; intermediate meshes are not recorded.')
    (path/'viewer.html').write_text(template.replace('__SIMULATION_DATA__', payload))
    report = {'case': name, 'config': asdict(config), 'history': history,
              'elapsed_seconds': time.monotonic()-started,
              'max_boundary_occupancy': max(row['boundary_occupancy'] for row in history),
              'max_cell_volume_error': max(row['max_cell_volume_error'] for row in history),
              'max_clipped_fraction': max(row['clipped_fraction'] for row in history),
              'min_radius_grid_cells': min(row['min_radius_grid_cells'] for row in history)}
    write_json(path/'analysis.json', report)
    return report


def compare(output):
    output = Path(output)
    protocol = json.loads((output/'protocol.json').read_text())
    reports = {name: json.loads((output/name/'analysis.json').read_text()) for name in protocol['cases']}
    pairs = {}
    for family in ('space', 'time'):
        names = protocol[f'{family}_cases']
        rows = []
        for left, right in zip(names[:-1], names[1:]):
            a, b = (Simulation.restore(output/name/'final_state.npz') for name in (left, right))
            if not np.array_equal(a.ids, b.ids) or not np.isclose(a.time, b.time):
                raise ValueError('comparison requires matched identities and physical times')
            # No projection is needed for a temporal comparison on the same mesh.
            if a.config.grid == b.config.grid:
                pa, pb = occupancy(a.phi), occupancy(b.phi)
            else:
                pa, pb = project_occupancy(a, protocol['grids'][0]), project_occupancy(b, protocol['grids'][0])
            ra, rb = reports[left]['history'][-1], reports[right]['history'][-1]
            rows.append({'cases': [left, right], 'field_relative_l2': float(np.linalg.norm(pa-pb)/np.linalg.norm(pb)),
                         'activator_relative_l2': float(np.linalg.norm(a.activator-b.activator)/np.linalg.norm(b.activator)),
                         'inhibitor_relative_l2': float(np.linalg.norm(a.inhibitor-b.inhibitor)/np.linalg.norm(b.inhibitor)),
                         'axis_ratio_relative_difference': abs(ra['axis_ratio']/rb['axis_ratio']-1),
                         'max_cell_volume_relative_difference': float(np.max(abs(a.volumes()/b.volumes()-1)))})
        pairs[family] = rows
    checks = {}
    limits = protocol['criteria']
    for family, rows in pairs.items():
        fine = rows[-1]
        checks[family+'_field_difference_below_1_percent'] = fine['field_relative_l2'] < limits['finest_pair_field_relative_l2_max']
        checks[family+'_signals_difference_below_1_percent'] = max(fine['activator_relative_l2'],fine['inhibitor_relative_l2']) < limits['finest_pair_signal_relative_l2_max']
        checks[family+'_shape_difference_below_1_percent'] = fine['axis_ratio_relative_difference'] < limits['finest_pair_axis_ratio_relative_difference_max']
        checks[family+'_field_difference_decreases'] = fine['field_relative_l2'] < rows[0]['field_relative_l2']
    finest = reports[protocol['space_cases'][-1]]
    checks['finest_cells_resolved'] = finest['min_radius_grid_cells'] >= limits['minimum_radius_grid_cells']
    checks['all_domains_clear'] = all(r['max_boundary_occupancy'] < limits['boundary_occupancy_max'] for r in reports.values())
    checks['all_cell_volumes_within_5_percent'] = all(r['max_cell_volume_error'] < limits['cell_volume_error_max'] for r in reports.values())
    checks['no_phase_field_clipping'] = all(r['max_clipped_fraction'] == 0 for r in reports.values())
    result = {'protocol': protocol, 'pairs': pairs, 'checks': checks, 'passed': all(checks.values()),
              'interpretation': 'A finite-duration manufactured-state screen, not a convergence proof. Spatial field comparisons include trilinear projection error. Unequal mesh ratios do not support a log2 order estimate. Cleavage, long-time pattern selection, and interface-width convergence remain separate.'}
    write_json(output/'comparison.json', result)
    lines = ['# Coupled resolution screen', '', f"All declared checks pass: **{result['passed']}**", '',
             '| Family | Cases | Occupancy L2 difference | Activator L2 difference | Axis-ratio difference |',
             '|---|---|---:|---:|---:|']
    for family, rows in pairs.items():
        for row in rows:
            lines.append(f"| {family} | {' / '.join(row['cases'])} | {row['field_relative_l2']:.6g} | {row['activator_relative_l2']:.6g} | {row['axis_ratio_relative_difference']:.6g} |")
    lines += ['', *[f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in checks.items()], '', result['interpretation']]
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    return result


def run(output, after=None, poll_interval=30.):
    output = Path(output)
    protocol = json.loads((output/'protocol.json').read_text())
    status = json.loads((output/'status.json').read_text())
    if status['state'] != 'prepared':
        raise ValueError('study is already queued, running, or finished')
    status = {'state': 'queued' if after else 'running', 'completed_cases': [],
              'dependency': str(Path(after).resolve()) if after else None}
    write_json(output/'status.json', status)
    try:
        if after:
            dependency = Path(after)
            while True:
                try:
                    upstream = json.loads((dependency/'status.json').read_text())
                except (OSError, ValueError):
                    time.sleep(poll_interval)
                    continue
                if upstream['state'] == 'failed':
                    raise RuntimeError('domain study failed; resolution jobs were not started')
                if upstream['state'] == 'completed':
                    blob = (dependency/'comparison.json').read_bytes()
                    report = json.loads(blob)
                    status['dependency_comparison_sha256'] = hashlib.sha256(blob).hexdigest()
                    if not (report['largest_domain_boundary_screen_pass'] and report['largest_pair_shape_screen_pass']):
                        status.update(state='blocked', reason='domain screens failed; review domain size before proceeding')
                        write_json(output/'status.json', status)
                        return
                    break
                time.sleep(poll_interval)
        for name in protocol['cases']:
            status.update(state='running', current_case=name)
            write_json(output/'status.json', status)
            run_case(output, name)
            status['completed_cases'].append(name)
        result = compare(output)
        status.update(state='completed', all_checks_pass=result['passed'], report='RESULTS.md')
        status.pop('current_case', None)
        write_json(output/'status.json', status)
    except Exception as error:
        status.update(state='failed', error=str(error))
        write_json(output/'status.json', status)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prepare_parser = sub.add_parser('prepare')
    prepare_parser.add_argument('--output', type=Path, required=True)
    runner = sub.add_parser('run')
    runner.add_argument('--output', type=Path, required=True)
    runner.add_argument('--after', type=Path)
    args = vars(parser.parse_args())
    command = args.pop('command')
    if command == 'prepare': prepare(**args)
    else: run(**args)


if __name__ == '__main__':
    main()
