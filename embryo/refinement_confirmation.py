"""Prospectively specified finer-grid confirmation of the coupled local screen."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .model import Config, Simulation, occupancy
from .projection_audit import project, relative_error
from .resolution import manufactured, run_case, write_json


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source, output, grid=112, probes=(88, 112, 128)):
    source, output = Path(source).resolve(), Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh confirmation directory')
    old = json.loads((source/'protocol.json').read_text())
    names = old['space_cases'][-2:]
    if grid <= old['cases'][names[-1]]['grid']:
        raise ValueError('confirmation must use a finer grid')
    if len(set(probes)) != 3 or any(type(n) is not int or n < 12 for n in probes):
        raise ValueError('three distinct integer probe grids >= 12 are required')
    configs = {name: old['cases'][name] for name in names}
    # Reject comparisons that change physics or the integration schedule.
    a, b = (dict(configs[name]) for name in names)
    a.pop('grid'); b.pop('grid')
    if a != b:
        raise ValueError('archived spatial cases differ in more than grid size')
    name = f'space-{grid}'
    config = Config(**{**configs[names[-1]], 'grid': grid})
    config.validate()
    configs[name] = asdict(config)
    paths = [source/'comparison.json', source/'protocol.json']
    for previous in names:
        paths += [source/previous/'final_state.npz', source/previous/'analysis.json']
    protocol = {
        'source': str(source), 'source_sha256': {str(p): digest(p) for p in paths},
        'space_cases': names+[name], 'cases': configs, 'new_case': name,
        'probe_grids': list(probes), 'orders': [3, 5],
        'criteria': {'field_relative_l2_max': .01, 'initial_analytic_error_max': .0025,
                     'diagnostic_spread_max': .001, 'signal_relative_l2_max': .01,
                     'axis_ratio_relative_difference_max': .01,
                     'boundary_occupancy_max': .01, 'cell_volume_error_max': .05,
                     'minimum_radius_grid_cells': 4.},
        'scope': 'Independent new finer-grid state; archived coarse states reused. Same finite-duration manufactured aggregate, no division. All six diagnostic variants must pass, with strictly decreasing pair error. Analytic calibration is not an evolved-error bound. Does not validate long-time development, sharp-interface convergence, or biological mechanisms.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json', protocol)
    write_json(output/'status.json', {'state': 'prepared', 'protocol_sha256': digest(output/'protocol.json')})
    return protocol


def evaluate(protocol, pairs, calibration, reports, signals, shape, sources_unchanged):
    limits = protocol['criteria']
    fine = [r for r in pairs if r['cases'] == protocol['space_cases'][-2:]]
    coarse = {(r['probe_grid'], r['order']): r for r in pairs if r['cases'] == protocol['space_cases'][:2]}
    values = [r['final_relative_l2'] for r in fine]
    return {
        'all_six_field_comparisons_below_1_percent': len(values) == 6 and max(values) < limits['field_relative_l2_max'],
        'field_difference_decreases_for_every_diagnostic': len(values) == 6 and all(r['final_relative_l2'] < coarse[r['probe_grid'], r['order']]['final_relative_l2'] for r in fine),
        'initial_reconstruction_below_0_25_percent': bool(calibration) and max(r['relative_l2'] for r in calibration) < limits['initial_analytic_error_max'],
        'diagnostic_spread_below_0_1_percentage_point': bool(values) and max(values)-min(values) < limits['diagnostic_spread_max'],
        'signals_below_1_percent': max(signals.values()) < limits['signal_relative_l2_max'],
        'shape_below_1_percent': shape < limits['axis_ratio_relative_difference_max'],
        'all_domains_clear': all(r['max_boundary_occupancy'] < limits['boundary_occupancy_max'] for r in reports),
        'all_volumes_within_5_percent': all(r['max_cell_volume_error'] < limits['cell_volume_error_max'] for r in reports),
        'no_phase_field_clipping': all(r['max_clipped_fraction'] == 0 for r in reports),
        'new_grid_cells_resolved': reports[-1]['min_radius_grid_cells'] >= limits['minimum_radius_grid_cells'],
        'archived_sources_unchanged': sources_unchanged}


def run(output):
    output = Path(output)
    protocol = json.loads((output/'protocol.json').read_text())
    status = json.loads((output/'status.json').read_text())
    if status['state'] != 'prepared' or digest(output/'protocol.json') != status['protocol_sha256']:
        raise ValueError('requires an unchanged prepared protocol')
    def unchanged():
        return all(digest(Path(p)) == value for p, value in protocol['source_sha256'].items())
    if not unchanged():
        raise ValueError('archived sources changed since preparation')
    try:
        write_json(output/'status.json', {**status, 'state': 'running', 'stage': 'evolution'})
        run_case(output, protocol['new_case'])
        write_json(output/'status.json', {**status, 'state': 'running', 'stage': 'comparison'})
        names = protocol['space_cases']
        paths = [Path(protocol['source'])/n if n != protocol['new_case'] else output/n for n in names]
        finals = [Simulation.restore(p/'final_state.npz') for p in paths]
        initials = [manufactured(Config(**protocol['cases'][n])) for n in names]
        if any(not np.array_equal(s.ids, finals[-1].ids) or not np.isclose(s.time, finals[-1].time) for s in finals):
            raise ValueError('unmatched identities or final times')
        reports = [json.loads((p/'analysis.json').read_text()) for p in paths]
        rows, calibration = [], []
        for probe in protocol['probe_grids']:
            exact = occupancy(manufactured(Config(**{**protocol['cases'][names[0]], 'grid': probe})).phi)
            for order in protocol['orders']:
                previous = None
                for name, initial, final in zip(names, initials, finals):
                    initial_field = project(initial, probe, order)
                    calibration.append({'case': name, 'probe_grid': probe, 'order': order,
                                        'relative_l2': relative_error(initial_field, exact)})
                    del initial_field
                    field = project(final, probe, order)
                    if previous is not None:
                        rows.append({'cases': [previous_name, name], 'probe_grid': probe, 'order': order,
                                     'final_relative_l2': relative_error(previous, field),
                                     'reconstructed_range': [float(field.min()), float(field.max())]})
                    previous, previous_name = field, name
                del previous, field
                print(f'probe={probe}, order={order}, finest difference={rows[-1]["final_relative_l2"]:.4%}', flush=True)
        a, b = finals[-2:]
        signals = {key: relative_error(getattr(a, key), getattr(b, key)) for key in ('activator', 'inhibitor')}
        shape = abs(reports[-2]['history'][-1]['axis_ratio']/reports[-1]['history'][-1]['axis_ratio']-1)
        checks = evaluate(protocol, rows, calibration, reports, signals, shape, unchanged())
        result = {'protocol': protocol, 'pairs': rows, 'calibration': calibration,
                  'signals': signals, 'axis_ratio_relative_difference': shape,
                  'checks': checks, 'passed': all(checks.values())}
        write_json(output/'comparison.json', result)
        lines = ['# Independent finer-grid confirmation', '', f'All declared checks pass: **{result["passed"]}**', '',
                 '| Probe | Order | Grid pair | Final occupancy difference |', '|---:|---:|---|---:|']
        lines += [f'| {r["probe_grid"]} | {r["order"]} | {" / ".join(r["cases"])} | {r["final_relative_l2"]:.4%} |' for r in rows]
        lines += ['', *[f'- {key}: {"PASS" if value else "FAIL"}' for key, value in checks.items()], '', protocol['scope']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json', {**status, 'state': 'completed', 'all_checks_pass': result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json', {**status, 'state': 'failed', 'error': str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run'])
    parser.add_argument('--source', type=Path, default=Path('outputs/coupled-resolution'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.source, args.output) if args.action == 'prepare' else run(args.output)


if __name__ == '__main__':
    main()
