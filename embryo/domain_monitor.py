"""Finish an independently running domain study and publish its status/report.

The monitor never changes simulation state. Worker PIDs are optional; when
provided, a vanished worker without final outputs records an explicit failure.
"""
import argparse
import json
import os
from pathlib import Path
import time

from .domain import compare, write_json


def result_markdown(report):
    protocol = report['protocol']
    lines = ['# Conservative domain continuation results', '',
             f"Completed through model time {protocol['until']:g} at fixed voxel spacing {protocol['dx']:g}.", '',
             '| Grid | Half-width | Peak boundary occupancy | Final axis ratio | First boundary crossing |',
             '|---:|---:|---:|---:|---:|']
    for row in report['runs']:
        crossing = row.get('first_boundary_threshold_crossing')
        lines.append(f"| {row['grid']}³ | {row['extent']:.2f} | {row['max_boundary_occupancy']:.6g} | {row['final']['axis_ratio']:.8f} | {crossing if crossing is not None else 'None'} |")
    pair = report['pairs'][-1]
    lines += ['', f"Largest-domain boundary screen: **{'PASS' if report['largest_domain_boundary_screen_pass'] else 'FAIL'}**.",
              f"Largest-pair shape screen: **{'PASS' if report['largest_pair_shape_screen_pass'] else 'FAIL'}**.", '',
              f"Maximum sampled relative axis-ratio difference, largest pair: {pair['max_axis_ratio_relative_difference']:.6g}.",
              f"Final shared-box phase-field relative L2 difference: {pair['final_common_box_phase_field_relative_l2']}.",
              f"Final maximum absolute activator difference: {pair['final_activator_max_absolute_difference']}.", '',
              'See comparison.json for full measurements and comparison.png for the histories. Each grid folder contains its full mesh playback and checkpoint.', '',
              'These are single-seed, finite-window domain tests. Passing the shape screen does not certify signaling or full-field convergence. Voxel spacing and the previously identified small-cell-resolution limitation are unchanged. No feedback-specific shape claim follows from this comparison.']
    return '\n'.join(lines) + '\n'


def monitor(output, workers=None, interval=10.):
    output = Path(output)
    protocol = json.loads((output / 'protocol.json').read_text())
    workers = workers or {}
    while True:
        completed, progress, missing = [], {}, []
        for grid in protocol['grids']:
            path = output / f'grid-{grid}'
            if (path / 'analysis.json').exists():
                completed.append(grid)
            else:
                pid = workers.get(grid)
                if pid:
                    try:
                        os.kill(pid, 0)
                        state = Path(f'/proc/{pid}/stat').read_text().split(') ', 1)[1][0]
                        if state == 'Z':
                            missing.append(grid)
                    except (ProcessLookupError, FileNotFoundError):
                        missing.append(grid)
            try:
                rows = json.loads((path / 'history.json').read_text())
                progress[str(grid)] = {key: rows[-1][key] for key in
                                      ('time', 'boundary_occupancy', 'axis_ratio', 'activator_std')}
            except (OSError, ValueError, IndexError):
                pass  # A worker may be replacing a progress file.
        status = {'state': 'running', 'target_time': protocol['until'],
                  'completed_grids': completed, 'progress': progress}
        if missing:
            status.update(state='failed', missing_workers=missing,
                          message='A worker exited without final analysis. Preserve its progress checkpoint and inspect the worker output.')
            write_json(output / 'status.json', status)
            raise RuntimeError(status['message'])
        if len(completed) == len(protocol['grids']):
            try:
                report = compare(output)
                (output / 'RESULTS.md').write_text(result_markdown(report))
            except Exception as error:
                status.update(state='failed', message=f'Final comparison failed: {error}')
                write_json(output / 'status.json', status)
                raise
            status.update(state='completed', report='RESULTS.md', comparison='comparison.json')
            write_json(output / 'status.json', status)
            return report
        write_json(output / 'status.json', status)
        time.sleep(interval)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', nargs='*', default=[], help='Optional grid:pid pairs')
    args = parser.parse_args()
    workers = {int(item.split(':')[0]): int(item.split(':')[1]) for item in args.workers}
    monitor(args.output, workers)


if __name__ == '__main__':
    main()
