"""Audit interpolation contamination in an archived voxel-refinement comparison.

No solver states, thresholds, or archived results are changed. Known analytic
initial fields calibrate diagnostic projection; evolved fields have no exact
solution, so interpolation agreement is evidence rather than an error bound.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from .model import Config, Simulation, occupancy
from .resolution import manufactured, write_json


def norm(field):
    return float(np.sqrt(np.sum(np.square(field), dtype=np.float64)))


def relative_error(a, b):
    denominator = norm(b)
    if denominator <= 0:
        raise ValueError('reference field must have positive L2 norm')
    return norm(a-b) / denominator


def project(simulation, probe_grid, order):
    """No clipping: expose spline overshoot instead of hiding diagnostic distortion."""
    if order not in (1, 3, 5):
        raise ValueError('supported interpolation orders are 1, 3, 5')
    if isinstance(probe_grid, bool) or not isinstance(probe_grid, int) or probe_grid < 2:
        raise ValueError('probe grid must be an integer >= 2')
    if probe_grid == simulation.config.grid:
        return occupancy(simulation.phi)
    extent = simulation.config.extent
    axis = (np.arange(probe_grid)+.5) * 2*extent/probe_grid - extent
    xyz = np.stack(np.meshgrid(*([axis]*3), indexing='ij'))
    indices = (xyz+extent)/simulation.dx - .5
    return np.stack([map_coordinates(occupancy(field), indices, order=order,
                                    mode='nearest', prefilter=order>1)
                     for field in simulation.phi])


def run(source, output, probes=(56,88,112), orders=(1,3,5)):
    source, output = Path(source).resolve(), Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh audit directory')
    original_bytes = (source/'comparison.json').read_bytes()
    comparison = json.loads(original_bytes)
    protocol = json.loads((source/'protocol.json').read_text())
    names = protocol['space_cases']
    if len(names) < 2 or len(set(probes)) != len(probes) or len(set(orders)) != len(orders):
        raise ValueError('need distinct probes/orders and at least two spatial cases')
    if any(order not in (1,3,5) for order in orders) or any(
            isinstance(n,bool) or not isinstance(n,int) or n<12 for n in probes):
        raise ValueError('invalid interpolation order or probe grid')
    states = {}
    provenance = {}
    for name in names:
        config = Config(**protocol['cases'][name])
        path = source/name/'final_state.npz'
        final = Simulation.restore(path)
        if final.config.grid != config.grid or not np.isclose(final.config.extent, protocol['extent']):
            raise ValueError('checkpoint geometry does not match source protocol')
        states[name] = (manufactured(config), final)
        provenance[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    output.mkdir(parents=True)
    audit_protocol = {'source':str(source), 'original_comparison_sha256':hashlib.sha256(original_bytes).hexdigest(),
                      'checkpoint_sha256':provenance, 'probe_grids':list(probes), 'orders':list(orders),
                      'original_field_tolerance':protocol['criteria']['finest_pair_field_relative_l2_max'],
                      'scope':__doc__}
    write_json(output/'protocol.json', audit_protocol)
    write_json(output/'status.json', {'state':'running'})
    rows, exact_errors = [], []
    for probe in probes:
        exact_config = Config(**{**protocol['cases'][names[0]], 'grid':probe})
        exact = occupancy(manufactured(exact_config).phi)
        for order in orders:
            projected = {}
            for name, (initial, final) in states.items():
                projected[name] = (project(initial,probe,order), project(final,probe,order))
                first, last = projected[name]
                exact_errors.append({'case':name, 'probe_grid':probe, 'order':order,
                                     'initial_relative_error_to_analytic':relative_error(first,exact),
                                     'initial_range':[float(first.min()),float(first.max())],
                                     'final_range':[float(last.min()),float(last.max())]})
            for left, right in zip(names[:-1], names[1:]):
                a0,a1 = projected[left]
                b0,b1 = projected[right]
                initial_error, final_error = relative_error(a0,b0), relative_error(a1,b1)
                difference0, difference1 = a0-b0, a1-b1
                denominator = norm(difference0)*norm(difference1)
                correlation = float(np.sum(difference0*difference1, dtype=np.float64)/denominator) if denominator>0 else None
                rows.append({'cases':[left,right], 'probe_grid':probe, 'order':order,
                             'initial_pair_relative_l2':initial_error, 'final_pair_relative_l2':final_error,
                             'initial_final_error_alignment':correlation,
                             'evolution_increment_difference_over_final_norm':norm(difference1-difference0)/norm(b1)})
            print(f'projection audit: probe={probe}, order={order}, finest initial={rows[-1]["initial_pair_relative_l2"]:.4%}, final={rows[-1]["final_pair_relative_l2"]:.4%}', flush=True)
            del projected
    finest_names = names[-2:]
    finest = [row for row in rows if row['cases']==finest_names]
    high = [row for row in finest if row['order']>1]
    tolerance = audit_protocol['original_field_tolerance']
    # Descriptive diagnostics, not replacement acceptance of the original test.
    summary = {'original_screen_passed':comparison['passed'],
               'original_finest_field_relative_l2':comparison['pairs']['space'][-1]['field_relative_l2'],
               'higher_order_final_range': [min(row['final_pair_relative_l2'] for row in high), max(row['final_pair_relative_l2'] for row in high)] if high else None,
               'all_tested_higher_order_pairs_below_original_tolerance': all(row['final_pair_relative_l2']<tolerance for row in high) if high else None,
               'largest_higher_order_initial_analytic_error_finest_pair':max((row['initial_relative_error_to_analytic'] for row in exact_errors if row['case'] in finest_names and row['order']>1),default=None),
               'interpretation':'Exploratory measurement audit. Initial interpolation error cannot be subtracted as a scalar correction to final error. Agreement among higher-order probes is not a rigorous bound on evolved-field error. Original pass/fail decisions are retained; an independently specified follow-up is needed.'}
    result = {'protocol':audit_protocol, 'summary':summary, 'pairs':rows, 'analytic_calibration':exact_errors}
    write_json(output/'analysis.json',result)
    lines=['# Projection-error audit', '', 'The archived resolution result is unchanged.', '',
           '| Probe grid | Order | Initial finest-pair L2 difference | Final finest-pair L2 difference |',
           '|---:|---:|---:|---:|']
    for row in finest:
        lines.append(f"| {row['probe_grid']}³ | {row['order']} | {row['initial_pair_relative_l2']:.4%} | {row['final_pair_relative_l2']:.4%} |")
    lines += ['', summary['interpretation']]
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for order in orders:
        selected=[row for row in finest if row['order']==order]
        for ax, key in zip(axes,('initial_pair_relative_l2','final_pair_relative_l2')):
            ax.plot([row['probe_grid'] for row in selected], [100*row[key] for row in selected], 'o-', label=f'order {order}')
            ax.axhline(100*tolerance, color='gray', ls=':', lw=.5)
            ax.set(xlabel='Probe grid per axis',ylabel='Finest-pair occupancy difference (%)')
    axes[0].set_title('Identical analytic initial geometry')
    axes[1].set_title('Archived evolved geometry')
    axes[0].legend()
    fig.suptitle('Diagnostic interpolation changes the measured discrepancy')
    fig.savefig(output/'projection.png',dpi=150)
    plt.close(fig)
    write_json(output/'status.json',{'state':'completed','report':'RESULTS.md','original_comparison_unchanged':hashlib.sha256((source/'comparison.json').read_bytes()).hexdigest()==audit_protocol['original_comparison_sha256']})
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.source,args.output)


if __name__=='__main__':main()
