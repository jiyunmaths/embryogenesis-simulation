"""Assess and plot the completed frozen neighbor-context screen.

Extra stability checks on unrecovered pulses are explicitly exploratory. The
original protocol, trajectories, and aggregate comparison remain unchanged.
"""
import argparse
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np

from .attribute_persistence import rhs, distance
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .neighbor_context import assess, job_name, read
from .resolution import write_json


def summarize(root, figure=None):
    root = Path(root).resolve(); p = read(root/'protocol.json'); saved = read(root/'results.json')
    if not saved.get('completed') or not saved.get('numerical_pass'):
        raise ValueError('Requires completed numerical acceptance')
    # Recompute all original classifications/metrics in a temporary path.
    with tempfile.TemporaryDirectory(prefix='neighbor-context-summary-') as name:
        check = Path(name); shutil.copy2(root/'protocol.json', check/'protocol.json')
        for context in p['contexts']:
            folder = check/context['key']; folder.mkdir()
            for job in context['jobs']:
                (folder/job_name(job)).symlink_to(root/context['key']/job_name(job), target_is_directory=True)
        if assess(check) != saved:
            raise ValueError('Reassessed outcomes differ from saved evidence')
    rows = [r for c in saved['contexts'] for r in c['comparisons']]
    records, unrecovered, baselines = [], [], []
    inputs = [root/'protocol.json', root/'results.json', root/'status.json']
    for context in p['contexts']:
        with np.load(context['source']) as source:
            delta, masses, ids = [source[k].copy() for k in ('delta', 'masses', 'ids')]
        fun = rhs(delta, context['beta'], context['da'], context['db'])
        for job in context['jobs']:
            folder = root/context['key']/job_name(job); record = read(folder/'result.json')
            records.append(record); inputs += [folder/'result.json', folder/'trajectories.npz']
            if job['arm'] == 'untouched': baselines += record['trials']
            if job['arm'] not in p['arms'][2:]: continue
            for f, trial in enumerate(record['trials']):
                if trial['target_recovery_time'] is not None and trial['network_recovery_time'] is not None: continue
                with np.load(folder/'trajectories.npz') as trajectories:
                    end = trajectories['pulse_paths'][f, -1]; control = trajectories['control'][-1]
                growth = float(np.linalg.eigvals(chemical_jacobian(end, delta, context['beta'], context['da'], context['db'])).real.max())
                residual = float(abs(fun(p['times'][-1], end.ravel())).max())
                i = list(ids).index(job['target'])
                unrecovered.append(dict(context=context['key'], **job, factor=trial['factor'],
                    target_recovery_time=trial['target_recovery_time'], network_recovery_time=trial['network_recovery_time'],
                    endpoint_rhs_max=residual, endpoint_jacobian_max_real=growth,
                    endpoint_stationary_stable=bool(residual < p['criteria']['stationarity_rhs_max'] and growth < 0),
                    network_endpoint_log_rms=float(distance(end, control, masses)),
                    target_control_endpoint=control[:, i].tolist(), target_pulse_endpoint=end[:, i].tolist()))
    challenged = [record for record in records if record['job']['arm'] in p['arms'][2:]]
    counts = dict(arm_jobs=len(records), trajectories=3*len(records), independent_histories=3,
        target_background_contexts=12, challenges=len(rows), informative=sum(r['informative'] for r in rows),
        state_shifted=sum(r['state_outcome']=='shifted' for r in rows),
        state_near_baseline=sum(r['state_outcome']=='near_baseline' for r in rows),
        pulse_comparisons=2*len(rows), response_shifted=sum(t['response_shift'] for r in rows for t in r['responses']),
        nearest_own=sum(t['nearest_baseline_response']=='own' for r in rows for t in r['responses']),
        target_unrecovered=sum(t['target_recovery_time'] is None for record in challenged for t in record['trials']),
        network_unrecovered=sum(t['network_recovery_time'] is None for record in challenged for t in record['trials']),
        baseline_target_unrecovered=sum(t['target_recovery_time'] is None for t in baselines),
        baseline_network_unrecovered=sum(t['network_recovery_time'] is None for t in baselines))
    by_arm = {arm: dict(challenges=sum(r['arm']==arm for r in rows),
        state_shifted=sum(r['arm']==arm and r['state_outcome']=='shifted' for r in rows),
        response_shifted=sum(t['response_shift'] for r in rows if r['arm']==arm for t in r['responses'])) for arm in p['arms'][2:]}
    report = dict(completed=True, numerical_pass=True, counts=counts, by_arm=by_arm,
        max_solver_log_error=max(record['max_solver_log_error'] for record in records),
        max_sham_error=max(c['sham_max_error'] for c in saved['contexts']),
        control_stationary_stable=sum(record['stationary_stable'] for record in records),
        exploratory_unrecovered_endpoints=unrecovered,
        exploratory_scope='Pulse-endpoint stability was inspected after observing nonrecovery. Pulses occur during neighborhood reorganization, not at settled equilibria. Distinct stable chemical endpoints support basin selection on the frozen graph, not instability of a settled endpoint, permanent memory, or cell identity.',
        protocol_sha256=digest(root/'protocol.json'), source_sha256={str(Path(__file__).resolve()):digest(__file__)},
        evidence_sha256={str(f):digest(f) for f in inputs}, limits=p['limits'])
    write_json(root/'summary.json', report)
    if figure is not None:
        plot(root, p, saved, Path(figure))
    return report


def plot(root, p, saved, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    state, response, labels = [], [], []
    for context in saved['contexts']:
        for target in sorted({r['target'] for r in context['comparisons']}):
            rows = [next(r for r in context['comparisons'] if r['target']==target and r['arm']==arm) for arm in p['arms'][2:]]
            state.append([r['late_state_ratio'] for r in rows])
            response.append([max(t['own_response_distance'] for t in r['responses']) for r in rows])
            family = 'fresh' if context['background']=='fresh_exchange' else 'untouched'
            labels.append(f"Seed {context['seed']} / {family} / cell {target}")
    specifications = [(np.array(state), 'A  Late target chemistry', 'Distance / initial pair separation', .1),
                      (np.array(response), 'B  Change in pulse response', 'RMS / initial log pulse', .01)]
    for ax, (values, title, colorlabel, threshold) in zip(axes[0], specifications):
        im = ax.imshow(values, cmap='magma', aspect='auto', norm=LogNorm(vmin=1e-4, vmax=max(values.max(), .01)))
        ax.set_yticks(range(len(labels)), labels, fontsize=8)
        ax.set_xticks(range(3), ['Conservative mix', 'Low reset', 'High reset'], fontsize=9)
        ax.set_title(title+f'\nDeclared change threshold: {threshold:g}', fontsize=11)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                value = values[i, j]; label = '<1e-4' if value < 1e-4 else f'{value:.3g}'
                ax.text(j, i, label, ha='center', va='center', fontsize=8,
                        color='black' if value > max(values.max(), .01)*.15 else 'white')
        fig.colorbar(im, ax=ax, shrink=.85, label=colorlabel)
    context = next(c for c in p['contexts'] if c['key']=='seed-8_unexchanged')
    with np.load(context['source']) as source: ids = source['ids']
    for ax, target, arm, title in [(axes[1, 0], 19, 'neighbor_low', 'C  Low-neighbor reset changes basin selection'),
                                  (axes[1, 1], 20, 'neighbor_mix', 'D  Conservative mixing changes recovery')]:
        index = list(ids).index(target)
        with np.load(root/context['key']/f'cell-{target}_untouched/trajectories.npz') as baseline:
            ax.plot(baseline['times'], baseline['control'][:, 0, index], color='#2864b3', label='Untouched control')
        with np.load(root/context['key']/f'cell-{target}_{arm}/trajectories.npz') as data:
            ax.plot(data['times'], data['control'][:, 0, index], color='#dd8c15', label='Changed neighbors: control')
            ax.plot(data['times'], data['pulse_paths'][0, :, 0, index], color='#b7314c', linestyle='--', label='Changed neighbors: -10% pulse')
            ax.plot(data['times'], data['pulse_paths'][1, :, 0, index], color='#348757', linestyle=':', label='Changed neighbors: +10% pulse')
        ax.set(title=title+f'\nSeed 8, untouched start, target {target}', xlabel='Elapsed chemical time', ylabel='Target activator')
        if target==19: ax.set_yscale('log')
        ax.legend(fontsize=8); ax.grid(alpha=.15)
    fig.suptitle('Neighbor chemistry changes state and response on fixed geometry\nTarget initial chemistry is held constant across neighbor arms; all cells then evolve freely', fontsize=13)
    filename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, dpi=170); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('outputs/neighbor-context'))
    parser.add_argument('--figure', type=Path, default=Path('docs/images/neighbor-context.png'))
    args = parser.parse_args()
    result = summarize(args.output, args.figure)
    print(json.dumps(dict(counts=result['counts'], by_arm=result['by_arm'],
                         unrecovered=result['exploratory_unrecovered_endpoints']), indent=2))
