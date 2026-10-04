"""Read-only three-timestep assessment, with post hoc timing diagnostics.

Keep the original raw-error gates. Time alignment is descriptive only.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify, frozen_spectrum, metrics
from .attribute_persistence import rhs
from .feedback_endpoint_bistability import chemical_jacobian
from .polarity_robustness import first_crossing
from .polarity_robustness_assessment import exceeded_intervals
from .polarity_robustness_refinement import parent_evidence, compare, check_halving
from .attribute_development import AttributeSimulation
from .resolution import write_json


def derive(root):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    state = read(root/'status.json')
    if state['state'] != 'completed_with_unresolved_checks' or state['completed'] != 1:
        raise ValueError('Requires the finished targeted refinement with retained failure')
    original = AttributeSimulation.restore(p['reference_job']['source'])
    refined = AttributeSimulation.restore(p['jobs'][0]['source']); check_halving(original, refined)
    parent = Path(p['parent']); old = read(parent/'protocol.json'); verify(old)
    jobs = [next(j for j in old['jobs'] if j['seed'] == 9 and j['family'] == 'uniform' and
        j['point']['chi'] == 0. and j['level'] == level) for level in ('coarse', 'fine')]+p['jobs']
    histories = []; evidence = {}; spectra = 0; spectrum_error = 0.
    runs = []
    for folder, protocol, job in zip((parent, parent, root), (old, old, p), jobs):
        result, history, files = parent_evidence(folder, protocol, job)
        evidence.update({str(f.resolve()): digest(f) for f in files})
        for row in history:
            s = frozen_spectrum(row['delta'], row['volumes'], 2., .02, 27.5)
            error = max(abs(s['uniform_jacobian_max_real']-row['uniform_growth_max']),
                float(abs(np.array(s['lambdas'])-row['laplacian_lambdas']).max()))
            if error > 1e-12 or s['unstable_modes'] != row['unstable_modes']:
                raise ValueError('Changed stored spectrum')
            if abs(float(np.std(np.log(np.array(row['chemistry'])[0])))-row['log_activator_sd']) > 1e-12:
                raise ValueError('Changed chemical contrast')
            spectrum_error = max(spectrum_error, error); spectra += 1
        t = np.array([r['elapsed'] for r in history]); sd = np.array([r['log_activator_sd'] for r in history])
        late = float(sd[t >= 216.-1e-9].min())
        if late != result['late_min_log_activator_sd'] or (late > .1) != result['persistent_contrast']:
            raise ValueError('Changed sustained-contrast classification')
        runs.append(dict(level=job['level'], dt=job['dt'], final_sd=float(sd[-1]), late_min_sd=late,
            onset=first_crossing(t, sd, .1, 'up'), growth_crossing=first_crossing(t,
                [r['uniform_growth_max'] for r in history]), audit=result['audit']))
        histories.append(history)
    x = [np.array([r['chemistry'] for r in h]) for h in histories]; pairs = []
    for i, j in ((0, 1), (1, 2), (0, 2)):
        e = abs(np.log(x[i]/x[j])); q = np.unravel_index(e.argmax(), e.shape)
        maximum = e.max(axis=(1, 2)); window = (t >= 75.) & (t <= 125.)
        logs = np.log(x[i]).reshape(len(t), -1); target = np.log(x[j][window]).reshape(window.sum(), -1)
        def shifted(shift):
            return np.stack([np.interp(t[window]+shift, t, col) for col in logs.T], axis=1)
        fit = minimize_scalar(lambda shift: np.mean((shifted(shift)-target)**2),
            bounds=(-.3, .3), method='bounded', options={'xatol': 1e-10})
        if not fit.success: raise ValueError('Timing diagnostic fit failed')
        aligned = abs(shifted(fit.x)-target)
        pairs.append(dict(levels=[runs[i]['level'], runs[j]['level']], timesteps=[runs[i]['dt'], runs[j]['dt']],
            maximum_log_error=float(e[q]), peak_elapsed=float(t[q[0]]), species=('activator', 'inhibitor')[q[1]],
            cell_id=histories[i][q[0]]['ids'][q[2]], concentrations=[float(x[i][q]), float(x[j][q])],
            concentration_ratio_discrepancy=float(np.expm1(e[q])), late_error=float(maximum[t >= 216.-1e-9].max()),
            samples_above_tolerance=int(np.sum(maximum > .01)), exceeded_intervals=exceeded_intervals(t, maximum, .01),
            error_series=maximum.tolist(), alignment=dict(window=[75., 125.], interpolated_log_states=True,
                fitted_shift=float(fit.x), raw_log_rms=float(np.sqrt(np.mean((logs[window]-target)**2))),
                aligned_log_rms=float(np.sqrt(fit.fun)), aligned_log_max=float(aligned.max()),
                aligned_max_series=aligned.max(axis=1).tolist(), times=t[window].tolist(),
                caveat='Post hoc fitted temporal offset; neither a changed acceptance criterion nor proof of the numerical cause.')))
    pilot = compare(root, p, p['pilot_duration']); long = compare(root, p, p['duration'])
    if not pilot['passed'] or long['passed'] or pairs[1]['maximum_log_error'] != long['errors']['chemical_log_max']:
        raise ValueError('Recorded refinement decision changed')
    endpoint_folder = root/jobs[-1]['key']/'endpoint'; ep = read(endpoint_folder/'protocol.json')
    endpoint = read(endpoint_folder/'assay/result.json')
    with np.load(endpoint_folder/'assay/paths.npz') as z:
        delta = z['delta']; masses = z['masses']; fun = rhs(delta, 2., .02, .55)
        jac = lambda time, y: chemical_jacobian(y.reshape(2, -1), delta, 2., .02, .55)
        for k, trial in enumerate(endpoint['trials']):
            computed = metrics(z[f'path_{k}'][-1], fun, jac, ep['criteria'])
            for field, value in computed.items():
                if isinstance(value, str):
                    if value != trial[field]: raise ValueError('Changed endpoint classification')
                elif not np.isclose(value, trial[field], rtol=1e-12, atol=1e-12):
                    raise ValueError('Changed endpoint diagnostic')
    for f in ('protocol.json', 'status.json', 'summary.json', 'pilot-refinement.json', 'long-refinement.json'):
        evidence[str(root/f)] = digest(root/f)
    summary = dict(assessed_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=digest(root/'protocol.json'),
        state=state, existing_history=9, new_histories=0, observations_recomputed=spectra,
        spectrum_recompute_max=spectrum_error, runs=runs, pairs=pairs, times=t.tolist(),
        pilot=pilot, long=long, endpoint=endpoint, evidence_sha256=evidence,
        source_sha256={**p['source_sha256'], str(Path(__file__).resolve()): digest(__file__),
            str(Path(__file__).with_name('polarity_robustness_assessment.py').resolve()):
                digest(Path(__file__).with_name('polarity_robustness_assessment.py'))},
        interpretation='All three steps initiate and maintain contrast, but raw formation-transient differences increase on further halving. Quantitative convergence is unresolved. A post hoc timing shift explains much of the mismatch without changing the gate. Float32 field/contact storage and split coupling are hypotheses, not established causes.')
    write_json(root/'assessment.json', summary)
    return summary, histories


def plot(summary, histories, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t = np.array(summary['times']); colors = ('#2878b5', '#e18d25', '#b84349')
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8), constrained_layout=True)
    for run, h, color in zip(summary['runs'], histories, colors):
        label = rf'$\Delta t={run["dt"]:g}$'
        axes[0, 0].plot(t, [r['log_activator_sd'] for r in h], color=color, label=label)
        cell = h[0]['ids'].index(27)
        axes[1, 0].plot(t, [r['chemistry'][0][cell] for r in h], color=color, label=label)
    axes[0, 0].set(title='All three timesteps form sustained contrast', yscale='log',
        ylabel='Across-cell SD of log activator', xlabel='Elapsed model time')
    axes[0, 0].axhline(.1, color='.4', ls=':'); axes[0, 0].legend(fontsize=8)
    for pair, color in zip(summary['pairs'][:2], colors):
        axes[0, 1].plot(t, pair['error_series'], color=color, label=' / '.join(pair['levels']))
    axes[0, 1].axhline(.01, color='black', ls='--', label='Original tolerance')
    axes[0, 1].set(title='Raw chemical error increases on further halving', yscale='log',
        ylim=(1e-8, .08), ylabel='Maximum absolute chemical log error', xlabel='Elapsed model time')
    axes[0, 1].legend(fontsize=8)
    axes[1, 0].set(title='Cell 27: same transition with shifted timing', xlim=(112, 123),
        yscale='log', ylabel='Activator activity', xlabel='Elapsed model time')
    axes[1, 0].legend(fontsize=8)
    pair = summary['pairs'][1]; fit = pair['alignment']; window = (t >= 75) & (t <= 125)
    axes[1, 1].plot(t[window], np.array(pair['error_series'])[window], label='Raw fine/finer discrepancy')
    axes[1, 1].plot(fit['times'], fit['aligned_max_series'], label=f'After fitted {fit["fitted_shift"]:.4f}-unit shift')
    axes[1, 1].set(title='Post hoc timing diagnostic; acceptance unchanged', yscale='log',
        ylabel='Maximum absolute chemical log error', xlabel='Elapsed model time')
    axes[1, 1].legend(fontsize=8)
    for ax in axes.flat: ax.grid(alpha=.15)
    fig.suptitle('History 9, zero directional tension: outcome repeats; quantitative transient remains unresolved', fontsize=13)
    destination = Path(destination); destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('outputs/polarity-robustness-refined'))
    parser.add_argument('--figure', type=Path, default=Path('docs/images/polarity-refinement-completed.png'))
    args = parser.parse_args(); result, histories = derive(args.output); plot(result, histories, args.figure)
    print(dict(full_refinement_pass=result['long']['passed'], endpoint_pass=result['endpoint']['numerical_pass'],
        observations_recomputed=result['observations_recomputed'], fine_finer_log_max=result['pairs'][1]['maximum_log_error']))
