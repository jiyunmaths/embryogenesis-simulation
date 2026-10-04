"""Delayed pulses at settled endpoints of the two neighbor-screen nonrecoveries.

Case follow-up on a frozen graph, not an independent developmental replication.
Compare fractional pulses and original absolute doses without reapplying the
neighbor intervention or restoring any cell's earlier concentrations.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np

from .attribute_persistence import rhs, distance
from .cell_exchange_response_moving import verify
from .cell_response import pulse, response_metrics
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .neighbor_context import assess as assess_original, job_name, load_result, read, validate_graph
from .resolution import write_json


CASES = (('seed-8_unexchanged', 19, 'neighbor_low', .9),
         ('seed-8_unexchanged', 20, 'neighbor_mix', .9))
KINDS = ('untouched', 'challenged', 'alternative')
CRITERIA = dict(solver_log_max=1e-5, stationarity_rhs_max=1e-6,
                control_drift_log_max=1e-7, endpoint_return_log_rms=1e-4,
                metric_recompute_rtol=1e-12, metric_recompute_atol=1e-14)


def pulse_conditions(current_a, original_a, volume):
    """Same fractional pulse OR same absolute original amount; record both."""
    if not all(np.isfinite(x) and x > 0 for x in (current_a, original_a, volume)):
        raise ValueError('Positive finite concentrations and volume required')
    conditions = []
    for dose in ('fractional', 'matched_amount'):
        for nominal in (.9, 1.1):
            factor = nominal if dose == 'fractional' else 1+(nominal-1)*original_a/current_a
            if not np.isfinite(factor) or factor <= 0 or factor == 1:
                raise ValueError('Nonzero positive delayed pulse required')
            conditions.append(dict(dose=dose, nominal_factor=nominal, factor=float(factor),
                injected_activator_amount=float((factor-1)*current_a*volume),
                original_activator_amount=float((nominal-1)*original_a*volume)))
    return conditions


def dynamics(job, arrays):
    fun = rhs(arrays['delta'], job['beta'], job['da'], job['db'])
    jac = lambda t, y: chemical_jacobian(y.reshape(2, -1), arrays['delta'], job['beta'], job['da'], job['db'])
    return fun, jac


def endpoint_metrics(state, fun, jac, criteria):
    residual = float(abs(fun(0., state.ravel())).max())
    growth = float(np.linalg.eigvals(jac(0., state.ravel())).real.max())
    if not np.isfinite(residual) or not np.isfinite(growth):
        raise ValueError('Nonfinite endpoint derivative or spectrum')
    return dict(rhs_max=residual, jacobian_max_real=growth,
                stationary_stable=bool(residual < criteria['stationarity_rhs_max'] and growth < 0))


def checked_snapshot(job):
    with np.load(job['source']) as d:
        arrays = {k: d[k].copy() for k in d.files}
    validate_graph(arrays['delta'], arrays['masses'])
    n = len(arrays['masses']); ids = arrays['ids']
    if (ids.shape != (n,) or len(set(ids.tolist())) != n or job['target'] not in ids or
            job['target_index'] != list(ids).index(job['target'])):
        raise ValueError('Invalid target/ID alignment')
    for key in ('initial', 'original_initial'):
        if arrays[key].shape != (2, n) or not np.isfinite(arrays[key]).all() or np.any(arrays[key] <= 0):
            raise ValueError('Invalid endpoint chemistry')
    return arrays


def prepare(root, source=Path('outputs/neighbor-context')):
    root, source = Path(root).resolve(), Path(source).resolve()
    if root.exists():
        raise FileExistsError(root)
    prior = read(source/'protocol.json'); verify(prior)
    status = read(source/'status.json'); saved = read(source/'results.json')
    if (status.get('state') != 'completed' or status.get('completed') != 60 or
            not status.get('numerical_pass') or not saved.get('numerical_pass')):
        raise ValueError('Requires the complete accepted neighbor screen')
    # Reconstruct all original outcomes in a temporary directory. Do not write
    # to the accepted screen or rely on its exploratory summary for selection.
    with tempfile.TemporaryDirectory(prefix='delayed-source-check-') as name:
        check = Path(name); shutil.copy2(source/'protocol.json', check/'protocol.json')
        for context in prior['contexts']:
            folder = check/context['key']; folder.mkdir()
            for job in context['jobs']:
                (folder/job_name(job)).symlink_to(source/context['key']/job_name(job), target_is_directory=True)
        if assess_original(check) != saved:
            raise ValueError('Original reassessment differs from saved outcomes')
    selected = []
    for context in prior['contexts']:
        for job in context['jobs']:
            if job['arm'] not in prior['arms'][2:]:
                continue
            record = read(source/context['key']/job_name(job)/'result.json')
            selected += [(context['key'], job['target'], job['arm'], trial['factor'])
                         for trial in record['trials']
                         if trial['target_recovery_time'] is None or trial['network_recovery_time'] is None]
    if set(selected) != set(CASES) or len(selected) != len(CASES):
        raise ValueError('Expected all and only the two original nonrecoveries')
    inputs = [source/f for f in ('protocol.json', 'status.json', 'results.json')]
    jobs, snapshots, cases = [], [], []
    for key, target, arm, factor in CASES:
        context = next(c for c in prior['contexts'] if c['key'] == key)
        challenged = source/key/job_name(dict(target=target, arm=arm))
        untouched = source/key/job_name(dict(target=target, arm='untouched'))
        record, paths = load_result(challenged, digest(source/'protocol.json'), prior['criteria'])
        _, base = load_result(untouched, digest(source/'protocol.json'), prior['criteria'])
        with np.load(context['source']) as d:
            graph = {k: d[k].copy() for k in ('delta', 'masses', 'ids')}
        i = list(graph['ids']).index(target); f = list(prior['factors']).index(factor)
        original = paths['initial'].copy()
        if not np.array_equal(original[:, i], base['initial'][:, i]):
            raise ValueError('Original intervention changed target chemistry')
        endpoints = dict(untouched=base['control'][-1], challenged=paths['control'][-1],
                         alternative=paths['pulse_paths'][f, -1])
        case = f'cell-{target}_{arm}'
        cases.append(dict(key=case, context=key, target=target, arm=arm,
            original_trial=record['trials'][f], original_factor=factor,
            original_target=original[:, i].tolist(), settled_target=endpoints['challenged'][:, i].tolist(),
            original_endpoint_separation=float(distance(endpoints['alternative'], endpoints['challenged'], graph['masses'])),
            immediate_folder=str(challenged)))
        for kind in KINDS:
            snapshot = dict(**graph, initial=endpoints[kind].copy(), original_initial=original.copy())
            job = dict(key=f'{case}_{kind}', case=case, kind=kind, target=target, target_index=i,
                       beta=context['beta'], da=context['da'], db=context['db'])
            fun, jac = dynamics(job, snapshot)
            stable = endpoint_metrics(snapshot['initial'], fun, jac, CRITERIA)
            if not stable['stationary_stable']:
                raise ValueError('Delayed pulses require stationary stable source endpoints')
            job['initial_stability'] = stable
            pulse_conditions(snapshot['initial'][0, i], original[0, i], graph['masses'][i])
            jobs.append(job); snapshots.append(snapshot)
        inputs.append(Path(context['source']))
        for folder in (challenged, untouched):
            inputs += [folder/'result.json', folder/'trajectories.npz']
    root.mkdir(parents=True)
    for job, snapshot in zip(jobs, snapshots):
        folder = root/job['key']; folder.mkdir()
        filename = folder/'source.npz'; np.savez_compressed(filename, **snapshot)
        job['source'] = str(filename); inputs.append(filename)
    sources = {**prior['source_sha256'], str(Path(__file__).resolve()): digest(__file__)}
    p = dict(cases=cases, jobs=jobs, times=prior['times'], criteria=CRITERIA,
        total_jobs=len(jobs), trajectories=5*len(jobs), independent_histories=1,
        frozen_settling_duration=prior['times'][-1], post_pulse_duration=prior['times'][-1],
        design='All two seed-8 immediate-pulse nonrecoveries, selected after the original screen. Continue each full-network untouched/control/negative-pulse endpoint at 240 frozen chemical units. No further neighbor intervention, cell reset, or geometry change. Paired unperturbed control and four target activator pulses: +/-10% of current concentration, and +/- the original 10% absolute amount. All cells evolve freely for a further 240 chemical units.',
        recovery='Existing sustained 10%-of-initial-displacement rule, at least 24 further observed units. Normalize by the actual log pulse, not its nominal percentage. Nonrecovery is censored at 240. Endpoint return requires stationary/local chemical stability and volume-weighted log RMS <1e-4.',
        interpretation='Primary: robustness of challenged settled equilibria, with percentage and amount-matched dose controls. Untouched endpoints are matched controls. Alternative endpoints are a planned secondary robustness test. Contrasting immediate and delayed outcomes tests dependence on network state at perturbation; it does not hold the full chemical state fixed or identify a single changing concentration.',
        limits='Case-selected follow-up in one history, not cross-history replication or cell-autonomous identity. Frozen delay is not a developmental age. Same-percent pulses have different amounts; matched-amount pulses have different percentages. Chemistry stability excludes moving mechanics, division, inheritance, transport-closure and full developmental-convergence claims.',
        source_sha256=sources, input_sha256={str(f.resolve()):digest(f) for f in inputs})
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=p['total_jobs']))
    return p


def trial_metrics(path, control, times, arrays, job, condition, criteria):
    fun, jac = dynamics(job, arrays)
    stability = endpoint_metrics(path[-1], fun, jac, criteria)
    metrics = response_metrics(path, control, times, arrays['masses'], job['target_index'], condition['factor'])
    classification = ('not_settled_by_horizon' if not stability['stationary_stable'] else
                      'returned' if metrics['final_log_rms_from_control'] < criteria['endpoint_return_log_rms'] else
                      'different_stable_endpoint')
    return dict(**condition, **metrics, **stability, classification=classification)


def equivalent(saved, expected, criteria):
    if set(saved) != set(expected):
        return False
    for key, value in expected.items():
        if value is None or isinstance(value, (str, bool)) or key.endswith('recovery_time'):
            if saved[key] != value:
                return False
        elif (not np.isfinite(saved[key]) or not np.isclose(saved[key], value,
                rtol=criteria['metric_recompute_rtol'], atol=criteria['metric_recompute_atol'])):
            return False
    return True


def load_checked(root, p, job):
    arrays = checked_snapshot(job); folder = root/job['key']; record = read(folder/'result.json')
    if (record.get('protocol_sha256') != digest(root/'protocol.json') or record.get('job') != job or
            record.get('trajectory_sha256') != digest(folder/'trajectories.npz') or not record.get('quality_pass') or
            not np.isfinite(record['max_solver_log_error']) or record['max_solver_log_error'] >= p['criteria']['solver_log_max']):
        raise ValueError('Invalid saved delayed-pulse evidence')
    with np.load(folder/'trajectories.npz') as d:
        paths = {k: d[k].copy() for k in d.files}
    times = np.asarray(p['times']); shape = (len(times), *arrays['initial'].shape)
    if (not np.array_equal(paths['times'], times) or paths['control'].shape != shape or
            paths['pulse_paths'].shape != (4, *shape) or paths['responses'].shape != (4, len(times), 2) or
            not np.array_equal(paths['initial'], arrays['initial'])):
        raise ValueError('Incomplete or misaligned delayed trajectories')
    for name in ('control', 'pulse_paths'):
        if not np.isfinite(paths[name]).all() or np.any(paths[name] <= 0):
            raise ValueError('Nonpositive or nonfinite delayed chemistry')
    if not np.array_equal(paths['control'][0], arrays['initial']):
        raise ValueError('Control did not start at the full settled endpoint')
    i = job['target_index']; initial = arrays['initial']
    conditions = pulse_conditions(initial[0, i], arrays['original_initial'][0, i], arrays['masses'][i])
    if len(record['trials']) != len(conditions):
        raise ValueError('Incomplete delayed pulse metrics')
    for f, condition in enumerate(conditions):
        if not np.array_equal(paths['pulse_paths'][f, 0], pulse(initial, i, condition['factor'])):
            raise ValueError('Delayed pulse changed its target, dose, or starting state')
        wave = np.log(paths['pulse_paths'][f, :, :, i]/paths['control'][:, :, i])/abs(np.log(condition['factor']))
        if not np.array_equal(paths['responses'][f], wave):
            raise ValueError('Invalid saved delayed waveform')
        expected = trial_metrics(paths['pulse_paths'][f], paths['control'], times, arrays, job, condition, p['criteria'])
        if not equivalent(record['trials'][f], expected, p['criteria']):
            raise ValueError('Saved recovery/endpoint metrics differ from trajectories')
    fun, jac = dynamics(job, arrays)
    stability = endpoint_metrics(initial, fun, jac, p['criteria'])
    drift = float(distance(paths['control'], initial, arrays['masses']).max())
    if (not equivalent(record['initial_stability'], stability, p['criteria']) or
            not stability['stationary_stable'] or not equivalent(job['initial_stability'], stability, p['criteria']) or
            drift > p['criteria']['control_drift_log_max'] or not np.isclose(drift, record['control_drift_log_max'],
                rtol=p['criteria']['metric_recompute_rtol'], atol=p['criteria']['metric_recompute_atol'])):
        raise ValueError('Unsettled or drifting delayed control')
    return record, paths


def worker(args):
    root, job = args; root = Path(root); p = read(root/'protocol.json'); verify(p)
    folder = root/job['key']
    if (folder/'result.json').exists():
        return load_checked(root, p, job)[0]
    arrays = checked_snapshot(job); initial = arrays['initial']; i = job['target_index']; times = np.asarray(p['times'])
    fun, jac = dynamics(job, arrays); started = time.monotonic()
    stability = endpoint_metrics(initial, fun, jac, p['criteria'])
    if not stability['stationary_stable']:
        raise ValueError('Delayed-pulse source is not a stable equilibrium')
    write_json(folder/'status.json', dict(state='running', phase='control'))
    control, error = checked_solve(fun, jac, initial, times)
    drift = float(distance(control, initial, arrays['masses']).max())
    if drift > p['criteria']['control_drift_log_max']:
        raise ValueError('Delayed control drift exceeds numerical criterion')
    paths, responses, trials, errors = [], [], [], [error]
    conditions = pulse_conditions(initial[0, i], arrays['original_initial'][0, i], arrays['masses'][i])
    for condition in conditions:
        write_json(folder/'status.json', dict(state='running', phase='pulse', **condition))
        path, error = checked_solve(fun, jac, pulse(initial, i, condition['factor']), times)
        paths.append(path); errors.append(error)
        responses.append(np.log(path[:, :, i]/control[:, :, i])/abs(np.log(condition['factor'])))
        trials.append(trial_metrics(path, control, times, arrays, job, condition, p['criteria']))
    filename = folder/'trajectories.npz'; temporary = folder/'trajectories.tmp.npz'
    np.savez_compressed(temporary, initial=initial, control=control, pulse_paths=np.asarray(paths),
                        responses=np.asarray(responses), times=times)
    temporary.replace(filename)
    record = dict(protocol_sha256=digest(root/'protocol.json'), trajectory_sha256=digest(filename),
        job=job, quality_pass=True, max_solver_log_error=max(errors), initial_stability=stability,
        control_drift_log_max=drift, trials=trials, wall_seconds=time.monotonic()-started)
    write_json(folder/'result.json', record)
    load_checked(root, p, job)
    write_json(folder/'status.json', dict(state='completed', quality_pass=True))
    return record


def assess(root, write=True):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    records = [load_checked(root, p, job)[0] for job in p['jobs']]
    counts = {}
    for kind in KINDS:
        trials = [t for r in records if r['job']['kind'] == kind for t in r['trials']]
        counts[kind] = dict(pulses=len(trials), target_recovered=sum(t['target_recovery_time'] is not None for t in trials),
            network_recovered=sum(t['network_recovery_time'] is not None for t in trials),
            returned=sum(t['classification'] == 'returned' for t in trials),
            different_stable_endpoint=sum(t['classification'] == 'different_stable_endpoint' for t in trials))
    result = dict(completed=True, numerical_pass=True, protocol_sha256=digest(root/'protocol.json'),
        arm_jobs=len(records), trajectories=p['trajectories'], independent_histories=p['independent_histories'],
        max_solver_log_error=max(r['max_solver_log_error'] for r in records),
        max_control_drift_log=max(r['control_drift_log_max'] for r in records), counts=counts,
        cases=p['cases'], records=records, interpretation=p['interpretation'], limits=p['limits'])
    if write:
        write_json(root/'results.json', result)
    return result


def run(root, workers=4):
    if not isinstance(workers, int) or workers < 1:
        raise ValueError('Positive integer worker count required')
    root = Path(root).resolve()
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); done = 0
        try:
            write_json(root/'status.json', dict(state='running', completed=0, total=p['total_jobs'], workers=workers))
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(worker, (root, job)):job['key'] for job in p['jobs']}
                for future in as_completed(futures):
                    future.result(); done += 1
                    write_json(root/'status.json', dict(state='running', completed=done, total=p['total_jobs'], workers=workers))
                    print(futures[future], 'completed', done, '/', p['total_jobs'], flush=True)
            result = assess(root)
            write_json(root/'status.json', dict(state='completed', completed=done, total=p['total_jobs'],
                                               numerical_pass=result['numerical_pass']))
        except Exception as exc:
            write_json(root/'status.json', dict(state='failed', completed=done, total=p['total_jobs'], error=str(exc)))
            raise


def plot(root, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    root = Path(root).resolve(); p = read(root/'protocol.json')
    if assess(root, write=False) != read(root/'results.json'):
        raise ValueError('Saved delayed assessment changed')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for row, case in enumerate(p['cases']):
        with np.load(Path(case['immediate_folder'])/'trajectories.npz') as d:
            original = {k:d[k].copy() for k in d.files}
        job = next(j for j in p['jobs'] if j['case'] == case['key'] and j['kind'] == 'challenged')
        record, delayed = load_checked(root, p, job); i = job['target_index']; arrays = checked_snapshot(job)
        immediate = distance(original['pulse_paths'][0], original['control'], arrays['masses'])
        ax = axes[row, 0]; ax.plot(original['times'], immediate, label='Immediate -10%', color='#bd3b32', lw=2)
        for f, label, color in ((0, 'Settled -10%', '#277da8'), (2, 'Settled, same removed amount', '#37955a')):
            ax.plot(delayed['times'], distance(delayed['pulse_paths'][f], delayed['control'], arrays['masses']),
                    label=label, color=color, lw=2, ls='--' if f == 2 else '-')
        ax.set_title(f"Cell {case['target']}: "+('low-neighbor reset' if case['arm'] == 'neighbor_low' else 'conservative neighbor mixing'))
        ax.set_ylabel('Network log RMS from paired control'); ax.set_xlabel('Chemical time since pulse')
        ax.legend(fontsize=8); ax.grid(alpha=.2)
        ax = axes[row, 1]
        ax.plot(original['times'], original['control'][:, 0, i], color='#30343b', label='Immediate control')
        ax.plot(original['times'], original['pulse_paths'][0, :, 0, i], color='#bd3b32', label='Immediate -10%')
        ax.plot(delayed['times'], delayed['control'][:, 0, i], color='#30343b', ls=':', label='Settled control')
        ax.plot(delayed['times'], delayed['pulse_paths'][0, :, 0, i], color='#277da8', label='Settled -10%')
        ax.plot(delayed['times'], delayed['pulse_paths'][2, :, 0, i], color='#37955a', ls='--', label='Settled, same amount')
        ax.set_yscale('log'); ax.set_ylabel('Target activator'); ax.set_xlabel('Chemical time since pulse')
        ax.set_title('Target chemical trajectories'); ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Selected seed-8 cases: immediate versus settled-state perturbations\nFrozen geometry; percentage and absolute-dose controls', fontsize=14)
    filename = Path(filename); filename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, dpi=180); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'assess', 'plot'))
    parser.add_argument('--output', type=Path, default=Path('outputs/neighbor-context-delayed'))
    parser.add_argument('--source', type=Path, default=Path('outputs/neighbor-context'))
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--figure', type=Path, default=Path('docs/images/neighbor-context-delayed.png'))
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output, args.source)
    elif args.command == 'run': run(args.output, args.workers)
    elif args.command == 'assess': assess(args.output)
    else: plot(args.output, args.figure)
