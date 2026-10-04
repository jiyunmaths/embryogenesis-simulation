"""Isolate production chemical stepping from moving geometry and dilution.

Freeze three measured history-9 contexts. Call the unchanged NumPy/SciPy and
PyTorch chemical routines, comparing to independent DOP853/Radau references.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from multiprocessing import get_context
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from .attribute_persistence import rhs
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .polarity_robustness import first_crossing, checked_history
from .resolution import write_json, _steps
from .transport import Transport, integrate_gm


CONTEXTS = ((0., 240.), (105., 60.), (117.75, 24.))
DTS = (.00375, .001875, .0009375)
CRITERIA = dict(reference_log_max=1e-8, production_log_max=.01,
    cpu_torch_log_max=1e-8, order_min=1.7, order_max=2.3,
    reference_noise_factor=100.)


def fixed_transport(delta, masses):
    delta, masses = validate_graph(delta, masses)
    g = masses[:, None]*delta; np.fill_diagonal(g, 0.)
    symmetric = -np.sqrt(masses)[:, None]*delta/np.sqrt(masses)[None, :]
    # Preserve the exact saved operator, rather than reconstructing/averaging it.
    return Transport(masses.copy(), sparse.csr_matrix(g), sparse.csr_matrix(delta), sparse.csr_matrix(symmetric))


def prepare(root, parent=Path('outputs/polarity-robustness-refined')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists(): raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    if read(parent/'status.json')['state'] != 'completed_with_unresolved_checks':
        raise ValueError('Requires the assessed targeted timestep failure')
    source = Path(old['reference_folder']); history = read(source/'history.json'); result = read(source/'result.json')
    checked_history(history, old, old['reference_job'], old['duration'])
    if result['history_sha256'] != digest(source/'history.json') or not result['quality_pass']:
        raise ValueError('Changed fine moving reference')
    root.mkdir(parents=True); contexts = []; jobs = []
    inputs = [parent/'protocol.json', parent/'status.json', parent/'long-refinement.json',
        source/'history.json', source/'result.json']
    for elapsed, duration in CONTEXTS:
        row = next(r for r in history if abs(r['elapsed']-elapsed) < 1e-10)
        delta, masses = validate_graph(row['delta'], row['volumes'])
        key = f'elapsed-{elapsed:g}'; file = root/(key+'.npz')
        np.savez_compressed(file, delta=delta, masses=masses, initial=np.array(row['chemistry']), ids=row['ids'])
        contexts.append(dict(key=key, elapsed=elapsed, duration=duration, source=str(file),
            description='Exact fine moving chemistry and measured graph/volumes at the stated elapsed time; thereafter freeze geometry with no dilution.'))
        inputs.append(file)
        for backend in ('native_chemistry', 'torch_cpu_chemistry'):
            for dt in DTS:
                jobs.append(dict(key=f'{key}_{backend}_dt-{dt:g}', context=key, backend=backend, dt=dt))
    files = [Path(__file__), *[Path(__file__).with_name(f) for f in ('transport.py', 'signaling.py',
        'gpu_backend.py', 'attribute_persistence.py', 'cell_response_exchange.py',
        'feedback_endpoint_bistability.py', 'neighbor_context.py', 'feedback_long.py',
        'polarity_robustness.py', 'parameter_robustness.py', 'resolution.py')]]
    p = dict(parent=str(parent), contexts=contexts, jobs=jobs, dts=list(DTS), interval=.15,
        beta=2., da=.02, db=.55, criteria=CRITERIA, existing_history=9, new_histories=0,
        source_sha256={str(f.resolve()): digest(f) for f in files},
        input_sha256={str(f.resolve()): digest(f) for f in inputs},
        design='Initial near-uniform formation on the exact saved elapsed-0 graph to 240; nonlinear contexts at elapsed 105 and 117.75 to 60 and 24. All geometry, masses and transport are fixed within each case, with no dilution or mechanical updates. Actual production SSP-RK2 functions, no surrogate reimplementation.',
        backend_scope='32 chemical ODEs: run the unchanged integrate_gm and gpu_backend.gm_step on CPU float64 arrays/tensors. This isolates their mathematics cheaply; it is not GPU arithmetic or mechanics validation. Spatial precision controls retain the resident PyTorch/custom-CUDA architecture separately.',
        interpretation='Compare absolute log errors to exact-context independent references, including nonlinear transients. Estimate order only when the finest production error exceeds 100 times the dual-solver discrepancy. Require the original .01 error limit and second-order reduction in resolvable cases. Passing frozen chemistry does not prove accuracy of split moving coupling or identify its error source.')
    write_json(root/'protocol.json', p); write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs)))
    return p


def references(root, p):
    root = Path(root); folder = root/'references'; folder.mkdir(exist_ok=True); target = folder/'results.json'
    ph = digest(root/'protocol.json')
    if target.exists():
        report = read(target)
        if report['protocol_sha256'] != ph: raise ValueError('Changed reference protocol')
        for f, h in report['paths_sha256'].items():
            if digest(f) != h: raise ValueError('Changed independent reference')
        return report
    rows = []; hashes = {}
    for context in p['contexts']:
        with np.load(context['source']) as z:
            delta, masses = validate_graph(z['delta'], z['masses']); initial = z['initial']
        times = np.arange(_steps(context['duration'], p['interval'])+1)*p['interval']
        fun = rhs(delta, p['beta'], p['da'], p['db'])
        jac = lambda t, y: chemical_jacobian(y.reshape(initial.shape), delta, p['beta'], p['da'], p['db'])
        path, error = checked_solve(fun, jac, initial, times)
        if error > p['criteria']['reference_log_max']: raise RuntimeError('Reference discrepancy above diagnostic limit')
        file = folder/(context['key']+'.npz'); np.savez_compressed(file, times=times, trajectory=path)
        hashes[str(file.resolve())] = digest(file)
        rows.append(dict(context=context['key'], solver_log_error=error, path=str(file.resolve())))
    report = dict(protocol_sha256=ph, results=rows, paths_sha256=hashes)
    write_json(target, report); return report


def worker(args):
    root, job = args; root = Path(root); p = read(root/'protocol.json'); verify(p)
    folder = root/job['key']; folder.mkdir(exist_ok=True); ph = digest(root/'protocol.json')
    if (folder/'result.json').exists():
        r = read(folder/'result.json')
        if r['protocol_sha256'] != ph or r['paths_sha256'] != digest(folder/'paths.npz'):
            raise ValueError('Changed completed chemistry trajectory')
        return r
    context = next(c for c in p['contexts'] if c['key'] == job['context'])
    with np.load(context['source']) as z:
        initial = z['initial'].copy(); delta, masses = validate_graph(z['delta'], z['masses'])
    transport = fixed_transport(delta, masses); a, b = initial.copy()
    if job['backend'] == 'torch_cpu_chemistry':
        import torch
        from .gpu_backend import gm_step
        torch.set_num_threads(1)
        a, b, matrix = [torch.from_numpy(x.copy()) for x in (a, b, delta)]
        stepper = lambda a, b: gm_step(a, b, matrix, job['dt'], p['beta'], p['da'], p['db'])
        snapshot = lambda a, b: np.stack((a.numpy(), b.numpy()))
    elif job['backend'] == 'native_chemistry':
        stepper = lambda a, b: integrate_gm(a, b, transport, job['dt'], p['beta'], p['da'], p['db'])
        snapshot = lambda a, b: np.stack((a, b))
    else: raise ValueError('Unknown chemistry backend')
    steps = _steps(context['duration'], job['dt']); every = _steps(p['interval'], job['dt'])
    path = [initial]; started = time.perf_counter(); minimum = float(initial.min())
    for step in range(1, steps+1):
        a, b = stepper(a, b); state = snapshot(a, b)
        if not np.isfinite(state).all() or np.any(state <= 0): raise FloatingPointError('Invalid chemical state')
        minimum = min(minimum, float(state.min()))
        if step % every == 0:
            path.append(state)
            if step % _steps(6., job['dt']) == 0:
                write_json(folder/'status.json', dict(state='running', elapsed=step*job['dt']))
    path = np.array(path); times = np.arange(len(path))*p['interval']
    reference = root/'references'/(context['key']+'.npz')
    with np.load(reference) as z:
        if not np.array_equal(times, z['times']): raise ValueError('Reference/sample clocks disagree')
        error = abs(np.log(path/z['trajectory']))
    np.savez_compressed(folder/'paths.npz', times=times, trajectory=path)
    result = dict(job=job, protocol_sha256=ph, paths_sha256=digest(folder/'paths.npz'),
        reference_sha256=digest(reference), max_log_error=float(error.max()),
        final_log_error=float(error[-1].max()), minimum_concentration=minimum,
        accuracy_pass=bool(error.max() <= p['criteria']['production_log_max']), wall_seconds=time.perf_counter()-started,
        onset=first_crossing(times, np.std(np.log(path[:, 0]), axis=1), .1, 'up'))
    write_json(folder/'result.json', result); write_json(folder/'status.json', dict(state='completed'))
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p); refs = references(root, p)
    rows = []; paths = {}
    for job in p['jobs']:
        folder = root/job['key']; r = read(folder/'result.json')
        if r['job'] != job or r['protocol_sha256'] != digest(root/'protocol.json') or r['paths_sha256'] != digest(folder/'paths.npz'):
            raise ValueError('Changed chemistry evidence')
        reference = root/'references'/(job['context']+'.npz')
        if digest(reference) != r['reference_sha256']: raise ValueError('Changed reference evidence')
        with np.load(folder/'paths.npz') as z, np.load(reference) as y:
            path = z['trajectory']; error = float(abs(np.log(path/y['trajectory'])).max())
            if not np.array_equal(z['times'], y['times']) or error != r['max_log_error']:
                raise ValueError('Changed error or physical clock')
            paths[job['context'], job['backend'], job['dt']] = path.copy()
        rows.append(r)
    comparisons = []
    for context in p['contexts']:
        referror = next(r['solver_log_error'] for r in refs['results'] if r['context'] == context['key'])
        for backend in ('native_chemistry', 'torch_cpu_chemistry'):
            selected = [next(r for r in rows if r['job']['context'] == context['key'] and r['job']['backend'] == backend and
                r['job']['dt'] == dt) for dt in p['dts']]
            errors = [r['max_log_error'] for r in selected]
            resolved = errors[-1] > p['criteria']['reference_noise_factor']*referror
            orders = np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved else None
            passed = all(r['accuracy_pass'] for r in selected) and (not resolved or all(
                p['criteria']['order_min'] <= order <= p['criteria']['order_max'] for order in orders))
            comparisons.append(dict(context=context['key'], backend=backend, errors=errors,
                reference_log_error=referror, order_resolved=resolved, orders=orders, passed=bool(passed)))
    pair_error = max(float(abs(np.log(paths[c['key'], 'native_chemistry', dt]/
        paths[c['key'], 'torch_cpu_chemistry', dt])).max()) for c in p['contexts'] for dt in p['dts'])
    summary = dict(completed=len(rows), total=len(p['jobs']), existing_history=9, new_histories=0,
        passed=bool(all(r['passed'] for r in comparisons) and pair_error <= p['criteria']['cpu_torch_log_max']),
        protocol_sha256=digest(root/'protocol.json'), results=rows, comparisons=comparisons,
        cpu_torch_log_max=pair_error, scope=p['interpretation'])
    write_json(root/'summary.json', summary); return summary


def run(root, workers=4):
    root = Path(root).resolve()
    if not 1 <= workers <= 4: raise ValueError('Use 1–4 small-ODE workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p)
        try:
            write_json(root/'status.json', dict(state='running', stage='independent_references', completed=0, total=len(p['jobs'])))
            references(root, p)
            done = 0
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                pending = {pool.submit(worker, (root, j)): j for j in p['jobs']}
                for future in as_completed(pending):
                    r = future.result(); done += 1
                    write_json(root/'status.json', dict(state='running', stage='production_stepping', completed=done, total=len(p['jobs'])))
                    print(f'{r["job"]["key"]}: max log error={r["max_log_error"]:.6g}', flush=True)
            summary = assess(root)
            write_json(root/'status.json', dict(state='completed' if summary['passed'] else 'completed_with_unresolved_checks',
                completed=done, total=len(p['jobs']), passed=summary['passed']))
            return summary
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/chemistry-accuracy'))
    parser.add_argument('--parent', type=Path, default=Path('outputs/polarity-robustness-refined'))
    parser.add_argument('--workers', type=int, default=4); args = parser.parse_args()
    report = prepare(args.output, args.parent) if args.action == 'prepare' else run(args.output, args.workers) if args.action == 'run' else assess(args.output)
    print({k: v for k, v in report.items() if k in ('completed', 'total', 'passed', 'existing_history', 'new_histories', 'cpu_torch_log_max', 'comparisons')})
