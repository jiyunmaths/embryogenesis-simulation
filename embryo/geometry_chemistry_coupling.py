"""Chemical time stepping on prescribed conservative changing geometry.

Use the actual production SSP-RK2 routine. Reference compartment amounts with
independent segmented DOP853/Radau; geometry is prescribed, not co-evolved.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import json
from multiprocessing import get_context
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
import torch

from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .gpu_backend import gm_step
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .polarity_robustness import checked_history
from .resolution import write_json, _steps

CONTEXTS = ((0., 240.), (90., 60.))
DTS = (.00375, .001875, .0009375)
METHODS = ('beginning', 'midpoint')
CRITERIA = dict(reference_log_max=1e-8, production_log_max=.01,
    dilution_error_max=2e-14, midpoint_order_min=1.7, midpoint_order_max=2.3,
    reference_noise_factor=100., spacing_sensitivity_log_max=.001)


def runtime():
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, torch=str(torch.__version__))


class Schedule:
    def __init__(self, times, masses, conductance):
        self.times = np.asarray(times, dtype=np.float64)
        self.masses = np.asarray(masses, dtype=np.float64)
        self.conductance = np.asarray(conductance, dtype=np.float64)
        frames, cells = self.masses.shape
        if (self.times.shape != (frames,) or frames < 2 or self.conductance.shape != (frames, cells, cells)
                or not all(np.isfinite(x).all() for x in (self.times, self.masses, self.conductance))
                or np.any(self.masses <= 0) or np.any(np.diff(self.times) <= 0)
                or not np.allclose(self.conductance, self.conductance.transpose(0, 2, 1), rtol=0, atol=1e-12)
                or np.any(self.conductance < 0) or np.any(np.diagonal(self.conductance, axis1=1, axis2=2) != 0)):
            raise ValueError('Finite positive masses and symmetric nonnegative conductances required')
        self.interval = float(self.times[1]-self.times[0])
        if not np.allclose(np.diff(self.times), self.interval, rtol=0, atol=1e-10):
            raise ValueError('Uniform source observation spacing required')
        self.exchange = self.conductance.copy()
        self.exchange[:, np.arange(cells), np.arange(cells)] = -self.conductance.sum(2)

    def at(self, t):
        if t < self.times[0]-1e-9 or t > self.times[-1]+1e-9: raise ValueError('Time outside prescribed geometry')
        coordinate = (np.clip(t, self.times[0], self.times[-1])-self.times[0])/self.interval
        index = min(int(coordinate), len(self.times)-2); fraction = coordinate-index
        volume = self.masses[index]+fraction*(self.masses[index+1]-self.masses[index])
        exchange = self.exchange[index]+fraction*(self.exchange[index+1]-self.exchange[index])
        return volume, exchange

    def delta(self, t):
        volume, exchange = self.at(t); return volume, exchange/volume[:, None]


def amount_rhs(schedule, beta, da, db):
    def fun(t, amount):
        volume, exchange = schedule.at(t); a, b = amount.reshape(2, -1)/volume
        reaction = np.stack((a*a/b-a, beta*(a*a-b)))
        flux = np.stack((da*(exchange@a), db*(exchange@b)))
        return (volume*reaction+flux).ravel()
    return fun


def amount_jacobian(schedule, beta, da, db):
    def jac(t, amount):
        volume, delta = schedule.delta(t); state = amount.reshape(2, -1)/volume
        matrix = chemical_jacobian(state, delta, beta, da, db)
        scale = np.tile(volume, 2)
        return scale[:, None]*matrix/scale[None, :]
    return jac


def production_step(a, b, schedule, t, dt, method, beta=2., da=.02, db=.55, stepper=gm_step):
    """Two time-coupling schemes, with the identical actual chemical step."""
    old_volume = schedule.at(t)[0]; new_volume = schedule.at(t+dt)[0]
    if method == 'beginning':
        sample_volume, delta = schedule.delta(t)
    elif method == 'midpoint':
        sample_volume, delta = schedule.delta(t+dt/2)
    else: raise ValueError('Unknown time-coupling method')
    before = torch.stack([torch.from_numpy(old_volume)@a, torch.from_numpy(old_volume)@b])
    if method == 'midpoint':
        ratio = torch.from_numpy(old_volume/sample_volume); a, b = a*ratio, b*ratio
    sample = torch.from_numpy(sample_volume)
    middle = torch.stack([sample@a, sample@b]); error = float(((middle/before)-1).abs().max())
    a, b = stepper(a, b, torch.from_numpy(delta), dt, beta, da, db)
    reacted = torch.stack([sample@a, sample@b])
    ratio = torch.from_numpy(sample_volume/new_volume); a, b = a*ratio, b*ratio
    end = torch.from_numpy(new_volume)
    error = max(error, float(((torch.stack([end@a, end@b])/reacted)-1).abs().max()))
    return a, b, error


def load_context(p, context, spacing):
    with np.load(context['source']) as z:
        stride = 1 if spacing == 'fine' else 2 if spacing == 'coarse' else None
        if stride is None: raise ValueError('Unknown prescribed observation spacing')
        return Schedule(z['times'][::stride], z['masses'][::stride], z['conductance'][::stride]), z['initial'].copy()


def prepare(root, parent=Path('outputs/polarity-robustness-refined')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists(): raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    folder = Path(old['reference_folder']); history = read(folder/'history.json'); r = read(folder/'result.json')
    checked_history(history, old, old['reference_job'], 240.)
    if not r['quality_pass'] or digest(folder/'history.json') != r['history_sha256']:
        raise ValueError('Changed measured geometry source')
    root.mkdir(parents=True); contexts = []; jobs = []; references = []; projection = 0.
    inputs = [parent/'protocol.json', folder/'history.json', folder/'result.json']
    for elapsed, duration in CONTEXTS:
        rows = [row for row in history if elapsed-1e-9 <= row['elapsed'] <= elapsed+duration+1e-9]
        times = np.arange(len(rows))*.15
        if not np.allclose([r['elapsed'] for r in rows], times+elapsed, rtol=0, atol=1e-9):
            raise ValueError('Geometry source clocks differ')
        masses = []; conductances = []
        for row in rows:
            delta, volume = validate_graph(row['delta'], row['volumes'])
            g = volume[:, None]*delta; np.fill_diagonal(g, 0.)
            symmetric = .5*(g+g.T)
            projection = max(projection, float(abs(symmetric-g).max()))
            masses.append(volume); conductances.append(symmetric)
        key = f'elapsed-{elapsed:g}'; file = root/(key+'.npz')
        np.savez_compressed(file, times=times, masses=masses, conductance=conductances,
                            initial=rows[0]['chemistry'], ids=rows[0]['ids'])
        context = dict(key=key, elapsed=elapsed, duration=duration, source=str(file.resolve()))
        contexts.append(context); inputs.append(file)
        for spacing in ('fine','coarse'):
            load_context({}, context, spacing)
            references.append(dict(key=f'{key}_{spacing}', context=key, spacing=spacing))
            for method in METHODS:
                for dt in DTS:
                    jobs.append(dict(key=f'{key}_{spacing}_{method}_dt-{dt:g}', context=key,
                                     spacing=spacing, method=method, dt=dt))
    sources = {**old['source_sha256'], **{str(Path(__file__).with_name(f).resolve()):digest(Path(__file__).with_name(f))
        for f in ('geometry_chemistry_coupling.py','gpu_backend.py','feedback_endpoint_bistability.py',
                  'neighbor_context.py','parameter_robustness.py','resolution.py')}}
    p = dict(contexts=contexts, jobs=jobs, references=references, dts=list(DTS), methods=list(METHODS),
        interval=.15, beta=2., da=.02, db=.55, criteria=CRITERIA, conductance_symmetrization_abs_max=projection,
        independent_histories=1, existing_history=9, new_histories=0, environment=runtime(),
        source_sha256=sources, input_sha256={str(f.resolve()):digest(f) for f in inputs},
        design='Prescribe fine history-9 geometry from elapsed 0 to 240 and from elapsed 90 to 150 (60 further units). Same chemical start within each context; interpolate positive volumes and symmetric nonnegative conductances at source spacing .15 and coarsened .3. Reconstruct Delta=M^-1(G-diag(sum G)). Never directly interpolate Delta while changing M. No mechanics/feedback updates.',
        reference_scope='DOP853/Radau solve compartment amounts Q=M c using analytic transformed Jacobian, stopping at every geometry knot with rtol=1e-12, atol=1e-14. Independent solver disagreement <1e-8. Geometry spacing sensitivity is screened separately at .001, not hidden in a solver tolerance.',
        methods_scope='Actual unchanged gpu_backend.gm_step in CPU float64 PyTorch for this 32-ODE diagnostic. Beginning matches current time order: chemistry on beginning geometry, then volume dilution. Midpoint converts starting amounts to midpoint concentrations, advances chemistry with midpoint geometry, then converts to ending concentrations. Both preserve the amount during each volume conversion. Moving confirmation uses resident GPU PyTorch plus custom CUDA separately.',
        interpretation='Measure absolute chemical errors and order on identical prescribed inputs, not backend acceptance or live-trajectory convergence. Quantitative production limit remains .01. Only midpoint is tested for second order when error is resolvable relative to independent-reference noise; the beginning scheme may be first order. .15-versus-.3 geometry spacing comparison is a sensitivity screen, not proof of faithful interpolation or source refinement. Passing cannot establish that splitting explains the original .0375 moving discrepancy. Preserve all old failures.')
    write_json(root/'protocol.json', p); write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs)))
    return p


def reference_worker(args):
    root, job = args; root = Path(root); p = read(root/'protocol.json'); verify(p)
    folder = root/'references'/job['key']; folder.mkdir(exist_ok=True, parents=True); ph = digest(root/'protocol.json')
    if (folder/'result.json').exists():
        r = read(folder/'result.json')
        if r['protocol_sha256'] != ph or r['job'] != job or r['paths_sha256'] != digest(folder/'paths.npz'):
            raise ValueError('Changed independent coupling reference')
        return r
    context = next(c for c in p['contexts'] if c['key'] == job['context'])
    schedule, initial = load_context(p, context, job['spacing'])
    fun = amount_rhs(schedule, p['beta'], p['da'], p['db']); jac = amount_jacobian(schedule, p['beta'], p['da'], p['db'])
    times = np.arange(_steps(context['duration'], p['interval'])+1)*p['interval']; paths = []; amount_paths = []
    started = time.perf_counter()
    for method in ('DOP853','Radau'):
        current = (initial*schedule.masses[0]).ravel(); path = [initial.copy()]; amounts = [current.reshape(initial.shape).copy()]
        for left, right in zip(schedule.times[:-1], schedule.times[1:]):
            selected = times[(times > left+1e-9)&(times <= right+1e-9)]
            kwargs = dict(jac=jac) if method == 'Radau' else {}
            sol = solve_ivp(fun, (left, right), current, method=method, t_eval=selected, rtol=1e-12, atol=1e-14, **kwargs)
            if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y <= 0):
                raise RuntimeError('Independent amount solver failed')
            for t, amount in zip(selected, sol.y.T):
                state = amount.reshape(initial.shape)
                amounts.append(state.copy()); path.append(state/schedule.at(t)[0])
            current = sol.y[:, -1].copy()
        paths.append(np.array(path)); amount_paths.append(np.array(amounts))
        write_json(folder/'status.json', dict(state='running', solver=method, completed_solvers=len(paths)))
    error = float(abs(np.log(paths[0]/paths[1])).max())
    if error > p['criteria']['reference_log_max']: raise RuntimeError(f'Independent reference gap {error} exceeds limit')
    np.savez_compressed(folder/'paths.npz', times=times, trajectory=paths[0], radau_trajectory=paths[1],
                        dop853_amounts=amount_paths[0], radau_amounts=amount_paths[1])
    result = dict(job=job, protocol_sha256=ph, paths_sha256=digest(folder/'paths.npz'),
                  solver_log_error=error, passed=True, wall_seconds=time.perf_counter()-started)
    write_json(folder/'result.json', result); write_json(folder/'status.json', dict(state='completed', passed=True))
    return result


def worker(args):
    root, job = args; root = Path(root); p = read(root/'protocol.json'); verify(p); torch.set_num_threads(1)
    folder = root/job['key']; folder.mkdir(exist_ok=True); ph = digest(root/'protocol.json')
    if (folder/'result.json').exists():
        r = read(folder/'result.json')
        if r['protocol_sha256'] != ph or r['job'] != job or r['paths_sha256'] != digest(folder/'paths.npz'):
            raise ValueError('Changed completed coupling trajectory')
        return r
    context = next(c for c in p['contexts'] if c['key'] == job['context'])
    schedule, initial = load_context(p, context, job['spacing']); dt = job['dt']
    steps = _steps(context['duration'], dt); every = _steps(p['interval'], dt)
    checkpoint_every = max(every, _steps(min(6., context['duration']), dt)); checkpoint = folder/'latest.npz'
    if checkpoint.exists():
        with np.load(checkpoint) as z:
            if str(z['protocol_sha256']) != ph or json.loads(str(z['job'])) != job:
                raise ValueError('Changed coupling checkpoint')
            step = int(z['step']); path = list(z['trajectory'].copy()); a, b = [torch.from_numpy(x.copy()) for x in z['state']]
            minimum = float(z['minimum']); dilution = float(z['dilution_error']); base_wall = float(z['wall_seconds'])
        if step % every or len(path) != step//every+1 or not np.array_equal(path[-1], np.stack((a.numpy(), b.numpy()))):
            raise ValueError('Coupling checkpoint and sampled clock disagree')
    else:
        step = 0; a, b = [torch.from_numpy(x.copy()) for x in initial]; path = [initial.copy()]
        minimum = float(initial.min()); dilution = 0.; base_wall = 0.
    started = time.perf_counter()
    while step < steps:
        a, b, error = production_step(a, b, schedule, step*dt, dt, job['method'], p['beta'], p['da'], p['db'])
        step += 1; state = np.stack((a.numpy(), b.numpy()))
        if not np.isfinite(state).all() or np.any(state <= 0): raise FloatingPointError('Invalid replay chemistry')
        dilution = max(dilution, error); minimum = min(minimum, float(state.min()))
        if dilution > p['criteria']['dilution_error_max']: raise RuntimeError('Amount conversion not conservative')
        if step % every == 0: path.append(state.copy())
        if step % checkpoint_every == 0 or step == steps:
            wall = base_wall+time.perf_counter()-started; temporary = folder/'latest.tmp.npz'
            np.savez_compressed(temporary, protocol_sha256=np.array(ph), job=np.array(json.dumps(job)), step=step,
                state=state, trajectory=np.array(path), minimum=minimum, dilution_error=dilution, wall_seconds=wall)
            temporary.replace(checkpoint)
            write_json(folder/'status.json', dict(state='running', elapsed=step*dt, until=context['duration'], wall_seconds=wall))
    path = np.array(path); times = np.arange(len(path))*p['interval']
    reference = root/'references'/f'{job["context"]}_{job["spacing"]}'/'paths.npz'
    with np.load(reference) as z:
        if not np.array_equal(times, z['times']): raise ValueError('Coupling reference clock differs')
        error = abs(np.log(path/z['trajectory']))
    np.savez_compressed(folder/'paths.npz', times=times, trajectory=path)
    result = dict(job=job, protocol_sha256=ph, paths_sha256=digest(folder/'paths.npz'), reference_sha256=digest(reference),
        max_log_error=float(error.max()), final_log_error=float(error[-1].max()), minimum_concentration=minimum,
        dilution_error_max=dilution, accuracy_pass=bool(error.max() <= p['criteria']['production_log_max']),
        wall_seconds=base_wall+time.perf_counter()-started)
    write_json(folder/'result.json', result); write_json(folder/'status.json', dict(state='completed'))
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p); rows = []; references = {}; paths = {}
    for job in p['references']:
        folder = root/'references'/job['key']; r = read(folder/'result.json')
        if (r['job'] != job or r['protocol_sha256'] != digest(root/'protocol.json') or
                r['paths_sha256'] != digest(folder/'paths.npz')): raise ValueError('Changed coupling reference')
        with np.load(folder/'paths.npz') as z:
            error = float(abs(np.log(z['trajectory']/z['radau_trajectory'])).max())
            if error != r['solver_log_error']: raise ValueError('Changed independent solver gap')
            references[job['context'], job['spacing']] = dict(result=r, trajectory=z['trajectory'].copy())
    for job in p['jobs']:
        folder = root/job['key']; r = read(folder/'result.json')
        reference = root/'references'/f'{job["context"]}_{job["spacing"]}'/'paths.npz'
        if (r['job'] != job or r['protocol_sha256'] != digest(root/'protocol.json') or
                digest(folder/'paths.npz') != r['paths_sha256'] or digest(reference) != r['reference_sha256']):
            raise ValueError('Changed coupling evidence')
        with np.load(folder/'paths.npz') as z, np.load(reference) as ref:
            error = float(abs(np.log(z['trajectory']/ref['trajectory'])).max())
            if not np.array_equal(z['times'], ref['times']) or error != r['max_log_error']:
                raise ValueError('Changed replay clock or error')
            paths[job['context'], job['spacing'], job['method'], job['dt']] = z['trajectory'].copy()
        rows.append(r)
    comparisons = []; spacing_checks = []
    for context in p['contexts']:
        for spacing in ('fine','coarse'):
            noise = references[context['key'], spacing]['result']['solver_log_error']
            for method in METHODS:
                selected = [next(r for r in rows if r['job']['context'] == context['key'] and
                    r['job']['spacing'] == spacing and r['job']['method'] == method and r['job']['dt'] == dt) for dt in p['dts']]
                errors = [r['max_log_error'] for r in selected]; resolved = errors[-1] > p['criteria']['reference_noise_factor']*noise
                orders = np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved else None
                order_pass = method != 'midpoint' or not resolved or all(p['criteria']['midpoint_order_min'] <= o <=
                    p['criteria']['midpoint_order_max'] for o in orders)
                comparisons.append(dict(context=context['key'], spacing=spacing, method=method, errors=errors,
                    reference_log_error=noise, order_resolved=resolved, orders=orders, order_pass=bool(order_pass),
                    passed=bool(all(r['accuracy_pass'] for r in selected) and order_pass)))
        sensitivity = float(abs(np.log(references[context['key'],'fine']['trajectory']/
            references[context['key'],'coarse']['trajectory'])).max())
        spacing_checks.append(dict(context=context['key'], reference_spacing_log_difference=sensitivity,
            passed=bool(sensitivity <= p['criteria']['spacing_sensitivity_log_max'])))
    summary = dict(protocol_sha256=digest(root/'protocol.json'), completed=len(rows), total=len(p['jobs']),
        independent_histories=1, new_histories=0, comparisons=comparisons, spacing_checks=spacing_checks,
        reference_pass=all(r['result']['passed'] for r in references.values()),
        time_stepping_pass=all(r['passed'] for r in comparisons), spacing_screen_pass=all(r['passed'] for r in spacing_checks),
        results=rows, scope=p['interpretation'])
    write_json(root/'summary.json', summary); return summary


def run(root, workers=4):
    root = Path(root).resolve()
    if not 1 <= workers <= 4: raise ValueError('Use 1–4 independent small-ODE workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); done = 0
        if p['environment'] != runtime(): raise ValueError('Changed CPU numerical runtime')
        try:
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                write_json(root/'status.json', dict(state='running', stage='independent_references', completed=0, total=len(p['jobs'])))
                for future in as_completed([pool.submit(reference_worker, (root, j)) for j in p['references']]):
                    r = future.result(); print(f'Reference {r["job"]["key"]}: gap={r["solver_log_error"]:.6g}', flush=True)
                write_json(root/'status.json', dict(state='running', stage='production_replay', completed=0, total=len(p['jobs'])))
                for future in as_completed([pool.submit(worker, (root, j)) for j in p['jobs']]):
                    r = future.result(); done += 1
                    write_json(root/'status.json', dict(state='running', stage='production_replay', completed=done, total=len(p['jobs'])))
                    print(f'{r["job"]["key"]}: error={r["max_log_error"]:.6g}', flush=True)
            summary = assess(root); passed = summary['time_stepping_pass'] and summary['spacing_screen_pass']
            write_json(root/'status.json', dict(state='completed' if passed else 'completed_with_unresolved_checks',
                completed=done, total=len(p['jobs']), time_stepping_pass=summary['time_stepping_pass'],
                spacing_screen_pass=summary['spacing_screen_pass']))
            return summary
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=len(p['jobs']), error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','run','assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/geometry-chemistry-coupling'))
    parser.add_argument('--workers', type=int, default=4); args = parser.parse_args()
    report = prepare(args.output) if args.action == 'prepare' else run(args.output, args.workers) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in report.items() if k in ('completed','total','time_stepping_pass','spacing_screen_pass','comparisons')})
