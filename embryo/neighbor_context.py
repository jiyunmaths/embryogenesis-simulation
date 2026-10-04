"""One-time neighbor-chemistry challenges on frozen mature contact graphs.

Target chemistry and all geometry are identical at intervention. Neighbors then
evolve freely; there are no clamped cells, reservoir, fate labels, or mechanics.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import json
from pathlib import Path
import time

import numpy as np

from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs
from .cell_exchange_response_moving import verify
from .cell_response import TIMES, pulse, response_metrics
from .cell_response_exchange import checked_solve, waveform_distance
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .resolution import write_json


ARMS = ('untouched', 'sham', 'neighbor_mix', 'neighbor_low', 'neighbor_high')
FACTORS = (.9, 1.1)
CRITERIA = dict(solver_log_max=1e-5, sham_log_max=1e-10,
                amount_relative_max=2e-14, neighbor_challenge_log_min=.05,
                initial_pair_log_min=.1, late_state_ratio_max=.1,
                response_shift_rms_min=.01, response_reference_min=.01,
                stationarity_rhs_max=1e-6, late_window=24.,
                metric_recompute_rtol=1e-12, metric_recompute_atol=1e-14)


def read(path):
    return json.loads(Path(path).read_text())


def validate_graph(delta, masses):
    delta, masses = np.asarray(delta, float), np.asarray(masses, float)
    n = len(masses)
    if (masses.ndim != 1 or n < 2 or delta.shape != (n, n) or
            not np.isfinite(delta).all() or not np.isfinite(masses).all() or np.any(masses <= 0)):
        raise ValueError('Expected finite contact operator and positive matching volumes')
    off = delta.copy(); np.fill_diagonal(off, 0.)
    scale = max(float(abs(delta).max()), 1.)
    if (np.any(off < 0) or np.any(np.diag(delta) > 0) or
            abs(delta.sum(axis=1)).max() > 1e-12*scale or
            abs(masses@delta).max() > 1e-12*scale*masses.sum() or
            not np.allclose(masses[:, None]*delta, (masses[:, None]*delta).T,
                            rtol=1e-12, atol=1e-14*scale)):
        raise ValueError('Expected conservative reciprocal contact transport')
    return delta, masses


def challenge(initial, masses, delta, ids, target, arm):
    """Only direct neighbors change; averaging conserves each species' amount."""
    delta, masses = validate_graph(delta, masses)
    initial, ids = np.asarray(initial, float), np.asarray(ids)
    n = len(masses)
    if (initial.shape != (2, n) or not np.isfinite(initial).all() or np.any(initial <= 0) or
            ids.shape != (n,) or len(set(ids.tolist())) != n or target not in ids or arm not in ARMS):
        raise ValueError('Invalid chemistry, cell IDs, target, or intervention')
    i = list(ids).index(target)
    neighbors = np.flatnonzero(delta[i] > 0)
    changed = initial.copy()
    order = sorted(range(n), key=lambda j: (initial[0, j], int(ids[j])))
    templates = dict(neighbor_low=order[0], neighbor_high=order[-1])
    if arm == 'neighbor_mix' and len(neighbors):
        changed[:, neighbors] = ((initial[:, neighbors]@masses[neighbors])/masses[neighbors].sum())[:, None]
    elif arm in templates:
        changed[:, neighbors] = initial[:, templates[arm], None]
    elif arm == 'sham':
        changed[:, neighbors] = initial[:, neighbors].copy()
    untouched = np.ones(n, bool); untouched[neighbors] = False
    if not np.array_equal(changed[:, untouched], initial[:, untouched]):
        raise ValueError('Intervention changed the target or a nonneighbor')
    before, after = initial@masses, changed@masses
    relative = float(abs((after-before)/before).max())
    if arm in ('untouched', 'sham', 'neighbor_mix') and relative > CRITERIA['amount_relative_max']:
        raise ValueError('Amount-preserving intervention failed conservation')
    size = (float(np.sqrt(np.sum(np.log(changed[:, neighbors]/initial[:, neighbors])**2*masses[neighbors]) /
                          (2*masses[neighbors].sum()))) if len(neighbors) else 0.)
    metadata = dict(target=int(target), target_index=i, neighbors=ids[neighbors].tolist(),
        template_cell=int(ids[templates[arm]]) if arm in templates else None,
        neighbor_log_rms=size, amount_change=(after-before).tolist(),
        relative_amount_change=relative, target_unchanged=True, nonneighbors_unchanged=True,
        target_transport_change=(delta[i]@(changed-initial).T).tolist())
    return changed, metadata


def job_name(job):
    return f"cell-{job['target']}_{job['arm']}"


def require_source(source, seed, family):
    source = Path(source); cp = read(source/'protocol.json'); verify(cp)
    status_file = source/'status.json'
    # The original seed-7 native coordinator records completion only at its
    # root. Later GPU adapters additionally provide a per-background status.
    if status_file.exists():
        if read(status_file).get('state') != 'completed':
            raise ValueError('Incomplete source background')
    elif seed != 7:
        raise ValueError('Missing new-history source status')
    result = read(source/f'{family}_control/result.json')
    if not result.get('quality_pass') or result.get('protocol_sha256') != digest(source/'protocol.json'):
        raise ValueError('Invalid source control evidence')
    return cp


def prepare(root, histories=Path('outputs/exchange-response-histories'),
            refined=Path('outputs/exchange-response-histories-refined'),
            historical=Path('outputs/cell-exchange-response-moving')):
    root, histories, refined, historical = [Path(p).resolve() for p in (root, histories, refined, historical)]
    if root.exists():
        raise FileExistsError(root)
    status = read(refined/'status.json'); fine = read(refined/'protocol.json'); verify(fine)
    if (status.get('state') != 'completed' or not status.get('passed') or status.get('completed') != 18 or
            not read(refined/'refinement.json').get('passed')):
        raise ValueError('Requires completed passing formation/response refinement')
    old = read(histories/'protocol.json'); verify(old)
    if read(histories/'status.json').get('completed') != 36 or not read(histories/'comparison.json').get('completed'):
        raise ValueError('Requires completed history replication')
    verify(read(historical/'protocol.json'))
    if read(historical/'status.json').get('state') != 'completed':
        raise ValueError('Historical seed-7 response study is incomplete')
    inputs = [refined/f for f in ('protocol.json', 'status.json', 'refinement.json')]
    inputs += [histories/f for f in ('protocol.json', 'status.json', 'comparison.json')]
    inputs += [historical/f for f in ('protocol.json', 'status.json', 'comparison.json')]
    contexts, snapshots = [], []
    roots = {7: historical, 8: histories/'seed-8/response', 9: histories/'seed-9/response'}
    for seed, parent in roots.items():
        for family in ('unexchanged', 'fresh_exchange'):
            source = parent/family; cp = require_source(source, seed, family)
            control = source/f'{family}_control'
            sim = AttributeSimulation.restore(cp['checkpoint']); graph = sim.signaling_graph()
            if (sim.time != 210. or len(sim.ids) != 16 or sim.divisions or
                    sim.attribute_mode != 'direct' or sim.config.signal_transport != 'conservative'):
                raise ValueError('Requires mature direct-feedback conservative t=210 source')
            with np.load(source/'initial_states.npz') as d:
                initial, ids, masses = d[family].copy(), d['ids'].copy(), d['masses'].copy()
            if (not np.array_equal(ids, sim.ids) or not np.array_equal(initial, [sim.activator, sim.inhibitor]) or
                    not np.allclose(masses, graph.masses, rtol=1e-12, atol=1e-14)):
                raise ValueError('Prepared chemistry/IDs/volumes differ from physical source')
            delta = graph.delta; validate_graph(delta, masses)
            first = read(control/'history.json')[0]
            if (first['time'] != 210. or first['ids'] != ids.tolist() or
                    not np.array_equal(first['chemistry'], initial) or
                    not np.allclose(first['delta'], delta, rtol=1e-5, atol=1e-8)):
                raise ValueError('Original control and frozen source disagree')
            targets = sorted(set(j['target'] for j in cp['jobs'] if j['target'] is not None))
            if len(targets) != 2 or any(target not in ids for target in targets):
                raise ValueError('Requires the two originally prespecified response targets')
            key = f'seed-{seed}_{family}'
            jobs = [dict(target=int(target), arm=arm) for target in targets for arm in ARMS]
            for job in jobs: challenge(initial, masses, delta, ids, **job)
            contexts.append(dict(key=key, seed=seed, background=family, targets=targets,
                checkpoint=str(Path(cp['checkpoint']).resolve()), start=210., jobs=jobs,
                beta=sim.config.signal_beta, da=sim.config.signal_da, db=sim.config.signal_dh))
            snapshots.append(dict(initial=initial, ids=ids, masses=masses, delta=delta, centers=sim.centers()))
            inputs += [source/'protocol.json', source/'initial_states.npz',
                       Path(cp['checkpoint']), control/'result.json', control/'history.json']
            if (source/'status.json').exists(): inputs.append(source/'status.json')
    sources = [Path(__file__), *[Path(__file__).with_name(f) for f in (
        'attribute_development.py', 'model.py', 'polarity.py', 'signaling.py', 'transport.py',
        'attribute_persistence.py', 'cell_response.py', 'cell_response_exchange.py',
        'cell_exchange_response_moving.py', 'feedback_endpoint_bistability.py', 'feedback_long.py', 'resolution.py')]]
    root.mkdir(parents=True)
    for context, snapshot in zip(contexts, snapshots):
        folder = root/context['key']; folder.mkdir()
        filename = folder/'source.npz'; np.savez_compressed(filename, **snapshot)
        context['source'] = str(filename); inputs.append(filename)
    p = dict(contexts=contexts, arms=ARMS, factors=FACTORS, times=TIMES.tolist(),
        criteria=CRITERIA, total_jobs=60, independent_histories=3,
        design='Six original t=210 untouched/fresh moving endpoints; two prespecified targets each. Five one-time neighbor interventions; a free frozen-graph control and immediate +/-10% target activator pulses per intervention, for 240 units. Target and nonneighbor initial chemistry stay identical. All cells then evolve without clamping.',
        templates='Min/max initial activator in each physical background, ties by cell ID; both species copied together. Neighbor mixing uses the measured-volume weighted mean and conserves amounts. Low/high resets externally add/remove recorded amounts.',
        interpretation='Continuous target-state changes and matched pulse responses relative to each intervention own control. Neighbor-only change is causal within the frozen model; it does not isolate a particular neighbor or establish cell autonomy, inheritance, moving mechanics, or biological function.',
        limits='Three histories; backgrounds, targets, arms and pulse signs are nested interventions. Frozen chemistry uses DOP853/Radau CPU references, not the GPU mechanics backend. Weak challenges and unresolved references remain explicit. Original t=210 starts match the completed same-state refinement; refined formation endpoints are not substituted. General transport closure and developmental convergence remain unresolved.',
        input_sha256={str(f.resolve()): digest(f) for f in inputs},
        source_sha256={str(f.resolve()): digest(f) for f in sources})
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=p['total_jobs']))
    return p


def load_result(folder, ph, criteria):
    record = read(folder/'result.json')
    if (record.get('protocol_sha256') != ph or not record.get('quality_pass') or
            record.get('trajectory_sha256') != digest(folder/'trajectories.npz') or
            not np.isfinite(record['max_solver_log_error']) or record['max_solver_log_error'] >= criteria['solver_log_max']):
        raise ValueError('Invalid saved neighbor-context evidence: '+str(folder))
    with np.load(folder/'trajectories.npz') as data:
        arrays = {k: data[k].copy() for k in data.files}
    for key in ('initial', 'control', 'pulse_paths'):
        if not np.isfinite(arrays[key]).all() or np.any(arrays[key] <= 0):
            raise ValueError('Nonpositive or nonfinite saved chemistry')
    if (arrays['control'].shape != (len(arrays['times']), *arrays['initial'].shape) or
            arrays['pulse_paths'].shape != (len(FACTORS), *arrays['control'].shape) or
            arrays['responses'].shape != (len(FACTORS), len(arrays['times']), 2) or
            not np.isfinite(arrays['responses']).all()):
        raise ValueError('Malformed or incomplete trajectories')
    return record, arrays


def worker(args):
    root, context, job = args
    root = Path(root); p = read(root/'protocol.json'); ph = digest(root/'protocol.json')
    source = Path(context['source'])
    if digest(source) != p['input_sha256'][str(source.resolve())]:
        raise ValueError('Frozen graph snapshot changed')
    folder = root/context['key']/job_name(job); folder.mkdir(exist_ok=True)
    if (folder/'result.json').exists():
        return load_result(folder, ph, p['criteria'])[0]
    with np.load(source) as d:
        initial, masses, delta, ids = [d[k].copy() for k in ('initial', 'masses', 'delta', 'ids')]
    changed, intervention = challenge(initial, masses, delta, ids, **job)
    i = intervention['target_index']; times = np.asarray(p['times'])
    fun = rhs(delta, context['beta'], context['da'], context['db'])
    jac = lambda t, y: chemical_jacobian(y.reshape(2, -1), delta, context['beta'], context['da'], context['db'])
    started = time.monotonic()
    write_json(folder/'status.json', dict(state='running', phase='control'))
    control, error = checked_solve(fun, jac, changed, times)
    errors, paths, responses, trials = [error], [], [], []
    for factor in p['factors']:
        write_json(folder/'status.json', dict(state='running', phase='pulse', factor=factor))
        path, error = checked_solve(fun, jac, pulse(changed, i, factor), times)
        errors.append(error); paths.append(path)
        responses.append(np.log(path[:, :, i]/control[:, :, i])/abs(np.log(factor)))
        trials.append(dict(factor=factor, injected_activator_amount=float((factor-1)*changed[0, i]*masses[i]),
                           **response_metrics(path, control, times, masses, i, factor)))
    end = control[-1]
    residual = float(abs(fun(times[-1], end.ravel())).max())
    growth = float(np.linalg.eigvals(jac(times[-1], end.ravel())).real.max())
    filename = folder/'trajectories.npz'; temporary = folder/'trajectories.tmp.npz'
    np.savez_compressed(temporary, initial=changed, control=control, pulse_paths=np.asarray(paths),
                       responses=np.asarray(responses), times=times)
    temporary.replace(filename)
    record = dict(protocol_sha256=ph, trajectory_sha256=digest(filename), job=job,
        context=context['key'], intervention=intervention, max_solver_log_error=max(errors),
        target_initial_rhs_change=(np.array([context['da'], context['db']])*
                                   intervention['target_transport_change']).tolist(),
        quality_pass=bool(max(errors) < p['criteria']['solver_log_max']), trials=trials,
        endpoint_rhs_max=residual, endpoint_jacobian_max_real=growth,
        stationary_stable=bool(residual < p['criteria']['stationarity_rhs_max'] and growth < 0),
        final_log_activator_sd=float(np.std(np.log(end[0]))), wall_seconds=time.monotonic()-started)
    write_json(folder/'result.json', record)
    write_json(folder/'status.json', dict(state='completed', quality_pass=record['quality_pass']))
    return record


def compare_context(root, context, p):
    root = Path(root); ph = digest(root/'protocol.json'); rows = []; sham_errors = []
    with np.load(context['source']) as d:
        initial, ids, masses, delta = [d[k].copy() for k in ('initial', 'ids', 'masses', 'delta')]
    loaded = {}
    for job in context['jobs']:
        record, arrays = load_result(root/context['key']/job_name(job), ph, p['criteria'])
        expected, metadata = challenge(initial, masses, delta, ids, **job)
        if (not np.array_equal(arrays['times'], p['times']) or not np.array_equal(arrays['initial'], expected) or
                not np.array_equal(arrays['control'][0], expected) or record['intervention'] != metadata or record['job'] != job):
            raise ValueError('Saved intervention or observation alignment changed')
        i = metadata['target_index']
        for f, factor in enumerate(p['factors']):
            if not np.array_equal(arrays['pulse_paths'][f, 0], pulse(expected, i, factor)):
                raise ValueError('Pulse changed target definition or starting state')
            wave = np.log(arrays['pulse_paths'][f, :, :, i]/arrays['control'][:, :, i])/abs(np.log(factor))
            if not np.array_equal(arrays['responses'][f], wave):
                raise ValueError('Saved response differs from paired trajectories')
            metrics = dict(factor=factor, injected_activator_amount=float((factor-1)*expected[0, i]*masses[i]),
                **response_metrics(arrays['pulse_paths'][f], arrays['control'], arrays['times'], masses, i, factor))
            saved = record['trials'][f]
            # Saving/reloading arrays changes reduction memory order. Preserve
            # exact recovery classifications/times; allow roundoff in gains.
            if set(saved) != set(metrics) or any(
                    (saved[k] != value if value is None or k.endswith('recovery_time') else
                     not np.isfinite(saved[k]) or not np.isclose(saved[k], value,
                         rtol=p['criteria']['metric_recompute_rtol'], atol=p['criteria']['metric_recompute_atol']))
                    for k, value in metrics.items()):
                raise ValueError('Saved recovery/response metrics differ from trajectories')
        loaded[(job['target'], job['arm'])] = (record, arrays)
    late = np.asarray(p['times']) >= p['times'][-1]-p['criteria']['late_window']
    for target in context['targets']:
        index = list(ids).index(target)
        other = next(x for x in context['targets'] if x != target)
        other_index = list(ids).index(other)
        base = loaded[(target, 'untouched')][1]; opposite = loaded[(other, 'untouched')][1]
        if not np.array_equal(base['control'], opposite['control']):
            raise ValueError('Untouched targets do not share the same physical background')
        separation = float(np.sqrt(np.mean(np.log(initial[:, index]/initial[:, other_index])**2)))
        sham = loaded[(target, 'sham')][1]
        sham_errors += [float(abs(np.log(sham['control']/base['control'])).max()),
                        float(abs(sham['responses']-base['responses']).max())]
        for arm in p['arms'][2:]:
            record, arrays = loaded[(target, arm)]
            change = np.sqrt(np.mean(np.log(arrays['control'][:, :, index]/base['control'][:, :, index])**2, axis=1))
            maximum = float(change[late].max())
            informative = (record['intervention']['neighbor_log_rms'] >= p['criteria']['neighbor_challenge_log_min'] and
                           separation > p['criteria']['initial_pair_log_min'])
            response_rows = []
            for f, factor in enumerate(p['factors']):
                own = waveform_distance(arrays['responses'][f], base['responses'][f], p['times'])
                alternate = waveform_distance(arrays['responses'][f], opposite['responses'][f], p['times'])
                reference = waveform_distance(base['responses'][f], opposite['responses'][f], p['times'])
                response_rows.append(dict(factor=factor, own_response_distance=own,
                    other_target_response_distance=alternate, reference_separation=reference,
                    response_shift=bool(informative and own > p['criteria']['response_shift_rms_min']),
                    nearest_baseline_response=('own' if own <= alternate else 'other_target')
                        if informative and reference > p['criteria']['response_reference_min'] else 'unresolved',
                    recovery_time=record['trials'][f]['target_recovery_time']))
            rows.append(dict(target=target, arm=arm, intervention=record['intervention'],
                informative=bool(informative), initial_target_pair_separation=separation,
                target_late_log_rms_max=maximum, late_state_ratio=maximum/separation if separation > 0 else None,
                state_outcome=('near_baseline' if maximum/separation <= p['criteria']['late_state_ratio_max'] else 'shifted')
                    if informative else 'uninformative',
                stationary_stable=record['stationary_stable'], responses=response_rows))
    report = dict(context=context['key'], seed=context['seed'], background=context['background'],
        numerical_pass=bool(max(sham_errors) <= p['criteria']['sham_log_max']),
        sham_max_error=max(sham_errors), comparisons=rows)
    write_json(root/context['key']/'comparison.json', report)
    return report


def assess(root):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    reports = [compare_context(root, c, p) for c in p['contexts']]
    result = dict(completed=True, numerical_pass=all(r['numerical_pass'] for r in reports),
                  contexts=reports, independent_histories=p['independent_histories'],
                  protocol_sha256=digest(root/'protocol.json'), interpretation=p['interpretation'], limits=p['limits'])
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
            write_json(root/'status.json', dict(state='running', completed=done, total=p['total_jobs'], workers=workers))
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(worker, (root, c, job)): (c['key'], job_name(job))
                           for c in p['contexts'] for job in c['jobs']}
                for future in as_completed(futures):
                    record = future.result()
                    if not record['quality_pass']: raise RuntimeError('Independent solver disagreement')
                    done += 1
                    write_json(root/'status.json', dict(state='running', completed=done,
                        total=p['total_jobs'], workers=workers, last_completed='/'.join(futures[future])))
                    print('/'.join(futures[future]), 'completed', done, '/', p['total_jobs'], flush=True)
            result = assess(root)
            write_json(root/'status.json', dict(state='completed', completed=done, total=p['total_jobs'],
                                               numerical_pass=result['numerical_pass']))
        except Exception as exc:
            write_json(root/'status.json', dict(state='failed', completed=done, total=p['total_jobs'], error=str(exc)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/neighbor-context'))
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output)
    elif args.command == 'run': run(args.output, args.workers)
    else: assess(args.output)
