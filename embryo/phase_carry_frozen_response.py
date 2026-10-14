"""Frozen t=510 chemical challenges on exact accepted measured operators.

Production calls the unchanged PyTorch SSP-RK2 kernel. Nine fixed contexts are
batched without communicating; every actual initial state also has independent
DOP853/Radau references and NumPy/SciPy production verification.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import json
import math
from multiprocessing import get_context
import os
from pathlib import Path
import time

import numpy as np
import scipy
from scipy import sparse
from scipy.integrate import solve_ivp
import torch

from .attribute_persistence import rhs
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .gpu_backend import gm_step
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .resolution import _steps
from .transport import Transport, integrate_gm

HISTORIES = (7, 8, 9)
LEVELS = (('coarse', .00375), ('fine', .001875), ('diagnostic', .0009375))
DOSES = ('fractional', 'amount')
CRITERIA = dict(raw_log_max=.01, normalized_response_max=.01,
    relative_auc_max=.02, recovery_time_max=.3, effect_error_max=.01,
    reference_log_max=1e-8, reference_normalized_max=1e-8,
    cpu_torch_log_max=1e-8, order_min=1.7, order_max=2.3,
    reference_noise_factor=100., pulse_amount_relative_max=2e-13)


def write_json(path, value):
    """Commit a report atomically; interrupted writes cannot authorize reuse."""
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def save_npz(path, **arrays):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp.npz')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def pin(paths):
    return {str(Path(f).resolve()): digest(f) for f in paths}


def accepted_parent(root, review_file):
    """Handle the compact public review plus its immutable archive manifest."""
    # Keep logical parent paths: the approved outputs alias resolves elsewhere,
    # while the independently accepted registry uses its original logical names.
    root, review_file = Path(root).absolute(), Path(review_file).absolute()
    p = read(root/'protocol.json'); review = read(review_file)
    originals = {key: root/file for key, file in (
        ('protocol_sha256', 'protocol.json'), ('original_summary_sha256', 'summary.json'),
        ('original_status_sha256', 'status.json'))}
    if (not review.get('full_study_accepted') or read(root/'status.json')['state'] != 'completed'
            or not read(root/'summary.json')['passed']):
        raise ValueError('Independently accepted completed moving parent required')
    for key, file in originals.items():
        if review[key] != digest(file):
            raise ValueError('Changed accepted parent: '+str(file))
    archive = Path(review['review_archive']).resolve()
    provenance_file = archive/'review-provenance.json'
    provenance = read(provenance_file)
    if not provenance['passed'] or any(provenance[key] != review[key] for key in originals):
        raise ValueError('Parent review provenance disagrees')
    inputs = [*originals.values(), review_file, provenance_file]
    # README/plan are historical reporting artifacts, never frozen model inputs.
    reporting_only = {Path('README.md').resolve(), Path('docs/plan.md').resolve()}
    for filename, expected in provenance['files_sha256'].items():
        file = Path(filename).resolve()
        if file in reporting_only:
            continue
        if digest(file) != expected:
            raise ValueError('Changed independent parent review artifact: '+filename)
        inputs.append(file)
    evidence_file = archive/'provenance_verification.json'
    evidence = read(evidence_file)
    if (not evidence['passed'] or evidence['protocol_sha256'] != review['protocol_sha256']
            or evidence.get('all_mismatches')):
        raise ValueError('Accepted immutable-evidence verification required')
    for filename, expected in p['source_sha256'].items():
        if digest(filename) != expected:
            raise ValueError('Changed inherited scientific source: '+filename)
    return p, review, evidence['evidence_sha256'], inputs


def parent_folder(root, p, job):
    item = p['reused'].get(job['key'])
    return Path(p['reuse_parent'])/item['job']['key'] if item else Path(root)/job['key']


def validate_context(delta, volumes, initial, ids):
    delta, volumes = validate_graph(delta, volumes)
    initial, ids = np.asarray(initial), np.asarray(ids)
    if (initial.dtype != np.float64 or initial.shape != (2, len(volumes))
            or not np.isfinite(initial).all() or np.any(initial <= 0)
            or ids.shape != volumes.shape or len(set(ids.tolist())) != len(ids)):
        raise ValueError('Positive float64 chemistry and unique matching IDs required')
    return delta, volumes, initial, ids


def trial_design(contexts):
    trials = []
    for context in contexts:
        for target, dose in [(None, None), *[(t, d) for t in context['targets'] for d in DOSES]]:
            suffix = 'control' if target is None else f'cell-{target}_{dose}'
            trials.append(dict(key=f'{context["key"]}_{suffix}', context=context['key'],
                seed=context['seed'], background=context['background'], target=target, dose=dose))
    return trials


def pulse_accounting(initial, changed, volumes, ids, target, record, common_amount):
    if target is None:
        if not np.array_equal(initial, changed):
            raise ValueError('Unperturbed chemistry changed')
        return 0.
    i = ids.tolist().index(target); expected = initial.copy()
    expected[0, i] *= record['factor']
    if not np.array_equal(changed, expected) or not 0 < record['factor'] < 1:
        raise ValueError('Pulse differs from exact accepted recipient-only preparation')
    before = float(initial[0, i]*volumes[i])
    removal = float((initial[0, i]-changed[0, i])*volumes[i])
    requested = .1*before if record['dose'] == 'fractional' else common_amount
    error = abs(removal/requested-1.)
    if error > CRITERIA['pulse_amount_relative_max']:
        raise ValueError('Pulse removal accounting failed')
    if (record['requested_removed_amount'] != requested or
            record['removed_activator_amount'] != removal or
            record['initial_target_activator_amount'] != before or
            record['initial_log_amplitude'] != float(abs(np.log(record['factor'])))):
        raise ValueError('Accepted pulse records disagree with actual source')
    return error


def substeps(delta, dt, beta, da, db):
    loss = max(1+da*float(-np.diag(delta).min()), beta+db*float(-np.diag(delta).min()))
    return max(1, math.ceil(dt*loss/.2))


def fixed_transport(delta, volumes):
    delta, volumes = validate_graph(delta, volumes)
    conductance = volumes[:, None]*delta; np.fill_diagonal(conductance, 0.)
    symmetric = -np.sqrt(volumes)[:, None]*delta/np.sqrt(volumes)[None, :]
    return Transport(volumes.copy(), sparse.csr_matrix(conductance),
        sparse.csr_matrix(delta), sparse.csr_matrix(symmetric))


def prepare(root, parent=Path('outputs/phase-carry-redistribution-response'),
            review_file=Path('docs/phase_carry_redistribution_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_frozen_response.md'), device='cuda:0'):
    root, parent = Path(root).resolve(), Path(parent).absolute()
    if root.exists():
        raise FileExistsError(root)
    old, review, evidence, inputs = accepted_parent(parent, review_file)
    contexts = []; old_starts = {}
    for seed in HISTORIES:
        targets = old['selections'][str(seed)]['selected_ids']
        for recipient in (None, *targets):
            background = 'sham' if recipient is None else f'redistributed-cell-{recipient}'
            prior = next(j for j in old['jobs'] if (j['seed'], j['level'], j['background'], j['target']) == (seed, 'fine', background, None))
            folder = parent_folder(parent, old, prior)
            consumed = [folder/'history.json', folder/'result.json', Path(prior['source'])]
            for file in consumed:
                if evidence.get(str(file), evidence.get(str(file.resolve()))) != digest(file):
                    raise ValueError('Source not present unchanged in accepted evidence: '+str(file))
            history, result = read(consumed[0]), read(consumed[1]); row = history[0]
            if (row['elapsed'] != 0. or row['time'] != 510. or not result['quality_pass']
                    or result['history_sha256'] != digest(consumed[0])
                    or result['source_sha256'] != digest(prior['source'])):
                raise ValueError('Actual challenge-time parent source/history required')
            initial = np.asarray(row['chemistry'], dtype=np.float64)
            delta, volumes, initial, ids = validate_context(row['delta'], row['volumes'], initial, row['ids'])
            with np.load(prior['source']) as z:
                meta = json.loads(str(z['metadata']))
                if (meta['time'] != 510. or not np.array_equal(z['ids'], ids)
                        or not np.array_equal([z['activator'], z['inhibitor']], initial)):
                    raise ValueError('Source checkpoint disagrees with t=510 observation')
                if (meta['config'] != old['accepted_configs'][f'{seed}_fine']
                        or any(meta['config'][k] != v for k, v in (
                            ('signal_beta', 2.), ('signal_da', .02), ('signal_dh', .55), ('polarity_tension', .35)))
                        or prior['point'] != dict(ratio=27.5, chi=.35)):
                    raise ValueError('Changed accepted frozen chemical/model parameters')
            context = dict(key=f'seed-{seed}_{background}', seed=seed, background=background,
                targets=targets if recipient is None else [recipient], source=str(root/'contexts'/f'seed-{seed}_{background}.npz'),
                parent_job=prior, physical_history=str(consumed[0].resolve()),
                original_source=str(Path(prior['source']).resolve()))
            contexts.append(context); old_starts[context['key']] = (delta, volumes, initial, ids)
            inputs += consumed
    trials = trial_design(contexts)
    if len(contexts) != 9 or len(trials) != 33:
        raise ValueError('Expected nine contexts and33 distinct starts')
    root.mkdir(parents=True)
    try:
        unique_problems = set()
        for context in contexts:
            delta, volumes, initial, ids = old_starts[context['key']]
            save_npz(context['source'], delta=delta, volumes=volumes, initial=initial, ids=ids)
            inputs.append(Path(context['source']))
            for trial in [t for t in trials if t['context'] == context['key']]:
                prior = next(j for j in old['jobs'] if (j['seed'], j['level'], j['background'], j['target'], j['dose']) ==
                    (trial['seed'], 'fine', trial['background'], trial['target'], trial['dose']))
                file = Path(prior['source'])
                if evidence.get(str(file), evidence.get(str(file.resolve()))) != digest(file):
                    raise ValueError('Unverified original pulse source')
                with np.load(file) as z:
                    changed = np.asarray([z['activator'], z['inhibitor']])
                    if not np.array_equal(z['ids'], ids) or json.loads(str(z['metadata']))['time'] != 510.:
                        raise ValueError('Pulse source IDs/clock changed')
                amount = 0. if trial['target'] is None else old['common_doses'][f'{trial["seed"]}_{trial["target"]}']['removed_amount']
                error = pulse_accounting(initial, changed, volumes, ids, trial['target'], prior, amount)
                trial.update({key:prior[key] for key in ('factor', 'initial_log_amplitude', 'removed_activator_amount',
                    'requested_removed_amount', 'initial_target_activator_amount')})
                trial.update(source=str(root/'starts'/(trial['key']+'.npz')), parent_job=prior,
                    original_source=str(file.resolve()), pulse_amount_relative_error=error)
                unique_problems.add((context['key'], delta.dtype.str, delta.shape, delta.tobytes(),
                    volumes.dtype.str, volumes.shape, volumes.tobytes(),
                    ids.dtype.str, ids.shape, ids.tobytes(),
                    changed.dtype.str, changed.shape, changed.tobytes()))
                save_npz(trial['source'], initial=changed); inputs += [file, Path(trial['source'])]
        beta, da, db = 2., .02, .55
        p = dict(parent=str(parent), parent_protocol_sha256=review['protocol_sha256'],
            review_file=str(Path(review_file).resolve()), contexts=contexts, trials=trials,
            histories=list(HISTORIES), selections=old['selections'], levels=[dict(level=l, dt=dt) for l, dt in LEVELS],
            beta=beta, da=da, db=db, start=510., duration=60., interval=.15, late_window=24., effect_min=.01,
            checkpoint_interval=3., independent_histories=3, new_histories=0,
            primary_paths=66, diagnostic_paths=33, total_production_paths=99,
            logical_reference_cases=33, unique_initial_value_problems=len(unique_problems),
            reference_trajectories=66, moving_paths_reused=0,
            device=device, runtime=dict(torch=old['device']['torch'], torch_cuda=old['device']['torch_cuda'],
                expected_gpu_name=old['device']['name'], numpy=np.__version__, scipy=scipy.__version__), criteria=CRITERIA,
            reference_settings=dict(rtol=1e-13, atol=1e-15, max_step=.15),
            scope='Exact nine t510 chemical/graph/volume backgrounds frozen through t570. No dilution or mechanical/polarity updates. Inherited backgrounds remain distinct; no shared-operator or matched-recipient intervention. Three existing histories,66primary+33diagnostic paths,33logical reference cases and66independent solver integrations. Identical actual starts across dose designs are retained and unique initial-value problems counted explicitly.')
        _, matrix, _, _ = batch_layout(p)
        for level, dt in LEVELS:
            counts = [substeps(old_starts[c['key']][0], dt, beta, da, db) for c in contexts]
            if len(set(counts)) != 1 or counts[0] != substeps(matrix, dt, beta, da, db):
                raise ValueError('Batch would change constituent positivity substeps')
        files = [Path(__file__), Path('tests/test_phase_carry_frozen_response.py'),
            *[Path(__file__).with_name(name) for name in ('attribute_persistence.py', 'cell_response.py',
                'cell_response_exchange.py', 'feedback_endpoint_bistability.py', 'feedback_long.py',
                'gpu_backend.py', 'neighbor_context.py', 'parameter_robustness.py', 'resolution.py',
                'transport.py', 'signaling.py')]]
        p['source_sha256'] = {**old['source_sha256'], **pin(files)}
        p['input_sha256'] = pin([*inputs, Path(protocol_doc)])
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', primary_paths=66, diagnostic_paths=33,
            completed_production_paths=0, independent_histories=3, new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error)))
        raise


def load_context(p, trial):
    context = next(c for c in p['contexts'] if c['key'] == trial['context'])
    if 'input_sha256' in p:
        for file in (context['source'], trial['source']):
            if p['input_sha256'][str(Path(file).resolve())] != digest(file):
                raise ValueError('Changed frozen preparation')
    with np.load(context['source']) as z, np.load(trial['source']) as q:
        delta, volumes, initial, ids = validate_context(z['delta'], z['volumes'], q['initial'], z['ids'])
        return context, delta.copy(), volumes.copy(), initial.copy(), ids.copy()


def unique_problem_count(p):
    """Count exact physical initial values within the retained context labels."""
    seen = set()
    for trial in p['trials']:
        _, delta, volumes, initial, ids = load_context(p, trial)
        seen.add((trial['context'], *[(x.dtype.str, x.shape, x.tobytes()) for x in (delta, volumes, ids, initial)]))
    return len(seen)


def batch_layout(p):
    contexts = p['contexts']; states = []; matrices = []; mapping = []
    width = max(sum(t['context'] == c['key'] for t in p['trials']) for c in contexts)
    for ci, context in enumerate(contexts):
        group = [t for t in p['trials'] if t['context'] == context['key']]
        _, delta, _, control, _ = load_context(p, group[0]); n = delta.shape[0]
        if ci and n != matrices[0].shape[0]:
            raise ValueError('Batch requires equal context cell counts')
        matrices.append(delta)
        columns = [load_context(p, t)[3] for t in group]
        columns += [control.copy() for _ in range(width-len(group))]
        states.append(np.stack(columns, axis=-1))
        mapping += [dict(trial=t['key'], context_index=ci, column=i) for i, t in enumerate(group)]
    matrix = sparse.block_diag(matrices).toarray()
    initial = np.concatenate(states, axis=1)
    return initial, matrix, mapping, matrices[0].shape[0]


def unpack(state, mapping, count):
    return np.stack([state[:, m['context_index']*count:(m['context_index']+1)*count, m['column']] for m in mapping])


def runtime_check(p):
    if str(torch.__version__) != p['runtime']['torch'] or torch.version.cuda != p['runtime']['torch_cuda']:
        raise ValueError('Changed PyTorch/CUDA runtime')
    if (p['runtime'].get('numpy', np.__version__) != np.__version__
            or p['runtime'].get('scipy', scipy.__version__) != scipy.__version__):
        raise ValueError('Changed NumPy/SciPy reference runtime')
    device = torch.device(p['device'])
    if device.type == 'cuda':
        if not torch.cuda.is_available() or torch.cuda.get_device_name(device) != p['runtime']['expected_gpu_name']:
            raise ValueError('Requested qualified CUDA device unavailable; no CPU fallback')
    elif device.type != 'cpu':
        raise ValueError('Only explicit CUDA or CPU devices supported')
    return device


def valid_state(state):
    if not np.isfinite(state).all() or np.any(state <= 0):
        raise FloatingPointError('Chemical stages require positive finite concentrations; no clipping')


def checked_torch_step(a, b, matrix, dt, p, count=None):
    # Independently inspect every forward-Euler stage in the unchanged SSP map.
    # The returned production value comes only from the original gm_step.
    if count is None:
        count = substeps(matrix.detach().cpu().numpy(), dt, p['beta'], p['da'], p['db'])
    step = dt/count; x, y = a, b
    for _ in range(count):
        fa = x*x/y-x+p['da']*(matrix@x)
        fb = p['beta']*(x*x-y)+p['db']*(matrix@y)
        first_a, first_b = x+step*fa, y+step*fb
        if not bool(torch.isfinite(first_a).all() & torch.isfinite(first_b).all() & (first_a > 0).all() & (first_b > 0).all()):
            raise FloatingPointError('SSP Euler stage lost positivity/finiteness')
        fa = first_a*first_a/first_b-first_a+p['da']*(matrix@first_a)
        fb = p['beta']*(first_a*first_a-first_b)+p['db']*(matrix@first_b)
        second_a, second_b = first_a+step*fa, first_b+step*fb
        if not bool(torch.isfinite(second_a).all() & torch.isfinite(second_b).all() & (second_a > 0).all() & (second_b > 0).all()):
            raise FloatingPointError('SSP second Euler stage lost positivity/finiteness')
        x, y = .5*x+.5*second_a, .5*y+.5*second_b
        if not bool(torch.isfinite(x).all() & torch.isfinite(y).all() & (x > 0).all() & (y > 0).all()):
            raise FloatingPointError('SSP final stage lost positivity/finiteness')
    out = gm_step(a, b, matrix, dt, p['beta'], p['da'], p['db'])
    if not bool(torch.isfinite(out[0]).all() & torch.isfinite(out[1]).all() & (out[0] > 0).all() & (out[1] > 0).all()):
        raise FloatingPointError('Production SSP state lost positivity/finiteness')
    return out


def reference_worker(args):
    root, trial = args; root = Path(root); p = read(root/'protocol.json'); ph = digest(root/'protocol.json')
    folder = root/'references'/trial['key']; report = folder/'result.json'; file = folder/'paths.npz'
    if report.exists():
        r = read(report)
        if r['protocol_sha256'] != ph or r['trial'] != trial or r['paths_sha256'] != digest(file):
            raise ValueError('Changed reference evidence')
        return r
    _, delta, _, initial, _ = load_context(p, trial)
    times = np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']
    fun = rhs(delta, p['beta'], p['da'], p['db'])
    jac = lambda t, y: chemical_jacobian(y.reshape(initial.shape), delta, p['beta'], p['da'], p['db'])
    started = time.perf_counter(); paths = []
    for method in ('DOP853', 'Radau'):
        opts = dict(jac=jac) if method == 'Radau' else {}
        sol = solve_ivp(fun, (0., p['duration']), initial.ravel(), method=method,
            t_eval=times, **p['reference_settings'], **opts)
        if not sol.success or sol.y.shape != (initial.size, len(times)):
            raise RuntimeError('Independent '+method+' solve failed')
        path = sol.y.T.reshape(len(times), *initial.shape); valid_state(path)
        if not np.array_equal(path[0], initial):
            raise ValueError('Independent reference changed initial state')
        paths.append(path)
    error = float(abs(np.log(paths[0]/paths[1])).max())
    save_npz(file, times=times, dop853=paths[0], radau=paths[1])
    r = dict(trial=trial, protocol_sha256=ph, paths_sha256=digest(file),
        raw_solver_log_error=error, passed=bool(error <= p['criteria']['reference_log_max']),
        wall_seconds=time.perf_counter()-started)
    write_json(report, r); return r


def cpu_worker(args):
    root, trial, level = args; root = Path(root); p = read(root/'protocol.json'); ph = digest(root/'protocol.json')
    folder = root/'cpu-verification'/level['level']/trial['key']; file = folder/'paths.npz'; report = folder/'result.json'
    if report.exists():
        r = read(report)
        if r['protocol_sha256'] != ph or r['trial'] != trial or r['level'] != level or r['paths_sha256'] != digest(file):
            raise ValueError('Changed CPU production verification')
        return r
    _, delta, volumes, initial, _ = load_context(p, trial)
    transport = fixed_transport(delta, volumes); a, b = initial.copy()
    steps, every = _steps(p['duration'], level['dt']), _steps(p['interval'], level['dt'])
    started = time.perf_counter(); paths = [initial]
    for index in range(1, steps+1):
        a, b = integrate_gm(a, b, transport, level['dt'], p['beta'], p['da'], p['db'])
        valid_state(np.array([a, b]))
        if index % every == 0:
            paths.append(np.array([a, b]))
    save_npz(file, times=np.arange(len(paths))*p['interval'], trajectory=np.array(paths))
    r = dict(trial=trial, level=level, protocol_sha256=ph, paths_sha256=digest(file),
        wall_seconds=time.perf_counter()-started)
    write_json(report, r); return r


def production(root, p, level):
    root = Path(root); ph = digest(root/'protocol.json'); folder = root/'production'/level['level']
    file, report, checkpoint = folder/'paths.npz', folder/'result.json', folder/'checkpoint.npz'
    if report.exists():
        r = read(report)
        if (r['protocol_sha256'] != ph or r['level'] != level or r['paths_sha256'] != digest(file)
                or r['checkpoint_sha256'] != digest(checkpoint)):
            raise ValueError('Changed completed production evidence')
        return r
    initial, delta, mapping, count = batch_layout(p); device = runtime_check(p)
    dt = level['dt']; total = _steps(p['duration'], dt); every = _steps(p['interval'], dt)
    save_every = _steps(p['checkpoint_interval'], dt)
    for c in p['contexts']:
        with np.load(c['source']) as z:
            if substeps(z['delta'], dt, p['beta'], p['da'], p['db']) != substeps(delta, dt, p['beta'], p['da'], p['db']):
                raise ValueError('Batched step changes positivity substeps')
    completed, wall, history = 0, 0., [unpack(initial, mapping, count)]
    state = initial
    if checkpoint.exists():
        commit = read(folder/'checkpoint-commit.json')
        if commit['protocol_sha256'] != ph or commit['checkpoint_sha256'] != digest(checkpoint):
            raise ValueError('Changed or incompletely committed production checkpoint; retain for review')
        with np.load(checkpoint) as z:
            meta = json.loads(str(z['metadata'])); completed = meta['step']; wall = meta['wall_seconds']
            state, history = z['state'].copy(), list(z['history'].copy())
            if (meta['protocol_sha256'] != ph or meta['level'] != level
                    or completed % save_every != 0 or not 0 <= completed <= total
                    or state.shape != initial.shape or np.asarray(history).shape != (completed//every+1, len(mapping), 2, count)
                    or len(history) != completed//every+1 or not np.array_equal(history[0], unpack(initial, mapping, count))
                    or not np.array_equal(history[-1], unpack(state, mapping, count))):
                raise ValueError('Checkpoint clock/preparation/history mismatch')
            valid_state(state); valid_state(np.array(history))
    a, b, matrix = [torch.as_tensor(x, dtype=torch.float64, device=device) for x in (state[0], state[1], delta)]
    batch_substeps = substeps(delta, dt, p['beta'], p['da'], p['db'])
    started = time.perf_counter()
    for step in range(completed, total+1):
        if step % save_every == 0:
            current = np.stack((a.detach().cpu().numpy(), b.detach().cpu().numpy()))
            meta = dict(protocol_sha256=ph, level=level, step=step, wall_seconds=wall+time.perf_counter()-started)
            save_npz(checkpoint, state=current, history=np.array(history), metadata=np.array(json.dumps(meta)))
            write_json(folder/'checkpoint-commit.json', dict(protocol_sha256=ph,
                checkpoint_sha256=digest(checkpoint), level=level, step=step))
            write_json(folder/'status.json', dict(state='running', elapsed=step*dt, horizon=p['duration']))
            print('Frozen production:', level['level'], step*dt, flush=True)
        if step == total:
            break
        a, b = checked_torch_step(a, b, matrix, dt, p, batch_substeps)
        if (step+1) % every == 0:
            history.append(unpack(np.stack((a.detach().cpu().numpy(), b.detach().cpu().numpy())), mapping, count))
    paths = np.array(history); times = np.arange(len(paths))*p['interval']
    if paths.shape != (len(times), len(mapping), 2, count) or len(times) != _steps(p['duration'], p['interval'])+1:
        raise ValueError('Incomplete batched production')
    save_npz(file, times=times, trajectory=paths, trial_keys=np.array([m['trial'] for m in mapping]))
    r = dict(protocol_sha256=ph, level=level, paths_sha256=digest(file), checkpoint_sha256=digest(checkpoint),
        actual_paths=len(mapping), padding_paths=len(p['contexts'])*initial.shape[-1]-len(mapping),
        wall_seconds=wall+time.perf_counter()-started)
    write_json(report, r); write_json(folder/'status.json', dict(state='completed', elapsed=p['duration']))
    return r


def reference_path_sets(root, p):
    root = Path(root); ph = digest(root/'protocol.json'); paths = {}; reference_errors = {}; evidence = {}
    for trial in p['trials']:
        folder = root/'references'/trial['key']; r = read(folder/'result.json'); file = folder/'paths.npz'
        if r['protocol_sha256'] != ph or r['trial'] != trial or r['paths_sha256'] != digest(file):
            raise ValueError('Changed independent reference')
        with np.load(file) as z:
            for label in ('dop853', 'radau'):
                paths[label, trial['key']] = z[label].copy()
            times = z['times'].copy()
            if not np.array_equal(times, np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']):
                raise ValueError('Changed individual reference clock')
        reference_errors[trial['key']] = float(abs(np.log(paths['dop853', trial['key']]/paths['radau', trial['key']])).max())
        evidence.update(pin([file, folder/'result.json']))
    expected_times = np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']
    if not np.array_equal(times, expected_times):
        raise ValueError('Changed reference clock')
    return paths, times, reference_errors, evidence


def reference_qualification(root, p):
    paths, times, errors, evidence = reference_path_sets(root, p)
    summaries = {label:summarize_responses(paths, label, p, times) for label in ('dop853', 'radau')}
    responses = compare_responses(paths, 'dop853', 'radau', p, times, summaries, p['criteria']['reference_normalized_max'])
    result = dict(protocol_sha256=digest(Path(root)/'protocol.json'),
        passed=bool(all(v <= p['criteria']['reference_log_max'] for v in errors.values()) and responses['passed']),
        raw_errors=errors, responses=responses, evidence_sha256=evidence)
    write_json(Path(root)/'reference-qualification.json', result)
    return result


def path_sets(root, p):
    root = Path(root); ph = digest(root/'protocol.json')
    paths, times, reference_errors, evidence = reference_path_sets(root, p)
    for level in p['levels']:
        folder = root/'production'/level['level']; file = folder/'paths.npz'; r = read(folder/'result.json')
        if (r['protocol_sha256'] != ph or r['level'] != level or r['paths_sha256'] != digest(file)
                or r['checkpoint_sha256'] != digest(folder/'checkpoint.npz')):
            raise ValueError('Changed production evidence')
        with np.load(file) as z:
            if not np.array_equal(z['times'], times) or set(z['trial_keys']) != {t['key'] for t in p['trials']}:
                raise ValueError('Changed production clock/keys')
            for key, path in zip(z['trial_keys'], z['trajectory'][:, :].transpose(1, 0, 2, 3)):
                paths[level['level'], str(key)] = path.copy()
        evidence.update(pin([file, folder/'result.json', folder/'checkpoint.npz', folder/'checkpoint-commit.json']))
        for trial in p['trials']:
            folder = root/'cpu-verification'/level['level']/trial['key']; file = folder/'paths.npz'; r = read(folder/'result.json')
            if r['protocol_sha256'] != ph or r['trial'] != trial or r['level'] != level or r['paths_sha256'] != digest(file):
                raise ValueError('Changed CPU verification evidence')
            with np.load(file) as z:
                if not np.array_equal(z['times'], times):
                    raise ValueError('Changed CPU verification clock')
                paths['cpu-'+level['level'], trial['key']] = z['trajectory'].copy()
            evidence.update(pin([file, folder/'result.json']))
    for (_, key), path in paths.items():
        valid_state(path)
        with np.load(next(t['source'] for t in p['trials'] if t['key'] == key)) as z:
            if path.shape != (len(times), *z['initial'].shape) or not np.array_equal(path[0], z['initial']):
                raise ValueError('Changed initial chemistry')
    return paths, times, reference_errors, evidence


def summarize_responses(paths, label, p, times):
    metrics, waves = {}, {}
    for trial in p['trials']:
        if trial['target'] is None:
            continue
        control = next(t for t in p['trials'] if t['context'] == trial['context'] and t['target'] is None)
        _, _, volumes, _, ids = load_context(p, trial); i = ids.tolist().index(trial['target'])
        values, baseline = paths[label, trial['key']], paths[label, control['key']]
        wave = np.log(values[:, :, i]/baseline[:, :, i])/abs(np.log(trial['factor']))
        waves[trial['key']] = wave
        metrics[trial['key']] = response_metrics(values, baseline, times, volumes, i, trial['factor'])
    effects = []
    for seed in p['histories']:
        for target in p['selections'][str(seed)]['selected_ids']:
            for dose in DOSES:
                sham, red = [next(t for t in p['trials'] if (t['seed'], t['background'], t['target'], t['dose']) ==
                    (seed, background, target, dose)) for background in ('sham', f'redistributed-cell-{target}')]
                effect = waveform_distance(waves[red['key']], waves[sham['key']], times)
                effects.append(dict(seed=seed, cell=target, dose=dose, normalized_response_shift_rms=effect,
                    response_shift_detected=bool(effect > p['effect_min'])))
    return dict(metrics=metrics, effects=effects)


def compare_responses(paths, a, b, p, times, summaries, limit=None):
    criteria = p['criteria']; rows = []
    for trial in p['trials']:
        if trial['target'] is None:
            continue
        control = next(t for t in p['trials'] if t['context'] == trial['context'] and t['target'] is None)
        deviations = [np.log(paths[label, trial['key']]/paths[label, control['key']]) for label in (a, b)]
        error = float(abs(deviations[0]-deviations[1]).max()/abs(np.log(trial['factor'])))
        m, n = [summaries[label]['metrics'][trial['key']] for label in (a, b)]
        auc = abs(n['target_activator_log_auc_per_log_pulse']/m['target_activator_log_auc_per_log_pulse']-1.)
        recovery = {}
        for metric in ('target_recovery_time', 'network_recovery_time'):
            x, y = m[metric], n[metric]
            recovery[metric] = 0. if x is None and y is None else None if x is None or y is None else abs(x-y)
        passed = (error <= (limit if limit is not None else criteria['normalized_response_max'])
            and auc <= criteria['relative_auc_max'] and all(v is not None and v <= criteria['recovery_time_max']+1e-9 for v in recovery.values()))
        rows.append(dict(trial=trial['key'], normalized_response_error=error, relative_auc_error=float(auc),
            recovery_time_errors=recovery, passed=bool(passed)))
    effects = []
    for x, y in zip(summaries[a]['effects'], summaries[b]['effects']):
        same = x['response_shift_detected'] == y['response_shift_detected']
        error = abs(x['normalized_response_shift_rms']-y['normalized_response_shift_rms'])
        effects.append(dict(seed=x['seed'], cell=x['cell'], dose=x['dose'], same_decision=same,
            rms_error=error, passed=bool(same and error <= criteria['effect_error_max'])))
    return dict(first=a, second=b, passed=bool(all(r['passed'] for r in rows+effects)), pulse_rows=rows, effects=effects)


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    paths, times, reference_errors, evidence = path_sets(root, p)
    labels = ['dop853', 'radau', *[l['level'] for l in p['levels']]]
    summaries = {label:summarize_responses(paths, label, p, times) for label in labels}
    comparisons = [compare_responses(paths, 'dop853', 'radau', p, times, summaries, p['criteria']['reference_normalized_max'])]
    comparisons += [compare_responses(paths, 'coarse', 'fine', p, times, summaries)]
    comparisons += [compare_responses(paths, 'fine', 'diagnostic', p, times, summaries)]
    comparisons += [compare_responses(paths, l['level'], method, p, times, summaries) for l in p['levels'] for method in ('dop853', 'radau')]
    errors = []; order_rows = []; raw_pairs = []
    for trial in p['trials']:
        for first, second in (('coarse', 'fine'), ('fine', 'diagnostic')):
            pair_error = float(abs(np.log(paths[first, trial['key']]/paths[second, trial['key']])).max())
            raw_pairs.append(dict(trial=trial['key'], first=first, second=second,
                raw_log_error=pair_error, passed=bool(pair_error <= p['criteria']['raw_log_max'])))
        values = []
        for level in p['levels']:
            name = level['level']
            error = max(float(abs(np.log(paths[name, trial['key']]/paths[ref, trial['key']])).max()) for ref in ('dop853', 'radau'))
            cpu = float(abs(np.log(paths[name, trial['key']]/paths['cpu-'+name, trial['key']])).max())
            values.append(error)
            errors.append(dict(trial=trial['key'], level=name, raw_reference_log_error=error, cpu_torch_log_error=cpu,
                passed=bool(error <= p['criteria']['raw_log_max'] and cpu <= p['criteria']['cpu_torch_log_max'])))
        noise = reference_errors[trial['key']]
        resolved = values[-1] > p['criteria']['reference_noise_factor']*noise
        orders = np.log2(np.array(values[:-1])/values[1:]).tolist() if resolved and all(v > 0 for v in values) else None
        order_pass = not resolved or (orders is not None and all(p['criteria']['order_min'] <= v <= p['criteria']['order_max'] for v in orders))
        order_rows.append(dict(trial=trial['key'], errors=values, reference_noise=noise,
            order_resolved=bool(resolved), orders=orders, passed=bool(order_pass)))
    history_outcomes = []
    for seed in p['histories']:
        by_cell = []
        for cell in p['selections'][str(seed)]['selected_ids']:
            row = [r for r in summaries['fine']['effects'] if r['seed'] == seed and r['cell'] == cell]
            by_cell.append(dict(cell=cell, shifted_at_both_doses=bool(all(r['response_shift_detected'] for r in row))))
        history_outcomes.append(dict(seed=seed, recipients=2, both_doses_above_screen=sum(r['shifted_at_both_doses'] for r in by_cell), cells=by_cell))
    reference_pass = all(e <= p['criteria']['reference_log_max'] for e in reference_errors.values())
    result = dict(protocol_sha256=digest(root/'protocol.json'),
        passed=bool(reference_pass and all(r['passed'] for r in errors+order_rows+comparisons+raw_pairs)),
        primary_paths=66, diagnostic_paths=33, total_production_paths=99,
        logical_reference_cases=33, unique_initial_value_problems=p.get('unique_initial_value_problems'),
        reference_trajectories=66, independent_histories=3, new_histories=0,
        references_raw_pass=reference_pass, reference_errors=reference_errors,
        production_accuracy=errors, orders=order_rows, raw_timestep_comparisons=raw_pairs, response_comparisons=comparisons,
        response_summaries=summaries, history_outcomes=history_outcomes, evidence_sha256=evidence, scope=p['scope'])
    write_json(root/'summary.json', result); return result


def run(root, workers=4):
    root = Path(root).resolve()
    if not 1 <= workers <= 4:
        raise ValueError('Use one to four independent CPU verification workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); runtime_check(p); torch.set_num_threads(1)
        # Spawned workers import NumPy/SciPy afresh with one native thread.
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ[name] = '1'
        try:
            write_json(root/'status.json', dict(state='running', stage='independent_references', completed=0, total=33))
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                pending = [pool.submit(reference_worker, (root, t)) for t in p['trials']]
                for i, future in enumerate(as_completed(pending), 1):
                    r = future.result()
                    if not r['passed']:
                        raise RuntimeError('Dual independent reference gate failed')
                    write_json(root/'status.json', dict(state='running', stage='independent_references', completed=i, total=33))
            if not reference_qualification(root, p)['passed']:
                raise RuntimeError('Normalized independent reference qualification failed; production blocked')
            for level in p['levels']:
                write_json(root/'status.json', dict(state='running', stage='production', level=level['level'], total_production_paths=99))
                production(root, p, level)
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                pending = [pool.submit(cpu_worker, (root, t, l)) for l in p['levels'] for t in p['trials']]
                for i, future in enumerate(as_completed(pending), 1):
                    future.result()
                    write_json(root/'status.json', dict(state='running', stage='cpu_verification', completed=i, total=99))
            result = assess(root)
            write_json(root/'status.json', dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'], total_production_paths=99, primary_paths=66, diagnostic_paths=33,
                independent_histories=3, new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('output/studies/phase-carry-frozen-response'))
    parser.add_argument('--parent', type=Path, default=Path('outputs/phase-carry-redistribution-response'))
    parser.add_argument('--review', type=Path, default=Path('docs/phase_carry_redistribution_response_assessment.json'))
    parser.add_argument('--protocol-doc', type=Path, default=Path('docs/phase_carry_frozen_response.md'))
    parser.add_argument('--device', default='cuda:0'); parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    result = prepare(args.output, args.parent, args.review, args.protocol_doc, args.device) if args.action == 'prepare' else run(args.output, args.workers) if args.action == 'run' else assess(args.output)
    print({key:value for key, value in result.items() if key in ('passed', 'primary_paths', 'diagnostic_paths', 'total_production_paths', 'independent_histories')})
