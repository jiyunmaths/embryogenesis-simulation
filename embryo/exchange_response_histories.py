"""Moving exchange/response replication on developmental histories 8 and 9.

Preserves the seed-7 assay, adds short native/GPU checks at each new starting
context, and executes the same mature protocol through the accepted GPU adapter.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import time

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .attribute_exchange import exchange
from .attribute_persistence import rhs, distance
from .benchmark_gpu_backend import cpu_audit
from .cell_exchange_response_moving import BACKGROUNDS, verify, assess as assess_responses
from .cell_response_exchange import checked_solve
from .cell_response_moving import initialize, name, observe
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .feedback_survival_validation import retime
from .gpu_backend import GpuSimulation
from .gpu_response_runner import require_validation, compatible, prepare as gpu_prepare, run as gpu_run
from .native_mechanics import NativeSimulation
from .resolution import write_json, _steps
from .validate_gpu_backend import CRITERIA, discrepancies


def read(path):
    return json.loads(Path(path).read_text())


def exchange_states(base, relaxed, ids, masses, expected_pair, expected_ids):
    """Select by initial activity only, then conserve both species' amounts."""
    base, relaxed, ids, masses = map(np.asarray, (base, relaxed, ids, masses))
    if (ids.ndim != 1 or len(ids) < 2 or len(np.unique(ids)) != len(ids) or
            masses.shape != ids.shape or base.shape != (2, len(ids)) or relaxed.shape != base.shape or
            not all(np.isfinite(x).all() for x in (base, relaxed, masses)) or
            any(np.any(x <= 0) for x in (base, relaxed, masses))):
        raise ValueError('Invalid initial exchange compartments')
    order = sorted(range(len(ids)), key=lambda i: (base[0, i], int(ids[i])))
    pair = [order[0], order[-1]]
    if pair != expected_pair or ids[pair].tolist() != expected_ids:
        raise ValueError('Prespecified target selection changed')
    fresh = exchange(base, masses, *pair)
    error = float(abs((fresh@masses)/(base@masses)-1).max())
    if error > 2e-14:
        raise ValueError('Exchange failed amount conservation')
    return fresh, pair, error


def prepare(root, validation=Path('outputs/gpu-backend-validation')):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError(root)
    accepted, evidence = require_validation(validation)
    exchange_root = Path('outputs/cell-response-exchange').resolve()
    ep = read(exchange_root/'protocol.json'); verify(ep)
    if read(exchange_root/'status.json')['state'] != 'completed':
        raise ValueError('Frozen exchange study must be complete')
    records, inputs = [], [exchange_root/'protocol.json', exchange_root/'status.json', *evidence]
    for seed in (8, 9):
        endpoint = Path(f'outputs/feedback-endpoint-bistability-seed-{seed}').resolve()
        survival = Path(f'outputs/feedback-survival-validation/seed-{seed}/survival').resolve()
        old = read(endpoint/'protocol.json')
        for f, h in old['sha256'].items():
            if digest(f) != h: raise ValueError('Endpoint evidence changed: '+f)
        if not read(endpoint/'results.json')['switch_on']['chemical_bistability_supported']:
            raise ValueError('Requires the previously validated patterned basin')
        survival_protocol = read(survival/'protocol.json')
        verify(survival_protocol)
        survival_status = read(survival/'status.json')
        if (survival_status['state'] != 'completed' or
                survival_status['protocol_sha256'] != digest(survival/'protocol.json') or
                not read(survival/'switch_on/result.json')['quality_pass']):
            raise ValueError('Source survival failed numerical quality')
        checkpoint = survival/'switch_on/latest_state.npz'
        host = AttributeSimulation.restore(checkpoint)
        if host.config.seed != seed or host.time != 150. or len(host.ids) != 16 or host.divisions:
            raise ValueError('Developmental source mismatch')
        comparison_config = type(host.config)(**vars(host.config))
        comparison_config.dt = accepted.dt
        compatible(comparison_config, accepted)
        graph = next(g for g in ep['graphs'] if g['key'] == f'seed-{seed}_switch_on')
        folder = exchange_root/graph['key']
        result = read(folder/'conservative.json')
        if (not result['settled'] or result['protocol_sha256'] != digest(exchange_root/'protocol.json') or
                digest(folder/'conservative.npz') != result['trajectory_sha256']):
            raise ValueError('Unvalidated pre-relaxed exchange')
        with np.load(endpoint/'switch_on.npz') as data:
            base, ids, masses, delta = data['trajectories'][0, -1].copy(), data['ids'], data['volumes'], data['delta']
        with np.load(folder/'conservative.npz') as data:
            relaxed = data['relaxation'][-1].copy()
            if not np.array_equal(ids, data['ids']) or not np.allclose(masses, data['masses'], rtol=1e-12, atol=1e-12):
                raise ValueError('Frozen relaxed state mismatch')
        live = host.signaling_graph()
        if (not np.array_equal(ids, host.ids) or not np.allclose(live.delta, delta, rtol=1e-12, atol=1e-12) or
                not np.allclose(live.masses, masses, rtol=1e-12, atol=1e-12)):
            raise ValueError('Frozen and mechanical source graphs differ')
        fresh, pair, amount_error = exchange_states(base, relaxed, ids, masses, graph['pair'], graph['ids'])
        records.append(dict(seed=seed, checkpoint=str(checkpoint), base=base, relaxed=relaxed,
                            fresh=fresh, ids=ids, masses=masses, delta=delta, pair=pair,
                            targets=ids[pair].tolist(), exchange_amount_error=amount_error,
                            raw_source_log_sd=float(np.std(np.log(host.activator))),
                            conditioned_log_sd=float(np.std(np.log(base[0]))), config=host.config))
        inputs.extend([checkpoint, survival/'protocol.json', survival/'status.json', survival/'switch_on/result.json',
                       endpoint/'protocol.json', endpoint/'results.json', endpoint/'switch_on.npz',
                       folder/'conservative.json', folder/'conservative.npz', folder/'comparison.json'])
    historical = Path('outputs/cell-exchange-response-moving').resolve()
    verify(read(historical/'protocol.json'))
    if read(historical/'status.json')['state'] != 'completed': raise ValueError('Seed-7 anchor incomplete')
    inputs.extend(historical/f for f in ('protocol.json', 'comparison.json', 'status.json'))
    sources = [Path(__file__), *[Path(__file__).with_name(f) for f in (
        'gpu_response_runner.py', 'gpu_backend.py', 'gpu_kernels.cu', 'cuda_mechanics.cu',
        'cuda_spatial_bench.cu', 'native_mechanics.py', 'native_mechanics.cpp', 'fast_mechanics.py',
        'fast_mechanics.c', 'attribute_exchange.py', 'attribute_persistence.py', 'attribute_development.py',
        'feedback_survival_validation.py', 'feedback_endpoint_bistability.py', 'benchmark_gpu_backend.py',
        'cell_response_exchange.py', 'cell_response_moving.py', 'cell_response.py',
        'cell_response_moving_refinement.py', 'cell_exchange_response_moving.py', 'validate_gpu_backend.py',
        'model.py', 'signaling.py', 'transport.py', 'polarity.py', 'resolution.py', 'feedback_long.py')]]
    hashes = {str(f.resolve()): digest(f) for f in sources}
    root.mkdir(parents=True)
    summaries = []
    for row in records:
        child = root/f"seed-{row['seed']}"/'formation-source'; child.mkdir(parents=True)
        sim = retime(row['checkpoint'], child/'source.npz', accepted.dt)
        if sim.time != 150. or not np.array_equal(sim.ids, row['ids']): raise ValueError('Retiming changed physical identity')
        np.savez_compressed(child/'initial_states.npz', unexchanged=row['base'], fresh_exchange=row['fresh'],
                            relaxed_exchange=row['relaxed'], ids=row['ids'], masses=row['masses'])
        frozen, errors = {}, {}
        c = row['config']; fun = rhs(row['delta'], c.signal_beta, c.signal_da, c.signal_dh)
        jac = lambda t, y: chemical_jacobian(y.reshape(2, -1), row['delta'], c.signal_beta, c.signal_da, c.signal_dh)
        times = np.arange(401)*.15
        for key, state in zip(BACKGROUNDS, (row['base'], row['fresh'], row['relaxed'])):
            frozen[key], errors[key] = checked_solve(fun, jac, state, times)
        np.savez_compressed(child/'frozen_reference.npz', times=times, **frozen)
        cp = dict(checkpoint=str(child/'source.npz'), start=150., dt=accepted.dt, duration=60., interval=.15,
                  jobs=[dict(family=k, target=None, factor=1.) for k in BACKGROUNDS],
                  scope='Matched moving formation/retention of untouched, fresh-conservative-exchange and pre-relaxed chemistry on one mature t=150 geometry. No division or prescribed identities.',
                  source_sha256=hashes, input_sha256={str(f): digest(f) for f in inputs+[child/'source.npz', child/'initial_states.npz', child/'frozen_reference.npz']})
        write_json(child/'protocol.json', cp)
        inputs.extend(child/f for f in ('protocol.json', 'source.npz', 'initial_states.npz', 'frozen_reference.npz'))
        summaries.append({k: row[k] for k in ('seed', 'targets', 'pair', 'exchange_amount_error', 'raw_source_log_sd', 'conditioned_log_sd')})
        summaries[-1]['frozen_reference_errors'] = errors
    protocol = dict(seeds=[8, 9], histories=summaries, backgrounds=BACKGROUNDS, factors=[.9, 1.1],
        validation=str(Path(validation).resolve()), historical_seed7=str(historical), dt=accepted.dt,
        formation_start=150., formation_duration=60., response_start=210., response_duration=60., interval=.15,
        prefix_duration=.6, prefix_native_threads=4, prefix_criteria=CRITERIA,
        criteria=dict(initial_pair_informative_min=.1, pair_return_ratio=.1, late_duration=24.,
                      response_reference_separation_min=.01, log_activator_sd_min=.1),
        scientific_jobs=36, prefix_checks=6,
        scope='Two additional pre-existing independent developmental histories. Same seed-7 assay: patterned equilibrium on the t=150 graph, untouched/fresh-conservative/pre-relaxed chemistry, moving formation to t=210; matched controls and +/-10% activator pulses in both prespecified cells through t=270. Six formation and thirty response continuations.',
        limits='Conditional on prior validated mature patterned basins; chemistry is equilibrated/transplanted at t=150, so this is not fresh zygote-to-identity formation. Three histories including the historical seed-7 anchor are a pilot, not population inference. Interventions/cells nested within histories. No assigned cell types, autonomous identity, full-developmental or spatial convergence claim. New exchange-formation and response timestep checks remain necessary.',
        source_sha256=hashes, input_sha256={str(f): digest(f) for f in inputs})
    write_json(root/'protocol.json', protocol)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=36, prefix_completed=0, prefix_total=6))
    return protocol


def prefix_check(root, seed, family, p):
    source = root/f'seed-{seed}'/'formation-source'
    folder = root/f'seed-{seed}'/'prefix'/family; folder.mkdir(parents=True, exist_ok=True)
    ph = digest(source/'protocol.json')
    if (folder/'comparison.json').exists():
        result = read(folder/'comparison.json')
        if not result['passed'] or result['source_protocol_sha256'] != ph: raise ValueError('Invalid saved prefix')
        for f, h in result['evidence_sha256'].items():
            if digest(f) != h: raise ValueError('Prefix evidence changed')
        return result
    cp = read(source/'protocol.json'); verify(cp)
    with np.load(source/'initial_states.npz') as data: state = data[family].copy()
    cpu = initialize(cp['checkpoint'], state, None, 1.); cpu.__class__ = NativeSimulation
    cpu.native_threads = p['prefix_native_threads']
    gpu = GpuSimulation(initialize(cp['checkpoint'], state, None, 1.))
    steps, every = _steps(p['prefix_duration'], cp['dt']), _steps(cp['interval'], cp['dt'])
    histories = [[], []]; errors = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                                               'relative_volume_max', 'relative_transport_max'), 0.)
    started = time.perf_counter()
    for step in range(steps+1):
        cpu_audit(cpu); gpu.audit()
        if step % every == 0:
            a, b = gpu.observe(step*cp['dt']), observe(cpu, step*cp['dt'])
            if max(a['boundary_occupancy'], b['boundary_occupancy']) >= .01: raise RuntimeError('Prefix boundary failure')
            for key, value in discrepancies(a, b).items(): errors[key] = max(errors[key], value)
            histories[0].append(a); histories[1].append(b)
            write_json(folder/'status.json', dict(state='running', elapsed=step*cp['dt'], errors=errors))
        if step == steps: break
        cpu.step(); gpu.step()
    phi_error = float(abs(gpu.phi.cpu().numpy()-cpu.phi).max())
    passed = all(v <= p['prefix_criteria'][k] for k, v in errors.items()) and phi_error <= p['prefix_criteria']['final_phi_abs_max']
    files = [folder/'gpu-history.json', folder/'cpu-history.json', folder/'gpu-end.npz', folder/'cpu-end.npz']
    write_json(files[0], histories[0]); write_json(files[1], histories[1]); gpu.checkpoint(files[2]); cpu.checkpoint(files[3])
    result = dict(passed=bool(passed), source_protocol_sha256=ph, errors=errors, phi_abs_max=phi_error,
                  wall_seconds=time.perf_counter()-started, evidence_sha256={str(f): digest(f) for f in files},
                  scope='Short same-state backend check at a new history/chemical context; not full-horizon equivalence or numerical convergence.')
    write_json(folder/'comparison.json', result); write_json(folder/'status.json', dict(state='completed', passed=bool(passed)))
    if not passed: raise RuntimeError('New-history native/GPU prefix disagreement')
    return result


def completed_history(child, job, times):
    child = Path(child); cp = read(child/'protocol.json'); verify(cp)
    folder = child/name(job); result = read(folder/'result.json')
    if not result['quality_pass'] or result['protocol_sha256'] != digest(child/'protocol.json'): raise ValueError('Invalid completed trajectory')
    h = read(folder/'history.json')
    with np.load(child/'initial_states.npz') as data: ids = data['ids']
    if (len(h) != len(times) or not np.allclose([r['elapsed'] for r in h], times, rtol=0, atol=1e-9) or
            any(not np.array_equal(r['ids'], ids) for r in h)): raise ValueError('Unaligned completed trajectory')
    path = np.array([r['chemistry'] for r in h])
    if path.shape != (len(times), 2, len(ids)) or not np.isfinite(path).all() or np.any(path <= 0):
        raise ValueError('Invalid chemistry')
    if not np.allclose([r['time'] for r in h], cp['start']+times, rtol=0, atol=1e-9):
        raise ValueError('Unaligned absolute trajectory clock')
    for row in h:
        if (np.shape(row['volumes']) != (len(ids),) or np.shape(row['polarity']) != (len(ids), 3) or
                np.shape(row['delta']) != (len(ids), len(ids)) or
                not np.isfinite(row['boundary_occupancy']) or not 0 <= row['boundary_occupancy'] < .01):
            raise ValueError('Invalid completed geometry')
        discrepancies(row, row)  # Independent finiteness/positivity check, despite worker quality flag.
    return h, path


def formation_assessment(root, seed, p):
    record = next(h for h in p['histories'] if h['seed'] == seed)
    source, child = root/f'seed-{seed}'/'formation-source', root/f'seed-{seed}'/'formation'
    times = np.arange(round(p['formation_duration']/p['interval'])+1)*p['interval']
    late = times >= p['formation_duration']-p['criteria']['late_duration']
    paths, histories = {}, {}
    for family in BACKGROUNDS:
        histories[family], paths[family] = completed_history(child, dict(family=family, target=None, factor=1.), times)
    with np.load(source/'initial_states.npz') as data:
        masses, initial = data['masses'], {k: data[k] for k in BACKGROUNDS}
    pair = record['pair']; base = paths['unexchanged']
    sep = float(distance(initial['fresh_exchange'][:, pair], initial['unexchanged'][:, pair], masses[pair]))
    informative = sep > p['criteria']['initial_pair_informative_min']
    donor = np.array([exchange(x, masses, *pair) for x in base])
    rows = []
    with np.load(source/'frozen_reference.npz') as frozen:
        if not np.allclose(frozen['times'], times, rtol=0, atol=1e-9):
            raise ValueError('Frozen reference clock mismatch')
        for family in BACKGROUNDS:
            path, end = paths[family], histories[family][-1]
            if (frozen[family].shape != path.shape or not np.isfinite(frozen[family]).all() or
                    np.any(frozen[family] <= 0) or not np.allclose(path[0], initial[family], rtol=1e-12, atol=0)):
                raise ValueError('Initial/frozen formation evidence mismatch')
            destination = float(distance(path[:, :, pair], base[:, :, pair], masses[pair])[late].max()/sep) if informative else None
            transferred = float(distance(path[:, :, pair], donor[:, :, pair], masses[pair])[late].max()/sep) if informative else None
            contrast = np.std(np.log(path[:, 0]), axis=1)
            rows.append(dict(family=family, informative_pair=informative, initial_pair_separation=sep,
                late_destination_ratio=destination, late_transferred_ratio=transferred,
                destination_like=bool(informative and destination < p['criteria']['pair_return_ratio']),
                transferred_like=bool(informative and transferred < p['criteria']['pair_return_ratio']),
                late_pair_order_reversed=bool(np.all(path[late, 0, pair[0]] > path[late, 0, pair[1]])),
                late_min_log_activator_sd=float(contrast[late].min()), final_log_activator_sd=float(contrast[-1]),
                contrast_retained=bool(contrast[late].min() > p['criteria']['log_activator_sd_min']),
                final_global_distance_from_unexchanged=float(distance(path[-1], base[-1], masses)),
                final_global_distance_from_frozen=float(distance(path[-1], frozen[family][-1], masses)),
                final_relative_operator_difference=float(np.linalg.norm(np.array(end['delta'])-histories['unexchanged'][-1]['delta'])/max(np.linalg.norm(histories['unexchanged'][-1]['delta']), 1e-30)),
                final_relative_axis_ratio_difference=float(abs(end['axis_ratio']/histories['unexchanged'][-1]['axis_ratio']-1))))
    report = dict(seed=seed, trials=rows, scope=p['scope'], limits=p['limits'])
    write_json(root/f'seed-{seed}'/'formation_comparison.json', report)
    return report


def materialize_response(root, seed, family, p):
    child = root/f'seed-{seed}'/'response-source'/family
    if (child/'protocol.json').exists(): verify(read(child/'protocol.json')); return child
    formation = root/f'seed-{seed}'/'formation'
    verify(read(formation/'protocol.json'))
    folder = formation/name(dict(family=family, target=None, factor=1.))
    result = read(folder/'result.json')
    if not result['quality_pass'] or result['protocol_sha256'] != digest(formation/'protocol.json'): raise ValueError('Incomplete formation')
    sim = AttributeSimulation.restore(folder/'latest_state.npz')
    history = read(folder/'history.json')
    if sim.time != p['response_start'] or history[-1]['ids'] != sim.ids.tolist() or not np.allclose(history[-1]['chemistry'], [sim.activator, sim.inhibitor], rtol=1e-12, atol=0):
        raise ValueError('Checkpoint/source history mismatch')
    targets = next(h['targets'] for h in p['histories'] if h['seed'] == seed)
    if child.exists(): raise ValueError('Incomplete prior response preparation')
    child.mkdir(parents=True); shutil.copy2(folder/'latest_state.npz', child/'source.npz')
    np.savez_compressed(child/'initial_states.npz', **{family: np.array([sim.activator, sim.inhibitor])}, ids=sim.ids, masses=sim.volumes())
    jobs = [dict(family=family, target=None, factor=1.)]+[dict(family=family, target=cell, factor=factor) for cell in targets for factor in p['factors']]
    files = [root/'protocol.json', formation/'protocol.json', folder/'result.json', folder/'history.json',
             folder/'latest_state.npz', child/'source.npz', child/'initial_states.npz']
    cp = dict(checkpoint=str(child/'source.npz'), start=p['response_start'], dt=p['dt'], duration=p['response_duration'],
              interval=p['interval'], jobs=jobs, scope=p['scope'], source_sha256=p['source_sha256'],
              input_sha256={str(f): digest(f) for f in files})
    write_json(child/'protocol.json', cp)
    return child


def assess(root):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    reports = []
    for seed in p['seeds']:
        formation = formation_assessment(root, seed, p)
        targets = next(h['targets'] for h in p['histories'] if h['seed'] == seed)
        rp = dict(backgrounds=BACKGROUNDS, targets=targets, factors=p['factors'], duration=p['response_duration'],
                  interval=p['interval'], scope=p['scope'], interpretation='Same-age donor/destination signed-response distances; descriptive nearest reference, not identity equivalence.', limitations=p['limits'])
        times = np.arange(round(p['response_duration']/p['interval'])+1)*p['interval']
        for family in BACKGROUNDS:
            child = root/f'seed-{seed}'/'response'/family
            for job in read(child/'protocol.json')['jobs']:
                completed_history(child, job, times)
        response = assess_responses(root/f'seed-{seed}'/'response', rp)
        counts = {key: sum(r['nearest_reference'] == key for r in response['comparisons']) for key in ('donor', 'destination', 'unresolved')}
        reports.append(dict(seed=seed, formation=formation, responses=response, descriptive_counts=counts))
    anchor = read(Path(p['historical_seed7'])/'comparison.json')
    report = dict(completed=True, new_histories=reports, historical_seed7_comparisons=anchor['comparisons'],
                  independent_histories=3, scope=p['scope'], limits=p['limits'], protocol_sha256=digest(root/'protocol.json'))
    write_json(root/'comparison.json', report)
    return report


def run(root):
    root = Path(root).resolve()
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); require_validation(p['validation']); torch.set_num_threads(1)
        prefixes = done = 0
        def status(stage, **extra):
            write_json(root/'status.json', dict(state='running', stage=stage, completed=done, total=36,
                prefix_completed=prefixes, prefix_total=6, backend='resident-gpu', **extra))
        try:
            for seed in p['seeds']:
                for family in BACKGROUNDS:
                    status('starting_state_backend_check', seed=seed, background=family)
                    prefix_check(root, seed, family, p); prefixes += 1
            for seed in p['seeds']:
                raw, output = root/f'seed-{seed}'/'formation-source', root/f'seed-{seed}'/'formation'
                if not (output/'protocol.json').exists(): gpu_prepare(output, raw, p['validation'])
                status('formation', seed=seed)
                gpu_run(output); done += 3
                formation_assessment(root, seed, p)
            for family in BACKGROUNDS:
                for seed in p['seeds']:
                    raw = materialize_response(root, seed, family, p)
                    output = root/f'seed-{seed}'/'response'/family
                    if not (output/'protocol.json').exists(): gpu_prepare(output, raw, p['validation'])
                    status('response', seed=seed, background=family)
                    gpu_run(output); done += 5
            assess(root)
            write_json(root/'status.json', dict(state='completed', completed=done, total=36,
                prefix_completed=prefixes, prefix_total=6, backend='resident-gpu'))
        except Exception as exc:
            write_json(root/'status.json', dict(state='failed', completed=done, total=36,
                prefix_completed=prefixes, prefix_total=6, error=str(exc)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/exchange-response-histories'))
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output)
    elif args.command == 'run': run(args.output)
    else: assess(args.output)
