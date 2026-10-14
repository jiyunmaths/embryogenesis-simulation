"""Full-horizon formation from inherited fresh developmental chemistry.

No solver/kernel changes. Keep original short qualification and failed studies
immutable; test all 24 mature paths through physical t78 at matched clocks.
"""
import argparse
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import os
import time

import numpy as np
import torch

from . import fresh_zygote_mature_qualification as f
from .fresh_zygote_mature_clock_fix import elapsed_steps
from .qualify_development_gpu import digest, read, write, pin, row, exact

PARENT = Path('output/studies/fresh-zygote-mature-clock-fix')
PARENT_SHA = '057f18713e37c562cc70743dd0dbff37be6e76b53040455a4535452752c19b04'
ASSESSMENT = Path('docs/fresh_zygote_mature_clock_fix_assessment.json')
DESIGN = Path('docs/fresh_zygote_formation.md')
LAUNCH = Path('archive/study-launches/2026-10-12-fresh-zygote-formation')
START = 18.
DURATION = 60.
SNAPSHOTS = (0, 12, 24, 36, 48, 60)
RESTARTS = (12, 36)
THRESHOLDS = (.05, .1, .2)


def accepted_parent():
    p = read(PARENT / 'protocol.json')
    a = read(ASSESSMENT)
    if digest(PARENT / 'protocol.json') != PARENT_SHA:
        raise ValueError('Immutable accepted mature-clock protocol required')
    record = a['independent_review']
    if (not a['review_completed'] or not a['scoped_mature_timestep_accepted'] or
            digest(record['path']) != record['sha256'] or
            not read(record['path'])['scoped_mature_timestep_accepted']):
        raise ValueError('Independent mature-clock acceptance required')
    return p, a


def design(root, parent):
    jobs = []
    for old in parent['jobs']:
        cfg = dict(old['config'], steps=elapsed_steps(START + DURATION, old['config']['dt']))
        jobs.append(dict(old, source_checkpoint=old['source'],
                         source_case=old['key'], source=str(Path(root) / old['key'] / 'source.npz'),
                         config=cfg))
    return jobs


def checked_transfer(source, destination, cfg):
    """Change only the declared horizon, preserving every inherited array byte."""
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError('Never overwrite a registered start')
    with np.load(source, allow_pickle=False) as z:
        data = {name: z[name].copy() for name in z.files}
    metadata = json.loads(str(data['metadata']))
    old = metadata['config']
    if (metadata['time'] != START or metadata['step_number'] * old['dt'] != START or
            metadata['divisions'] or len(data['ids']) != 16 or
            len(json.loads(str(data['event_audits']))) != 15 or
            str(data['chemical_method']) != f.td.gpu.METHOD or
            int(data['development_schema']) != f.td.gpu.SCHEMA):
        raise ValueError('Complete accepted t18 checkpoint required')
    if cfg != dict(old, steps=elapsed_steps(START + DURATION, old['dt'])):
        raise ValueError('Only the declared final step horizon may change')
    metadata['config'] = cfg
    data['metadata'] = np.array(json.dumps(metadata))
    np.savez_compressed(destination, **data)
    with np.load(destination, allow_pickle=False) as z:
        current = {name: z[name] for name in z.files}
    if set(data) != set(current):
        raise ValueError('Transferred checkpoint schema changed')
    for name in data:
        if (data[name].dtype != current[name].dtype or data[name].shape != current[name].shape or
                data[name].tobytes() != current[name].tobytes()):
            raise ValueError('Transferred checkpoint changed array: ' + name)
    return len(data) - 1


def verify(p):
    f.q.verify(p)
    parent, _ = accepted_parent()
    if p['parent_protocol_sha256'] != PARENT_SHA or p['jobs'] != design(p['study_root'], parent):
        raise ValueError('Registered inherited endpoint design changed')
    inherited = ('histories', 'branches', 'levels', 'time_pairs', 'physical_limits', 'time_limits',
                 'homologous_limits', 'dense_time_spectrum_limit', 'prefix_steps', 'prefix_limits',
                 'method', 'gpu_uuid', 'runtime_binary_sha256', 'interval', 'diagnostic_interval')
    if any(p[k] != parent[k] for k in inherited):
        raise ValueError('Inherited numerical/model contract changed')
    for mapping in ('source_sha256', 'input_sha256'):
        if any(p[mapping].get(k) != v for k, v in parent[mapping].items()):
            raise ValueError('Original evidence pin changed')
    if (p['start_time'] != START or p['duration'] != DURATION or p['snapshots'] != list(SNAPSHOTS) or
            p['restart_elapsed'] != list(RESTARTS) or p['restart_steps'] != 8 or
            p['thresholds'] != list(THRESHOLDS) or p['sustained_window'] != 24. or
            p['late_window'] != 24. or p['onset_time_limit'] != .3 or p['growth_abs_limit'] != .001 or
            p['new_scientific_histories'] != 0 or p['new_numerical_paths'] != 24):
        raise ValueError('Unregistered formation scope or decision rule')


def prepare(root):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('Fresh study directory required')
    parent, a = accepted_parent()
    f.q.verify(parent)
    root.mkdir(parents=True)
    jobs = design(root, parent)
    arrays = 0
    for job in jobs:
        (root / job['key']).mkdir()
        arrays += checked_transfer(job['source_checkpoint'], job['source'], job['config'])
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__), DESIGN, Path('tests/test_fresh_zygote_formation.py'),
                        LAUNCH / 'launch.py', LAUNCH / 'audit_prelaunch.py']))
    inputs = dict(parent['input_sha256'])
    inputs.update(pin([ASSESSMENT, a['independent_review']['path'], PARENT / 'protocol.json',
                       PARENT / 'status.json', PARENT / 'summary.json',
                       *[job['source'] for job in jobs]]))
    p = deepcopy(parent)
    for key in ('reused_prefix_case', 'scheduling_fix'):
        p.pop(key, None)
    p.update(jobs=jobs, study_root=str(root), parent_protocol_sha256=PARENT_SHA,
             duration=DURATION, snapshots=list(SNAPSHOTS), restart_elapsed=list(RESTARTS),
             thresholds=list(THRESHOLDS), sustained_window=24., late_window=24.,
             onset_time_limit=.3, growth_abs_limit=.001, source_sha256=sources, input_sha256=inputs,
             stop_policy='Stop on any prefix, physical, restart or full-horizon pair failure. '
                         'Preserve failed/unexecuted cases. No automatic timestep fallback, threshold '
                         'relaxation, history replacement, rerun or 240-unit extension.',
             scope='Three existing fresh developmental histories, four matched inherited contexts, '
                   'two mature timesteps, physical t18–78. Test sustained chemical formation while '
                   'mechanics/transport co-evolve, with no chemistry reset. Numerical acceptance '
                   'requires independent review. No biological identity or long maintenance claim.')
    write(root / 'protocol.json', p)
    write(root / 'transfer-verification.json', dict(passed=True, checked_arrays=arrays,
        changed_metadata=['config.steps'], all24_starts_lossless=True,
        protocol_sha256=digest(root / 'protocol.json')))
    verify(p)
    write(root / 'status.json', dict(state='prepared', protocol_sha256=digest(root / 'protocol.json')))
    return p


def growth(values, cfg):
    values = np.asarray(values)
    beta = cfg['signal_beta']
    jac = np.array([[1., -1.], [2 * beta, -beta]])
    rates = [np.linalg.eigvals(jac - ell * np.diag([cfg['signal_da'], cfg['signal_dh']])).real.max()
             for ell in values if ell > 1e-10]
    return dict(frozen_uniform_max_growth=float(max(rates)) if rates else None,
                frozen_uniform_unstable_modes=sum(bool(rate > 1e-10) for rate in rates))


def formation(history, p):
    times = np.asarray([x['elapsed'] for x in history])
    sd = np.asarray([x['log_activator_sd'] for x in history])
    if (len(times) != 401 or abs(times[0]) > 1e-9 or abs(times[-1] - DURATION) > 1e-9 or
            not np.allclose(times, np.arange(401) * .15, rtol=0., atol=1e-9) or
            not np.isfinite(sd).all() or sd.min() < 0):
        raise ValueError('Complete finite formation history required')
    outcomes = []
    for threshold in p['thresholds']:
        above = sd > threshold
        crossed = np.flatnonzero(above)
        onset = None
        for i, stamp in enumerate(times):
            if stamp + p['sustained_window'] > times[-1] + 1e-9:
                break
            interval = (times >= stamp - 1e-9) & (times <= stamp + p['sustained_window'] + 1e-9)
            if above[interval].all():
                onset = float(stamp)
                break
        late = times >= DURATION - p['late_window'] - 1e-9
        transitions = np.flatnonzero(above[1:] != above[:-1]) + 1
        outcomes.append(dict(threshold=threshold,
            first_crossing_elapsed=float(times[crossed[0]]) if len(crossed) else None,
            first_sustained_onset_elapsed=onset, sustained_formation=onset is not None,
            late_persistent=bool(above[late].all()), late_min_log_sd=float(sd[late].min()),
            final_log_sd=float(sd[-1]), sampled_transitions=[dict(elapsed=float(times[i]),
                direction='up' if above[i] else 'down') for i in transitions]))
    return outcomes


def run_path(root, job, p):
    folder = Path(root) / job['key']
    output = folder / 'gpu-result.json'
    if output.exists():
        raise FileExistsError('Never rerun a recorded formation path')
    sim = f.td.gpu.DevelopmentGpu.restore(job['source'])
    if dict(sim.config.__dict__) != job['config']:
        raise ValueError('Restored start config does not match frozen job')
    start_step = sim.step_number
    initial = f.invariant(sim)
    initial_target = float(sim.target.sum())
    history, checkpoints = [], []
    failure = None
    maxima = dict(max_volume_error=0., min_radius=1e100, max_clipping=0., dilution_amount_error=0.)
    clone, remaining, active_restart = None, 0, None
    restarted = {}
    trace_rows = 0
    started = time.perf_counter()
    dt = job['config']['dt']
    obs_steps = elapsed_steps(p['interval'], dt)
    trace_steps = elapsed_steps(p['diagnostic_interval'], dt)
    snapshot_steps = {elapsed_steps(x, dt): x for x in p['snapshots']}
    restart_steps = {elapsed_steps(x, dt): x for x in p['restart_elapsed']}
    try:
        with (folder / 'network.jsonl').open('x') as network:
            for k in range(elapsed_steps(DURATION, dt) + 1):
                if sim.step_number != start_step + k or sim.time != (start_step + k) * dt:
                    raise ValueError('Formation continuation clock drift')
                now = sim.audit()
                for name in maxima:
                    maxima[name] = (min if name == 'min_radius' else max)(maxima[name], now[name])
                if not f.q.path_quality(maxima, p['physical_limits']):
                    raise FloatingPointError('Every-step physical gate failed')
                u = torch.stack((sim.activator, sim.inhibitor))
                if not bool(torch.isfinite(u).all()) or float(u.min()) <= 0:
                    raise FloatingPointError('Chemical positivity failed')
                if sim.host.divisions or len(sim.ids) != 16 or len(sim.event_audits) != 15:
                    raise ValueError('Unexpected division during formation')
                if abs(float(sim.target.sum()) / initial_target - 1) > 2e-14:
                    raise ValueError('Inherited target drift')
                if k % trace_steps == 0:
                    f.td.json_line(network, f.td.network_trace(sim))
                    trace_rows += 1
                if k % obs_steps == 0:
                    if f.invariant(sim) != initial:
                        raise ValueError('Lineage/cycle/random state changed without division')
                    observed = row(sim, True)
                    observed.update(elapsed=k * dt, parents=sim.host.parents.tolist())
                    observed['spectrum'] = f.q.spectrum(observed['volumes'], observed['delta']).tolist()
                    observed.update(growth(observed['spectrum'], job['config']))
                    history.append(observed)
                    if (observed['boundary_occupancy'] >= p['physical_limits']['boundary_max'] or
                            any(v != 1 for v in observed['nondividing_components'])):
                        raise FloatingPointError('Boundary/connectivity gate failed')
                    write(folder / 'gpu-history.json', history)
                    network.flush()
                    write(folder / 'gpu-progress.json', dict(state='running', time=sim.time,
                        elapsed=k * dt, cells=16, seconds=time.perf_counter() - started, audit=maxima))
                if k in snapshot_steps:
                    cp = folder / f'gpu-e{int(snapshot_steps[k]):02d}.npz'
                    sim.checkpoint(cp)
                    checkpoints.append(cp)
                if k in restart_steps:
                    if clone is not None:
                        raise ValueError('Overlapping registered restart windows')
                    active_restart = int(restart_steps[k])
                    cp = folder / f'restart-e{active_restart:02d}.npz'
                    sim.checkpoint(cp)
                    clone = f.td.gpu.DevelopmentGpu.restore(cp)
                    exact(sim, clone)
                    remaining = p['restart_steps']
                if k == elapsed_steps(DURATION, dt):
                    break
                sim.step()
                if clone is not None:
                    clone.step()
                    exact(sim, clone)
                    remaining -= 1
                    if remaining == 0:
                        restarted[str(active_restart)] = p['restart_steps']
                        clone = None
        if (restarted != {str(int(x)): p['restart_steps'] for x in p['restart_elapsed']} or
                sim.time != START + DURATION or f.invariant(sim) != initial):
            raise ValueError('Formation horizon or exact restarts incomplete')
    except Exception as exc:
        failure = dict(type=type(exc).__name__, message=str(exc), time=sim.time, elapsed=sim.time - START)
        sim.checkpoint(folder / 'failure-state.npz')
    write(folder / 'gpu-history.json', history)
    report = dict(passed=failure is None, job=job,
        protocol_sha256=digest(Path(root) / 'protocol.json'), failure=failure,
        final_time=sim.time, final_cells=len(sim.ids), events=sim.event_audits,
        audit=maxima, exact_restart_steps=restarted, lineage_cycle_rng_unchanged=f.invariant(sim) == initial,
        network_trace_rows=trace_rows, seconds=time.perf_counter() - started,
        formation=formation(history, p) if failure is None else None,
        evidence_sha256=pin([folder / 'gpu-history.json', folder / 'network.jsonl', *checkpoints,
                            *folder.glob('*state.npz'), *folder.glob('restart-e*.npz')]))
    write(output, report)
    write(folder / 'gpu-progress.json', dict(state='completed' if report['passed'] else 'failed',
        time=sim.time, elapsed=sim.time - START, cells=len(sim.ids), failure=failure, seconds=report['seconds']))
    del sim, clone
    torch.cuda.empty_cache()
    return report


def compare_pair(root, pair, p):
    root = Path(root)
    output = root / f"time-h{pair['seed']}-{pair['branch']}.json"
    if output.exists():
        raise FileExistsError('Never overwrite a formation comparison')
    a, b = pair['cases']
    histories = [read(root / key / 'gpu-history.json') for key in (a, b)]
    errors, unmatched = f.q.compare_observations(*histories)
    homologous, missing = f.q.homologous_errors(*histories)
    if unmatched or missing or any(x['ids'] != y['ids'] or x['parents'] != y['parents']
                                 for x, y in zip(*histories)):
        raise ValueError('Fixed ancestry/cell correspondence required')
    fields = []
    for elapsed in p['snapshots']:
        densities = [f.q.checkpoint_density(root / key / f'gpu-e{elapsed:02d}.npz', None, 3) for key in (a, b)]
        fields.append(dict(elapsed=elapsed, time=START + elapsed,
            errors=f.q.density_errors(densities[0][0], densities[1][0]), phase_ranges=[x[1] for x in densities]))
    field_max = max(v for item in fields for v in item['errors'].values())
    streams = [[json.loads(x) for x in (root / key / 'network.jsonl').read_text().splitlines()] for key in (a, b)]
    if len(streams[0]) != 4001 or len(streams[1]) != 4001 or len(histories[0]) != 401 or len(histories[1]) != 401:
        raise ValueError('Incomplete formation observations')
    dense, topology = [], []
    for x, y in zip(*streams):
        if abs(x['time'] - y['time']) > 1e-9 or x['ids'] != y['ids'] or x['parents'] != y['parents']:
            raise ValueError('Dense formation correspondence changed')
        sa, sb = [f.q.spectrum(z['volumes'], z['delta']) for z in (x, y)]
        dense.append(dict(time=x['time'], elapsed=x['time'] - START,
            relative_spectrum=float(np.linalg.norm(sa - sb) / (np.linalg.norm(sb) or 1.)), edges=[x['edges'], y['edges']]))
        if not np.array_equal(np.array(x['conductance']) > 0, np.array(y['conductance']) > 0):
            topology.append(x['time'])
    peak = max(dense, key=lambda x: x['relative_spectrum'])
    failed = {k: v for k, v in errors.items() if v > p['time_limits'][k]}
    failed.update({'homologous_' + k: v for k, v in homologous.items() if v > p['homologous_limits'][k]})
    if field_max > p['time_limits']['field_max']:
        failed['field_max'] = field_max
    if peak['relative_spectrum'] > p['dense_time_spectrum_limit']:
        failed['dense_spectrum'] = peak['relative_spectrum']
    grow = 0.
    for x, y in zip(*histories):
        gx, gy = x['frozen_uniform_max_growth'], y['frozen_uniform_max_growth']
        if (gx is None) != (gy is None):
            failed['growth_support'] = True
        if gx is not None and gy is not None:
            grow = max(grow, abs(gx - gy))
    if grow > p['growth_abs_limit']:
        failed['growth_abs_max'] = grow
    outcomes = [formation(h, p) for h in histories]
    outcome_checks = []
    for x, y in zip(*outcomes):
        concordant = x['sustained_formation'] == y['sustained_formation'] and x['late_persistent'] == y['late_persistent']
        onset_error = abs(x['first_sustained_onset_elapsed'] - y['first_sustained_onset_elapsed']) if x['sustained_formation'] and y['sustained_formation'] else None
        passed = concordant and (onset_error is None or onset_error <= p['onset_time_limit'] + 1e-9)
        outcome_checks.append(dict(threshold=x['threshold'], classifications_agree=concordant,
            sustained_onset_time_error=onset_error, passed=passed))
        if not passed:
            failed['formation_threshold_' + str(x['threshold'])] = outcome_checks[-1]
    report = dict(pair, passed=not failed, errors=errors, homologous_errors=homologous,
        fields=fields, field_max=field_max, dense_peak=peak, topology_mismatch_times=topology,
        growth_abs_max=grow, formation=outcomes, formation_checks=outcome_checks,
        failed_gates=failed, protocol_sha256=digest(root / 'protocol.json'),
        evidence_sha256=pin(root / key / 'gpu-result.json' for key in (a, b)))
    write(output, report)
    return report


def run(root):
    root = Path(root).resolve()
    with (root / 'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root / 'protocol.json')
        verify(p)
        ph = digest(root / 'protocol.json')
        if read(root / 'status.json') != dict(state='prepared', protocol_sha256=ph):
            raise ValueError('Unchanged prepared study required')
        if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']:
            raise ValueError('Dedicated GPU selection required')
        torch.set_num_threads(1)
        lib = f.td.gpu.library()
        if digest(lib.binary_path) != p['runtime_binary_sha256']:
            raise ValueError('Accepted runtime binary changed')
        write(root / 'runtime.json', dict(gpu=torch.cuda.get_device_name(0), torch=torch.__version__,
            binary=lib.binary_path, binary_sha256=digest(lib.binary_path)))
        started = time.perf_counter()
        paths, pairs, prefixes = [], [], []
        active = None
        def status(state, **extra):
            write(root / 'status.json', dict(state=state, active=active, completed=len(paths), total=24,
                pairs_completed=len(pairs), pairs_total=12, independent_histories=3, new_scientific_histories=0,
                qualification_accepted=False, scientific_accepted=False, maintenance_extension_started=False,
                seconds=time.perf_counter() - started, protocol_sha256=ph, **extra))
        try:
            for pair in p['time_pairs']:
                for key in pair['cases']:
                    active = key
                    job = next(j for j in p['jobs'] if j['key'] == key)
                    status('running_prefix')
                    prefix = f.pilot.prefix(root, job, p)
                    prefixes.append(prefix)
                    if not prefix['passed']:
                        raise RuntimeError('Formation actual-start prefix failed: ' + key)
                    status('running_formation')
                    result = run_path(root, job, p)
                    paths.append(dict(case=key, passed=result['passed'], seconds=result['seconds']))
                    if not result['passed']:
                        raise RuntimeError('Formation path failed: ' + key)
                status('comparing')
                result = compare_pair(root, pair, p)
                pairs.append(result)
                if not result['passed']:
                    raise RuntimeError('Formation full-horizon gate failed: ' + pair['branch'])
            verify(p)
            write(root / 'summary.json', dict(execution_completed=True, formation_screens_passed=True,
                qualification_accepted=False, scientific_accepted=False, independent_review_pending=True,
                independent_histories=3, new_scientific_histories=0, new_numerical_paths=24,
                maintenance_extension_started=False, paths=paths, pairs=pairs, prefixes=prefixes,
                seconds=time.perf_counter() - started, protocol_sha256=ph, scope=p['scope']))
            status('completed_provisional', formation_screens_passed=True, independent_review_pending=True)
        except Exception as exc:
            status('failed', failure=dict(type=type(exc).__name__, message=str(exc)), paths=paths, pairs=pairs, prefixes=prefixes)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    (prepare if args.action == 'prepare' else run)(args.root)


if __name__ == '__main__':
    main()
