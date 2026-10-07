"""Moving chemical exchange and response from qualified carry maintenance states.

Keep the scientific kernels and earlier evidence unchanged. Both new timesteps
start from the same fine-level mature state, including its phase carry.
"""
import argparse
from dataclasses import asdict
import fcntl
from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .attribute_exchange import exchange
from .attribute_persistence import distance
from .benchmark_gpu_backend import cpu_audit
from .cell_exchange_response_moving import local_response
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .cell_response_moving import observe
from .cell_response_moving_refinement import CRITERIA as RESPONSE_CRITERIA, response_error
from .feedback_long import digest, restore_checkpoint, save_checkpoint
from .gpu_backend import library as resident_library
from .gpu_precision_control import PrecisionSimulation, library
from .moving_initiation_controls import cpu_chemical_step, exact_state
from .native_mechanics import NativeSimulation
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .parameter_robustness_moving import retain
from .phase_carry_maintenance import current_evidence as maintenance_evidence
from .polarity_robustness import first_crossing
from .resolution import _steps, write_json
from .validate_gpu_backend import CRITERIA as GPU_CRITERIA, discrepancies

HISTORIES = (7, 8, 9)
LEVELS = (('coarse', .00375), ('fine', .001875))
BACKGROUNDS = ('unexchanged', 'fresh_exchange')
MECHANICAL_LIMITS = dict(chemical_log_max=.01, polarity_abs_max=.01,
    relative_axis_max=.01, relative_volume_max=.005, relative_transport_max=.01,
    final_phi_abs_max=.02)


def payload(path):
    with np.load(path, allow_pickle=False) as z:
        out = {key: z[key].copy() for key in z.files}
    if (str(out.get('precision_arm')) != 'phase_carry' or
            out['phase_carry'].dtype != np.float64 or
            out['phase_carry'].shape != out['phi'].shape or
            not np.isfinite(out['phase_carry']).all() or
            out['rounding'].dtype != np.int32 or out['rounding'].shape != (2,) or
            np.any(out['rounding'] < 0)):
        raise ValueError('Requires a valid carry checkpoint')
    return out


def prepared_checkpoint(source, destination, dt, chemistry=None):
    """Retiming/chemical editing preserves all other fields and carry bytes."""
    original = payload(source)
    out = {key: value.copy() for key, value in original.items() if key != 'long_experiment'}
    meta = json.loads(str(out['metadata'])); old_dt = meta['config']['dt']
    meta['step_number'] = _steps(meta['time'], dt)
    for key in ('steps', 'save_every'):
        meta['config'][key] = _steps(meta['config'][key]*old_dt, dt)
    meta['config']['dt'] = dt
    out['metadata'] = np.array(json.dumps(meta))
    if chemistry is not None:
        chemistry = np.asarray(chemistry, dtype=np.float64)
        if (chemistry.shape != (2, len(out['ids'])) or
                not np.isfinite(chemistry).all() or np.any(chemistry <= 0)):
            raise ValueError('Positive finite two-species chemistry required')
        out['activator'], out['inhibitor'] = chemistry.copy()
    changed = {'metadata', 'activator', 'inhibitor', 'long_experiment'}
    if any(not np.array_equal(original[key], value) for key, value in out.items() if key not in changed):
        raise ValueError('Nonchemical checkpoint field changed')
    destination = Path(destination); destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(destination)
    np.savez_compressed(destination, **out)
    host = AttributeSimulation.restore(destination)
    if (host.divisions or host.attribute_mode != 'direct' or not host.config.feedback or
            len(host.ids) != host.config.max_cells or abs(host.time-host.step_number*dt) > 1e-9):
        raise ValueError('Invalid mature physical start')
    return host


def select_pair(chemistry, ids, context):
    """Maximal chemical contrast among above-median context contrasts."""
    chemistry, ids, context = np.asarray(chemistry), np.asarray(ids), np.asarray(context, float)
    if (chemistry.shape != (2, len(ids)) or context.shape != (len(ids), 3) or
            len(set(ids.tolist())) != len(ids) or not np.isfinite(context).all() or
            not np.isfinite(chemistry).all() or np.any(chemistry <= 0)):
        raise ValueError('Invalid selection inputs')
    standardized = (context-context.mean(axis=0))/np.maximum(context.std(axis=0), 1e-12)
    pairs = list(combinations(range(len(ids)), 2))
    differences = [float(np.linalg.norm(standardized[i]-standardized[j])) for i, j in pairs]
    cutoff = float(np.median(differences))
    if cutoff <= 1e-12:
        raise ValueError('No resolved geometric-context contrast')
    eligible = [(i, j, cd) for (i, j), cd in zip(pairs, differences) if cd >= cutoff]
    i, j, cd = min(eligible, key=lambda item: (
        -abs(np.log(chemistry[0, item[0]]/chemistry[0, item[1]])),
        min(int(ids[item[0]]), int(ids[item[1]])), max(int(ids[item[0]]), int(ids[item[1]]))))
    pair = sorted((i, j), key=lambda index: (chemistry[0, index], int(ids[index])))
    if abs(np.log(chemistry[0, pair[0]]/chemistry[0, pair[1]])) <= .1:
        raise ValueError('Selected pair lacks chemical contrast')
    return dict(pair=pair, selected_ids=[int(ids[k]) for k in pair],
        context=context.tolist(), standardized_context=standardized.tolist(),
        context_features=['exposure', 'conservative_exit_rate', 'measured_volume'],
        context_distance=cd, median_context_distance=cutoff,
        log_activator_separation=float(abs(np.log(chemistry[0, pair[0]]/chemistry[0, pair[1]]))))


def job_design(root, selections):
    root = Path(root)
    jobs = []
    for seed in HISTORIES:
        selected = selections[str(seed)]
        for stage in ('formation', 'response'):
            for background in BACKGROUNDS:
                for level, dt in LEVELS:
                    trials = [(None, 1.)] if stage == 'formation' else [(None, 1.), *[(cell, .9) for cell in selected['selected_ids']]]
                    for target, factor in trials:
                        suffix = 'control' if target is None else f'cell-{target}_neg10'
                        key = f'seed-{seed}_{level}_{background}_{stage}_{suffix}'
                        source = root/f'seed-{seed}'/f'{level}-{background}.npz' if stage == 'formation' else root/key/'source.npz'
                        jobs.append(dict(key=key, seed=seed, stage=stage, background=background,
                            level=level, dt=dt, target=target, factor=factor, source=str(source),
                            start=390. if stage == 'formation' else 450., duration=60.,
                            point=dict(ratio=27.5, chi=.35)))
    return jobs


def configuration(host, job, p):
    if (asdict(host.config) != p['accepted_configs'][f'{job["seed"]}_{job["level"]}'] or host.divisions or
            host.attribute_mode != 'direct' or len(host.ids) != host.config.max_cells or
            abs(host.time-host.step_number*job['dt']) > 1e-9):
        raise ValueError('Changed mature configuration, mode or clock')


def quality_pass(audit, p):
    return bool(audit['max_volume_error'] < .05 and audit['min_radius'] >= 4 and
        audit['max_clipping'] == 0 and audit['boundary_max'] < .01 and
        audit['dilution_error_max'] <= p['amount_error_max'])


def validate_history(history, job, p, horizon):
    times = np.arange(_steps(horizon, p['interval'])+1)*p['interval']
    if len(history) != len(times) or not np.allclose([r['elapsed'] for r in history], times, rtol=0, atol=1e-9):
        raise ValueError('Incomplete or misaligned observation history')
    host = AttributeSimulation.restore(job['source']); configuration(host, job, p)
    expected = np.array([host.activator, host.inhibitor])
    if not np.array_equal(history[0]['chemistry'], expected):
        raise ValueError('Changed initial chemistry, exchange or pulse')
    for row in history:
        values = np.asarray(row['chemistry'])
        if (row['ids'] != host.ids.tolist() or abs(row['time']-job['start']-row['elapsed']) > 1e-9 or
                values.shape != expected.shape or not np.isfinite(values).all() or np.any(values <= 0)):
            raise ValueError('Changed cell IDs, physical clock or invalid chemistry')
        if abs(float(np.std(np.log(values[0])))-row['log_activator_sd']) > 1e-12:
            raise ValueError('Contrast inconsistent with chemical state')
        validate_graph(row['delta'], row['volumes'])
    return history


def contrast_summary(history, p, horizon):
    time_values = np.array([r['elapsed'] for r in history])
    values = np.array([r['log_activator_sd'] for r in history])
    late = time_values >= horizon-min(p['late_window'], horizon)-1e-9
    threshold = p['contrast_min']
    return dict(final_log_sd=float(values[-1]), late_min_log_sd=float(values[late].min()),
        persistent_contrast=bool(values[late].min() > threshold),
        whole_window_min_log_sd=float(values.min()), whole_window_retained=bool(values.min() > threshold),
        first_loss=first_crossing(time_values, values, threshold, 'down'))


def current_evidence(root, job, p):
    folder = Path(root)/job['key']; ph = digest(Path(root)/'protocol.json')
    result, history = read(folder/'result.json'), read(folder/'history.json')
    if (result['job'] != job or result['protocol_sha256'] != ph or not result['quality_pass'] or
            result['history_sha256'] != digest(folder/'history.json') or
            result['checkpoint_sha256'] != digest(folder/'latest_state.npz') or
            result['source_sha256'] != digest(job['source'])):
        raise ValueError('Changed completed scientific evidence')
    host, audit, saved = restore_checkpoint(folder/'latest_state.npz', ph, job)
    configuration(host, job, p); payload(folder/'latest_state.npz')
    validate_history(history, job, p, job['duration'])
    if (saved != history or audit != result['audit'] or not quality_pass(audit, p) or
            abs(host.time-job['start']-job['duration']) > 1e-9 or
            not np.array_equal(history[-1]['chemistry'], [host.activator, host.inhibitor])):
        raise ValueError('Checkpoint/history mismatch or failed quality')
    return result, history


def context_gate(root, job, p):
    """Independent chemistry, exact residual restart and strict native prefix."""
    folder = Path(root)/job['key']/'prefix'; file = folder/'result.json'
    ph = digest(Path(root)/'protocol.json') if (Path(root)/'protocol.json').exists() else None
    source_sha = digest(job['source'])
    if file.exists():
        result = read(file)
        if (not result['passed'] or result['job'] != job or result['source_sha256'] != source_sha or
                result['protocol_sha256'] not in (None, ph)):
            raise ValueError('Changed context gate')
        if result['protocol_sha256'] is None and ph is not None and p['input_sha256'].get(str(file.resolve())) != digest(file):
            raise ValueError('Unpinned prelaunch gate')
        for path, sha in result['evidence_sha256'].items():
            if digest(path) != sha:
                raise ValueError('Changed gate evidence')
        return result
    folder.mkdir(parents=True, exist_ok=True)
    sim = PrecisionSimulation.restore(job['source']); configuration(sim.host, job, p)
    chemical_error = accounting = 0.
    for _ in range(8):
        delta = sim.matrices()[1].cpu().numpy(); old_volume = sim.geometry[:, 0].cpu().numpy().copy()
        before = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
        c = sim.config
        expected = cpu_chemical_step(before, delta, job['dt'], c.signal_beta, c.signal_da, c.signal_dh)
        sim.step(); expected *= old_volume/sim.geometry[:, 0].cpu().numpy()
        actual = np.array([sim.activator.cpu().numpy(), sim.inhibitor.cpu().numpy()])
        chemical_error = max(chemical_error, float(abs(np.log(actual/expected)).max()))
        accounting = max(accounting, sim.audit()['dilution_amount_error'])
    if chemical_error > 1e-11 or accounting > p['amount_error_max']:
        write_json(folder/'chemical-check-failure.json', dict(job=job, source_sha256=source_sha,
            cpu_chemical_log_max=chemical_error, amount_accounting_max=accounting,
            chemical_limit=1e-11, amount_limit=p['amount_error_max']))
        raise RuntimeError(f'Independent chemistry/accounting gate failed: log={chemical_error:g}, amount={accounting:g}')
    checkpoint = folder/'restart.npz'; sim.checkpoint(checkpoint)
    restart = PrecisionSimulation.restore(checkpoint); exact_state(sim, restart)
    for _ in range(4):
        sim.step(); restart.step(); exact_state(sim, restart)
    del sim, restart
    cpu = AttributeSimulation.restore(job['source']); configuration(cpu, job, p)
    cpu.__class__ = NativeSimulation; cpu.native_threads = p['prefix_native_threads']
    gpu = PrecisionSimulation.restore(job['source'])
    histories = [[], []]
    errors = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                           'relative_volume_max', 'relative_transport_max'), 0.)
    steps, every = _steps(p['prefix_duration'], job['dt']), _steps(p['interval'], job['dt'])
    for step in range(steps+1):
        cpu_audit(cpu); q = gpu.audit()
        if q['dilution_amount_error'] > p['amount_error_max']:
            raise RuntimeError('Prefix dilution failed')
        if step % every == 0:
            a, b = gpu.observe(step*job['dt']), observe(cpu, step*job['dt'])
            if max(a['boundary_occupancy'], b['boundary_occupancy']) >= .01:
                raise RuntimeError('Prefix boundary failed')
            for key, value in discrepancies(a, b).items():
                errors[key] = max(errors[key], value)
            histories[0].append(a); histories[1].append(b)
        if step < steps:
            cpu.step(); gpu.step()
    phi_error = float(abs(gpu.phi.cpu().numpy()-cpu.phi).max())
    passed = all(value <= GPU_CRITERIA[key] for key, value in errors.items()) and phi_error <= GPU_CRITERIA['final_phi_abs_max']
    files = [folder/'gpu-history.json', folder/'cpu-history.json', folder/'gpu-end.npz', folder/'cpu-end.npz', checkpoint]
    write_json(files[0], histories[0]); write_json(files[1], histories[1]); gpu.checkpoint(files[2]); cpu.checkpoint(files[3])
    del gpu, cpu; torch.cuda.empty_cache()
    result = dict(passed=bool(passed), job=job, protocol_sha256=ph, source_sha256=source_sha,
        errors=errors, final_phi_abs_max=phi_error, independent_chemical_steps=8,
        cpu_chemical_log_max=chemical_error, amount_accounting_max=accounting, restart_exact_steps=4,
        evidence_sha256={str(f.resolve()): digest(f) for f in files},
        scope='Exact inherited-carry context, short native/carry prefix; not full-horizon backend equivalence.')
    write_json(file, result)
    if not passed:
        raise RuntimeError('Strict native/carry context gate failed: '+job['key'])
    print('Context gate passed:', job['key'], flush=True)
    return result


def materialize_response(root, job, p):
    folder = Path(root)/job['key']; file = folder/'handoff.json'
    formation = next(j for j in p['jobs'] if j['seed'] == job['seed'] and j['level'] == job['level'] and
        j['background'] == job['background'] and j['stage'] == 'formation')
    prior = Path(root)/formation['key']
    if file.exists():
        record = read(file)
        if record['job'] != job or record['protocol_sha256'] != digest(Path(root)/'protocol.json'):
            raise ValueError('Changed response handoff')
        for path, sha in record['evidence_sha256'].items():
            if digest(path) != sha:
                raise ValueError('Changed response source/handoff evidence')
        return record
    _, history = current_evidence(root, formation, p)
    values = np.array(history[-1]['chemistry']); ids = history[-1]['ids']
    injected = 0.
    if job['target'] is not None:
        index = ids.index(job['target'])
        injected = float((job['factor']-1)*values[0, index]*history[-1]['volumes'][index])
        values[0, index] *= job['factor']
    prepared_checkpoint(prior/'latest_state.npz', job['source'], job['dt'], values)
    files = [prior/f for f in ('result.json', 'history.json', 'latest_state.npz')] + [Path(job['source'])]
    record = dict(job=job, protocol_sha256=digest(Path(root)/'protocol.json'),
        formation_job=formation['key'], injected_activator_amount=injected,
        evidence_sha256={str(f.resolve()): digest(f) for f in files},
        preservation='Own completed moving context and phase carry retained; pulse applied once before initial checkpoint.')
    write_json(file, record)
    return record


def advance(root, job, p, horizon):
    if horizon not in (p['pilot_duration'], job['duration']):
        raise ValueError('Undeclared horizon')
    folder = Path(root)/job['key']; folder.mkdir(parents=True, exist_ok=True)
    if job['stage'] == 'response':
        materialize_response(root, job, p)
    context_gate(root, job, p)
    ph = digest(Path(root)/'protocol.json'); cp = folder/'latest_state.npz'
    if (folder/'result.json').exists():
        _, history = current_evidence(root, job, p)
        return [r for r in history if r['elapsed'] <= horizon+1e-9]
    if cp.exists():
        host, audit, history = restore_checkpoint(cp, ph, job)
        configuration(host, job, p); validate_history(history, job, p, host.time-job['start'])
        if not np.array_equal(history[-1]['chemistry'], [host.activator, host.inhibitor]):
            raise ValueError('Restart chemistry/history mismatch')
        sim = PrecisionSimulation.restore(cp)
    else:
        sim = PrecisionSimulation.restore(job['source']); history = []
        audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0.,
                     boundary_max=0., dilution_error_max=0., wall_seconds=0.)
    start, stop = _steps(job['start'], job['dt']), _steps(job['start']+horizon, job['dt'])
    if sim.step_number > stop:
        del sim; torch.cuda.empty_cache()
        return [r for r in history if r['elapsed'] <= horizon+1e-9]
    if not start <= sim.step_number <= stop:
        raise ValueError('Invalid continuation clock')
    every, save_every = _steps(p['interval'], job['dt']), _steps(p['checkpoint_interval'], job['dt'])
    clock, base_wall = time.perf_counter(), audit['wall_seconds']
    while sim.step_number <= stop:
        q = sim.audit()
        for key in ('max_volume_error', 'max_clipping'):
            audit[key] = max(audit[key], q[key])
        audit['min_radius'] = min(audit['min_radius'], q['min_radius'])
        audit['dilution_error_max'] = max(audit['dilution_error_max'], q['dilution_amount_error'])
        if audit['dilution_error_max'] > p['amount_error_max']:
            raise RuntimeError('Dilution accounting failed')
        offset = sim.step_number-start
        if offset % every == 0:
            elapsed = offset*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed) > 1e-10:
                row = retain(sim.observe(elapsed), job); row['precision'] = sim.precision_diagnostics(); history.append(row)
            audit['boundary_max'] = max(audit['boundary_max'], history[-1]['boundary_occupancy'])
            if audit['boundary_max'] >= .01:
                raise RuntimeError('Boundary screen failed')
            audit['wall_seconds'] = base_wall+time.perf_counter()-clock
            write_json(folder/'status.json', dict(state='running', elapsed=elapsed,
                target_horizon=horizon, audit=audit, log_activator_sd=history[-1]['log_activator_sd']))
            if offset % save_every == 0 or sim.step_number == stop:
                save_checkpoint(sim, cp, audit, history, ph, job); write_json(folder/'history.json', history)
                print(f'{job["key"]}: elapsed={elapsed:g}, target={horizon:g}, wall={audit["wall_seconds"]:.1f}s', flush=True)
        if sim.step_number == stop:
            break
        sim.step()
    validate_history(history, job, p, horizon)
    if horizon == job['duration']:
        if not quality_pass(audit, p):
            raise RuntimeError('Final quality gate failed')
        write_json(folder/'result.json', dict(job=job, protocol_sha256=ph, source_sha256=digest(job['source']),
            history_sha256=digest(folder/'history.json'), checkpoint_sha256=digest(cp),
            quality_pass=True, audit=audit, **contrast_summary(history, p, horizon)))
        write_json(folder/'status.json', dict(state='completed', elapsed=horizon))
    else:
        write_json(folder/'status.json', dict(state='awaiting_pair_refinement', elapsed=horizon, audit=audit))
    del sim; torch.cuda.empty_cache()
    return history


def stage_histories(root, p, seed, stage, level, horizon):
    out = {}
    for job in p['jobs']:
        if (job['seed'], job['stage'], job['level']) != (seed, stage, level):
            continue
        history = [r for r in read(Path(root)/job['key']/'history.json') if r['elapsed'] <= horizon+1e-9]
        validate_history(history, job, p, horizon)
        out[(job['background'], job['target'])] = history
    return out


def state_comparison(histories, selection, p, horizon):
    pair, m = selection['pair'], np.asarray(selection['masses'])
    baseline = np.array([r['chemistry'] for r in histories[('unexchanged', None)]])
    transferred = np.array([exchange(row, m, *pair) for row in baseline])
    separation = float(distance(transferred[0][:, pair], baseline[0][:, pair], m[pair]))
    if separation <= .1:
        raise ValueError('Uninformative conservative exchange')
    late = np.array([r['elapsed'] for r in histories[('unexchanged', None)]]) >= horizon-min(p['late_window'], horizon)-1e-9
    rows = {}
    for background in BACKGROUNDS:
        history = histories[(background, None)]; values = np.array([r['chemistry'] for r in history])
        destination = float(distance(values[:, :, pair], baseline[:, :, pair], m[pair])[late].max()/separation)
        donor = float(distance(values[:, :, pair], transferred[:, :, pair], m[pair])[late].max()/separation)
        rows[background] = dict(**contrast_summary(history, p, horizon),
            initial_pair_separation=separation, destination_ratio=destination, transferred_ratio=donor,
            moving_state_outcome='destination_like' if destination < .1 else 'transferred_like' if donor < .1 else 'reorganized',
            global_destination_distance=float(distance(values[-1], baseline[-1], m)))
    return rows


def response_summary(histories, selection, p, horizon):
    ids = histories[('unexchanged', None)][0]['ids']
    times = np.arange(_steps(horizon, p['interval'])+1)*p['interval']
    datasets, metrics = {}, {}
    for background in BACKGROUNDS:
        control = np.array([r['chemistry'] for r in histories[(background, None)]])
        weights = np.array(histories[(background, None)][0]['volumes'])
        for target in selection['selected_ids']:
            cell = ids.index(target); path = np.array([r['chemistry'] for r in histories[(background, target)]])
            datasets[(background, target)] = local_response(path, control, cell, .9)
            metrics[(background, target)] = dict(injected_activator_amount=float(-.1*control[0, 0, cell]*weights[cell]),
                **response_metrics(path, control, times, weights, cell, .9))
    rows = []
    for target in selection['selected_ids']:
        donor = next(cell for cell in selection['selected_ids'] if cell != target)
        wave = datasets[('fresh_exchange', target)]; destination_wave = datasets[('unexchanged', target)]
        donor_wave = datasets[('unexchanged', donor)]
        sep = waveform_distance(destination_wave, donor_wave, times)
        dd, ds = waveform_distance(wave, destination_wave, times), waveform_distance(wave, donor_wave, times)
        rows.append(dict(cell=target, donor=donor, factor=.9, baseline_separation=sep,
            destination_distance=dd, donor_distance=ds, donor_minus_destination_distance=ds-dd,
            informative=bool(sep > .01), nearest_reference=('donor' if ds < dd else 'destination') if sep > .01 else 'unresolved'))
    return dict(metrics=[dict(background=b, cell=t, **value) for (b, t), value in metrics.items()],
        comparisons=rows, scope='Two negative-pulse targets within one mature history; nearest reference is descriptive, not donor equivalence.')


def mechanical_comparison(coarse, fine, coarse_phi, fine_phi, p, horizon):
    if len(coarse) != len(fine):
        raise ValueError('Unaligned timestep observations')
    errors = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                           'relative_volume_max', 'relative_transport_max'), 0.)
    for a, b in zip(coarse, fine):
        for key, value in discrepancies(a, b).items():
            errors[key] = max(errors[key], value)
    errors['final_phi_abs_max'] = float(abs(np.asarray(coarse_phi)-fine_phi).max())
    summaries = [contrast_summary(h, p, horizon) for h in (coarse, fine)]
    loss_times = [r['first_loss'] for r in summaries]
    loss_error = 0. if loss_times == [None, None] else None if None in loss_times else abs(loss_times[0]-loss_times[1])
    same = all(summaries[0][key] == summaries[1][key] for key in ('persistent_contrast', 'whole_window_retained'))
    passed = all(value <= p['mechanical_limits'][key] for key, value in errors.items()) and same and loss_error is not None and loss_error <= .3
    return dict(passed=bool(passed), errors=errors, loss_times=loss_times, loss_time_error=loss_error,
        contrast_summaries=summaries, same_contrast_outcomes=same)


def response_refinement(coarse_histories, fine_histories, selection, p, horizon):
    summaries = [response_summary(h, selection, p, horizon) for h in (coarse_histories, fine_histories)]
    metric = 'target_activator_log_auc_per_log_pulse'; rows = []
    for background in BACKGROUNDS:
        controls = [np.array([r['chemistry'] for r in h[(background, None)]]) for h in (coarse_histories, fine_histories)]
        for target in selection['selected_ids']:
            paths = [np.array([r['chemistry'] for r in h[(background, target)]]) for h in (coarse_histories, fine_histories)]
            error = response_error(paths[0], controls[0], paths[1], controls[1], .9)
            a, b = [next(r for r in summary['metrics'] if r['background'] == background and r['cell'] == target) for summary in summaries]
            auc_error = abs(b[metric]/a[metric]-1)
            recovery = {}
            for name in ('target_recovery_time', 'network_recovery_time'):
                x, y = a[name], b[name]
                recovery[name] = 0. if x is None and y is None else None if x is None or y is None else abs(x-y)
            c = p['response_limits']
            passed = error <= c['max_normalized_response_error'] and auc_error <= c['max_relative_auc_error'] and all(
                value is not None and value <= c['max_recovery_time_error']+1e-9 for value in recovery.values())
            rows.append(dict(background=background, cell=target, normalized_response_error=error,
                relative_auc_error=auc_error, recovery_time_errors=recovery, passed=bool(passed)))
    decisions = [r['nearest_reference'] for r in summaries[0]['comparisons']] == [r['nearest_reference'] for r in summaries[1]['comparisons']]
    return dict(passed=bool(all(r['passed'] for r in rows) and decisions), rows=rows,
        same_nearest_reference_decisions=decisions, coarse=summaries[0], fine=summaries[1])


def pair_report(root, p, jobs, horizon):
    first = jobs[0]; stage = first['stage']; seed = first['seed']
    selected = p['selections'][str(seed)]; histories = {}
    for level in ('coarse', 'fine'):
        histories[level] = stage_histories(root, p, seed, stage, level, horizon)
    folder = Path(root)/'pairs'/f'seed-{seed}_{stage}'; folder.mkdir(parents=True, exist_ok=True)
    label = 'pilot' if horizon == p['pilot_duration'] else 'full'; file = folder/f'{label}-refinement.json'
    old = read(file) if file.exists() else None
    if old:
        if old['protocol_sha256'] != digest(Path(root)/'protocol.json'):
            raise ValueError('Changed paired protocol')
        for path, sha in old['evidence_sha256'].items():
            if digest(path) != sha:
                raise ValueError('Changed paired history/fields')
    files, comparisons = [], []
    for background, target in sorted(histories['coarse'], key=lambda key: (key[0], -1 if key[1] is None else key[1])):
        arrays = {}; matched = []
        for level in ('coarse', 'fine'):
            job = next(j for j in jobs if (j['level'], j['background'], j['target']) == (level, background, target))
            stem = f'{label}_{background}_'+('control' if target is None else f'cell-{target}')+f'_{level}'
            hf, pf = folder/f'{stem}.json', folder/f'{stem}.npz'
            h = histories[level][(background, target)]
            if old:
                if read(hf) != h:
                    raise ValueError('Changed compared history')
            else:
                z = payload(Path(root)/job['key']/'latest_state.npz'); meta = json.loads(str(z['metadata']))
                if abs(meta['time']-job['start']-horizon) > 1e-9:
                    raise ValueError('Field snapshot does not match paired horizon')
                write_json(hf, h); np.savez_compressed(pf, phi=z['phi'])
            with np.load(pf) as z:
                arrays[level] = z['phi'].copy()
            matched.append(h); files.extend((hf, pf))
        comparisons.append(dict(background=background, target=target,
            **mechanical_comparison(*matched, arrays['coarse'], arrays['fine'], p, horizon)))
    result = dict(seed=seed, stage=stage, horizon=horizon,
        passed=bool(all(r['passed'] for r in comparisons)), mechanical=comparisons,
        protocol_sha256=digest(Path(root)/'protocol.json'),
        evidence_sha256={str(f.resolve()): digest(f) for f in files})
    if stage == 'formation':
        outcomes = {level: state_comparison(histories[level], selected, p, horizon) for level in ('coarse', 'fine')}
        a, b = [outcomes[level]['fresh_exchange'] for level in ('coarse', 'fine')]
        error = max(abs(a[key]-b[key]) for key in ('destination_ratio', 'transferred_ratio'))
        same = a['moving_state_outcome'] == b['moving_state_outcome']
        result.update(state=outcomes, pair_ratio_abs_error=error, same_state_outcome=same)
        result['passed'] = bool(result['passed'] and error <= .05 and same)
    else:
        response = response_refinement(histories['coarse'], histories['fine'], selected, p, horizon)
        result['response'] = response; result['passed'] = bool(result['passed'] and response['passed'])
    if old:
        if result != old:
            raise ValueError('Paired decision no longer reproduces')
        return result
    write_json(file, result)
    return result


def require_device(device):
    if (torch.cuda.get_device_name(0) != device['name'] or torch.__version__ != device['torch'] or
            torch.version.cuda != device['torch_cuda'] or digest(library()._name) != device['experimental_binary_sha256'] or
            digest(resident_library()._name) != device['accepted_binary_sha256']):
        raise ValueError('Changed qualified GPU/runtime/kernels')


def prepare(root, parent=Path('outputs/phase-carry-maintenance')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists():
        raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    review_file = Path('docs/phase_carry_maintenance_assessment.json').resolve(); review = read(review_file)
    if (read(parent/'status.json')['state'] != 'completed' or not read(parent/'summary.json')['passed'] or
            not review['accepted'] or review['protocol_sha256'] != digest(parent/'protocol.json') or
            review['original_summary_sha256'] != digest(parent/'summary.json') or
            review['full_review_sha256'] != digest(review['full_review'])):
        raise ValueError('Requires completed qualified and reviewed carry maintenance')
    require_device(old['device']); torch.set_num_threads(1)
    root.mkdir(parents=True)
    inputs = [review_file, Path(review['full_review']), Path('docs/phase_carry_exchange_response.md').resolve()]
    inputs += [parent/f for f in ('protocol.json', 'status.json', 'summary.json', 'implementation-verification.json', 'float64-reference/result.json')]
    selections, accepted = {}, {}
    estimate = 0.
    try:
        for seed in HISTORIES:
            source_job = next(j for j in old['jobs'] if j['seed'] == seed and j['family'] == 'pattern' and j['level'] == 'fine' and j['point']['chi'] == .35)
            result, history = maintenance_evidence(parent, old, source_job)
            folder = parent/source_job['key']; source = folder/'latest_state.npz'; host = AttributeSimulation.restore(source)
            if not result['whole_window_retained'] or abs(host.time-390.) > 1e-9:
                raise ValueError('Unqualified developed starting state')
            ep = read(folder/'endpoint/protocol.json'); verify(ep)
            endpoint = read(folder/'endpoint/assay/result.json')
            if (not endpoint['local_bistability_supported'] or not endpoint['numerical_pass'] or
                    endpoint['protocol_sha256'] != digest(folder/'endpoint/protocol.json') or
                    endpoint['paths_sha256'] != digest(folder/'endpoint/assay/paths.npz')):
                raise ValueError('Changed qualified chemical endpoint')
            initial = np.array(history[-1]['chemistry']); masses = np.array(history[-1]['volumes'])
            _, exposure = host.contacts(); delta = np.array(history[-1]['delta'])
            selection = select_pair(initial, host.ids, np.column_stack([exposure, -np.diag(delta), masses]))
            pair = selection['pair']; swapped = exchange(initial, masses, *pair)
            amount_error = float(abs((swapped@masses)/(initial@masses)-1).max())
            separation = float(distance(swapped[:, pair], initial[:, pair], masses[pair]))
            if amount_error > 2e-14 or separation <= .1:
                raise ValueError('Exchange amount/informativeness gate failed')
            sham = initial.copy(); sham[:, pair] = sham[:, pair].copy()
            if not np.array_equal(sham, initial):
                raise ValueError('Sham assignment changed chemistry')
            selections[str(seed)] = dict(**selection, source_job=source_job['key'], source_sha256=digest(source),
                ids=host.ids.tolist(), masses=masses.tolist(), initial_chemistry=initial.tolist(),
                exchanged_chemistry=swapped.tolist(), initial_exchange_separation=separation,
                amount_error_max=amount_error, sham_exact=True,
                pair_context=[selection['context'][i] for i in pair],
                initial_centers=history[-1]['centers'], initial_polarity=history[-1]['polarity'])
            for level, dt in LEVELS:
                for background, values in (('unexchanged', initial), ('fresh_exchange', swapped)):
                    copied = root/f'seed-{seed}'/f'{level}-{background}.npz'
                    prepared = prepared_checkpoint(source, copied, dt, values)
                    config = asdict(prepared.config)
                    if level in accepted and accepted[level] != config:
                        # Only the inherited history seed can differ across histories.
                        if {k:v for k,v in accepted[level].items() if k != 'seed'} != {k:v for k,v in config.items() if k != 'seed'}:
                            raise ValueError('Unexpected between-history physical parameters')
                    accepted.setdefault(level, config); inputs.append(copied)
            inputs += [folder/f for f in ('result.json', 'history.json', 'latest_state.npz', 'endpoint/protocol.json', 'endpoint/assay/result.json', 'endpoint/assay/paths.npz')]
            estimate += sum(r['audit']['wall_seconds']/240.*480. for r in read(review['full_review'])['paths'] if not r['reused'] and r['job']['seed'] == seed and r['job']['point']['chi'] == .35)
        jobs = job_design(root, selections)
        # Keep a separate configuration for each history: seed is part of checkpoint provenance.
        accepted = {f'{seed}_{level}':asdict(AttributeSimulation.restore(root/f'seed-{seed}'/f'{level}-unexchanged.npz').config) for seed in HISTORIES for level, _ in LEVELS}
        p = dict(jobs=jobs, selections=selections, accepted_configs=accepted,
            device=old['device'], parent=str(parent), histories=list(HISTORIES),
            independent_histories=3, new_histories=0, total_jobs=48, formation_jobs=12, response_jobs=36,
            start=390., formation_end=450., response_end=510., duration=60., interval=.15,
            pilot_duration=6., late_window=24., checkpoint_interval=3., prefix_duration=.6,
            prefix_native_threads=4, amount_error_max=2e-14, contrast_min=.1,
            mechanical_limits=MECHANICAL_LIMITS, response_limits=RESPONSE_CRITERIA,
            estimated_moving_seconds=estimate,
            design='Three qualified carry-maintenance fine positive-chi endpoints, retimed without resetting residuals. Above-median geometric-context pair with largest log-activator contrast. Fresh amount-preserving two-species exchange versus untouched. Sixty moving units then controls and -10% activator pulses in both targets for sixty units, at two timesteps.',
            scope='Three existing mature histories, 48 nested paths. Negative-pulse/fresh-exchange assay only. No autonomous/committed biological identity, fresh carry zygote development, new shape axis, spatial convergence, positive-pulse or broader ledger promotion.')
        # Context lookup is explicit and permits only the original history seed.
        gates = []
        for job in (j for j in jobs if j['stage'] == 'formation'):
            write_json(root/'status.json', dict(state='preparing', current_job=job['key'], contexts_completed=len(gates), total_contexts=12))
            gate = context_gate(root, job, p); gates.append(gate)
            inputs += [root/job['key']/'prefix/result.json', *[Path(f) for f in gate['evidence_sha256']]]
        write_json(root/'selection.json', selections); inputs.append(root/'selection.json')
        write_json(root/'implementation-verification.json', dict(passed=True, actual_formation_contexts=gates,
            inherited_carry_preserved=True, sham_exact=True, amount_error_max=max(s['amount_error_max'] for s in selections.values())))
        inputs.append(root/'implementation-verification.json')
        files = [Path(__file__), Path('tests/test_phase_carry_exchange_response.py')]
        p['source_sha256'] = {**old['source_sha256'], **{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256'] = {**old['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}}
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', completed=0, total=48, actual_formation_contexts=12,
            estimated_moving_seconds=estimate, independent_histories=3, new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error)))
        raise


def stage_jobs(p, seed, stage):
    return [j for j in p['jobs'] if j['seed'] == seed and j['stage'] == stage]


def schedule_stage(root, p, seed, stage, done):
    jobs = stage_jobs(p, seed, stage)
    for horizon in (p['pilot_duration'], p['duration']):
        for job in jobs:
            write_json(Path(root)/'status.json', dict(state='running', stage=stage,
                current_job=job['key'], target_horizon=horizon, completed=done, total=48))
            advance(root, job, p, horizon)
            if horizon == p['duration']:
                done += 1
        report = pair_report(root, p, jobs, horizon)
        if not report['passed']:
            return done, False
    return done, True


def assess(root, p=None):
    root = Path(root); p = read(root/'protocol.json') if p is None else p; verify(p)
    rows, histories, missing = [], {}, []
    for job in p['jobs']:
        if not (root/job['key']/'result.json').exists():
            missing.append(job['key']); continue
        result, history = current_evidence(root, job, p)
        rows.append(result); histories[job['key']] = history
    by_history = []
    for seed in p['histories']:
        reports = {}
        for stage in ('formation', 'response'):
            jobs = stage_jobs(p, seed, stage)
            if all(j['key'] in histories for j in jobs):
                reports[stage] = {label:pair_report(root, p, jobs, horizon) for label,horizon in (('pilot', p['pilot_duration']), ('full', p['duration']))}
        qualified = len(reports) == 2 and all(r['passed'] for stage in reports.values() for r in stage.values())
        row = dict(seed=seed, qualified=bool(qualified), selected_ids=p['selections'][str(seed)]['selected_ids'], stage_reports=reports)
        if 'response' in reports:
            comparisons = reports['response']['full']['response']['fine']['comparisons']
            row.update(donor_nearer=sum(r['nearest_reference'] == 'donor' for r in comparisons),
                destination_nearer=sum(r['nearest_reference'] == 'destination' for r in comparisons),
                unresolved=sum(r['nearest_reference'] == 'unresolved' for r in comparisons))
        by_history.append(row)
    result = dict(protocol_sha256=digest(root/'protocol.json'), completed=len(rows), total=48,
        passed=bool(not missing and all(r['qualified'] for r in by_history)),
        independent_histories=3, new_histories=0, results=rows, missing=missing,
        by_history=by_history, scope=p['scope'])
    write_json(root/'summary.json', result)
    return result


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); require_device(p['device']); torch.set_num_threads(1)
        done = 0; failures = []
        try:
            for seed in p['histories']:
                for stage in ('formation', 'response'):
                    done, passed = schedule_stage(root, p, seed, stage, done)
                    if not passed:
                        failures.append(dict(seed=seed, stage=stage, reason='prespecified_timestep_gate_failed'))
                        break
            result = assess(root, p)
            write_json(root/'status.json', dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'], completed=done, total=48, failed_stages=failures))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=48, error=str(error)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-exchange-response'))
    args = parser.parse_args()
    result = prepare(args.output) if args.action == 'prepare' else run(args.output) if args.action == 'run' else assess(args.output)
    print({key:value for key,value in result.items() if key in ('passed', 'completed', 'total', 'estimated_moving_seconds')})
