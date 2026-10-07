"""Delayed moving responses from independently qualified t=510 contexts.

Use the existing carry mechanics/chemistry without changing scientific kernels.
Fractional and matched-amount pulses each have their own continuing background.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path

import numpy as np
import torch

from . import phase_carry_exchange_response as moving
from . import phase_carry_network_context as context
from .cell_exchange_response_moving import local_response
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .cell_response_moving_refinement import response_error
from .feedback_long import digest, restore_checkpoint
from .neighbor_context import read
from .parameter_robustness import verify
from .resolution import write_json

HISTORIES = (7, 8, 9)
LEVELS = moving.LEVELS
DOSES = ('fractional', 'amount')


def matched_removal(first, second):
    """Ten percent of the smaller measured activator amount, always feasible."""
    amounts = np.asarray([first, second], dtype=np.float64)
    if not np.isfinite(amounts).all() or np.any(amounts <= 0):
        raise ValueError('Positive finite target amounts required')
    return float(.1*amounts.min())


def pulse(initial, masses, ids, target, dose, common_amount=None):
    values = np.asarray(initial, dtype=np.float64)
    masses, ids = np.asarray(masses), np.asarray(ids)
    if (values.shape != (2, len(masses)) or masses.ndim != 1 or ids.shape != masses.shape or
            len(set(ids.tolist())) != len(ids) or target not in ids or
            not np.isfinite(values).all() or np.any(values <= 0) or
            not np.isfinite(masses).all() or np.any(masses <= 0) or dose not in DOSES):
        raise ValueError('Invalid pulse state, mass, ID or dose')
    i = ids.tolist().index(target)
    before = float(values[0, i]*masses[i])
    removal = .1*before if dose == 'fractional' else common_amount
    if removal is None or not np.isfinite(removal) or not 0 < removal < before:
        raise ValueError('Removal must leave strictly positive activator')
    factor = .9 if dose == 'fractional' else float(1-removal/before)
    if not 0 < factor < 1:
        raise ValueError('Pulse is not representable as a nonzero depletion')
    changed = values.copy(); changed[0, i] *= factor
    actual = float((values[0, i]-changed[0, i])*masses[i])
    if not np.isclose(actual, removal, rtol=2e-13, atol=0):
        raise ValueError('Measured pulse amount disagrees with requested dose')
    return changed, dict(factor=factor, removed_activator_amount=actual,
        requested_removed_amount=float(removal), initial_target_activator_amount=before,
        initial_log_amplitude=float(abs(np.log(factor))))


def job_design(root, selections):
    root = Path(root); contexts, jobs = [], []
    for seed in HISTORIES:
        targets = selections[str(seed)]['selected_ids']
        for recipient in (None, *targets):
            name = 'sham' if recipient is None else f'reset-cell-{recipient}'
            contexts.append(dict(seed=seed, name=name, recipient=recipient,
                targets=targets if recipient is None else [recipient]))
            for level, dt in LEVELS:
                for target, dose in [(None, None), *[(t, d) for t in (targets if recipient is None else [recipient]) for d in DOSES]]:
                    suffix = 'control' if target is None else f'cell-{target}_{dose}'
                    key = f'seed-{seed}_{level}_{name}_{suffix}'
                    jobs.append(dict(key=key, seed=seed, level=level, dt=dt,
                        stage='delayed_response', background=name, recipient=recipient,
                        target=target, dose=dose, source=str(root/key/'source.npz'),
                        start=510., duration=60., point=dict(ratio=27.5, chi=.35)))
    return contexts, jobs


def parent_endpoint(parent, old, job):
    if job['key'] in old['reused']:
        item = old['reused'][job['key']]
        return Path(old['parent'])/item['job']['key']/'latest_state.npz'
    return Path(parent)/job['key']/'latest_state.npz'


def prepare(root, parent=Path('outputs/phase-carry-network-context')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists():
        raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    review_file = Path('docs/phase_carry_network_context_assessment.json').resolve()
    review = read(review_file)
    if (read(parent/'status.json')['state'] != 'completed' or not read(parent/'summary.json')['passed'] or
            not review['full_study_accepted'] or review['protocol_sha256'] != digest(parent/'protocol.json') or
            review['original_summary_sha256'] != digest(parent/'summary.json') or
            review['original_status_sha256'] != digest(parent/'status.json') or
            review['full_review_sha256'] != digest(review['full_review'])):
        raise ValueError('Requires independently reviewed completed moving contexts')
    root.mkdir(parents=True)
    contexts, jobs = job_design(root, old['selections'])
    inputs = [parent/f for f in ('protocol.json', 'status.json', 'summary.json')]
    inputs += [review_file, Path(review['full_review']), Path('docs/phase_carry_delayed_response.md').resolve()]
    accepted, preparations, common_doses = {}, {}, {}
    estimate = 0.
    p = dict(parent=str(parent), device=old['device'], selections=old['selections'],
        contexts=contexts, jobs=jobs, accepted_configs=accepted, histories=list(HISTORIES),
        independent_histories=3, new_histories=0, start=510., duration=60., pilot_duration=6.,
        interval=.15, late_window=24., checkpoint_interval=3., prefix_duration=.6,
        prefix_native_threads=4, amount_error_max=2e-14, contrast_min=.1, effect_min=.01,
        mechanical_limits=moving.MECHANICAL_LIMITS, response_limits=moving.RESPONSE_CRITERIA,
        total_jobs=len(jobs), new_jobs=len(jobs), reused_jobs=0,
        design='Delayed pulses from each actual fine unperturbed t=510 sham/reset endpoint. Same fine physical state per background at both dt levels. Shared own-background moving controls, -10% fractional and common feasible removed-amount pulses. Live carry mechanics, polarity, conservative transport, reactions and dilution.',
        scope='Three existing mature histories, six nested recipient comparisons; no new developmental histories. Backgrounds may differ in chemistry and physical context. Tested delayed moving response, not stationary equilibria, autonomous/inherited biological identity, neighbor-only context, spatial-versus-bulk-dose causation or general backend acceptance.')
    try:
        for seed in HISTORIES:
            states = {}
            for ctx in [c for c in contexts if c['seed'] == seed]:
                previous = next(j for j in old['jobs'] if j['seed'] == seed and j['level'] == 'fine'
                    and j['recipient'] == ctx['recipient'] and j['target'] is None)
                result, history = context.current_evidence(parent, previous, old)
                source = parent_endpoint(parent, old, previous)
                values, masses, ids = np.array(history[-1]['chemistry']), np.array(history[-1]['volumes']), history[-1]['ids']
                state = moving.payload(source); meta = json.loads(str(state['metadata']))
                if (abs(meta['time']-510.) > 1e-9 or not np.array_equal(state['ids'], ids) or
                        not np.array_equal([state['activator'], state['inhibitor']], values)):
                    raise ValueError('Parent endpoint clock/chemistry mismatch')
                ctx.update(parent_job=previous, source=str(source))
                states[ctx['name']] = (source, values, masses, ids, result['audit']['wall_seconds'])
                evidence_folder = source.parent
                inputs.extend(evidence_folder/f for f in ('latest_state.npz', 'history.json', 'result.json', 'prefix/result.json'))
                inputs.extend(Path(f) for f in read(evidence_folder/'prefix/result.json')['evidence_sha256'])
            for target in old['selections'][str(seed)]['selected_ids']:
                amounts = []
                for name in ('sham', f'reset-cell-{target}'):
                    _, values, masses, ids, _ = states[name]; i = ids.index(target)
                    amounts.append(float(values[0, i]*masses[i]))
                common_doses[f'{seed}_{target}'] = dict(background_target_amounts=amounts,
                    removed_amount=matched_removal(*amounts))
            for ctx in [c for c in contexts if c['seed'] == seed]:
                source, values, masses, ids, wall = states[ctx['name']]
                preparations[f'{seed}_{ctx["name"]}'] = dict(source=str(source), ids=ids,
                    chemistry=values.tolist(), measured_gpu_masses=masses.tolist())
                for j in [j for j in jobs if j['seed'] == seed and j['background'] == ctx['name']]:
                    changed = values.copy()
                    record = dict(factor=1., removed_activator_amount=0., requested_removed_amount=0.,
                        initial_target_activator_amount=None, initial_log_amplitude=0.)
                    if j['target'] is not None:
                        changed, record = pulse(values, masses, ids, j['target'], j['dose'],
                            common_doses[f'{seed}_{j["target"]}']['removed_amount'])
                    j.update(record)
                    host = moving.prepared_checkpoint(source, j['source'], j['dt'], changed)
                    key = f'{seed}_{j["level"]}'
                    if key in accepted and accepted[key] != asdict(host.config):
                        raise ValueError('Background configurations differ')
                    accepted[key] = asdict(host.config)
                    inputs.append(Path(j['source']))
                    estimate += wall*(.001875/j['dt'])
            print('Prepared delayed backgrounds and feasible doses:', seed, flush=True)
        p.update(preparations=preparations, common_doses=common_doses, estimated_moving_seconds=estimate)
        files = [Path(__file__), Path('tests/test_phase_carry_delayed_response.py')]
        p['source_sha256'] = {**old['source_sha256'], **{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256'] = {**old['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}}
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs),
            new_jobs=len(jobs), reused_jobs=0, estimated_moving_seconds=estimate,
            independent_histories=3, new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error))); raise


def archive_exact_zero_restart(root, job, p):
    """Handle the old positive-duration validator without editing pinned sources.

    Only an exactly source-equivalent committed time-zero state can be archived
    and restarted. Nonzero or changed states are never silently discarded.
    """
    folder = Path(root)/job['key']; cp = folder/'latest_state.npz'
    if not cp.exists() or (folder/'result.json').exists():
        return False
    host, audit, history = restore_checkpoint(cp, digest(Path(root)/'protocol.json'), job)
    if abs(host.time-job['start']) > 1e-9:
        return False
    moving.configuration(host, job, p)
    initial, saved = moving.payload(job['source']), moving.payload(cp)
    keys = set(initial)-{'long_experiment'}
    if keys != set(saved)-{'long_experiment'}:
        raise ValueError('Zero-time checkpoint schema differs')
    for key in keys:
        equal = (json.loads(str(initial[key])) == json.loads(str(saved[key]))) if key == 'metadata' else (
            initial[key].dtype == saved[key].dtype and np.array_equal(initial[key], saved[key]))
        if not equal:
            raise ValueError('Zero-time checkpoint differs from pinned source: '+key)
    if (len(history) != 1 or abs(history[0]['elapsed']) > 1e-9 or
            abs(history[0]['time']-job['start']) > 1e-9 or history[0]['ids'] != host.ids.tolist() or
            not np.array_equal(history[0]['chemistry'], [host.activator, host.inhibitor]) or
            not moving.quality_pass(audit, p)):
        raise ValueError('Invalid zero-time history or physical audit')
    if (folder/'history.json').exists() and read(folder/'history.json') != history:
        raise ValueError('Zero-time checkpoint/history differs')
    archive = folder/'restart-archives'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive.mkdir(parents=True)
    files = [f for f in (cp, folder/'history.json', folder/'status.json',
        folder/'latest_state.raw.npz', folder/'latest_state.tmp.npz') if f.exists()]
    manifest = dict(reason='Exact pinned-source time-zero restart; positive-duration validator bypassed by reinitialization',
        protocol_sha256=digest(Path(root)/'protocol.json'), source_sha256=digest(job['source']),
        files={f.name:digest(f) for f in files})
    for f in files:
        f.rename(archive/f.name)
    write_json(archive/'manifest.json', manifest)
    return True


def advance(root, job, p, horizon):
    archive_exact_zero_restart(root, job, p)
    return moving.advance(root, job, p, horizon)


def response_summary(histories, jobs, targets, p, horizon):
    times = np.arange(round(horizon/p['interval'])+1)*p['interval']
    rows = []
    for target in targets:
        contexts, waves, raw_waves = {}, {}, {}
        for name in ('sham', f'reset-cell-{target}'):
            h = histories[name, None, None]; ids = h[0]['ids']; i = ids.index(target)
            control = np.array([r['chemistry'] for r in h]); masses = np.array(h[0]['volumes'])
            doses = {}
            for dose in DOSES:
                job = next(j for j in jobs if j['background'] == name and j['target'] == target and j['dose'] == dose)
                path = np.array([r['chemistry'] for r in histories[name, target, dose]])
                if not np.array_equal(path[0, 0, i], control[0, 0, i]*job['factor']):
                    raise ValueError('Pulse factor disagrees with actual start')
                wave = local_response(path, control, i, job['factor'])
                waves[name, dose] = wave
                raw_waves[name, dose] = wave*abs(np.log(job['factor']))
                doses[dose] = dict(factor=job['factor'], removed_activator_amount=job['removed_activator_amount'],
                    initial_log_amplitude=job['initial_log_amplitude'],
                    metrics=response_metrics(path, control, times, masses, i, job['factor']))
            contexts[name] = dict(**moving.contrast_summary(h, p, horizon), doses=doses,
                dose_dependence_normalized_rms=waveform_distance(waves[name, 'fractional'], waves[name, 'amount'], times))
        differences = {}
        for dose in DOSES:
            shift = waveform_distance(waves['sham', dose], waves[f'reset-cell-{target}', dose], times)
            differences[dose] = dict(normalized_response_shift_rms=shift, response_shift_detected=bool(shift>p['effect_min']),
                raw_log_response_shift_rms=waveform_distance(raw_waves['sham', dose], raw_waves[f'reset-cell-{target}', dose], times))
        rows.append(dict(cell=target, contexts=contexts, differences=differences,
            shifted_at_both_doses=bool(all(r['response_shift_detected'] for r in differences.values()))))
    return rows


def response_comparison(coarse, fine, jobs, targets, p, horizon):
    summaries = [response_summary(h, [j for j in jobs if j['level'] == level], targets, p, horizon)
        for h, level in ((coarse, 'coarse'), (fine, 'fine'))]
    rows = []
    for target in targets:
        for name in ('sham', f'reset-cell-{target}'):
            for dose in DOSES:
                pair = [next(j for j in jobs if j['level'] == level and j['background'] == name and
                    j['target'] == target and j['dose'] == dose) for level in ('coarse', 'fine')]
                if pair[0]['factor'] != pair[1]['factor']:
                    raise ValueError('Paired pulse factors differ')
                controls = [np.array([r['chemistry'] for r in h[name, None, None]]) for h in (coarse, fine)]
                paths = [np.array([r['chemistry'] for r in h[name, target, dose]]) for h in (coarse, fine)]
                error = response_error(paths[0], controls[0], paths[1], controls[1], pair[0]['factor'])
                metrics = [next(r for r in s if r['cell'] == target)['contexts'][name]['doses'][dose]['metrics'] for s in summaries]
                metric = 'target_activator_log_auc_per_log_pulse'
                auc_error = abs(metrics[1][metric]/metrics[0][metric]-1)
                recovery = {}
                for key in ('target_recovery_time', 'network_recovery_time'):
                    a, b = [r[key] for r in metrics]
                    recovery[key] = 0. if a is None and b is None else None if None in (a, b) else abs(a-b)
                limits = p['response_limits']
                passed = error <= limits['max_normalized_response_error'] and auc_error <= limits['max_relative_auc_error'] and all(
                    value is not None and value <= limits['max_recovery_time_error']+1e-9 for value in recovery.values())
                rows.append(dict(context=name, cell=target, dose=dose, normalized_response_error=error,
                    relative_auc_error=auc_error, recovery_time_errors=recovery, passed=bool(passed)))
    same = True; shift_error = 0.
    for a, b in zip(*summaries):
        for dose in DOSES:
            same &= a['differences'][dose]['response_shift_detected'] == b['differences'][dose]['response_shift_detected']
            shift_error = max(shift_error, abs(a['differences'][dose]['normalized_response_shift_rms']-
                b['differences'][dose]['normalized_response_shift_rms']))
    return dict(passed=bool(all(r['passed'] for r in rows) and same and shift_error <= .01), rows=rows,
        same_response_effect_decisions=bool(same), response_shift_error_max=shift_error,
        coarse=summaries[0], fine=summaries[1])


def pair_report(root, p, seed, horizon):
    root = Path(root); label = 'pilot' if horizon == p['pilot_duration'] else 'full'
    folder = root/'pairs'/f'seed-{seed}'; folder.mkdir(parents=True, exist_ok=True)
    file = folder/f'{label}-refinement.json'; old = read(file) if file.exists() else None
    ph = digest(root/'protocol.json')
    if old:
        if old['protocol_sha256'] != ph:
            raise ValueError('Changed paired protocol')
        for path, sha in old['evidence_sha256'].items():
            if digest(path) != sha:
                raise ValueError('Changed paired evidence')
    jobs = [j for j in p['jobs'] if j['seed'] == seed]; histories = {}; files = []
    for level in ('coarse', 'fine'):
        histories[level] = {}
        for job in [j for j in jobs if j['level'] == level]:
            history = read(root/job['key']/'history.json')
            history = [r for r in history if r['elapsed'] <= horizon+1e-9]
            moving.validate_history(history, job, p, horizon)
            histories[level][job['background'], job['target'], job['dose']] = history
            stem = folder/f'{label}_{job["key"]}'; hf, pf = stem.with_suffix('.json'), stem.with_suffix('.npz')
            if old:
                if read(hf) != history:
                    raise ValueError('Changed paired history')
            else:
                state = moving.payload(root/job['key']/'latest_state.npz')
                if abs(json.loads(str(state['metadata']))['time']-job['start']-horizon) > 1e-9:
                    raise ValueError('Field snapshot clock differs')
                write_json(hf, history); np.savez_compressed(pf, phi=state['phi'])
            files.extend((hf, pf))
    mechanical = []
    for key in histories['coarse']:
        fields = []
        for level in ('coarse', 'fine'):
            job = next(j for j in jobs if j['level'] == level and (j['background'], j['target'], j['dose']) == key)
            with np.load(folder/f'{label}_{job["key"]}.npz') as z:
                fields.append(z['phi'].copy())
        mechanical.append(dict(context=key[0], target=key[1], dose=key[2], **moving.mechanical_comparison(
            histories['coarse'][key], histories['fine'][key], fields[0], fields[1], p, horizon)))
    responses = response_comparison(histories['coarse'], histories['fine'], jobs,
        p['selections'][str(seed)]['selected_ids'], p, horizon)
    result = dict(seed=seed, horizon=horizon, protocol_sha256=ph,
        passed=bool(all(r['passed'] for r in mechanical) and responses['passed']),
        mechanical=mechanical, responses=responses, evidence_sha256={str(f.resolve()):digest(f) for f in files})
    if old and old != result:
        raise ValueError('Changed paired decision')
    write_json(file, result)
    print('Delayed paired gate:', seed, horizon, result['passed'], flush=True)
    return result


def schedule_history(root, p, seed, done):
    jobs = [j for j in p['jobs'] if j['seed'] == seed]
    for horizon in (p['pilot_duration'], p['duration']):
        for job in jobs:
            write_json(Path(root)/'status.json', dict(state='running', stage='pilot' if horizon == p['pilot_duration'] else 'full',
                current_job=job['key'], target_horizon=horizon, completed=done, total=p['total_jobs'],
                independent_histories=3, new_histories=0))
            advance(root, job, p, horizon)
            if horizon == p['duration']:
                done += 1
        if not pair_report(root, p, seed, horizon)['passed']:
            return done, False
    return done, True


def preflight(root):
    """Test the smallest matched-amount pulses at full size and both dt levels."""
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    moving.require_device(p['device']); torch.set_num_threads(1)
    jobs = [j for j in p['jobs'] if j['seed'] in (8, 9) and j['background'] == 'sham' and
        j['target'] == p['selections'][str(j['seed'])]['selected_ids'][-1] and j['dose'] == 'amount']
    records = []
    for job in jobs:
        write_json(root/'status.json', dict(state='preflighting', current_job=job['key'],
            prefix_contexts_completed=len(records), prefix_contexts_total=len(jobs), completed=0, total=p['total_jobs']))
        records.append(moving.context_gate(root, job, p))
    result = dict(passed=bool(all(r['passed'] for r in records)), protocol_sha256=digest(root/'protocol.json'),
        full_size_contexts=len(records), contexts=records,
        evidence_sha256={str((root/j['key']/'prefix/result.json').resolve()):digest(root/j['key']/'prefix/result.json') for j in jobs},
        scope='Four 72^3 low-amplitude matched-dose sham contexts; remaining contexts require individual gates and every history requires its paired pilot before long continuation.')
    write_json(root/'preflight.json', result)
    write_json(root/'status.json', dict(state='prepared_gpu_preflight_passed', completed=0, total=p['total_jobs'],
        actual_preflight_contexts=len(records), independent_histories=3, new_histories=0))
    return result


def assess(root, p=None):
    root = Path(root); p = read(root/'protocol.json') if p is None else p; verify(p)
    results, missing, by_history = [], [], []
    for job in p['jobs']:
        if not (root/job['key']/'result.json').exists():
            missing.append(job['key']); continue
        result, _ = moving.current_evidence(root, job, p)
        results.append(dict(job=job, evidence=result))
    for seed in p['histories']:
        reports = {}
        if all(j['key'] not in missing for j in p['jobs'] if j['seed'] == seed):
            reports = {label:pair_report(root, p, seed, horizon) for label, horizon in
                (('pilot', p['pilot_duration']), ('full', p['duration']))}
        by_history.append(dict(seed=seed, qualified=bool(len(reports) == 2 and all(r['passed'] for r in reports.values())), reports=reports))
    result = dict(protocol_sha256=digest(root/'protocol.json'), completed=len(results), total=p['total_jobs'],
        passed=bool(not missing and all(r['qualified'] for r in by_history)), missing=missing,
        independent_histories=3, new_histories=0, results=results, by_history=by_history, scope=p['scope'])
    write_json(root/'summary.json', result)
    return result


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); moving.require_device(p['device']); torch.set_num_threads(1)
        preflight_record = read(root/'preflight.json')
        if not preflight_record['passed'] or preflight_record['protocol_sha256'] != digest(root/'protocol.json'):
            raise ValueError('Full-size preflight required')
        for path, sha in preflight_record['evidence_sha256'].items():
            if digest(path) != sha:
                raise ValueError('Changed preflight evidence')
        done, failed = 0, []
        try:
            for seed in p['histories']:
                done, passed = schedule_history(root, p, seed, done)
                if not passed:
                    failed.append(dict(seed=seed, reason='prespecified_timestep_gate_failed')); break
            result = assess(root, p)
            write_json(root/'status.json', dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'], completed=done, total=p['total_jobs'], failed_histories=failed,
                independent_histories=3, new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=p['total_jobs'], error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'preflight', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-delayed-response'))
    args = parser.parse_args()
    result = dict(prepare=prepare, preflight=preflight, run=run, assess=assess)[args.action](args.output)
    print({k:v for k,v in result.items() if k in ('passed', 'completed', 'total', 'estimated_moving_seconds')})
