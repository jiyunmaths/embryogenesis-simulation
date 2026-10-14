"""Delayed responses after amount-preserving surrounding redistribution.

The existing carry, transport and response kernels remain unchanged. Continue
the actual independently reviewed fine t=510 backgrounds at two timesteps;
reuse a historical path only after exact regenerated source verification.
"""
import argparse
from dataclasses import asdict
import fcntl
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import torch

from . import phase_carry_context_dosage as dosage
from . import phase_carry_delayed_response as delayed
from . import phase_carry_exchange_response as moving
from . import phase_carry_network_context as context
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify
from .resolution import write_json

HISTORIES = (7, 8, 9)
LEVELS = moving.LEVELS
DOSES = delayed.DOSES
pulse = delayed.pulse
matched_removal = delayed.matched_removal
source_equivalence = context.source_equivalence
RESERVE_BYTES = 2*2**30


def job_design(root, selections):
    root = Path(root); contexts, jobs = [], []
    for seed in HISTORIES:
        targets = selections[str(seed)]['selected_ids']
        for recipient in (None, *targets):
            name = 'sham' if recipient is None else f'redistributed-cell-{recipient}'
            contexts.append(dict(seed=seed, name=name, recipient=recipient,
                targets=targets if recipient is None else [recipient]))
            for level, dt in LEVELS:
                for target, dose in [(None, None), *[(t, d) for t in (targets if recipient is None else [recipient]) for d in DOSES]]:
                    suffix = 'control' if target is None else f'cell-{target}_{dose}'
                    key = f'seed-{seed}_{level}_{name}_{suffix}'
                    jobs.append(dict(key=key, seed=seed, level=level, dt=dt,
                        stage='redistribution_response', background=name, recipient=recipient,
                        target=target, dose=dose, source=str(root/key/'source.npz'),
                        start=510., duration=60., point=dict(ratio=27.5, chi=.35)))
    return contexts, jobs


def reviewed_parent(root, review_file):
    """Require a completed parent and its immutable independent assessment."""
    root, review_file = Path(root).resolve(), Path(review_file).resolve()
    old = read(root/'protocol.json'); verify(old); review = read(review_file)
    assessment = review.get('full_assessment', review.get('full_review'))
    assessment_sha = review.get('full_assessment_sha256', review.get('full_review_sha256'))
    if not assessment or not assessment_sha:
        raise ValueError('Independent full assessment and hash required')
    if (read(root/'status.json')['state'] != 'completed' or not read(root/'summary.json')['passed'] or
            not review['full_study_accepted'] or review['protocol_sha256'] != digest(root/'protocol.json') or
            review['original_summary_sha256'] != digest(root/'summary.json') or
            review['original_status_sha256'] != digest(root/'status.json') or
            assessment_sha != digest(assessment)):
        raise ValueError('Requires independently reviewed completed parent: '+str(root))
    full = read(assessment)
    for key in ('full_study_accepted', 'protocol_sha256', 'original_summary_sha256', 'original_status_sha256'):
        if full.get(key) != review[key]:
            raise ValueError('Independent assessment disagrees with public review: '+key)
    evidence = full.get('evidence_sha256', {})
    for path, sha in evidence.items():
        if digest(path) != sha:
            raise ValueError('Changed independently assessed evidence: '+path)
    inputs = [root/f for f in ('protocol.json', 'status.json', 'summary.json')]
    inputs += [review_file, Path(assessment), *[Path(f) for f in evidence]]
    return old, inputs


def parent_endpoint(root, p, job):
    """Resolve the actual checkpoint through the dosage/network reuse chain."""
    root = Path(root)
    if job['key'] not in p.get('reused', {}):
        return root/job['key']/'latest_state.npz'
    prior = p['reused'][job['key']]['job']; parent = Path(p['parent'])
    return parent_endpoint(parent, read(parent/'protocol.json'), prior)


def storage_estimate(root, p, new_paths=48):
    """Observed complete delayed path plus paired evidence, with 12% padding."""
    root = Path(root); measured = {}
    for job in p['jobs']:
        folder = root/job['key']; pairs = root/'pairs'/f'seed-{job["seed"]}'
        files = [f for f in folder.rglob('*') if f.is_file()]
        files += [pairs/f'{label}_{job["key"]}{suffix}'
            for label in ('pilot', 'full') for suffix in ('.npz', '.json')]
        if not files or any(not f.is_file() for f in files):
            raise ValueError('Complete historical delayed footprint required')
        measured[job['key']] = sum(f.stat().st_size for f in files)
    if not measured:
        raise ValueError('Nonempty completed delayed paths required')
    maximum = max(measured.values())
    return dict(observed_path_bytes=measured, maximum_observed_path_bytes=maximum,
        padding_fraction=.12, new_paths=int(new_paths),
        estimated_additional_bytes=int(np.ceil(new_paths*maximum*1.12)),
        reserve_bytes=RESERVE_BYTES)


def prepare(root, parent=Path('outputs/phase-carry-context-dosage'),
            reuse_parent=Path('outputs/phase-carry-delayed-response')):
    root, parent, reuse_parent = [Path(f).resolve() for f in (root, parent, reuse_parent)]
    if root.exists():
        raise FileExistsError(root)
    old, inputs = reviewed_parent(parent, Path('docs/phase_carry_context_dosage_assessment.json'))
    last, previous_inputs = reviewed_parent(reuse_parent, Path('docs/phase_carry_delayed_response_assessment.json'))
    if (old['latest_review_parent'] != str(reuse_parent) or old['device'] != last['device'] or
            old['selections'] != last['selections'] or old['parent'] != last['parent']):
        raise ValueError('Historical delayed assay and dosage parent are incompatible')
    inputs += previous_inputs + [Path('docs/phase_carry_redistribution_response.md').resolve()]
    footprint = storage_estimate(reuse_parent, last)
    expected_bytes = footprint['estimated_additional_bytes']
    if shutil.disk_usage(root.parent).free < expected_bytes+RESERVE_BYTES:
        raise OSError('Insufficient free space for 48 new paths plus 2 GiB reserve')
    root.mkdir(parents=True)
    contexts, jobs = job_design(root, old['selections'])
    reused, preparations, accepted, common_doses = {}, {}, {}, {}
    p = dict(parent=str(parent), reuse_parent=str(reuse_parent), device=old['device'],
        selections=old['selections'], contexts=contexts, jobs=jobs, reused=reused,
        preparations=preparations, accepted_configs=accepted, common_doses=common_doses,
        histories=list(HISTORIES), independent_histories=3, new_histories=0,
        start=510., duration=60., pilot_duration=6., interval=.15, late_window=24.,
        checkpoint_interval=3., prefix_duration=.6, prefix_native_threads=4,
        amount_error_max=2e-14, contrast_min=.1, effect_min=.01,
        mechanical_limits=moving.MECHANICAL_LIMITS, response_limits=moving.RESPONSE_CRITERIA,
        total_jobs=66, new_jobs=48, reused_jobs=18, estimated_additional_bytes=expected_bytes,
        storage_reserve_bytes=RESERVE_BYTES, storage_estimate=footprint,
        design='Delayed -10% and common feasible removed-amount activator pulses from each actual fine t=510 sham/amount-preserving redistribution endpoint. Identical physical state per background at both timesteps; own live moving controls. Exact sham control/fractional paths reused; new common-amount paths and redistributed backgrounds individually qualified.',
        scope='Three existing mature histories and six nested recipient comparisons, no new developmental histories. Initial amounts were matched at t=450, not necessarily at t=510 or later. Backgrounds may differ in chemistry and physical context. Delayed response behavior, not geometry mediation, neighbor-only topology, autonomous/inherited biological identity, spatial convergence or population inference.')
    estimate = 0.
    try:
        for seed in HISTORIES:
            states = {}
            for ctx in [c for c in contexts if c['seed'] == seed]:
                arm = 'sham' if ctx['recipient'] is None else 'redistributed'
                prior = next(j for j in old['jobs'] if (j['seed'], j['level'], j['arm'], j['recipient']) ==
                    (seed, 'fine', arm, ctx['recipient']))
                result, history = dosage.current_evidence(parent, prior, old)
                source = parent_endpoint(parent, old, prior)
                state = moving.payload(source); meta = json.loads(str(state['metadata']))
                values, masses, ids = np.asarray(history[-1]['chemistry']), np.asarray(history[-1]['volumes']), history[-1]['ids']
                if (abs(meta['time']-510.) > 1e-9 or not np.array_equal(state['ids'], ids) or
                        not np.array_equal([state['activator'], state['inhibitor']], values)):
                    raise ValueError('Fine parent endpoint clock/IDs/chemistry mismatch')
                ctx.update(parent_job=prior, source=str(source))
                states[ctx['name']] = (source, values, masses, ids, result['audit']['wall_seconds'])
                preparations[f'{seed}_{ctx["name"]}'] = dict(source=str(source), ids=ids,
                    chemistry=values.tolist(), measured_gpu_masses=masses.tolist())
                inputs.extend(source.parent/f for f in ('latest_state.npz', 'history.json', 'result.json', 'prefix/result.json'))
                inputs.extend(Path(f) for f in read(source.parent/'prefix/result.json')['evidence_sha256'])
            for target in old['selections'][str(seed)]['selected_ids']:
                amounts = []
                for name in ('sham', f'redistributed-cell-{target}'):
                    _, values, masses, ids, _ = states[name]; i = ids.index(target)
                    amounts.append(float(values[0, i]*masses[i]))
                common_doses[f'{seed}_{target}'] = dict(background_target_amounts=amounts,
                    removed_amount=matched_removal(*amounts))
            for j in [j for j in jobs if j['seed'] == seed]:
                source, values, masses, ids, wall = states[j['background']]
                changed = values.copy()
                record = dict(factor=1., removed_activator_amount=0., requested_removed_amount=0.,
                    initial_target_activator_amount=None, initial_log_amplitude=0.)
                if j['target'] is not None:
                    changed, record = pulse(values, masses, ids, j['target'], j['dose'],
                        common_doses[f'{seed}_{j["target"]}']['removed_amount'])
                j.update(record)
                if j['background'] == 'sham' and (j['target'] is None or j['dose'] == 'fractional'):
                    previous = next(a for a in last['jobs'] if (a['seed'], a['level'], a['background'], a['target'], a['dose']) ==
                        (seed, j['level'], 'sham', j['target'], j['dose']))
                    moving.current_evidence(reuse_parent, previous, last)
                    with tempfile.TemporaryDirectory(dir=root) as temporary:
                        expected = Path(temporary)/'expected.npz'
                        host = moving.prepared_checkpoint(source, expected, j['dt'], changed)
                        source_equivalence(expected, previous['source'])
                    # No reuse based on similar endpoints or equal parameters.
                    j['source'] = previous['source']
                    folder = reuse_parent/previous['key']; pair_folder = reuse_parent/'pairs'/f'seed-{seed}'
                    paths = [folder/f for f in ('latest_state.npz', 'history.json', 'result.json', 'prefix/result.json')]
                    paths += [Path(previous['source']), *[Path(f) for f in read(folder/'prefix/result.json')['evidence_sha256']]]
                    paths += [pair_folder/f'{label}_{previous["key"]}{suffix}' for label in ('pilot', 'full') for suffix in ('.npz', '.json')]
                    reused[j['key']] = dict(job=previous, protocol_sha256=digest(reuse_parent/'protocol.json'),
                        source_exact=True, evidence_sha256={str(f.resolve()):digest(f) for f in paths},
                        pilot_phi=str(pair_folder/f'pilot_{previous["key"]}.npz'),
                        full_phi=str(pair_folder/f'full_{previous["key"]}.npz'))
                    inputs += paths
                else:
                    host = moving.prepared_checkpoint(source, j['source'], j['dt'], changed)
                    inputs.append(Path(j['source']))
                    estimate += wall*(.001875/j['dt'])
                key = f'{seed}_{j["level"]}'
                if key in accepted and accepted[key] != asdict(host.config):
                    raise ValueError('Background/pulse configurations differ')
                accepted[key] = asdict(host.config)
            print('Prepared delayed amount-preserving redistribution contexts:', seed, flush=True)
        if len(reused) != 18 or len(jobs)-len(reused) != 48:
            raise ValueError('Unexpected exact-reuse design counts')
        files = [Path(__file__), Path('tests/test_phase_carry_redistribution_response.py')]
        p['source_sha256'] = {**old['source_sha256'], **{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256'] = {**old['input_sha256'], **last['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}}
        p['estimated_moving_seconds'] = estimate
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', completed=0, total=66, new_jobs=48,
            reused_jobs=18, independent_histories=3, new_histories=0, estimated_moving_seconds=estimate))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error))); raise


def current_evidence(root, job, p):
    if job['key'] not in p['reused']:
        return moving.current_evidence(root, job, p)
    item = p['reused'][job['key']]; parent = Path(p['reuse_parent'])
    if not item.get('source_exact') or digest(parent/'protocol.json') != item['protocol_sha256']:
        raise ValueError('Changed or nonexact reused protocol/source')
    source_equivalence(job['source'], item['job']['source'])
    result, history = moving.current_evidence(parent, item['job'], read(parent/'protocol.json'))
    moving.validate_history(history, job, p, job['duration'])
    for path, sha in item['evidence_sha256'].items():
        if digest(path) != sha:
            raise ValueError('Changed reused delayed evidence')
    return result, history


def advance(root, job, p, horizon):
    if job['key'] in p['reused']:
        _, history = current_evidence(root, job, p)
        return [r for r in history if r['elapsed'] <= horizon+1e-9]
    if shutil.disk_usage(root).free < p['storage_reserve_bytes']:
        raise OSError('Less than 2 GiB storage reserve; committed checkpoints retained')
    return delayed.advance(root, job, p, horizon)


def _delayed_view(histories, jobs):
    """Adapt background labels only; reuse unchanged response calculations."""
    def name(value):
        return value.replace('redistributed-cell-', 'reset-cell-')
    h = {(name(bg), target, dose): rows for (bg, target, dose), rows in histories.items()}
    j = [dict(row, background=name(row['background'])) for row in jobs]
    return h, j


def _redistributed_labels(report):
    if isinstance(report, dict):
        return {key.replace('reset-cell-', 'redistributed-cell-'): _redistributed_labels(value)
            for key, value in report.items()}
    if isinstance(report, list):
        return [_redistributed_labels(value) for value in report]
    if isinstance(report, str):
        return report.replace('reset-cell-', 'redistributed-cell-')
    return report


def response_summary(histories, jobs, targets, p, horizon):
    h, j = _delayed_view(histories, jobs)
    return _redistributed_labels(delayed.response_summary(h, j, targets, p, horizon))


def response_comparison(coarse, fine, jobs, targets, p, horizon):
    a, j = _delayed_view(coarse, jobs); b, _ = _delayed_view(fine, jobs)
    return _redistributed_labels(delayed.response_comparison(a, b, j, targets, p, horizon))


def pilot_quality(root, p, seed):
    """Block long continuation when any pilot already violates physical limits."""
    root = Path(root); file = root/'pairs'/f'seed-{seed}'/'pilot-quality.json'
    file.parent.mkdir(parents=True, exist_ok=True); ph = digest(root/'protocol.json')
    if file.exists():
        record = read(file)
        if record['protocol_sha256'] != ph or not record['passed']:
            raise ValueError('Changed or failed pilot physical gate')
        return record
    rows = []
    for j in [j for j in p['jobs'] if j['seed'] == seed]:
        if j['key'] in p['reused']:
            audit = current_evidence(root, j, p)[0]['audit']; evidence = p['reused'][j['key']]['evidence_sha256']
        else:
            cp = root/j['key']/'latest_state.npz'; state = moving.payload(cp)
            saved = json.loads(str(state['long_experiment'])); meta = json.loads(str(state['metadata']))
            if (saved['protocol_hash'] != ph or saved['job'] != j or
                    abs(meta['time']-j['start']-p['pilot_duration']) > 1e-9 or
                    abs(saved['history'][-1]['elapsed']-p['pilot_duration']) > 1e-9):
                raise ValueError('Pilot physical checkpoint does not match horizon/protocol')
            audit = saved['audit']; evidence = {str(cp.resolve()):digest(cp)}
        rows.append(dict(job=j['key'], audit=audit, passed=bool(moving.quality_pass(audit, p)), evidence_sha256=evidence))
    record = dict(protocol_sha256=ph, seed=seed, passed=bool(rows and all(r['passed'] for r in rows)), rows=rows,
        scope='New pilot audits and reused full-window audits; full audits bound their pilot interval.')
    write_json(file, record)
    if not record['passed']:
        raise RuntimeError('Pilot physical quality failed; long continuation blocked')
    return record


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
    jobs = [j for j in p['jobs'] if j['seed'] == seed]; histories = {}; fields = {}; paths = []
    if horizon == p['pilot_duration']:
        quality = pilot_quality(root, p, seed); paths.append(folder/'pilot-quality.json')
    for level in ('coarse', 'fine'):
        histories[level] = {}
        for j in [j for j in jobs if j['level'] == level]:
            reused = j['key'] in p['reused']
            raw = current_evidence(root, j, p)[1] if reused else read(root/j['key']/'history.json')
            history = [r for r in raw if r['elapsed'] <= horizon+1e-9]
            moving.validate_history(history, j, p, horizon)
            key = (j['background'], j['target'], j['dose']); histories[level][key] = history
            hf = folder/f'{label}_{j["key"]}.json'
            pf = Path(p['reused'][j['key']][label+'_phi']) if reused else folder/f'{label}_{j["key"]}.npz'
            if old:
                if read(hf) != history:
                    raise ValueError('Changed paired history')
            else:
                write_json(hf, history)
                if not reused:
                    state = moving.payload(root/j['key']/'latest_state.npz')
                    if abs(json.loads(str(state['metadata']))['time']-j['start']-horizon) > 1e-9:
                        raise ValueError('Field snapshot clock differs')
                    np.savez_compressed(pf, phi=state['phi'])
            with np.load(pf) as z:
                fields[level, key] = z['phi'].copy()
            paths += [hf, pf]
    mechanical = [dict(context=key[0], target=key[1], dose=key[2], **moving.mechanical_comparison(
        histories['coarse'][key], histories['fine'][key], fields['coarse', key], fields['fine', key], p, horizon))
        for key in histories['coarse']]
    responses = response_comparison(histories['coarse'], histories['fine'], jobs,
        p['selections'][str(seed)]['selected_ids'], p, horizon)
    result = dict(seed=seed, horizon=horizon, protocol_sha256=ph,
        passed=bool(all(r['passed'] for r in mechanical) and responses['passed']),
        mechanical=mechanical, responses=responses, evidence_sha256={str(f.resolve()):digest(f) for f in paths})
    if horizon == p['pilot_duration']:
        result['physical_quality_passed'] = quality['passed']
    if old and old != result:
        raise ValueError('Changed paired decision')
    if not old:
        write_json(file, result)
    print('Redistribution delayed paired gate:', seed, horizon, result['passed'], flush=True)
    return result


def schedule_history(root, p, seed, done):
    jobs = [j for j in p['jobs'] if j['seed'] == seed]
    for horizon in (p['pilot_duration'], p['duration']):
        for j in jobs:
            write_json(Path(root)/'status.json', dict(state='running', stage='pilot' if horizon==p['pilot_duration'] else 'full',
                current_job=j['key'], target_horizon=horizon, completed=done, total=66,
                new_jobs=48, reused_jobs=18, independent_histories=3, new_histories=0))
            advance(root, j, p, horizon)
            if horizon == p['duration']:
                done += 1
        # Explicitly check quality even when the pair implementation is replaced
        # in a unit test; no good refinement result can override bad physics.
        if horizon == p['pilot_duration']:
            pilot_quality(root, p, seed)
        if not pair_report(root, p, seed, horizon)['passed']:
            return done, False
    return done, True


def preflight(root):
    """Qualify four weak common-amount sham pulses before expensive work."""
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    moving.require_device(p['device']); torch.set_num_threads(1)
    jobs = [j for j in p['jobs'] if j['seed'] in (8, 9) and j['background'] == 'sham' and
        j['target'] == p['selections'][str(j['seed'])]['selected_ids'][-1] and j['dose'] == 'amount']
    if len(jobs) != 4 or any(j['key'] in p['reused'] for j in jobs):
        raise ValueError('Expected four newly prepared weak common-amount sham cases')
    records = []
    for j in jobs:
        if shutil.disk_usage(root).free < p['storage_reserve_bytes']:
            raise OSError('Insufficient preflight storage reserve')
        write_json(root/'status.json', dict(state='preflighting', current_job=j['key'], completed=0, total=66,
            prefix_contexts_completed=len(records), prefix_contexts_total=len(jobs)))
        records.append(moving.context_gate(root, j, p))
    result = dict(passed=bool(all(r['passed'] for r in records)), protocol_sha256=digest(root/'protocol.json'),
        full_size_contexts=len(records), contexts=records,
        evidence_sha256={str((root/j['key']/'prefix/result.json').resolve()):digest(root/j['key']/'prefix/result.json') for j in jobs},
        scope='Four 72^3 weak common-amount sham contexts; every remaining new context requires its own gate and each history requires physical and paired pilot qualification before full continuation.')
    write_json(root/'preflight.json', result)
    write_json(root/'status.json', dict(state='prepared_gpu_preflight_passed', completed=0, total=66,
        new_jobs=48, reused_jobs=18, independent_histories=3, new_histories=0))
    return result


def assess(root, p=None):
    root = Path(root); p = read(root/'protocol.json') if p is None else p; verify(p)
    results, missing, rows = [], [], []
    for j in p['jobs']:
        if j['key'] not in p['reused'] and not (root/j['key']/'result.json').exists():
            missing.append(j['key']); continue
        result, _ = current_evidence(root, j, p)
        results.append(dict(job=j, reused=j['key'] in p['reused'], evidence=result))
    for seed in p['histories']:
        reports = {}
        if all(j['key'] not in missing for j in p['jobs'] if j['seed'] == seed):
            reports = {label:pair_report(root, p, seed, horizon) for label, horizon in
                (('pilot', p['pilot_duration']), ('full', p['duration']))}
        rows.append(dict(seed=seed, qualified=bool(len(reports)==2 and all(r['passed'] for r in reports.values())), reports=reports))
    result = dict(protocol_sha256=digest(root/'protocol.json'), passed=bool(not missing and all(r['qualified'] for r in rows)),
        completed=len(results), total=66, new_jobs=48, reused_jobs=18, independent_histories=3,
        new_histories=0, missing=missing, results=results, by_history=rows, scope=p['scope'])
    write_json(root/'summary.json', result)
    return result


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); moving.require_device(p['device']); torch.set_num_threads(1)
        pre = read(root/'preflight.json')
        if not pre['passed'] or pre['protocol_sha256'] != digest(root/'protocol.json'):
            raise ValueError('Qualified full-size preflight required')
        for path, sha in pre['evidence_sha256'].items():
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
                passed=result['passed'], completed=done, total=66, new_jobs=48, reused_jobs=18,
                failed_histories=failed, independent_histories=3, new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=66,
                new_jobs=48, reused_jobs=18, error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'preflight', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-redistribution-response'))
    args = parser.parse_args()
    result = dict(prepare=prepare, preflight=preflight, run=run, assess=assess)[args.action](args.output)
    print({k:v for k,v in result.items() if k in ('passed', 'completed', 'total', 'estimated_moving_seconds')})
