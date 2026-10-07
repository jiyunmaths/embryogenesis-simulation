"""Separate surrounding spatial chemistry from initial bulk chemical dosage.

Use the qualified moving carry kernels unchanged. Exact sham/raw-reset paths
are reused; the two new factorial arms preserve the recipient at preparation.
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

from . import phase_carry_delayed_response as delayed
from . import phase_carry_exchange_response as moving
from . import phase_carry_network_context as context
from .attribute_development import AttributeSimulation
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify
from .resolution import write_json

HISTORIES = (7, 8, 9)
LEVELS = moving.LEVELS
ARMS = ('sham', 'reset', 'redistributed', 'bulk')
FACTORS = dict(sham=('original', 'original'), reset=('reference', 'reference'),
               redistributed=('reference', 'original'), bulk=('original', 'reference'))
CONTRASTS = dict(spatial_original_dose=('redistributed', 'sham'),
                 spatial_reference_dose=('reset', 'bulk'),
                 dosage_original_distribution=('bulk', 'sham'),
                 dosage_reference_distribution=('reset', 'redistributed'),
                 total_reset=('reset', 'sham'))


def surrounding_state(initial, reference, masses, ids, recipient, arm):
    """Factorial distribution/dosage intervention, separately for each species.

    Distribution denotes relative concentrations over nonrecipient cells.
    Dosage denotes their volume-weighted amount. The recipient is untouched.
    """
    x, ref, m, ids = [np.asarray(a) for a in (initial, reference, masses, ids)]
    if (m.ndim != 1 or len(m) < 2 or x.shape != (2, len(m)) or ref.shape != x.shape or
            ids.shape != m.shape or len(set(ids.tolist())) != len(m) or recipient not in ids or
            arm not in ARMS or any(not np.isfinite(a).all() or np.any(a <= 0) for a in (x, ref, m))):
        raise ValueError('Invalid positive chemical state, masses, IDs, recipient or arm')
    i = ids.tolist().index(recipient); selected = np.arange(len(m)) != i
    original = x[:, selected] @ m[selected]
    reference_amounts = ref[:, selected] @ m[selected]
    distribution, dosage = FACTORS[arm]
    base = ref if distribution == 'reference' else x
    desired = reference_amounts if dosage == 'reference' else original
    base_amounts = base[:, selected] @ m[selected]
    if any(not np.isfinite(a).all() or np.any(a <= 0) for a in (original, reference_amounts, base_amounts)):
        raise ValueError('Surrounding amount is not representable')
    factors = desired / base_amounts
    out = np.array(x, dtype=np.float64, copy=True)
    # Exact copies for the two historical arms, avoiding a needless roundoff.
    if arm == 'reset':
        out[:, selected] = ref[:, selected]
    elif arm != 'sham':
        out[:, selected] = base[:, selected] * factors[:, None]
    actual = out[:, selected] @ m[selected]
    relative_error = float(np.max(abs(actual-desired)/desired))
    if (not np.isfinite(out).all() or np.any(out <= 0) or
            relative_error > 2e-14 or not np.array_equal(out[:, i], x[:, i])):
        raise ValueError('Intervention positivity, amount or recipient preservation failed')
    log_shift = np.log(out[:, selected])-np.log(x[:, selected])
    rms = float(np.sqrt(np.sum(log_shift**2*m[selected])/ (2*m[selected].sum())))
    return out, dict(arm=arm, distribution=distribution, dosage=dosage, recipient=int(recipient),
        changed_ids=ids[selected].astype(int).tolist(), scale_factors=factors.tolist(),
        original_surrounding_amounts=original.tolist(), reference_surrounding_amounts=reference_amounts.tolist(),
        desired_surrounding_amounts=desired.tolist(), actual_surrounding_amounts=actual.tolist(),
        relative_amount_error=relative_error, added_amounts=((out-x)@m).tolist(),
        initial_total_amounts=(x@m).tolist(), final_total_amounts=(out@m).tolist(),
        surrounding_log_rms=rms, recipient_exact=True)


def job_design(root, selections):
    jobs = []
    for seed in HISTORIES:
        contexts = [('sham', None)] + [(arm, target) for target in selections[str(seed)]['selected_ids']
            for arm in ARMS if arm != 'sham']
        for arm, recipient in contexts:
            for level, dt in LEVELS:
                name = arm if recipient is None else f'{arm}-cell-{recipient}'
                key = f'seed-{seed}_{level}_{name}_control'
                jobs.append(dict(key=key, seed=seed, level=level, dt=dt, arm=arm,
                    stage='context_dosage', background=name, recipient=recipient, target=None, factor=1.,
                    source=str(Path(root)/key/'source.npz'), start=450., duration=60.,
                    point=dict(ratio=27.5, chi=.35)))
    return jobs


def reviewed_parent(root, review_file):
    old = read(root/'protocol.json'); verify(old); review = read(review_file)
    if (read(root/'status.json')['state'] != 'completed' or not read(root/'summary.json')['passed'] or
            not review['full_study_accepted'] or review['protocol_sha256'] != digest(root/'protocol.json') or
            review['original_summary_sha256'] != digest(root/'summary.json') or
            review['original_status_sha256'] != digest(root/'status.json') or
            review['full_review_sha256'] != digest(review['full_review'])):
        raise ValueError('Requires independently reviewed completed parent: '+str(root))
    return old, [root/f for f in ('protocol.json', 'status.json', 'summary.json')] + [review_file, Path(review['full_review'])]


def prepare(root, parent=Path('outputs/phase-carry-network-context')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists():
        raise FileExistsError(root)
    old, inputs = reviewed_parent(parent, Path('docs/phase_carry_network_context_assessment.json').resolve())
    last, latest_inputs = reviewed_parent(Path('outputs/phase-carry-delayed-response').resolve(),
        Path('docs/phase_carry_delayed_response_assessment.json').resolve())
    if last['parent'] != str(parent) or last['device'] != old['device']:
        raise ValueError('Latest reviewed study does not match the original reset parent')
    inputs += latest_inputs + [Path('docs/phase_carry_context_dosage.md').resolve()]
    # Reserve the observed historical path footprint plus assessment overhead.
    expected_bytes = int(24*.35*2**30)
    if shutil.disk_usage(root.parent).free < expected_bytes+2*2**30:
        raise OSError('Insufficient free space for 24 new paths plus 2 GiB reserve')
    root.mkdir(parents=True)
    jobs = job_design(root, old['selections']); reused, preparations, accepted = {}, {}, {}
    p = dict(parent=str(parent), latest_review_parent=str(Path('outputs/phase-carry-delayed-response').resolve()), device=old['device'], selections=old['selections'],
        jobs=jobs, reused=reused, preparations=preparations, accepted_configs=accepted,
        histories=list(HISTORIES), independent_histories=3, new_histories=0, start=450., duration=60.,
        pilot_duration=6., interval=.15, late_window=24., checkpoint_interval=3., prefix_duration=.6,
        prefix_native_threads=4, amount_error_max=2e-14, contrast_min=.1, effect_min=.01,
        state_metric_error_max=.01, mechanical_limits=moving.MECHANICAL_LIMITS,
        total_jobs=42, new_jobs=24, reused_jobs=18, estimated_additional_bytes=expected_bytes,
        design='Two-species distribution/dosage factorial over all nonrecipient cells. Identical fine t=450 physical start at both dt levels; recipient chemistry exact. Live carry chemistry, polarity, mechanics, conservative transport and dilution. Unperturbed state stage only.',
        scope='Three existing mature histories, six nested recipients; no fresh developmental histories. Initial species amounts or distribution are controlled, not their later values. No neighbor-only topology, geometry mediation, response behavior, autonomous or inherited biological identity, spatial convergence or population inference.')
    estimate = 0.
    try:
        for seed in HISTORIES:
            record = old['interventions'][f'seed-{seed}_sham']
            x, ref, masses = [np.array(record[k]) for k in ('initial_chemistry', 'reference_chemistry', 'masses')]
            ids = record['ids']; targets = old['selections'][str(seed)]['selected_ids']
            original_job = next(j for j in old['jobs'] if (j['seed'], j['level'], j['recipient'], j['target']) == (seed, 'fine', None, None))
            source = Path(original_job['source']); z = moving.payload(source)
            if (not np.array_equal([z['activator'], z['inhibitor']], x) or not np.array_equal(z['ids'], ids)):
                raise ValueError('Fine original physical start disagrees with reviewed preparation')
            inputs.append(source)
            for j in [j for j in jobs if j['seed'] == seed]:
                target = targets[0] if j['recipient'] is None else j['recipient']
                changed, accounting = surrounding_state(x, ref, masses, ids, target, j['arm'])
                if j['arm'] == 'sham':
                    accounting['recipient'] = None
                context_key = f'{seed}_{j["background"]}'
                preparations[context_key] = dict(source=str(source), ids=ids, masses=masses.tolist(),
                    initial_chemistry=x.tolist(), reference_chemistry=ref.tolist(), chemistry=changed.tolist(), accounting=accounting)
                if j['arm'] in ('sham', 'reset'):
                    prior = next(a for a in old['jobs'] if (a['seed'], a['level'], a['recipient'], a['target']) ==
                        (seed, j['level'], j['recipient'], None))
                    result, _ = context.current_evidence(parent, prior, old)
                    with tempfile.TemporaryDirectory(dir=root) as temporary:
                        host = moving.prepared_checkpoint(source, Path(temporary)/'expected.npz', j['dt'], changed)
                        context.source_equivalence(Path(temporary)/'expected.npz', prior['source'])
                    j['source'] = prior['source']
                    pair_folder = parent/'pairs'/f'seed-{seed}'
                    paths = [pair_folder/f'{label}_{prior["key"]}{suffix}' for label in ('pilot', 'full') for suffix in ('.npz', '.json')]
                    if prior['key'] in old['reused']:
                        nested = old['reused'][prior['key']]
                        paths += [Path(f) for f in nested['evidence_sha256']]
                    else:
                        folder = parent/prior['key']; paths += [folder/f for f in ('history.json', 'result.json', 'latest_state.npz', 'prefix/result.json')]
                        paths += [Path(f) for f in read(folder/'prefix/result.json')['evidence_sha256']]
                    paths += [Path(prior['source'])]
                    evidence = {str(f.resolve()):digest(f) for f in paths}
                    inputs += paths
                    reused[j['key']] = dict(job=prior, protocol_sha256=digest(parent/'protocol.json'), evidence_sha256=evidence,
                        pilot_phi=str(pair_folder/f'pilot_{prior["key"]}.npz'), full_phi=str(pair_folder/f'full_{prior["key"]}.npz'))
                else:
                    host = moving.prepared_checkpoint(source, j['source'], j['dt'], changed)
                    inputs.append(Path(j['source']))
                    prior = next(a for a in old['jobs'] if (a['seed'], a['level'], a['recipient'], a['target']) == (seed, j['level'], j['recipient'], None))
                    result, _ = context.current_evidence(parent, prior, old)
                    estimate += result['audit']['wall_seconds']
                key = f'{seed}_{j["level"]}'
                if key in accepted and accepted[key] != asdict(host.config):
                    raise ValueError('Factorial arms have different configurations')
                accepted[key] = asdict(host.config)
            print('Prepared distribution/dosage factorial:', seed, flush=True)
        files = [Path(__file__), Path('tests/test_phase_carry_context_dosage.py')]
        p['source_sha256'] = {**last['source_sha256'], **{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256'] = {**last['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}}
        p['estimated_moving_seconds'] = estimate
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', completed=0, total=42, new_jobs=24, reused_jobs=18,
            estimated_moving_seconds=estimate, independent_histories=3, new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error))); raise


def current_evidence(root, job, p):
    if job['key'] not in p['reused']:
        return moving.current_evidence(root, job, p)
    item = p['reused'][job['key']]; parent = Path(p['parent'])
    if digest(parent/'protocol.json') != item['protocol_sha256']:
        raise ValueError('Changed reused parent protocol')
    context.source_equivalence(job['source'], item['job']['source'])
    result, history = context.current_evidence(parent, item['job'], read(parent/'protocol.json'))
    moving.validate_history(history, job, p, job['duration'])
    for path, sha in item['evidence_sha256'].items():
        if digest(path) != sha:
            raise ValueError('Changed reused evidence')
    return result, history


def temporal_rms(values, times):
    values, times = np.asarray(values), np.asarray(times)
    if len(times) < 2 or np.any(np.diff(times) <= 0) or values.shape != (len(times), 2):
        raise ValueError('Aligned two-species waveforms require at least two increasing times')
    squares = np.mean(values**2, axis=1)
    integral = np.sum(.5*(squares[1:]+squares[:-1])*np.diff(times))
    return float(np.sqrt(integral/(times[-1]-times[0])))


def state_summary(histories, targets, p, horizon):
    times = np.arange(round(horizon/p['interval'])+1)*p['interval']
    late = times >= max(0., horizon-p['late_window'])-1e-9
    rows = []
    for target in targets:
        log_states = {}
        for arm in ARMS:
            history = histories[arm, None if arm == 'sham' else target]
            i = history[0]['ids'].index(target)
            if len(history) != len(times):
                raise ValueError('Unaligned factorial histories')
            log_states[arm] = np.log(np.asarray([r['chemistry'] for r in history])[:, :, i])
        differences = {key:log_states[a]-log_states[b] for key, (a, b) in CONTRASTS.items()}
        # Signed interaction is formed before taking a magnitude.
        differences['interaction'] = log_states['reset']-log_states['bulk']-log_states['redistributed']+log_states['sham']
        contrasts = {}
        for key, wave in differences.items():
            value = temporal_rms(wave[late], times[late])
            contrasts[key] = dict(late_log_rms=value, effect_detected=bool(value>p['effect_min']),
                full_log_rms=temporal_rms(wave, times), final_log_difference=wave[-1].tolist(),
                final_species_ratios=np.exp(wave[-1]).tolist())
        rows.append(dict(cell=target, contrasts=contrasts))
    return rows


def state_comparison(coarse, fine, targets, p, horizon):
    a, b = [state_summary(h, targets, p, horizon) for h in (coarse, fine)]
    error = 0.; same = True
    for ca, cb in zip(a, b):
        for key in ca['contrasts']:
            x, y = ca['contrasts'][key], cb['contrasts'][key]
            same &= x['effect_detected'] == y['effect_detected']
            error = max(error, abs(x['late_log_rms']-y['late_log_rms']), abs(x['full_log_rms']-y['full_log_rms']))
    return dict(passed=bool(same and error<=p['state_metric_error_max']), same_effect_decisions=bool(same),
        metric_error_max=error, coarse=a, fine=b)


def pilot_quality(root, p, seed):
    """Block long work if a pilot already fails the physical audit."""
    root = Path(root); file = root/'pairs'/f'seed-{seed}'/'pilot-quality.json'
    file.parent.mkdir(parents=True, exist_ok=True)
    ph = digest(root/'protocol.json')
    if file.exists():
        record = read(file)
        if record['protocol_sha256'] != ph or not record['passed']:
            raise ValueError('Changed or failed pilot physical gate')
        return record
    rows = []
    for j in [j for j in p['jobs'] if j['seed'] == seed]:
        if j['key'] in p['reused']:
            audit = current_evidence(root, j, p)[0]['audit']
            evidence = p['reused'][j['key']]['evidence_sha256']
        else:
            cp = root/j['key']/'latest_state.npz'; state = moving.payload(cp)
            saved = json.loads(str(state['long_experiment'])); meta = json.loads(str(state['metadata']))
            if (saved['protocol_hash'] != ph or saved['job'] != j or
                    abs(meta['time']-j['start']-p['pilot_duration']) > 1e-9 or
                    abs(saved['history'][-1]['elapsed']-p['pilot_duration']) > 1e-9):
                raise ValueError('Physical pilot checkpoint does not match its horizon')
            audit = saved['audit']; evidence = {str(cp.resolve()):digest(cp)}
        rows.append(dict(job=j['key'], audit=audit, passed=bool(moving.quality_pass(audit, p)), evidence_sha256=evidence))
    record = dict(protocol_sha256=ph, seed=seed, passed=bool(rows and all(r['passed'] for r in rows)), rows=rows,
        scope='New pilot audits and reused full-window audits; passing full audits also bound their pilot interval.')
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
        quality = pilot_quality(root, p, seed)
        paths.append(folder/'pilot-quality.json')
    for level in ('coarse', 'fine'):
        histories[level] = {}
        for j in [j for j in jobs if j['level'] == level]:
            reused = j['key'] in p['reused']
            raw = current_evidence(root, j, p)[1] if reused else read(root/j['key']/'history.json')
            history = [r for r in raw if r['elapsed'] <= horizon+1e-9]
            moving.validate_history(history, j, p, horizon)
            key = (j['arm'], j['recipient']); histories[level][key] = history
            hf = folder/f'{label}_{j["key"]}.json'
            if reused:
                pf = Path(p['reused'][j['key']][label+'_phi'])
            else:
                pf = folder/f'{label}_{j["key"]}.npz'
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
    mechanical = [dict(arm=key[0], recipient=key[1], **moving.mechanical_comparison(
        histories['coarse'][key], histories['fine'][key], fields['coarse', key], fields['fine', key], p, horizon))
        for key in histories['coarse']]
    states = state_comparison(histories['coarse'], histories['fine'], p['selections'][str(seed)]['selected_ids'], p, horizon)
    result = dict(seed=seed, horizon=horizon, protocol_sha256=ph,
        passed=bool(all(r['passed'] for r in mechanical) and states['passed']), mechanical=mechanical, states=states,
        evidence_sha256={str(f.resolve()):digest(f) for f in paths})
    if horizon == p['pilot_duration']:
        result['physical_quality_passed'] = quality['passed']
    if old and old != result:
        raise ValueError('Changed paired decision')
    if not old:
        write_json(file, result)
    print('Distribution/dosage paired gate:', seed, horizon, result['passed'], flush=True)
    return result


def preflight(root):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    moving.require_device(p['device']); torch.set_num_threads(1)
    jobs = [j for j in p['jobs'] if j['seed']==8 and j['recipient']==19 and j['arm'] in ('redistributed', 'bulk')]
    records = []
    for j in jobs:
        write_json(root/'status.json', dict(state='preflighting', current_job=j['key'], completed=0, total=42,
            prefix_contexts_completed=len(records), prefix_contexts_total=len(jobs)))
        records.append(moving.context_gate(root, j, p))
    result = dict(passed=bool(all(r['passed'] for r in records)), protocol_sha256=digest(root/'protocol.json'),
        full_size_contexts=len(records), contexts=records,
        evidence_sha256={str((root/j['key']/'prefix/result.json').resolve()):digest(root/j['key']/'prefix/result.json') for j in jobs})
    write_json(root/'preflight.json', result)
    write_json(root/'status.json', dict(state='prepared_gpu_preflight_passed', completed=0, total=42,
        new_jobs=24, reused_jobs=18, independent_histories=3, new_histories=0))
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
        completed=len(results), total=42, new_jobs=24, reused_jobs=18, independent_histories=3,
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
                for horizon in (p['pilot_duration'], p['duration']):
                    for j in [j for j in p['jobs'] if j['seed'] == seed]:
                        write_json(root/'status.json', dict(state='running', stage='pilot' if horizon==p['pilot_duration'] else 'full',
                            current_job=j['key'], target_horizon=horizon, completed=done, total=42, new_jobs=24, reused_jobs=18,
                            independent_histories=3, new_histories=0))
                        if j['key'] in p['reused']:
                            current_evidence(root, j, p)
                        else:
                            if shutil.disk_usage(root).free < 2*2**30:
                                raise OSError('Less than 2 GiB storage reserve; committed checkpoints retained')
                            delayed.advance(root, j, p, horizon)
                        if horizon == p['duration']:
                            done += 1
                    if not pair_report(root, p, seed, horizon)['passed']:
                        failed.append(dict(seed=seed, reason='prespecified_timestep_gate_failed')); break
                if failed:
                    break
            result = assess(root, p)
            write_json(root/'status.json', dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'], completed=done, total=42, new_jobs=24, reused_jobs=18,
                failed_histories=failed, independent_histories=3, new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=42, error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'preflight', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-context-dosage'))
    args = parser.parse_args()
    result = dict(prepare=prepare, preflight=preflight, run=run, assess=assess)[args.action](args.output)
    print({k:v for k,v in result.items() if k in ('passed', 'completed', 'total', 'estimated_moving_seconds')})
