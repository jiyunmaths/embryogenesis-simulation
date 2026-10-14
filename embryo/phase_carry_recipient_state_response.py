"""Match recipient chemistry while retaining distinct surrounding networks.

All scientific stepping remains in the accepted frozen helpers. This module
prepares local resets and compares their own-control responses; no identity
labels, ongoing clamp, mechanics or parent trajectory reuse.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from multiprocessing import get_context
import os
from pathlib import Path
import time

import numpy as np
import torch

from . import phase_carry_frozen_response as frozen
from . import phase_carry_shared_context_response as shared
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

HISTORIES = frozen.HISTORIES
LEVELS = frozen.LEVELS
DOSES = ('fractional',)
CRITERIA = frozen.CRITERIA
write_json = frozen.write_json
save_npz = frozen.save_npz


def match_recipient(data, target, pair):
    delta, volumes, initial, ids = frozen.validate_context(data['delta'], data['volumes'], data['initial'], data['ids'])
    pair = np.asarray(pair)
    if pair.dtype != np.float64 or pair.shape != (2,) or not np.isfinite(pair).all() or np.any(pair <= 0) or target not in ids:
        raise ValueError('A measured positive float64 two-species recipient pair is required')
    i = ids.tolist().index(target); matched = initial.copy(); matched[:, i] = pair
    before, after = initial*volumes, matched*volumes
    accounting = dict(pre_matching_per_cell_species_amounts=before.tolist(),
        post_matching_per_cell_species_amounts=after.tolist(),
        matching_per_cell_species_amount_change=(after-before).tolist(),
        pre_matching_global_species_amounts=before.sum(axis=1).tolist(),
        post_matching_global_species_amounts=after.sum(axis=1).tolist(),
        matching_global_species_amount_change=(after-before).sum(axis=1).tolist())
    return dict(delta=delta.copy(), volumes=volumes.copy(), initial=matched, ids=ids.copy()), accounting


def design(parent, saved, root):
    """Deduplicate unchanged native controls, retaining each recipient pulse."""
    root = Path(root); contexts = {}; arrays = {}; factorials = []
    for previous in parent['factorials']:
        seed, cell = previous['seed'], previous['cell']
        for reference in ('S', 'R'):
            donor = saved[previous['contexts'][reference+'_'+reference]]
            index = donor['ids'].tolist().index(cell); pair = donor['initial'][:, index]
            combination = {}
            for name, parent_key in previous['contexts'].items():
                source = saved[parent_key]
                data, accounting = match_recipient(source, cell, pair)
                native = np.array_equal(data['initial'], source['initial'])
                key = parent_key + ('_native' if native else f'_cell-{cell}_matched-{reference}')
                if key not in contexts:
                    origin = next(c for c in parent['contexts'] if c['key'] == parent_key)
                    contexts[key] = dict(key=key, seed=seed, background=key,
                        parent_context=parent_key, chemical_context=origin['chemical_context'],
                        physical_context=origin['physical_context'], matching_reference=reference,
                        matching_targets=[], targets=[], source=str(root/'contexts'/(key+'.npz')),
                        matching_accounting=accounting, native_preparation=native,
                        donor_amount_accounting=parent.get('assignments', {}).get(parent_key))
                    arrays[key] = data
                elif not all(np.array_equal(arrays[key][k], data[k]) for k in data):
                    raise ValueError('Context deduplication changed an initial-value problem')
                if cell not in contexts[key]['targets']:
                    contexts[key]['targets'].append(cell)
                    contexts[key]['matching_targets'].append(dict(cell=cell, reference=reference,
                        original_pair=source['initial'][:, index].tolist(), matched_pair=pair.tolist()))
                combination[name] = key
            factorials.append(dict(seed=seed, cell=cell, matching_reference=reference, contexts=combination))
    contexts = list(contexts.values()); trials = []
    for context in contexts:
        for target in (None, *context['targets']):
            suffix = 'control' if target is None else f'cell-{target}_fractional'
            trials.append(dict(key=context['key']+'_'+suffix, context=context['key'],
                seed=context['seed'], background=context['background'], target=target,
                dose=None if target is None else 'fractional', matching_reference=context['matching_reference']))
    return contexts, trials, factorials, arrays


def prepare(root, parent=Path('output/studies/phase-carry-shared-context-response'),
            review_file=Path('docs/phase_carry_shared_context_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_recipient_state_response.md')):
    root, parent = Path(root).resolve(), Path(parent).absolute()
    if root.exists():
        raise FileExistsError(root)
    old, audit, inputs = shared.accepted_parent(parent, review_file)
    if len(old['contexts']) != 21 or len(old['factorials']) != 6 or old['histories'] != list(HISTORIES):
        raise ValueError('The accepted reciprocal shared-context parent is required')
    saved = {}
    for context in old['contexts']:
        file = Path(context['source'])
        if audit['inputs_sha256'].get(str(file)) != digest(file) or old['input_sha256'].get(str(file.resolve())) != digest(file):
            raise ValueError('Changed actual t510 shared-context input')
        with np.load(file, allow_pickle=False) as z:
            saved[context['key']] = {k:z[k].copy() for k in ('delta', 'volumes', 'initial', 'ids')}
        inputs.append(file)
    contexts, trials, factorials, arrays = design(old, saved, root)
    if (len(contexts), len(trials), len(factorials)) != (45, 93, 12):
        raise ValueError('Expected45 backgrounds,93 starts and12 reference-specific factorials')
    root.mkdir(parents=True)
    try:
        for context in contexts:
            data = arrays[context['key']]; save_npz(context['source'], **data)
            inputs.append(Path(context['source']))
            for trial in [t for t in trials if t['context'] == context['key']]:
                initial = data['initial']; changed = initial.copy()
                record = dict(factor=1., initial_log_amplitude=0., requested_removed_amount=0.,
                    removed_activator_amount=0., initial_target_activator_amount=None,
                    initial_target_volume=None, pulse_amount_relative_error=0.)
                if trial['target'] is not None:
                    i = data['ids'].tolist().index(trial['target'])
                    changed, record = shared.shared_pulse(initial, data['volumes'], data['ids'],
                        trial['target'], 'fractional', [initial[0,i], initial[0,i]])
                trial.update(record); trial['source'] = str(root/'starts'/(trial['key']+'.npz'))
                save_npz(trial['source'], initial=changed); inputs.append(Path(trial['source']))
                if context['native_preparation']:
                    prior = next(t for t in old['trials'] if (t['context'],t['target'],t['dose']) ==
                        (context['parent_context'],trial['target'],trial['dose']))
                    file = Path(prior['source'])
                    if audit['inputs_sha256'].get(str(file)) != digest(file):
                        raise ValueError('Unaccepted unchanged comparison input')
                    with np.load(file, allow_pickle=False) as z:
                        if not np.array_equal(changed,z['initial']):
                            raise ValueError('Unchanged arm differs from its accepted parent')
                    trial['native_parent_trial'] = prior['key']; inputs.append(file)
        p = dict(parent=str(parent), parent_protocol_sha256=digest(parent/'protocol.json'),
            contexts=contexts, trials=trials, factorials=factorials, histories=list(HISTORIES),
            selections=old['selections'], levels=old['levels'], beta=old['beta'], da=old['da'], db=old['db'],
            start=510., duration=60., interval=.15, late_window=24., checkpoint_interval=3., effect_min=.01,
            criteria=dict(CRITERIA), reference_settings=old['reference_settings'], runtime=old['runtime'], device=old['device'],
            primary_paths=186, diagnostic_paths=93, total_production_paths=279, cpu_verification_paths=279,
            logical_reference_cases=93, reference_trajectories=186, independent_histories=3, new_histories=0,
            parent_paths_reused=0, pulse_design='fractional; common amount is identical after local matching',
            scope='Initially match both recipient species using two measured reference pairs on both shared physical contexts. Preserve all other concentrations. Immediate pulse/reset with own freely evolving controls; surrounding amounts remain unmatched. Three mature histories, no identity autonomy, neighbor specificity or new development.')
        p['unique_initial_value_problems'] = frozen.unique_problem_count(p)
        if p['unique_initial_value_problems'] != 93:
            raise ValueError('Unexpected exact initial-value count')
        _, matrix, _, _ = frozen.batch_layout(p)
        for level in p['levels']:
            count = frozen.substeps(matrix,level['dt'],p['beta'],p['da'],p['db'])
            for context in contexts:
                if count != frozen.substeps(arrays[context['key']]['delta'],level['dt'],p['beta'],p['da'],p['db']):
                    raise ValueError('Batch changes constituent positivity substeps')
        p['source_sha256'] = {**old['source_sha256'], **frozen.pin([Path(__file__),Path('tests/test_phase_carry_recipient_state_response.py')])}
        inputs += [parent/'production'/level/f for level in ('coarse','fine') for f in ('paths.npz','result.json')]
        for file in inputs[-4:]:
            if audit['inputs_sha256'].get(str(file)) != digest(file):
                raise ValueError('Changed parent comparison paths')
        p['input_sha256'] = frozen.pin([*inputs,Path(protocol_doc)])
        write_json(root/'protocol.json',p)
        write_json(root/'status.json',dict(state='prepared',primary_paths=186,diagnostic_paths=93,
            logical_reference_cases=93,unique_initial_value_problems=93))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def summarize(paths, label, p, times):
    metrics, waves = {}, {}
    for trial in p['trials']:
        if trial['target'] is None:
            continue
        control = next(t for t in p['trials'] if t['context'] == trial['context'] and t['target'] is None)
        _, _, volumes, _, ids = frozen.load_context(p, trial); i = ids.tolist().index(trial['target'])
        values, baseline = paths[label, trial['key']], paths[label, control['key']]
        waves[trial['key']] = np.log(values[:, :, i]/baseline[:, :, i])/trial['initial_log_amplitude']
        metrics[trial['key']] = frozen.response_metrics(values, baseline, times, volumes, i, trial['factor'])
    effects, contrast_waves, native_diagonal = [], {}, []
    for item in p['factorials']:
        for dose in DOSES:
            four = {}
            for name, context in item['contexts'].items():
                trial = next(t for t in p['trials'] if (t['context'], t['target'], t['dose']) == (context, item['cell'], dose))
                four[name] = waves[trial['key']]
            native_diagonal.append(dict(seed=item['seed'], cell=item['cell'], dose=dose, matching_reference=item['matching_reference'],
                normalized_response_shift_rms=frozen.waveform_distance(four['R_R'], four['S_S'], times),
                scope='Descriptive diagonal physical-context contrast after recipient matching; no scientific screen.'))
            contrasts = [('chemical', 'sham_physics', four['R_S']-four['S_S']),
                ('chemical', 'redistributed_physics', four['R_R']-four['S_R']),
                ('physical', 'sham_chemistry', four['S_R']-four['S_S']),
                ('physical', 'redistributed_chemistry', four['R_R']-four['R_S']),
                ('interaction', 'difference_of_differences', four['R_R']-four['S_R']-four['R_S']+four['S_S'])]
            for kind, condition, wave in contrasts:
                key = f'{item["seed"]}_{item["cell"]}_{item["matching_reference"]}_{dose}_{kind}_{condition}'
                value = frozen.waveform_distance(wave, np.zeros_like(wave), times)
                row = dict(key=key, seed=item['seed'], cell=item['cell'], dose=dose, matching_reference=item['matching_reference'],
                    kind=kind, condition=condition, normalized_response_shift_rms=value)
                # Only chemical contrasts have a prespecified scientific screen.
                if kind == 'chemical':
                    row['response_shift_detected'] = bool(value > p['effect_min'])
                effects.append(row); contrast_waves[key] = wave
    return dict(metrics=metrics, effects=effects, native_diagonal=native_diagonal), contrast_waves


def compare(paths, first, second, p, times, summaries, waves, reference=False):
    return shared.compare(paths,first,second,p,times,summaries,waves,reference)


def reference_qualification(root, p):
    paths, times, errors, evidence = frozen.reference_path_sets(root, p)
    unpacked = {label:summarize(paths, label, p, times) for label in ('dop853', 'radau')}
    summaries = {label:item[0] for label, item in unpacked.items()}; waves = {label:item[1] for label, item in unpacked.items()}
    responses = compare(paths, 'dop853', 'radau', p, times, summaries, waves, reference=True)
    result = dict(protocol_sha256=digest(Path(root)/'protocol.json'),
        passed=bool(all(v <= p['criteria']['reference_log_max'] for v in errors.values()) and responses['passed']),
        raw_errors=errors, responses=responses, evidence_sha256=evidence)
    write_json(Path(root)/'reference-qualification.json', result); return result


def baseline_descriptions(paths, p, times):
    controls = []; paired = []; late = times >= p['duration']-p['late_window']
    def network_rms(difference, observations, volumes):
        # Equal species weights; capacity weights over all cells.
        square = np.sum(np.mean(difference**2, axis=1)*volumes[None, :], axis=1)/volumes.sum()
        return float(np.sqrt(np.trapezoid(square, observations)/(observations[-1]-observations[0])))
    for context in p['contexts']:
        trial = next(t for t in p['trials'] if t['context'] == context['key'] and t['target'] is None)
        values = paths['fine', trial['key']]; _, _, volumes, _, ids = frozen.load_context(p, trial)
        contrast = np.std(np.log(values[:, 0]), axis=1)
        drift = np.log(values/values[0][None])
        drift_curve = np.sqrt(np.sum(np.mean(drift**2, axis=1)*volumes[None, :], axis=1)/volumes.sum())
        controls.append(dict(context=context['key'], final_chemistry=values[-1].tolist(), ids=ids.tolist(),
            final_global_species_amounts=(values[-1]@volumes).tolist(),
            final_per_cell_species_ratios_to_reassigned_start=(values[-1]/values[0]).tolist(),
            final_network_two_species_log_drift_rms=float(drift_curve[-1]),
            maximum_network_two_species_log_drift_rms=float(drift_curve.max()),
            late_network_two_species_log_drift_rms=network_rms(drift[late], times[late], volumes),
            final_log_activator_sd=float(contrast[-1]), whole_min_log_activator_sd=float(contrast.min()),
            late_min_log_activator_sd=float(contrast[late].min())))
    for item in p['factorials']:
        for physics in ('S', 'R'):
            trials = [next(t for t in p['trials'] if t['context'] == item['contexts'][chem+'_'+physics] and t['target'] is None) for chem in ('S', 'R')]
            a, b = [paths['fine', t['key']] for t in trials]
            _, _, volumes, _, ids = frozen.load_context(p, trials[0]); i = ids.tolist().index(item['cell'])
            difference = np.log(b[:, :, i]/a[:, :, i])
            paired.append(dict(seed=item['seed'], cell=item['cell'], matching_reference=item['matching_reference'], physics='sham' if physics == 'S' else 'redistributed',
                initial_species_ratios=(b[0, :, i]/a[0, :, i]).tolist(), final_species_ratios=(b[-1, :, i]/a[-1, :, i]).tolist(),
                late_two_species_state_log_rms=frozen.waveform_distance(difference[late], np.zeros_like(difference[late]), times[late]),
                late_network_two_species_state_log_rms=network_rms(np.log(b[late]/a[late]), times[late], volumes)))
    return dict(controls=controls, paired_chemical_states=paired)


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    paths, times, reference_errors, evidence = frozen.path_sets(root, p)
    labels = ['dop853', 'radau', *[l['level'] for l in p['levels']]]
    unpacked = {label:summarize(paths, label, p, times) for label in labels}
    summaries = {label:item[0] for label, item in unpacked.items()}; waves = {label:item[1] for label, item in unpacked.items()}
    comparisons = [compare(paths, 'dop853', 'radau', p, times, summaries, waves, reference=True)]
    comparisons += [compare(paths, a, b, p, times, summaries, waves) for a, b in (('coarse', 'fine'), ('fine', 'diagnostic'))]
    comparisons += [compare(paths, l['level'], method, p, times, summaries, waves) for l in p['levels'] for method in ('dop853', 'radau')]
    accuracy, raw_pairs, orders = [], [], []
    for trial in p['trials']:
        errors = []
        for level in p['levels']:
            name = level['level']
            error = max(float(abs(np.log(paths[name, trial['key']]/paths[ref, trial['key']])).max()) for ref in ('dop853', 'radau'))
            cpu = float(abs(np.log(paths[name, trial['key']]/paths['cpu-'+name, trial['key']])).max())
            accuracy.append(dict(trial=trial['key'], level=name, raw_reference_log_error=error,
                cpu_torch_log_error=cpu, passed=bool(error <= p['criteria']['raw_log_max'] and cpu <= p['criteria']['cpu_torch_log_max'])))
            errors.append(error)
        for a, b in (('coarse', 'fine'), ('fine', 'diagnostic')):
            error = float(abs(np.log(paths[a, trial['key']]/paths[b, trial['key']])).max())
            raw_pairs.append(dict(trial=trial['key'], first=a, second=b, raw_log_error=error,
                passed=bool(error <= p['criteria']['raw_log_max'])))
        noise = reference_errors[trial['key']]; resolved = errors[-1] > p['criteria']['reference_noise_factor']*noise
        values = np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved and all(x > 0 for x in errors) else None
        orders.append(dict(trial=trial['key'], errors=errors, reference_noise=noise, order_resolved=bool(resolved), orders=values,
            passed=bool(not resolved or (values is not None and all(p['criteria']['order_min'] <= x <= p['criteria']['order_max'] for x in values)))))
    native_checks = []
    for level in ('coarse', 'fine'):
        file = Path(p['parent'])/'production'/level/'paths.npz'
        with np.load(file) as z:
            old_paths = {str(k):v for k, v in zip(z['trial_keys'], z['trajectory'].transpose(1, 0, 2, 3))}
            if not np.array_equal(z['times'], times):
                raise ValueError('Native comparison clock differs')
            for trial in p['trials']:
                if 'native_parent_trial' in trial:
                    error = float(abs(np.log(paths[level, trial['key']]/old_paths[trial['native_parent_trial']])).max())
                    native_checks.append(dict(trial=trial['key'], level=level, raw_log_error=error,
                        passed=bool(error <= p['criteria']['cpu_torch_log_max'])))
    outcomes = []
    for seed in p['histories']:
        for reference in ('S','R'):
            for physics in ('sham_physics','redistributed_physics'):
                effects = [r for r in summaries['fine']['effects'] if (r['seed'],r['matching_reference'],r['kind'],r['condition']) == (seed,reference,'chemical',physics)]
                if len(effects) != 2:
                    raise ValueError('Missing recipient comparison')
                outcomes.append(dict(seed=seed,matching_reference=reference,physics=physics,recipients=2,
                    above_screen=sum(r['response_shift_detected'] for r in effects),
                    cells=[dict(cell=r['cell'],response_shift_detected=r['response_shift_detected']) for r in effects]))
    reciprocal_outcomes = []
    for seed in p['histories']:
        cells = []
        for cell in p['selections'][str(seed)]['selected_ids']:
            effects = [r for r in summaries['fine']['effects'] if (r['seed'],r['cell'],r['kind']) == (seed,cell,'chemical')]
            if len(effects) != 4:
                raise ValueError('Missing matching-reference/physical-context comparison')
            cells.append(dict(cell=cell,all_four_above_screen=bool(all(r['response_shift_detected'] for r in effects))))
        reciprocal_outcomes.append(dict(seed=seed,recipients=2,cells=cells,all_four_above_screen=sum(r['all_four_above_screen'] for r in cells)))
    result = dict(protocol_sha256=digest(root/'protocol.json'),
        passed=bool(all(e <= p['criteria']['reference_log_max'] for e in reference_errors.values())
            and all(r['passed'] for r in comparisons+accuracy+raw_pairs+orders+native_checks)),
        primary_paths=p['primary_paths'], diagnostic_paths=p['diagnostic_paths'], total_production_paths=p['total_production_paths'], cpu_verification_paths=p['cpu_verification_paths'],
        logical_reference_cases=p['logical_reference_cases'], unique_initial_value_problems=p['unique_initial_value_problems'], reference_trajectories=p['reference_trajectories'],
        independent_histories=3, new_histories=0, references_raw_errors=reference_errors,
        production_accuracy=accuracy, raw_timestep_comparisons=raw_pairs, orders=orders,
        response_comparisons=comparisons, response_summaries=summaries,
        baseline_descriptions=baseline_descriptions(paths, p, times), history_outcomes=outcomes,
        reciprocal_history_outcomes=reciprocal_outcomes,
        native_parent_comparisons=native_checks, evidence_sha256=evidence, scope=p['scope'])
    write_json(root/'summary.json', result); return result


def run(root, workers=4):
    root = Path(root).resolve()
    if not 1 <= workers <= 4:
        raise ValueError('Use one to four single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
            os.environ[name]='1'
        started=time.perf_counter()
        try:
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                write_json(root/'status.json',dict(state='running',stage='independent_references',completed=0,total=len(p['trials'])))
                pending=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,future in enumerate(as_completed(pending),1):
                    if not future.result()['passed']:
                        raise RuntimeError('Raw dual-reference gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=len(p['trials'])))
                if not reference_qualification(root,p)['passed']:
                    raise RuntimeError('Normalized reference gate failed; production blocked')
                # CPU verification and GPU production have disjoint output paths.
                pending=[pool.submit(frozen.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for level in p['levels']:
                    done=[f for f in pending if f.done()]
                    for future in done:future.result()
                    write_json(root/'status.json',dict(state='running',stage='production_and_cpu_verification',
                        level=level['level'],cpu_completed=len(done),cpu_total=len(pending),total_production_paths=p['total_production_paths']))
                    frozen.production(root,p,level)
                for i,future in enumerate(as_completed(pending),1):
                    future.result()
                    write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=len(pending)))
            result=assess(root)
            write_json(root/'timing.json',dict(coordinator_wall_seconds=time.perf_counter()-started,
                cpu_verification_overlaps_gpu=True,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'],primary_paths=p['primary_paths'],diagnostic_paths=p['diagnostic_paths'],
                total_production_paths=p['total_production_paths'],independent_histories=3,new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/phase-carry-recipient-state-response'))
    parser.add_argument('--parent',type=Path,default=Path('output/studies/phase-carry-shared-context-response'))
    parser.add_argument('--review',type=Path,default=Path('docs/phase_carry_shared_context_response_assessment.json'))
    parser.add_argument('--protocol-doc',type=Path,default=Path('docs/phase_carry_recipient_state_response.md'))
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    result=prepare(args.output,args.parent,args.review,args.protocol_doc) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('passed','primary_paths','diagnostic_paths','total_production_paths')})
