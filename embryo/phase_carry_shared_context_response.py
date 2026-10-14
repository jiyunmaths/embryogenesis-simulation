"""Reciprocal shared-physics challenges on accepted frozen t=510 states.

Reassign complete chemical concentration vectors, preserving their cell order,
to either measured operator/volume context. Reuse the unchanged frozen chemistry,
reference, CPU verification and checkpoint machinery; never reuse a trajectory
as a new shared-context result.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from multiprocessing import get_context
import os
from pathlib import Path

import numpy as np
import torch

from . import phase_carry_frozen_response as frozen
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

HISTORIES = frozen.HISTORIES
LEVELS = frozen.LEVELS
DOSES = frozen.DOSES
CRITERIA = frozen.CRITERIA
write_json = frozen.write_json
save_npz = frozen.save_npz


def accepted_parent(parent, review_file):
    parent, review_file = Path(parent).absolute(), Path(review_file).absolute()
    p, review = read(parent/'protocol.json'), read(review_file)
    core = {key:parent/file for key, file in (
        ('protocol_sha256', 'protocol.json'), ('original_summary_sha256', 'summary.json'),
        ('original_status_sha256', 'status.json'))}
    if (not review['full_study_accepted'] or read(parent/'status.json')['state'] != 'completed'
            or not read(parent/'summary.json')['passed']):
        raise ValueError('Accepted completed frozen parent required')
    for key, file in core.items():
        if digest(file) != review[key]:
            raise ValueError('Changed accepted frozen parent')
    independent = Path(review['independent_review'])
    if digest(independent) != review['independent_review_sha256']:
        raise ValueError('Changed independent frozen review')
    audit = read(independent)
    if not audit['passed'] or not audit['inputs_unchanged'] or audit['protocol_sha256'] != review['protocol_sha256']:
        raise ValueError('Independent frozen audit did not accept the parent')
    archive = Path(review['review_archive']).absolute()
    provenance_file = archive/'review-provenance.json'; provenance = read(provenance_file)
    if not provenance['passed'] or any(provenance[k] != review[k] for k in core):
        raise ValueError('Frozen review provenance disagrees')
    inputs = [*core.values(), review_file, independent, provenance_file]
    reporting_only = {Path('README.md').resolve(), Path('docs/plan.md').resolve()}
    for filename, expected in provenance['files_sha256'].items():
        file = Path(filename).resolve()
        if file in reporting_only:
            continue
        if digest(file) != expected:
            raise ValueError('Changed scientific frozen review asset: '+filename)
        inputs.append(file)
    for filename, expected in p['source_sha256'].items():
        if digest(filename) != expected:
            raise ValueError('Changed inherited scientific source')
    if (p['beta'], p['da'], p['db']) != (2., .02, .55):
        raise ValueError('Changed accepted chemical parameters')
    return p, audit, inputs


def design(selections, root):
    """Share only the identical S-chemistry/S-physics control across recipients."""
    root = Path(root); contexts = {}; factorials = []
    for seed in HISTORIES:
        for cell in selections[str(seed)]['selected_ids']:
            combinations = {}
            for chemistry, physics in (('S', 'S'), ('R', 'S'), ('S', 'R'), ('R', 'R')):
                chemical_name = 'sham' if chemistry == 'S' else f'redistributed-cell-{cell}'
                physical_name = 'sham' if physics == 'S' else f'redistributed-cell-{cell}'
                key = f'seed-{seed}_chemistry-{chemical_name}_physics-{physical_name}'
                if key not in contexts:
                    contexts[key] = dict(key=key, seed=seed,
                        background=f'chemistry-{chemical_name}_physics-{physical_name}', targets=[],
                        chemical_context=f'seed-{seed}_{chemical_name}',
                        physical_context=f'seed-{seed}_{physical_name}',
                        source=str(root/'contexts'/(key+'.npz')))
                if cell not in contexts[key]['targets']:
                    contexts[key]['targets'].append(cell)
                combinations[chemistry+'_'+physics] = key
            factorials.append(dict(seed=seed, cell=cell, contexts=combinations))
    contexts = list(contexts.values()); trials = frozen.trial_design(contexts)
    return contexts, trials, factorials


def transferred_context(chemical, physical):
    """Preserve concentrations; explicitly account for assignment to new volumes."""
    if not np.array_equal(chemical['ids'], physical['ids']):
        raise ValueError('Chemical and physical contexts must have identical cell IDs/order')
    delta, volumes, initial, ids = frozen.validate_context(physical['delta'], physical['volumes'],
        chemical['initial'], chemical['ids'])
    old_amounts = initial*chemical['volumes'][None, :]
    new_amounts = initial*volumes[None, :]
    return dict(delta=delta.copy(), volumes=volumes.copy(), initial=initial.copy(), ids=ids.copy()), dict(
        concentrations_exact=True, original_per_cell_species_amounts=old_amounts.tolist(),
        original_volumes=chemical['volumes'].tolist(), assigned_volumes=volumes.tolist(),
        volume_ratios=(volumes/chemical['volumes']).tolist(),
        assigned_per_cell_species_amounts=new_amounts.tolist(),
        added_per_cell_species_amounts=(new_amounts-old_amounts).tolist(),
        original_global_species_amounts=old_amounts.sum(axis=1).tolist(),
        assigned_global_species_amounts=new_amounts.sum(axis=1).tolist(),
        added_global_species_amounts=(new_amounts-old_amounts).sum(axis=1).tolist())


def shared_pulse(initial, volumes, ids, target, dose, paired_activators):
    initial, volumes, ids = np.asarray(initial), np.asarray(volumes), np.asarray(ids)
    paired_activators = np.asarray(paired_activators)
    if (initial.dtype != np.float64 or initial.shape != (2, len(volumes))
            or volumes.ndim != 1 or ids.shape != volumes.shape or len(set(ids.tolist())) != len(ids)
            or not np.isfinite(volumes).all() or np.any(volumes <= 0)
            or paired_activators.shape != (2,) or not np.isfinite(paired_activators).all()
            or np.any(paired_activators <= 0) or dose not in DOSES or target not in ids):
        raise ValueError('Invalid shared-context dose')
    frozen.valid_state(initial)
    i = ids.tolist().index(target); before = float(initial[0, i]*volumes[i])
    minimum = float(min(paired_activators))
    factor = .9 if dose == 'fractional' else float(1.-.1*minimum/float(initial[0, i]))
    requested = .1*before if dose == 'fractional' else float(.1*minimum*volumes[i])
    if not 0 < factor < 1 or requested <= 0:
        raise ValueError('Nonzero representable positive depletion required')
    changed = initial.copy(); changed[0, i] *= factor
    measured = float((initial[0, i]-changed[0, i])*volumes[i])
    error = abs(measured/requested-1.)
    if error > CRITERIA['pulse_amount_relative_max']:
        raise ValueError('Shared-context pulse amount accounting failed')
    frozen.valid_state(changed)
    return changed, dict(factor=factor, initial_log_amplitude=float(abs(np.log(factor))),
        requested_removed_amount=requested, removed_activator_amount=measured,
        initial_target_activator_amount=before, initial_target_volume=float(volumes[i]),
        pulse_amount_relative_error=error)


def prepare(root, parent=Path('output/studies/phase-carry-frozen-response'),
            review_file=Path('docs/phase_carry_frozen_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_shared_context_response.md')):
    root, parent = Path(root).resolve(), Path(parent).absolute()
    if root.exists():
        raise FileExistsError(root)
    old, audit, inputs = accepted_parent(parent, review_file)
    saved = {}
    for context in old['contexts']:
        file = Path(context['source'])
        if (audit['inputs_sha256'].get(str(file)) != digest(file)
                or old['input_sha256'].get(str(file.resolve())) != digest(file)):
            raise ValueError('Changed or unaccepted actual t=510 context')
        with np.load(file) as z:
            data = {key:z[key].copy() for key in ('delta', 'volumes', 'initial', 'ids')}
        frozen.validate_context(data['delta'], data['volumes'], data['initial'], data['ids'])
        saved[context['key']] = data; inputs.append(file)
    contexts, trials, factorials = design(old['selections'], root)
    if len(contexts) != 21 or len(trials) != 69 or len(factorials) != 6:
        raise ValueError('Unexpected reciprocal factorial counts')
    root.mkdir(parents=True)
    try:
        assignments = {}
        for context in contexts:
            data, accounting = transferred_context(saved[context['chemical_context']], saved[context['physical_context']])
            save_npz(context['source'], **data)
            inputs.append(Path(context['source'])); assignments[context['key']] = accounting
            native = context['chemical_context'] == context['physical_context']
            for trial in [t for t in trials if t['context'] == context['key']]:
                initial = data['initial']; record = dict(factor=1., initial_log_amplitude=0.,
                    requested_removed_amount=0., removed_activator_amount=0.,
                    initial_target_activator_amount=None, initial_target_volume=None, pulse_amount_relative_error=0.)
                changed = initial.copy()
                if trial['target'] is not None:
                    i = data['ids'].tolist().index(trial['target'])
                    amounts = [saved[f'seed-{trial["seed"]}_{name}']['initial'][0, i]
                        for name in ('sham', f'redistributed-cell-{trial["target"]}')]
                    changed, record = shared_pulse(initial, data['volumes'], data['ids'], trial['target'], trial['dose'], amounts)
                trial.update(record); trial['source'] = str(root/'starts'/(trial['key']+'.npz'))
                save_npz(trial['source'], initial=changed); inputs.append(Path(trial['source']))
                if native:
                    prior = next(t for t in old['trials'] if (t['context'], t['target'], t['dose']) ==
                        (context['chemical_context'], trial['target'], trial['dose']))
                    original = Path(prior['source'])
                    if audit['inputs_sha256'].get(str(original)) != digest(original):
                        raise ValueError('Unaccepted native comparison preparation')
                    with np.load(original) as z:
                        same = np.array_equal(z['initial'], changed)
                    if trial['dose'] in (None, 'fractional') and not same:
                        raise ValueError('Native control/fractional preparation changed')
                    trial['native_parent_preparation_exact'] = same
                    if same:
                        trial['native_parent_trial'] = prior['key']
                    inputs.append(original)
        p = dict(parent=str(parent), parent_protocol_sha256=digest(parent/'protocol.json'),
            contexts=contexts, trials=trials, factorials=factorials, assignments=assignments,
            histories=list(HISTORIES), selections=old['selections'], levels=old['levels'],
            beta=old['beta'], da=old['da'], db=old['db'], start=510., duration=60., interval=.15,
            late_window=24., checkpoint_interval=3., effect_min=.01, criteria=dict(CRITERIA),
            reference_settings=old['reference_settings'], runtime=old['runtime'], device=old['device'],
            primary_paths=138, diagnostic_paths=69, total_production_paths=207,
            cpu_verification_paths=207, logical_reference_cases=69, reference_trajectories=138,
            independent_histories=3, new_histories=0, parent_paths_reused=0,
            scope='Reciprocal fixed t510 whole-chemical-state×measured-operator/volume factorial. Concentrations preserved at immediate reassignment/pulse, amounts recorded and not matched. Chemical state includes levels and distribution. Three inherited histories. No recipient-state matching, equilibration, ongoing mechanics, identity autonomy/inheritance, new development or population inference.')
        p['unique_initial_value_problems'] = frozen.unique_problem_count(p)
        _, matrix, _, _ = frozen.batch_layout(p)
        for level in p['levels']:
            count = frozen.substeps(matrix, level['dt'], p['beta'], p['da'], p['db'])
            for context in contexts:
                with np.load(context['source']) as z:
                    if count != frozen.substeps(z['delta'], level['dt'], p['beta'], p['da'], p['db']):
                        raise ValueError('Batch changes constituent positivity substeps')
        p['source_sha256'] = {**old['source_sha256'], **frozen.pin([Path(__file__), Path('tests/test_phase_carry_shared_context_response.py')])}
        # Consumed old primary trajectories are comparison evidence, never reuse.
        inputs += [parent/'production'/level/f for level in ('coarse', 'fine') for f in ('paths.npz', 'result.json')]
        for file in inputs[-4:]:
            if audit['inputs_sha256'].get(str(file)) != digest(file):
                raise ValueError('Changed accepted native trajectory comparison')
        p['input_sha256'] = frozen.pin([*inputs, Path(protocol_doc)])
        write_json(root/'protocol.json', p)
        write_json(root/'status.json', dict(state='prepared', primary_paths=138, diagnostic_paths=69,
            logical_reference_cases=69, unique_initial_value_problems=p['unique_initial_value_problems']))
        return p
    except Exception as error:
        write_json(root/'status.json', dict(state='preparation_failed', error=str(error))); raise


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
            native_diagonal.append(dict(seed=item['seed'], cell=item['cell'], dose=dose,
                normalized_response_shift_rms=frozen.waveform_distance(four['R_R'], four['S_S'], times),
                scope='Descriptive native diagonal at the current shared-context dose convention; common-amount factors may differ from the parent.'))
            contrasts = [('chemical', 'sham_physics', four['R_S']-four['S_S']),
                ('chemical', 'redistributed_physics', four['R_R']-four['S_R']),
                ('physical', 'sham_chemistry', four['S_R']-four['S_S']),
                ('physical', 'redistributed_chemistry', four['R_R']-four['R_S']),
                ('interaction', 'difference_of_differences', four['R_R']-four['S_R']-four['R_S']+four['S_S'])]
            for kind, condition, wave in contrasts:
                key = f'{item["seed"]}_{item["cell"]}_{dose}_{kind}_{condition}'
                value = frozen.waveform_distance(wave, np.zeros_like(wave), times)
                row = dict(key=key, seed=item['seed'], cell=item['cell'], dose=dose,
                    kind=kind, condition=condition, normalized_response_shift_rms=value)
                # Only chemical contrasts have a prespecified scientific screen.
                if kind == 'chemical':
                    row['response_shift_detected'] = bool(value > p['effect_min'])
                effects.append(row); contrast_waves[key] = wave
    return dict(metrics=metrics, effects=effects, native_diagonal=native_diagonal), contrast_waves


def compare(paths, first, second, p, times, summaries, waves, reference=False):
    # The unchanged helper handles all pulse/own-control response/AUC/recovery
    # rows. Supply chemical effects only because descriptive rows have no labels.
    filtered = {label:dict(metrics=s['metrics'], effects=[r for r in s['effects'] if r['kind'] == 'chemical'])
        for label, s in summaries.items()}
    limit = p['criteria']['reference_normalized_max'] if reference else p['criteria']['normalized_response_max']
    report = frozen.compare_responses(paths, first, second, p, times, filtered, limit)
    rms_limit = p['criteria']['reference_normalized_max'] if reference else p['criteria']['effect_error_max']
    contrasts = []
    for x, y in zip(summaries[first]['effects'], summaries[second]['effects']):
        if x['key'] != y['key']:
            raise ValueError('Changed factorial contrast order')
        waveform_error = float(abs(waves[first][x['key']]-waves[second][y['key']]).max())
        rms_error = abs(x['normalized_response_shift_rms']-y['normalized_response_shift_rms'])
        same = None if x['kind'] != 'chemical' else x['response_shift_detected'] == y['response_shift_detected']
        passed = waveform_error <= limit and rms_error <= rms_limit and same is not False
        contrasts.append(dict(key=x['key'], kind=x['kind'], condition=x['condition'],
            direct_waveform_error=waveform_error, direct_waveform_limit=limit, rms_error=rms_error, rms_limit=rms_limit,
            same_scientific_decision=same, passed=bool(passed)))
    report['factorial_contrasts'] = contrasts
    report['passed'] = bool(report['passed'] and all(r['passed'] for r in contrasts))
    return report


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
            paired.append(dict(seed=item['seed'], cell=item['cell'], physics='sham' if physics == 'S' else 'redistributed',
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
        for physics in ('sham_physics', 'redistributed_physics'):
            effects = [r for r in summaries['fine']['effects'] if r['seed'] == seed and r['kind'] == 'chemical' and r['condition'] == physics]
            cells = [dict(cell=c, shifted_at_both_doses=bool(all(r['response_shift_detected'] for r in effects if r['cell'] == c)))
                for c in p['selections'][str(seed)]['selected_ids']]
            outcomes.append(dict(seed=seed, physics=physics, recipients=2,
                both_doses_above_screen=sum(r['shifted_at_both_doses'] for r in cells), cells=cells,
                dose_counts=[dict(dose=dose, above_screen=sum(r['response_shift_detected'] for r in effects if r['dose'] == dose)) for dose in DOSES]))
    reciprocal_outcomes = []
    for seed in p['histories']:
        cells = []
        for cell in p['selections'][str(seed)]['selected_ids']:
            effects = [r for r in summaries['fine']['effects'] if (r['seed'], r['cell'], r['kind']) == (seed, cell, 'chemical')]
            cells.append(dict(cell=cell, both_doses_both_physical_contexts_above_screen=bool(all(r['response_shift_detected'] for r in effects))))
        reciprocal_outcomes.append(dict(seed=seed, recipients=2, cells=cells,
            both_doses_both_physical_contexts_above_screen=sum(r['both_doses_both_physical_contexts_above_screen'] for r in cells)))
    result = dict(protocol_sha256=digest(root/'protocol.json'),
        passed=bool(all(e <= p['criteria']['reference_log_max'] for e in reference_errors.values())
            and all(r['passed'] for r in comparisons+accuracy+raw_pairs+orders+native_checks)),
        primary_paths=138, diagnostic_paths=69, total_production_paths=207, cpu_verification_paths=207,
        logical_reference_cases=69, unique_initial_value_problems=p['unique_initial_value_problems'], reference_trajectories=138,
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
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); frozen.runtime_check(p); torch.set_num_threads(1)
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ[name] = '1'
        try:
            write_json(root/'status.json', dict(state='running', stage='independent_references', completed=0, total=69))
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                pending = [pool.submit(frozen.reference_worker, (root, t)) for t in p['trials']]
                for i, future in enumerate(as_completed(pending), 1):
                    if not future.result()['passed']:
                        raise RuntimeError('Raw dual-reference gate failed')
                    write_json(root/'status.json', dict(state='running', stage='independent_references', completed=i, total=69))
            if not reference_qualification(root, p)['passed']:
                raise RuntimeError('Normalized factorial reference gate failed; production blocked')
            for level in p['levels']:
                write_json(root/'status.json', dict(state='running', stage='production', level=level['level'], total_production_paths=207))
                frozen.production(root, p, level)
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn')) as pool:
                pending = [pool.submit(frozen.cpu_worker, (root, t, l)) for l in p['levels'] for t in p['trials']]
                for i, future in enumerate(as_completed(pending), 1):
                    future.result()
                    write_json(root/'status.json', dict(state='running', stage='cpu_verification', completed=i, total=207))
            result = assess(root)
            write_json(root/'status.json', dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'], primary_paths=138, diagnostic_paths=69,
                total_production_paths=207, independent_histories=3, new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('output/studies/phase-carry-shared-context-response'))
    parser.add_argument('--parent', type=Path, default=Path('output/studies/phase-carry-frozen-response'))
    parser.add_argument('--review', type=Path, default=Path('docs/phase_carry_frozen_response_assessment.json'))
    parser.add_argument('--protocol-doc', type=Path, default=Path('docs/phase_carry_shared_context_response.md'))
    parser.add_argument('--workers', type=int, default=4); args = parser.parse_args()
    result = prepare(args.output, args.parent, args.review, args.protocol_doc) if args.action == 'prepare' else run(args.output, args.workers) if args.action == 'run' else assess(args.output)
    print({k:v for k, v in result.items() if k in ('passed', 'primary_paths', 'diagnostic_paths', 'total_production_paths')})
