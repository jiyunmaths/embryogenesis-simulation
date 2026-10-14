"""Qualify reciprocal transfer, dose invariance and signed factorial contrasts."""
from pathlib import Path

import numpy as np
import pytest
import torch

from embryo import phase_carry_shared_context_response as study


def fixture(tmp_path, duration=.6):
    """Small nonidentical chemical/physical contexts with the full logical design."""
    ids = np.array([11, 17]); selections = {str(s):dict(selected_ids=ids.tolist()) for s in study.HISTORIES}
    saved = {}
    for seed in study.HISTORIES:
        for recipient in (None, *ids):
            name = 'sham' if recipient is None else f'redistributed-cell-{recipient}'
            initial = np.array([[.8, 1.1], [.9, 1.2]], dtype=np.float64)
            initial *= 1.+.03*(seed-7)
            volumes = np.array([1., 2.])
            if recipient is not None:
                index = ids.tolist().index(recipient)
                initial[0, index] *= .8 if index == 0 else 1.2
                initial[1, 1-index] *= 1.05
                volumes *= np.array([1.3, .7]) if index == 0 else np.array([.8, 1.1])
            conductance = .1+.02*(seed-7)+(0. if recipient is None else .015*(index+1))
            delta = np.array([[-conductance, conductance], [conductance, -conductance]])/volumes[:, None]
            saved[f'seed-{seed}_{name}'] = dict(delta=delta, volumes=volumes, initial=initial, ids=ids)
    contexts, trials, factorials = study.design(selections, tmp_path)
    assignments = {}
    for context in contexts:
        data, assignments[context['key']] = study.transferred_context(saved[context['chemical_context']], saved[context['physical_context']])
        study.save_npz(context['source'], **data)
    for trial in trials:
        context = next(c for c in contexts if c['key'] == trial['context'])
        with np.load(context['source']) as z:
            data = {k:z[k].copy() for k in z.files}
        changed = data['initial'].copy(); factor = 1.
        if trial['target'] is not None:
            index = ids.tolist().index(trial['target'])
            pair = [saved[f'seed-{trial["seed"]}_{name}']['initial'][0, index]
                for name in ('sham', f'redistributed-cell-{trial["target"]}')]
            changed, record = study.shared_pulse(changed, data['volumes'], ids, trial['target'], trial['dose'], pair)
            factor = record['factor']; trial.update(record)
        trial.update(source=str(tmp_path/'starts'/(trial['key']+'.npz')), factor=factor,
            initial_log_amplitude=float(abs(np.log(factor))))
        study.save_npz(trial['source'], initial=changed)
    parent = tmp_path/'parent'
    times = np.arange(round(duration/.15)+1)*.15
    for level in ('coarse', 'fine'):
        study.save_npz(parent/'production'/level/'paths.npz', times=times,
            trial_keys=np.array([], dtype=str), trajectory=np.empty((len(times), 0, 2, 2)))
    p = dict(contexts=contexts, trials=trials, factorials=factorials, assignments=assignments,
        histories=list(study.HISTORIES), selections=selections,
        levels=[dict(level=l, dt=dt) for l, dt in study.LEVELS], beta=2., da=.02, db=.55,
        start=510., duration=duration, interval=.15, checkpoint_interval=.3, effect_min=.01,
        late_window=24., criteria=dict(study.CRITERIA), device='cpu', parent=str(parent),
        reference_settings=dict(rtol=1e-13, atol=1e-15, max_step=.15),
        runtime=dict(torch=str(torch.__version__), torch_cuda=torch.version.cuda, expected_gpu_name='unused'),
        source_sha256={}, input_sha256=study.frozen.pin([*[c['source'] for c in contexts], *[t['source'] for t in trials]]),
        scope='CPU numerical qualification fixture')
    p['unique_initial_value_problems'] = study.frozen.unique_problem_count(p)
    study.write_json(tmp_path/'protocol.json', p)
    return p


def synthetic_paths(p, times, label='fine'):
    paths = {}
    for context in p['contexts']:
        control = next(t for t in p['trials'] if t['context'] == context['key'] and t['target'] is None)
        _, _, _, initial, ids = study.frozen.load_context(p, control)
        drift = np.exp(times[:, None, None]*np.array([[.002, -.003], [-.001, .002]])[None])
        baseline = initial[None]*drift
        paths[label, control['key']] = baseline
        chemistry = .05 if 'redistributed' in context['chemical_context'] else 0.
        physics = .07 if 'redistributed' in context['physical_context'] else 0.
        interaction = .09 if chemistry and physics else 0.
        for trial in [t for t in p['trials'] if t['context'] == context['key'] and t['target'] is not None]:
            i = ids.tolist().index(trial['target'])
            response = np.column_stack((-np.exp(-times/2), .1*times*np.exp(-times/2)))
            response += (chemistry+physics+interaction)*times[:, None]*np.exp(-times[:, None]/4)
            values = baseline.copy(); values[:, :, i] *= np.exp(response*trial['initial_log_amplitude'])
            paths[label, trial['key']] = values
    return paths


def test_full_design_has_21_contexts_69_cases_and_57_exact_ivps(tmp_path):
    p = fixture(tmp_path)
    assert len(p['contexts']) == 21 and len(p['factorials']) == 6
    assert len(p['trials']) == len({t['key'] for t in p['trials']}) == 69
    assert sum(t['target'] is None for t in p['trials']) == 21
    assert sum(t['dose'] == 'fractional' for t in p['trials']) == 24
    assert sum(t['dose'] == 'amount' for t in p['trials']) == 24
    assert p['unique_initial_value_problems'] == 57
    for seed in study.HISTORIES:
        pair = [f for f in p['factorials'] if f['seed'] == seed]
        assert pair[0]['contexts']['S_S'] == pair[1]['contexts']['S_S']
        assert all(pair[0]['contexts'][c] != pair[1]['contexts'][c] for c in ('R_S', 'S_R', 'R_R'))


def test_batch_mapping_matches_trial_order_and_preserves_each_physical_block(tmp_path):
    p = fixture(tmp_path)
    initial, delta, mapping, count = study.frozen.batch_layout(p)
    assert count == 2 and initial.shape == (2, 42, 5)
    assert [m['trial'] for m in mapping] == [t['key'] for t in p['trials']]
    assert np.all(initial > 0)
    for index, context in enumerate(p['contexts']):
        with np.load(context['source']) as z:
            np.testing.assert_array_equal(delta[2*index:2*index+2, 2*index:2*index+2], z['delta'])
    for level in p['levels']:
        global_count = study.frozen.substeps(delta, level['dt'], p['beta'], p['da'], p['db'])
        assert global_count == 1
        for context in p['contexts']:
            with np.load(context['source']) as z:
                assert study.frozen.substeps(z['delta'], level['dt'], p['beta'], p['da'], p['db']) == global_count


def test_transfer_preserves_concentrations_and_exact_operator_capacity_pair():
    chemical = dict(ids=np.array([11, 17]), initial=np.array([[.2, 3.], [.8, 2.]]),
        volumes=np.array([1., 2.]))
    physical = dict(ids=chemical['ids'].copy(), volumes=np.array([4., 3.]),
        delta=np.array([[-.3, .3], [.4, -.4]]))
    result, accounting = study.transferred_context(chemical, physical)
    for key in ('delta', 'volumes'):
        np.testing.assert_array_equal(result[key], physical[key])
    np.testing.assert_array_equal(result['initial'], chemical['initial'])
    np.testing.assert_array_equal(accounting['volume_ratios'], [4., 1.5])
    old = chemical['initial']*chemical['volumes']; new = chemical['initial']*physical['volumes']
    np.testing.assert_array_equal(accounting['added_per_cell_species_amounts'], new-old)
    np.testing.assert_array_equal(accounting['assigned_global_species_amounts'], new.sum(axis=1))
    result['initial'][0, 0] = 9.
    assert chemical['initial'][0, 0] == .2
    physical['ids'] = physical['ids'][::-1]
    with pytest.raises(ValueError, match='identical cell IDs'):
        study.transferred_context(chemical, physical)


def test_invalid_operator_volume_pair_is_rejected(tmp_path):
    p = fixture(tmp_path); context = p['contexts'][0]
    with np.load(context['source']) as z:
        values = {k:z[k].copy() for k in z.files}
    changed = dict(values); changed['volumes'] = np.array([3., 2.])
    with pytest.raises(ValueError):
        study.transferred_context(values, changed)


def test_canonical_common_factors_are_identical_across_physics_and_amounts_rescale():
    initial = np.array([[2., .4], [.8, 1.3]])
    a, first = study.shared_pulse(initial, np.array([1., 2.]), [11, 17], 11, 'amount', [.3, 2.])
    b, second = study.shared_pulse(initial, np.array([7., 2.]), [11, 17], 11, 'amount', [.3, 2.])
    np.testing.assert_array_equal(a, b)
    assert first['factor'].hex() == second['factor'].hex()
    assert second['requested_removed_amount'] == first['requested_removed_amount']*7.
    assert second['removed_activator_amount'] == first['removed_activator_amount']*7.
    expected = initial.copy(); expected[0, 0] *= first['factor']
    np.testing.assert_array_equal(a, expected)
    assert first['pulse_amount_relative_error'] < 2e-13


def test_lower_common_and_fractional_preparations_are_identical():
    initial = np.array([[.3, 2.], [.8, 1.3]])
    a, first = study.shared_pulse(initial, [1., 2.], [11, 17], 11, 'amount', [.3, 2.])
    b, second = study.shared_pulse(initial, [1., 2.], [11, 17], 11, 'fractional', [.3, 2.])
    assert first['factor'] == second['factor'] == .9
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize('change', ['negative_volume', 'duplicate_id', 'wrong_pair', 'float32', 'zero_species'])
def test_invalid_pulse_preparation_cannot_be_silently_clipped(change):
    initial = np.array([[.3, 2.], [.8, 1.3]]); volumes = np.array([1., 2.]); ids = np.array([11, 17]); pair = [.3, 2.]
    if change == 'negative_volume': volumes[0] = -1.
    if change == 'duplicate_id': ids[:] = 11
    if change == 'wrong_pair': pair = [.3]
    if change == 'float32': initial = initial.astype(np.float32)
    if change == 'zero_species': initial[1, 0] = 0.
    with pytest.raises((ValueError, FloatingPointError)):
        study.shared_pulse(initial, volumes, ids, 11, 'amount', pair)


def test_signed_factorial_algebra_and_own_drifting_controls(tmp_path):
    p = fixture(tmp_path, duration=60.); times = np.arange(401)*.15
    paths = synthetic_paths(p, times)
    summary, waves = study.summarize(paths, 'fine', p, times)
    assert len(summary['metrics']) == 48 and len(summary['effects']) == 60
    assert len(summary['native_diagonal']) == 12
    assert sum('response_shift_detected' in e for e in summary['effects']) == 24
    item = p['factorials'][0]; cell = item['cell']; four = {}
    for name, context in item['contexts'].items():
        pulse = next(t for t in p['trials'] if (t['context'], t['target'], t['dose']) == (context, cell, 'fractional'))
        control = next(t for t in p['trials'] if t['context'] == context and t['target'] is None)
        _, _, _, _, ids = study.frozen.load_context(p, pulse); index = ids.tolist().index(cell)
        four[name] = np.log(paths['fine', pulse['key']][:, :, index]/paths['fine', control['key']][:, :, index])/pulse['initial_log_amplitude']
    key = f'{item["seed"]}_{cell}_fractional_'
    chemical_s = waves[key+'chemical_sham_physics']; chemical_r = waves[key+'chemical_redistributed_physics']
    physical_s = waves[key+'physical_sham_chemistry']; physical_r = waves[key+'physical_redistributed_chemistry']
    interaction = waves[key+'interaction_difference_of_differences']
    np.testing.assert_allclose(chemical_s+physical_r, four['R_R']-four['S_S'], atol=1e-15)
    np.testing.assert_allclose(chemical_r+physical_s, four['R_R']-four['S_S'], atol=1e-15)
    np.testing.assert_allclose(interaction, chemical_r-chemical_s, atol=1e-15)
    np.testing.assert_allclose(interaction, physical_r-physical_s, atol=1e-15)
    expected = study.frozen.waveform_distance(interaction, np.zeros_like(interaction), times)
    recorded = next(e for e in summary['effects'] if e['key'] == key+'interaction_difference_of_differences')
    assert recorded['normalized_response_shift_rms'] == expected
    assert expected > 0


def test_reference_contrast_gate_detects_four_wave_cancellation_error(tmp_path):
    p = fixture(tmp_path, duration=60.); times = np.arange(401)*.15
    paths = synthetic_paths(p, times, 'dop853')
    paths.update({('radau', key):values.copy() for (label, key), values in paths.items()})
    item = p['factorials'][0]; envelopes = 1.-np.exp(-times)
    for name, sign in (('R_R', 1.), ('S_R', -1.), ('R_S', -1.), ('S_S', 1.)):
        trial = next(t for t in p['trials'] if (t['context'], t['target'], t['dose']) == (item['contexts'][name], item['cell'], 'fractional'))
        _, _, _, _, ids = study.frozen.load_context(p, trial); index = ids.tolist().index(item['cell'])
        paths['radau', trial['key']][:, :, index] *= np.exp(sign*6e-9*envelopes[:, None]*trial['initial_log_amplitude'])
    values = {label:study.summarize(paths, label, p, times) for label in ('dop853', 'radau')}
    report = study.compare(paths, 'dop853', 'radau', p, times, {k:v[0] for k,v in values.items()}, {k:v[1] for k,v in values.items()}, reference=True)
    assert all(r['passed'] for r in report['pulse_rows'])
    assert not report['passed']
    interaction = next(r for r in report['factorial_contrasts'] if r['key'] == f'{item["seed"]}_{item["cell"]}_fractional_interaction_difference_of_differences')
    assert interaction['direct_waveform_error'] > 1e-8
    assert interaction['rms_limit'] == 1e-8


def test_primary_screen_is_strict_and_classification_disagreement_fails(tmp_path):
    p = fixture(tmp_path, duration=60.); times = np.arange(401)*.15
    paths = synthetic_paths(p, times, 'coarse')
    first, _ = study.summarize(paths, 'coarse', p, times)
    target = next(e for e in first['effects'] if e['kind'] == 'chemical')
    p['effect_min'] = target['normalized_response_shift_rms']
    paths.update({('fine', key):values.copy() for (label, key), values in paths.items()})
    values = {label:study.summarize(paths, label, p, times) for label in ('coarse', 'fine')}
    row = next(e for e in values['coarse'][0]['effects'] if e['key'] == target['key'])
    assert row['normalized_response_shift_rms'] == p['effect_min'] and not row['response_shift_detected']
    other = next(e for e in values['fine'][0]['effects'] if e['key'] == target['key'])
    other['response_shift_detected'] = True
    report = study.compare(paths, 'coarse', 'fine', p, times, {k:v[0] for k,v in values.items()}, {k:v[1] for k,v in values.items()})
    assert not report['passed']
    failed = next(r for r in report['factorial_contrasts'] if r['key'] == target['key'])
    assert failed['direct_waveform_error'] == 0. and failed['same_scientific_decision'] is False


def test_baseline_drift_and_network_separation_use_capacity_weights(tmp_path):
    p = fixture(tmp_path, duration=60.); times = np.arange(401)*.15
    paths = synthetic_paths(p, times)
    descriptions = study.baseline_descriptions(paths, p, times)
    assert len(descriptions['controls']) == 21 and len(descriptions['paired_chemical_states']) == 12
    first = descriptions['controls'][0]
    control = p['trials'][0]; path = paths['fine', control['key']]
    _, _, volumes, _, _ = study.frozen.load_context(p, control)
    change = np.log(path[-1]/path[0])
    expected = np.sqrt(np.sum(np.mean(change**2, axis=0)*volumes)/volumes.sum())
    assert first['final_network_two_species_log_drift_rms'] == pytest.approx(expected)
    assert first['late_network_two_species_log_drift_rms'] > 0
    assert all(r['late_network_two_species_state_log_rms'] > 0 for r in descriptions['paired_chemical_states'])


def test_changed_transferred_context_fails_before_chemical_integration(tmp_path):
    p = fixture(tmp_path); file = Path(p['contexts'][0]['source'])
    with np.load(file) as z: arrays = {k:z[k].copy() for k in z.files}
    arrays['volumes'] *= 2.; study.save_npz(file, **arrays)
    with pytest.raises(ValueError, match='preparation'):
        study.frozen.load_context(p, p['trials'][0])


def test_full_cpu_factorial_pipeline_checks_every_case_and_all_contrasts(tmp_path):
    p = fixture(tmp_path)
    for trial in p['trials']:
        assert study.frozen.reference_worker((tmp_path, trial))['passed']
    assert study.reference_qualification(tmp_path, p)['passed']
    for level in p['levels']:
        study.frozen.production(tmp_path, p, level)
        for trial in p['trials']:
            study.frozen.cpu_worker((tmp_path, trial, level))
    result = study.assess(tmp_path)
    assert result['passed']
    assert len(result['production_accuracy']) == 207 and len(result['orders']) == 69
    assert len(result['raw_timestep_comparisons']) == 138
    assert len(result['response_comparisons']) == 9
    assert all(len(r['pulse_rows']) == 48 and len(r['factorial_contrasts']) == 60 for r in result['response_comparisons'])
    assert result['unique_initial_value_problems'] == 57
    assert all('response_shift_detected' not in e for e in result['response_summaries']['fine']['effects'] if e['kind'] != 'chemical')
