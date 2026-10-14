"""Qualify delayed responses after amount-preserving context redistribution.

These CPU tests distinguish response behavior from background concentration,
require exact carry reuse, and keep independent units at developmental histories.
"""
import json

import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.optimize import brentq

from embryo import phase_carry_redistribution_response as study
from embryo import phase_carry_exchange_response as moving
from embryo.feedback_long import digest
from embryo.resolution import write_json
from test_phase_carry_exchange_response import source_fixture


def selections():
    return {str(seed): dict(selected_ids=ids) for seed, ids in zip(
        study.HISTORIES, ([22, 23], [19, 20], [29, 30]))}


def test_design_has_66_unique_paths_and_background_specific_controls(tmp_path):
    contexts, jobs = study.job_design(tmp_path, selections())
    assert study.HISTORIES == (7, 8, 9)
    assert len(contexts) == 9 and len(jobs) == len({j['key'] for j in jobs}) == 66
    assert sum(j['target'] is None for j in jobs) == 18
    assert sum(j['background'] == 'sham' for j in jobs) == 30
    assert sum(j['dose'] == 'amount' for j in jobs) == 24
    assert sum(j['dose'] == 'fractional' for j in jobs) == 24
    assert all(j['start'] == 510. and j['duration'] == 60. for j in jobs)
    for seed in study.HISTORIES:
        selected = selections()[str(seed)]['selected_ids']
        for level, _ in study.LEVELS:
            group = [j for j in jobs if j['seed'] == seed and j['level'] == level]
            assert len(group) == 11
            for name in ('sham', *[f'redistributed-cell-{t}' for t in selected]):
                background = [j for j in group if j['background'] == name]
                assert sum(j['target'] is None for j in background) == 1
                targets = selected if name == 'sham' else [int(name.split('-')[-1])]
                assert {(j['target'], j['dose']) for j in background if j['target'] is not None} == {
                    (t, dose) for t in targets for dose in study.DOSES}


def test_matched_actual_amount_is_feasible_with_unequal_endpoint_volumes():
    original = np.array([[2., .4], [.6, .3]])
    redistributed = np.array([[.05, .7], [.8, .4]])
    original_masses, redistributed_masses = [1., 3.], [2., 4.]
    amount = study.matched_removal(original[0, 0]*original_masses[0],
        redistributed[0, 0]*redistributed_masses[0])
    assert amount == pytest.approx(.01)
    factors = []
    for x, masses in ((original, original_masses), (redistributed, redistributed_masses)):
        changed, record = study.pulse(x, masses, [11, 17], 11, 'amount', amount)
        assert record['removed_activator_amount'] == pytest.approx(amount, rel=2e-13)
        assert record['initial_log_amplitude'] == pytest.approx(abs(np.log(record['factor'])))
        np.testing.assert_array_equal(changed[:, 1], x[:, 1])
        np.testing.assert_array_equal(changed[1], x[1])
        assert np.all(changed > 0)
        factors.append(record['factor'])
    assert factors == pytest.approx([.995, .9])
    _, fractional = study.pulse(original, original_masses, [11, 17], 11, 'fractional')
    assert fractional['removed_activator_amount'] == pytest.approx(.2)
    assert fractional['removed_activator_amount'] != pytest.approx(amount)


@pytest.mark.parametrize('amount', [None, 0., -1., 1., 2., np.nan, np.inf])
def test_amount_pulse_cannot_clip_empty_or_use_a_nonfinite_dose(amount):
    with pytest.raises(ValueError):
        study.pulse(np.ones((2, 2)), [1., 1.], [11, 17], 11, 'amount', amount)


def synthetic_histories(redistributed_tau=4.):
    choices = {str(seed): dict(selected_ids=[11, 17]) for seed in study.HISTORIES}
    _, jobs = study.job_design('/tmp/redistribution-response-fixture', choices)
    jobs = [j for j in jobs if j['seed'] == 7]
    times = np.arange(401)*.15
    ids, initial = [11, 17], np.array([[.1, 2.], [.5, 3.]])
    histories = {}
    for name in ('sham', 'redistributed-cell-11', 'redistributed-cell-17'):
        start = initial.copy()
        if name != 'sham':
            start[0, ids.index(int(name.split('-')[-1]))] *= .05
        # Each control drifts independently, so using another background's
        # absolute chemistry must not masquerade as response behavior.
        drift = .002 if name == 'sham' else -.003
        control = np.repeat(start[None], len(times), axis=0)*np.exp(drift*times[:, None, None])
        histories[name, None, None] = [dict(elapsed=float(t), time=510.+float(t), ids=ids,
            chemistry=x.tolist(), volumes=[1., 2.], log_activator_sd=float(np.std(np.log(x[0]))))
            for t, x in zip(times, control)]
        targets = ids if name == 'sham' else [int(name.split('-')[-1])]
        for target in targets:
            for dose in study.DOSES:
                factor = .995 if name == 'sham' and dose == 'amount' else .9
                i = ids.index(target)
                tau = 1. if name == 'sham' else redistributed_tau
                path = control.copy()
                path[:, 0, i] *= np.exp(np.log(factor)*np.exp(-times/tau))
                path[0, 0, i] = control[0, 0, i]*factor
                histories[name, target, dose] = [dict(elapsed=float(t), time=510.+float(t), ids=ids,
                    chemistry=x.tolist(), volumes=[1., 2.], log_activator_sd=float(np.std(np.log(x[0]))))
                    for t, x in zip(times, path)]
                for job in jobs:
                    if (job['background'], job['target'], job['dose']) == (name, target, dose):
                        job.update(factor=factor, initial_log_amplitude=float(abs(np.log(factor))),
                            removed_activator_amount=float((1-factor)*start[0, i]*[1., 2.][i]))
        for job in jobs:
            if job['background'] == name and job['target'] is None:
                job['factor'] = 1.
    p = dict(interval=.15, late_window=24., contrast_min=.1, effect_min=.01,
        response_limits=moving.RESPONSE_CRITERIA)
    return histories, jobs, p


def test_behavioral_shift_uses_each_own_control_and_actual_log_pulse_amplitude():
    histories, jobs, p = synthetic_histories()
    summary = study.response_summary(histories, [j for j in jobs if j['level'] == 'fine'], [11, 17], p, 60.)
    assert all(row['shifted_at_both_doses'] for row in summary)
    for row in summary:
        assert set(row['contexts']) == {'sham', f'redistributed-cell-{row["cell"]}'}
        for context in row['contexts'].values():
            assert context['dose_dependence_normalized_rms'] < 1e-13
            assert all(dose['metrics']['target_recovery_time'] is not None for dose in context['doses'].values())
        assert row['differences']['fractional']['normalized_response_shift_rms'] == pytest.approx(
            row['differences']['amount']['normalized_response_shift_rms'], abs=1e-13)
        assert row['differences']['fractional']['raw_log_response_shift_rms'] != pytest.approx(
            row['differences']['amount']['raw_log_response_shift_rms'])
    assert study.response_comparison(histories, histories, jobs, [11, 17], p, 60.)['passed']


def test_different_levels_and_drifting_controls_do_not_imply_different_behavior():
    histories, jobs, p = synthetic_histories(redistributed_tau=1.)
    summary = study.response_summary(histories, [j for j in jobs if j['level'] == 'fine'], [11, 17], p, 60.)
    assert all(not row['shifted_at_both_doses'] for row in summary)
    assert all(dose['normalized_response_shift_rms'] < 1e-13
        for row in summary for dose in row['differences'].values())
    # The target baselines still differ by twenty-fold at time zero.
    assert histories['sham', None, None][0]['chemistry'][0][0]/histories[
        'redistributed-cell-11', None, None][0]['chemistry'][0][0] == pytest.approx(20.)


def test_small_pulse_raw_agreement_cannot_mask_a_normalized_response_error():
    coarse, jobs, p = synthetic_histories()
    fine = {key: json.loads(json.dumps(rows)) for key, rows in coarse.items()}
    fine['sham', 11, 'amount'][5]['chemistry'][0][0] *= np.exp(.0002)
    report = study.response_comparison(coarse, fine, jobs, [11, 17], p, 60.)
    assert not report['passed']
    failed = next(row for row in report['rows'] if
        (row['context'], row['cell'], row['dose']) == ('sham', 11, 'amount'))
    assert failed['normalized_response_error'] > .01


def test_refinement_requires_the_same_response_effect_decision_near_the_screen():
    times = np.arange(401)*.15
    def difference(tau):
        squared = (np.exp(-times)-np.exp(-times/tau))**2/2
        return np.sqrt(trapezoid(squared, times)/60.)-.01
    threshold_tau = brentq(difference, 1., 2.)
    coarse, jobs, p = synthetic_histories(redistributed_tau=threshold_tau*(1-1e-4))
    fine, _, _ = synthetic_histories(redistributed_tau=threshold_tau*(1+1e-4))
    report = study.response_comparison(coarse, fine, jobs, [11, 17], p, 60.)
    assert all(row['passed'] for row in report['rows'])
    assert report['response_shift_error_max'] < .01
    assert not report['same_response_effect_decisions'] and not report['passed']


def test_wrong_saved_pulse_factor_is_rejected_instead_of_renormalized():
    histories, jobs, p = synthetic_histories()
    job = next(j for j in jobs if (j['level'], j['background'], j['target'], j['dose']) ==
        ('fine', 'sham', 11, 'amount'))
    job['factor'] = .99
    with pytest.raises(ValueError, match='Pulse factor'):
        study.response_summary(histories, [j for j in jobs if j['level'] == 'fine'], [11, 17], p, 60.)


@pytest.mark.parametrize('changed_key', ['phase_carry', 'activator', 'metadata', 'rounding'])
def test_reuse_rejects_any_change_in_full_initial_state(tmp_path, changed_key):
    host, source = source_fixture(tmp_path, start=510.)
    copied = tmp_path/'copied.npz'
    moving.prepared_checkpoint(source, copied, host.config.dt)
    study.source_equivalence(source, copied)
    payload = moving.payload(copied)
    if changed_key == 'metadata':
        meta = json.loads(str(payload['metadata']))
        meta['rng']['state']['state'] += 1
        payload['metadata'] = np.array(json.dumps(meta))
    elif changed_key == 'rounding':
        payload[changed_key][0] += 1
    else:
        payload[changed_key].flat[0] = np.nextafter(payload[changed_key].flat[0], np.inf)
    np.savez_compressed(copied, **payload)
    with pytest.raises(ValueError, match=changed_key):
        study.source_equivalence(source, copied)


def test_delayed_pulse_and_retiming_preserve_carry_geometry_and_rng(tmp_path):
    host, source = source_fixture(tmp_path, start=510.)
    before = moving.payload(source)
    chemistry = np.array([host.activator, host.inhibitor])
    values, _ = study.pulse(chemistry, host.volumes(), host.ids, int(host.ids[0]), 'fractional')
    copied = tmp_path/'pulsed.npz'
    new = moving.prepared_checkpoint(source, copied, host.config.dt/2, values)
    after = moving.payload(copied)
    for key in before:
        if key not in ('metadata', 'activator', 'inhibitor', 'long_experiment'):
            assert before[key].dtype == after[key].dtype
            np.testing.assert_array_equal(before[key], after[key])
    first, second = [json.loads(str(state['metadata'])) for state in (before, after)]
    for key in ('rng', 'fate_rng', 'signal_rng', 'time', 'lineage', 'divisions'):
        assert first[key] == second[key]
    assert second['time'] == 510. and new.step_number == 2*host.step_number
    np.testing.assert_array_equal([new.activator, new.inhibitor], values)


def accepted_parent(tmp_path, review_style):
    parent = tmp_path/'parent'
    parent.mkdir()
    old = dict(source_sha256={}, input_sha256={})
    write_json(parent/'protocol.json', old)
    write_json(parent/'summary.json', dict(passed=True))
    write_json(parent/'status.json', dict(state='completed', passed=True))
    full_review = tmp_path/'full-assessment.json'
    review = dict(full_study_accepted=True,
        protocol_sha256=digest(parent/'protocol.json'),
        original_summary_sha256=digest(parent/'summary.json'),
        original_status_sha256=digest(parent/'status.json'))
    write_json(full_review, {**review, 'evidence_sha256': {}})
    review[review_style] = str(full_review)
    review[review_style+'_sha256'] = digest(full_review)
    review_file = tmp_path/'assessment.json'
    write_json(review_file, review)
    return parent, review_file, full_review, old


@pytest.mark.parametrize('review_style', ['full_assessment', 'full_review'])
def test_accepted_parent_requires_the_completed_original_and_full_review_hash(tmp_path, review_style):
    parent, review_file, full_review, old = accepted_parent(tmp_path, review_style)
    restored, inputs = study.reviewed_parent(parent, review_file)
    assert restored == old
    assert {parent/'protocol.json', parent/'summary.json', parent/'status.json', review_file, full_review}.issubset(
        set(inputs))


@pytest.mark.parametrize('changed', ['protocol.json', 'summary.json', 'status.json', 'full-assessment.json'])
def test_stale_parent_review_cannot_authorize_new_scientific_runs(tmp_path, changed):
    parent, review_file, full_review, _ = accepted_parent(tmp_path, 'full_assessment')
    path = full_review if changed == 'full-assessment.json' else parent/changed
    document = json.loads(path.read_text())
    document['changed_after_review'] = True
    write_json(path, document)
    with pytest.raises(ValueError):
        study.reviewed_parent(parent, review_file)


def test_an_unaccepted_parent_does_not_become_accepted_by_matching_hashes(tmp_path):
    parent, review_file, _, _ = accepted_parent(tmp_path, 'full_assessment')
    review = json.loads(review_file.read_text())
    review['full_study_accepted'] = False
    write_json(review_file, review)
    with pytest.raises(ValueError):
        study.reviewed_parent(parent, review_file)


def test_matching_full_review_hash_does_not_override_a_disagreeing_assessment(tmp_path):
    parent, review_file, full_review, _ = accepted_parent(tmp_path, 'full_assessment')
    full = json.loads(full_review.read_text())
    full['full_study_accepted'] = False
    write_json(full_review, full)
    review = json.loads(review_file.read_text())
    review['full_assessment_sha256'] = digest(full_review)
    write_json(review_file, review)
    with pytest.raises(ValueError, match='disagrees'):
        study.reviewed_parent(parent, review_file)


def test_review_acceptance_rechecks_its_underlying_immutable_evidence(tmp_path):
    parent, review_file, full_review, _ = accepted_parent(tmp_path, 'full_assessment')
    evidence_file = tmp_path/'reviewed-evidence.json'
    write_json(evidence_file, {'value': 1})
    full = json.loads(full_review.read_text())
    full['evidence_sha256'] = {str(evidence_file): digest(evidence_file)}
    write_json(full_review, full)
    review = json.loads(review_file.read_text())
    review['full_assessment_sha256'] = digest(full_review)
    write_json(review_file, review)
    _, inputs = study.reviewed_parent(parent, review_file)
    assert evidence_file in inputs
    write_json(evidence_file, {'value': 2})
    with pytest.raises(ValueError, match='assessed evidence'):
        study.reviewed_parent(parent, review_file)


def test_parent_endpoint_resolves_nested_reuse_to_the_actual_checkpoint(tmp_path):
    newest, parent, original = [tmp_path/name for name in ('newest', 'parent', 'original')]
    for root in (newest, parent, original):
        root.mkdir()
    write_json(original/'protocol.json', {'reused': {}})
    write_json(parent/'protocol.json', dict(parent=str(original),
        reused={'prior': {'job': {'key': 'physical'}}}))
    p = dict(parent=str(parent), reused={'current': {'job': {'key': 'prior'}}})
    assert study.parent_endpoint(newest, p, {'key': 'current'}) == original/'physical'/'latest_state.npz'


def pilot_checkpoint(tmp_path, **audit_changes):
    host, source = source_fixture(tmp_path, start=510.)
    seed = host.config.seed
    job = dict(key='pilot', seed=seed, start=510.)
    p = dict(jobs=[job], reused={}, pilot_duration=6., amount_error_max=2e-14)
    write_json(tmp_path/'protocol.json', p)
    payload = moving.payload(source)
    meta = json.loads(str(payload['metadata']))
    meta.update(time=516., step_number=round(516./host.config.dt))
    payload['metadata'] = np.array(json.dumps(meta))
    audit = dict(max_volume_error=.01, min_radius=5., max_clipping=0., boundary_max=0., dilution_error_max=0.)
    audit.update(audit_changes)
    payload['long_experiment'] = np.array(json.dumps(dict(protocol_hash=digest(tmp_path/'protocol.json'),
        job=job, audit=audit, history=[dict(elapsed=6.)])))
    folder = tmp_path/'pilot'
    folder.mkdir()
    np.savez_compressed(folder/'latest_state.npz', **payload)
    return p, seed


def test_qualified_physical_pilot_allows_long_continuation(tmp_path):
    p, seed = pilot_checkpoint(tmp_path)
    assert study.pilot_quality(tmp_path, p, seed)['passed']


@pytest.mark.parametrize('bad_audit', [dict(max_volume_error=.06), dict(min_radius=3.9),
    dict(max_clipping=1.), dict(boundary_max=.02), dict(dilution_error_max=3e-14)])
def test_physical_pilot_failure_blocks_long_work(tmp_path, bad_audit):
    p, seed = pilot_checkpoint(tmp_path, **bad_audit)
    with pytest.raises(RuntimeError, match='long continuation blocked'):
        study.pilot_quality(tmp_path, p, seed)
    record = json.loads((tmp_path/'pairs'/f'seed-{seed}'/'pilot-quality.json').read_text())
    assert not record['passed']


def test_paired_pilot_failure_never_advances_to_long_duration(tmp_path, monkeypatch):
    _, jobs = study.job_design(tmp_path, selections())
    p = dict(jobs=jobs, reused={}, pilot_duration=6., duration=60., total_jobs=66,
        new_jobs=48, reused_jobs=18)
    calls = []
    monkeypatch.setattr(study, 'advance', lambda root, job, protocol, horizon: calls.append(horizon))
    monkeypatch.setattr(study, 'pair_report', lambda *args: dict(passed=False))
    monkeypatch.setattr(study, 'pilot_quality', lambda *args: dict(passed=True))
    done, passed = study.schedule_history(tmp_path, p, 7, 0)
    assert not passed and done == 0
    assert calls == [6.]*22


def test_physical_failure_also_blocks_long_duration_when_the_pair_passes(tmp_path, monkeypatch):
    _, jobs = study.job_design(tmp_path, selections())
    p = dict(jobs=jobs, reused={}, pilot_duration=6., duration=60., total_jobs=66,
        new_jobs=48, reused_jobs=18)
    calls = []
    monkeypatch.setattr(study, 'advance', lambda root, job, protocol, horizon: calls.append(horizon))
    monkeypatch.setattr(study, 'pair_report', lambda *args: dict(passed=True))
    def physical_failure(*args):
        raise RuntimeError('Pilot physical quality failed; long continuation blocked')
    monkeypatch.setattr(study, 'pilot_quality', physical_failure)
    with pytest.raises(RuntimeError, match='long continuation blocked'):
        study.schedule_history(tmp_path, p, 7, 0)
    assert calls == [6.]*22
