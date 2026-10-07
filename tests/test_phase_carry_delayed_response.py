"""Protect delayed dose normalization, own controls and safe exact restarts."""
from dataclasses import asdict
import json
import os

import numpy as np
import pytest

from embryo import phase_carry_delayed_response as study
from embryo import phase_carry_exchange_response as moving
from embryo.attribute_development import AttributeSimulation
from embryo.feedback_long import digest
from embryo.resolution import write_json
from test_phase_carry_exchange_response import source_fixture


def test_common_amount_is_feasible_for_severely_depleted_background():
    amount = study.matched_removal(2., .01)
    assert amount == .001
    assert study.matched_removal(.01, 2.) == amount
    for invalid in (0., -1., np.nan, np.inf):
        with pytest.raises(ValueError):
            study.matched_removal(invalid, 2.)


def test_equal_amount_across_different_volumes_is_not_equal_fraction():
    a = np.array([[2., .2], [.5, .3]]); b = np.array([[.05, .3], [.7, .4]])
    ma, mb, ids = np.array([1., 3.]), np.array([2., 4.]), [11, 17]
    amount = study.matched_removal(a[0, 0]*ma[0], b[0, 0]*mb[0])
    x, ra = study.pulse(a, ma, ids, 11, 'amount', amount)
    y, rb = study.pulse(b, mb, ids, 11, 'amount', amount)
    np.testing.assert_allclose([ra['removed_activator_amount'], rb['removed_activator_amount']], amount, rtol=2e-13)
    assert ra['factor'] == .995 and rb['factor'] == .9
    np.testing.assert_array_equal(x[1], a[1]); np.testing.assert_array_equal(x[:, 1], a[:, 1])
    np.testing.assert_array_equal(y[1], b[1]); np.testing.assert_array_equal(y[:, 1], b[:, 1])
    assert np.all(x > 0) and np.all(y > 0)
    _, fractional = study.pulse(a, ma, ids, 11, 'fractional', amount)
    assert fractional['factor'] == .9
    assert fractional['removed_activator_amount'] > ra['removed_activator_amount']


@pytest.mark.parametrize('amount', [None, 0., -1., 1., 2., np.nan])
def test_invalid_absolute_pulse_cannot_clip_or_empty_a_cell(amount):
    with pytest.raises(ValueError):
        study.pulse(np.ones((2, 2)), [1., 1.], [11, 17], 11, 'amount', amount)


def test_delayed_retiming_preserves_physical_state_carry_and_random_streams(tmp_path):
    host, source = source_fixture(tmp_path, start=510.)
    original = moving.payload(source); x = np.array([host.activator, host.inhibitor])
    values, record = study.pulse(x, host.volumes(), host.ids, int(host.ids[0]), 'fractional')
    dest = tmp_path/'pulse.npz'; new = moving.prepared_checkpoint(source, dest, host.config.dt/2, values)
    after = moving.payload(dest)
    for key in original:
        if key not in ('metadata', 'activator', 'inhibitor', 'long_experiment'):
            np.testing.assert_array_equal(after[key], original[key])
    a, b = [json.loads(str(z['metadata'])) for z in (original, after)]
    for key in ('rng', 'fate_rng', 'signal_rng', 'time', 'lineage', 'divisions'):
        assert a[key] == b[key]
    assert b['time'] == 510. and new.step_number == host.step_number*2
    assert record['factor'] == .9


def synthetic_histories():
    selections = {str(s):dict(selected_ids=[11, 17]) for s in study.HISTORIES}
    _, jobs = study.job_design('/tmp/delayed-fixture', selections)
    jobs = [j for j in jobs if j['seed'] == 7]
    times = np.arange(401)*.15; ids = [11, 17]; initial = np.array([[.1, 2.], [.5, 3.]])
    result = {}
    for name in ('sham', 'reset-cell-11', 'reset-cell-17'):
        start = initial.copy()
        if name != 'sham':
            start[0, ids.index(int(name.split('-')[-1]))] *= .05
        control = np.repeat(start[None], len(times), axis=0)*np.exp(.002*times[:, None, None])
        result[name, None, None] = [dict(elapsed=float(t), time=510.+float(t), ids=ids,
            chemistry=x.tolist(), volumes=[1., 2.], log_activator_sd=float(np.std(np.log(x[0])))) for t, x in zip(times, control)]
        for target in (ids if name == 'sham' else [int(name.split('-')[-1])]):
            for dose in study.DOSES:
                factor = .995 if name == 'sham' and dose == 'amount' else .9
                i = ids.index(target); path = control.copy(); tau = 1. if name == 'sham' else 4.
                path[:, 0, i] *= np.exp(np.log(factor)*np.exp(-times/tau))
                path[0, 0, i] = control[0, 0, i]*factor
                result[name, target, dose] = [dict(elapsed=float(t), time=510.+float(t), ids=ids,
                    chemistry=x.tolist(), volumes=[1., 2.], log_activator_sd=float(np.std(np.log(x[0])))) for t, x in zip(times, path)]
                for j in jobs:
                    if (j['background'], j['target'], j['dose']) == (name, target, dose):
                        j.update(factor=factor, initial_log_amplitude=float(abs(np.log(factor))),
                            removed_activator_amount=float((1-factor)*start[0, i]*[1., 2.][i]))
        for j in jobs:
            if j['background'] == name and j['target'] is None:
                j['factor'] = 1.
    p = dict(interval=.15, late_window=24., contrast_min=.1, effect_min=.01,
        response_limits=moving.RESPONSE_CRITERIA)
    return result, jobs, p


def test_each_dose_uses_actual_log_amplitude_and_its_own_drifting_control():
    h, jobs, p = synthetic_histories()
    fine = [j for j in jobs if j['level'] == 'fine']
    rows = study.response_summary(h, fine, [11, 17], p, 60.)
    assert all(r['shifted_at_both_doses'] for r in rows)
    for row in rows:
        for c in row['contexts'].values():
            assert c['dose_dependence_normalized_rms'] < 1e-13
            assert all(d['metrics']['target_recovery_time'] is not None for d in c['doses'].values())
        assert abs(row['differences']['fractional']['normalized_response_shift_rms']-
            row['differences']['amount']['normalized_response_shift_rms']) < 1e-13
    assert study.response_comparison(h, h, jobs, [11, 17], p, 60.)['passed']


def test_small_pulse_raw_agreement_cannot_hide_normalized_error():
    coarse, jobs, p = synthetic_histories(); fine = {k:json.loads(json.dumps(v)) for k, v in coarse.items()}
    fine['sham', 11, 'amount'][5]['chemistry'][0][0] *= np.exp(.0002)
    report = study.response_comparison(coarse, fine, jobs, [11, 17], p, 60.)
    assert not report['passed']
    bad = next(r for r in report['rows'] if (r['context'], r['cell'], r['dose']) == ('sham', 11, 'amount'))
    assert bad['normalized_response_error'] > .01


def test_design_shares_background_controls_but_not_controls_across_contexts(tmp_path):
    selections = {str(s):dict(selected_ids=[11, 17]) for s in study.HISTORIES}
    contexts, jobs = study.job_design(tmp_path, selections)
    assert len(contexts) == 9 and len(jobs) == len({j['key'] for j in jobs}) == 66
    assert sum(j['target'] is None for j in jobs) == 18
    assert sum(j['dose'] == 'amount' for j in jobs) == 24
    assert all(j['start'] == 510. and j['stage'] != 'response' for j in jobs)


def test_failed_pilot_prevents_long_continuation(tmp_path, monkeypatch):
    selections = {str(s):dict(selected_ids=[11, 17]) for s in study.HISTORIES}
    _, jobs = study.job_design(tmp_path, selections)
    p = dict(jobs=jobs, pilot_duration=6., duration=60., total_jobs=66)
    calls = []; monkeypatch.setattr(study, 'advance', lambda r, j, p, h:calls.append(h))
    monkeypatch.setattr(study, 'pair_report', lambda *args:dict(passed=False))
    done, passed = study.schedule_history(tmp_path, p, 7, 0)
    assert not passed and done == 0 and calls == [6.]*22


def zero_restart(tmp_path):
    host, source = source_fixture(tmp_path, start=510.)
    root = tmp_path/'study'; root.mkdir(); write_json(root/'protocol.json', dict(test=True))
    job = dict(key='test', seed=host.config.seed, level='fine', dt=host.config.dt,
        source=str(source), start=510., duration=60.)
    p = dict(accepted_configs={f'{host.config.seed}_fine':asdict(host.config)}, amount_error_max=2e-14)
    history = [dict(elapsed=0., time=510., ids=host.ids.tolist(), chemistry=[host.activator.tolist(), host.inhibitor.tolist()])]
    audit = dict(max_volume_error=0., min_radius=5., max_clipping=0., boundary_max=0., dilution_error_max=0.)
    folder = root/job['key']; folder.mkdir()
    out = moving.payload(source)
    out['long_experiment'] = np.array(json.dumps(dict(audit=audit, history=history,
        protocol_hash=digest(root/'protocol.json'), job=job)))
    np.savez_compressed(folder/'latest_state.npz', **out)
    write_json(folder/'history.json', history)
    return root, job, p, out


def test_zero_restart_archives_exact_source_and_leaves_source_unchanged(tmp_path):
    root, job, p, _ = zero_restart(tmp_path); sha = digest(job['source'])
    assert study.archive_exact_zero_restart(root, job, p)
    assert not (root/job['key']/'latest_state.npz').exists() and digest(job['source']) == sha
    manifests = list((root/job['key']/'restart-archives').glob('*/manifest.json'))
    assert len(manifests) == 1 and read_json(manifests[0])['source_sha256'] == sha


def read_json(path):
    return json.loads(path.read_text())


def test_changed_zero_checkpoint_is_rejected_without_archiving(tmp_path):
    root, job, p, out = zero_restart(tmp_path); out['phase_carry'].flat[0] += 1e-16
    cp = root/job['key']/'latest_state.npz'; np.savez_compressed(cp, **out); sha = digest(cp)
    with pytest.raises(ValueError, match='phase_carry'):
        study.archive_exact_zero_restart(root, job, p)
    assert digest(cp) == sha and not (cp.parent/'restart-archives').exists()


def test_nonzero_checkpoint_is_left_to_original_continuation_validator(tmp_path):
    root, job, p, out = zero_restart(tmp_path)
    meta = json.loads(str(out['metadata'])); meta['time'] += job['dt']; meta['step_number'] += 1
    out['metadata'] = np.array(json.dumps(meta)); cp = root/job['key']/'latest_state.npz'
    np.savez_compressed(cp, **out); sha = digest(cp)
    assert not study.archive_exact_zero_restart(root, job, p) and digest(cp) == sha


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')
@pytest.mark.parametrize('dt', [.00375, .001875])
def test_small_matched_pulse_preserves_nonzero_carry_and_passes_gpu_context(tmp_path, dt):
    host, source = source_fixture(tmp_path, start=510.)
    config = type(host.config)(**{**asdict(host.config), 'grid':24, 'extent':1.2})
    expanded = AttributeSimulation(config, mode='direct')
    partition = .5*(1+np.tanh(expanded.xyz[0]/config.interface_width))
    expanded.phi = np.array([expanded.phi[0]*partition, expanded.phi[0]*(1-partition)], dtype=np.float32)
    expanded.target = expanded.volumes()
    for key in ('ids', 'parents', 'due', 'fate', 'activator', 'inhibitor', 'polarity'):
        setattr(expanded, key, getattr(host, key).copy())
    expanded.time, expanded.step_number = host.time, host.step_number; expanded.checkpoint(source)
    with np.load(source) as z:
        out = {k:z[k].copy() for k in z.files}
    out.update(precision_arm=np.array('phase_carry'), phase_carry=np.full(expanded.phi.shape, 1e-8), rounding=np.array([5, 2], np.int32))
    np.savez_compressed(source, **out)
    x = np.array([expanded.activator, expanded.inhibitor]); m = expanded.volumes()
    values, record = study.pulse(x, m, expanded.ids, int(expanded.ids[0]), 'amount', .005*x[0, 0]*m[0])
    assert record['factor'] == .995
    dest = tmp_path/'small-pulse.npz'; new = moving.prepared_checkpoint(source, dest, dt, values)
    job = dict(key='small-matched-pulse', seed=host.config.seed, level='fine', dt=dt, source=str(dest))
    p = dict(accepted_configs={f'{host.config.seed}_fine':asdict(new.config)}, amount_error_max=2e-14,
        prefix_native_threads=2, prefix_duration=.015, interval=.00375)
    gate = moving.context_gate(tmp_path, job, p)
    assert gate['passed'] and gate['restart_exact_steps'] == 4 and gate['cpu_chemical_log_max'] <= 1e-11
