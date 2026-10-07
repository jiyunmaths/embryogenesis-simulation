"""Protect amount conservation, inherited carry and response-specific decisions."""
from dataclasses import asdict
import json
import os
from pathlib import Path

import numpy as np
import pytest

from embryo import phase_carry_exchange_response as study
from embryo.attribute_development import AttributeSimulation
from embryo.attribute_exchange import exchange
from embryo.resolution import write_json
from test_polarity_robustness import small_source


def source_fixture(tmp_path, carry=1e-8, start=390.):
    host, source, _ = small_source(tmp_path)
    with np.load(source) as z:
        out = {key:z[key].copy() for key in z.files}
    meta = json.loads(str(out['metadata'])); meta['time'] = start; meta['step_number'] = round(start/host.config.dt)
    out.update(metadata=np.array(json.dumps(meta)), precision_arm=np.array('phase_carry'),
        phase_carry=np.full(host.phi.shape, carry, dtype=np.float64), rounding=np.array([5, 2], np.int32))
    np.savez_compressed(source, **out)
    return AttributeSimulation.restore(source), source


def test_pair_selection_requires_distinct_context_and_is_deterministic():
    chemistry = np.array([[.1, 10., .3, 3.], [1., 2., 1.5, 1.8]])
    ids = np.array([20, 21, 22, 23])
    context = np.array([[0., 0., 0.], [0., 0., 0.], [1., 2., 3.], [3., 1., 2.]])
    selection = study.select_pair(chemistry, ids, context)
    assert set(selection['pair']) != {0, 1}
    assert selection['context_distance'] >= selection['median_context_distance']
    assert selection == study.select_pair(chemistry, ids, context)
    with pytest.raises(ValueError, match='geometric-context'):
        study.select_pair(chemistry, ids, np.ones((4, 3)))


def test_exchange_and_retiming_preserve_inherited_carry_and_random_streams(tmp_path):
    host, source = source_fixture(tmp_path)
    before = study.payload(source); values = np.array([host.activator, host.inhibitor])
    masses = np.array([1., 3.]); values[0] = [.1, 2.]; values[1] = [.3, 1.4]
    transplanted = exchange(values, masses, 0, 1)
    np.testing.assert_allclose(transplanted@masses, values@masses, rtol=1e-15)
    dest = tmp_path/'edited.npz'; new = study.prepared_checkpoint(source, dest, host.config.dt/2, transplanted)
    after = study.payload(dest)
    for key in before:
        if key not in ('metadata', 'activator', 'inhibitor', 'long_experiment'):
            np.testing.assert_array_equal(before[key], after[key])
    a, b = [json.loads(str(z['metadata'])) for z in (before, after)]
    for key in ('rng', 'fate_rng', 'signal_rng', 'time', 'lineage', 'divisions'):
        assert a[key] == b[key]
    assert new.step_number == host.step_number*2
    assert new.config.dt == host.config.dt/2
    np.testing.assert_array_equal([new.activator, new.inhibitor], transplanted)
    with pytest.raises(FileExistsError):
        study.prepared_checkpoint(source, dest, host.config.dt/2, transplanted)


def test_missing_or_invalid_carry_cannot_silently_fall_back(tmp_path):
    _, source = source_fixture(tmp_path)
    z = study.payload(source); z['phase_carry'] = z['phase_carry'].astype(np.float32)
    broken = tmp_path/'bad.npz'; np.savez_compressed(broken, **z)
    with pytest.raises(ValueError, match='carry'):
        study.payload(broken)
    z['phase_carry'] = np.full(z['phi'].shape, np.nan)
    np.savez_compressed(broken, **z)
    with pytest.raises(ValueError, match='carry'):
        study.payload(broken)


def response_fixture():
    times = np.arange(401)*.15
    ids = [10, 11]; masses = np.array([1., 2.])
    initial = np.array([[.1, 2.], [.3, 1.4]])
    transplanted = exchange(initial, masses, 0, 1)
    histories = {}
    for background, start, tau in (('unexchanged', initial, [1., 4.]), ('fresh_exchange', transplanted, [4., 1.])):
        control = np.repeat(start[None], len(times), axis=0)*np.exp(.002*times[:, None, None])
        for target in (None, *ids):
            path = control.copy()
            if target is not None:
                cell = ids.index(target); path[:, 0, cell] *= np.exp(np.log(.9)*np.exp(-times/tau[cell]))
            histories[(background, target)] = [dict(elapsed=float(t), time=450.+float(t), ids=ids,
                chemistry=x.tolist(), volumes=masses.tolist()) for t, x in zip(times, path)]
    return histories, dict(selected_ids=ids, pair=[0, 1], masses=masses.tolist()), dict(
        interval=.15, response_limits=study.RESPONSE_CRITERIA)


def test_donor_response_uses_own_drifting_control_and_signed_waveform():
    histories, selection, p = response_fixture()
    result = study.response_summary(histories, selection, p, 60.)
    assert all(r['nearest_reference'] == 'donor' and r['donor_distance'] < 1e-12 for r in result['comparisons'])
    assert all(r['baseline_separation'] > .01 for r in result['comparisons'])
    assert study.response_refinement(histories, histories, selection, p, 60.)['passed']
    assert all(r['injected_activator_amount'] < 0 for r in result['metrics'])


def test_response_gate_rejects_waveform_error_even_when_raw_error_is_small():
    coarse, selection, p = response_fixture()
    fine = {key:json.loads(json.dumps(rows)) for key, rows in coarse.items()}
    fine[('fresh_exchange', 10)][5]['chemistry'][0][0] *= np.exp(.005)
    result = study.response_refinement(coarse, fine, selection, p, 60.)
    assert not result['passed']
    assert max(r['normalized_response_error'] for r in result['rows']) > .01


def test_lost_contrast_cannot_be_hidden_by_small_raw_error():
    p = dict(late_window=.15, contrast_min=.1, mechanical_limits=study.MECHANICAL_LIMITS)
    coarse = []
    for t in (0., .15, .3):
        coarse.append(dict(elapsed=t, time=390.+t, ids=[1, 2], chemistry=[np.exp([.105, -.105]).tolist(), [1., 1.]],
            polarity=[[0., 0., 0.]]*2, volumes=[1., 1.], delta=[[-1., 1.], [1., -1.]],
            axis_ratio=1.2, log_activator_sd=.105))
    fine = json.loads(json.dumps(coarse))
    fine[-1]['chemistry'][0] = np.exp([.098, -.098]).tolist(); fine[-1]['log_activator_sd'] = .098
    result = study.mechanical_comparison(coarse, fine, np.zeros(2), np.zeros(2), p, .3)
    assert result['errors']['chemical_log_max'] < .01
    assert not result['passed'] and result['loss_time_error'] is None


def test_design_and_failed_pilot_block_long_stage(tmp_path, monkeypatch):
    selections = {str(s):dict(selected_ids=[20, 27]) for s in (7, 8, 9)}
    jobs = study.job_design(tmp_path, selections)
    assert len(jobs) == 48 and len(set(j['key'] for j in jobs)) == 48
    assert sum(j['stage'] == 'formation' for j in jobs) == 12
    assert sum(j['stage'] == 'response' for j in jobs) == 36
    assert all(j['factor'] == .9 for j in jobs if j['target'] is not None)
    p = dict(jobs=jobs, pilot_duration=6., duration=60.)
    calls = []
    monkeypatch.setattr(study, 'advance', lambda root, job, p, h:calls.append(h))
    monkeypatch.setattr(study, 'pair_report', lambda *args:dict(passed=False))
    done, passed = study.schedule_stage(tmp_path, p, 7, 'formation', 0)
    assert not passed and done == 0 and calls == [6.]*4


def test_response_handoff_preserves_carry_and_pulse_is_not_repeated(tmp_path, monkeypatch):
    host, source = source_fixture(tmp_path, start=450.)
    seed = host.config.seed; ids = host.ids.tolist(); pair = ids[:2]
    jobs = study.job_design(tmp_path, {str(s):dict(selected_ids=pair) for s in (7, 8, 9)})
    formation = next(j for j in jobs if j['seed'] == seed and j['stage'] == 'formation' and j['level'] == 'coarse' and j['background'] == 'fresh_exchange')
    response = next(j for j in jobs if j['seed'] == seed and j['stage'] == 'response' and j['level'] == 'coarse' and j['background'] == 'fresh_exchange' and j['target'] == pair[0])
    folder = tmp_path/formation['key']; folder.mkdir()
    folder.joinpath('latest_state.npz').write_bytes(source.read_bytes())
    values = np.array([host.activator, host.inhibitor]); volumes = host.volumes()
    history = [dict(chemistry=values.tolist(), ids=ids, volumes=volumes.tolist())]
    write_json(folder/'history.json', history); write_json(folder/'result.json', dict(quality_pass=True))
    write_json(tmp_path/'protocol.json', dict(test=True))
    monkeypatch.setattr(study, 'current_evidence', lambda *args:(dict(quality_pass=True), history))
    record = study.materialize_response(tmp_path, response, dict(jobs=jobs))
    edited = study.payload(response['source']); original = study.payload(source)
    np.testing.assert_array_equal(edited['phase_carry'], original['phase_carry'])
    np.testing.assert_array_equal(edited['phi'], original['phi'])
    assert edited['activator'][0] == original['activator'][0]*.9
    assert record['injected_activator_amount'] == pytest.approx(-.1*values[0, 0]*volumes[0])
    sha = study.digest(response['source'])
    assert study.materialize_response(tmp_path, response, dict(jobs=jobs)) == record
    assert study.digest(response['source']) == sha


def test_interrupted_continuation_preserves_initial_carry_and_saved_history(tmp_path, monkeypatch):
    host, source = source_fixture(tmp_path, carry=.25)
    seed = host.config.seed
    job = dict(key='continuation', seed=seed, level='coarse', dt=host.config.dt,
        stage='formation', background='unexchanged', source=str(source), target=None, factor=1.,
        start=390., duration=.3, point=dict(ratio=27.5, chi=.35))
    p = dict(accepted_configs={f'{seed}_coarse':asdict(host.config)}, interval=.15,
        checkpoint_interval=.15, pilot_duration=.15, late_window=.15, contrast_min=.1, amount_error_max=2e-14)
    write_json(tmp_path/'protocol.json', p); fail = [False]

    class FakeCarry:
        def __init__(self, path):
            self.host = AttributeSimulation.restore(path)
            out = study.payload(path); self.phase_carry = out['phase_carry']; self.rounding = out['rounding']
        def __getattr__(self, key):
            return getattr(self.host, key)
        @classmethod
        def restore(cls, path):
            return cls(path)
        def step(self):
            if fail[0] and self.host.step_number == round(390./job['dt'])+45:
                fail[0] = False; raise RuntimeError('interrupted')
            self.host.activator += .001*(1+float(self.phase_carry.flat[0]))
            self.phase_carry += .1; self.rounding[0] += 1
            self.host.step_number += 1; self.host.time = self.host.step_number*self.host.config.dt
        def audit(self):
            return dict(max_volume_error=0., min_radius=5., max_clipping=0., dilution_amount_error=0.)
        def precision_diagnostics(self):
            return dict(marker=float(self.phase_carry.flat[0]))
        def observe(self, elapsed):
            graph = self.host.signaling_graph()
            return dict(elapsed=elapsed, time=self.host.time, ids=self.host.ids.tolist(),
                chemistry=np.array([self.host.activator, self.host.inhibitor]).tolist(),
                volumes=graph.masses.tolist(), delta=graph.delta.tolist(), polarity=self.host.polarity.tolist(),
                axis_ratio=1.2, boundary_occupancy=0., log_activator_sd=float(np.std(np.log(self.host.activator))))
        def checkpoint(self, path):
            self.host.checkpoint(path)
            with np.load(path) as z:out={key:z[key].copy() for key in z.files}
            out.update(precision_arm=np.array('phase_carry'), phase_carry=self.phase_carry, rounding=self.rounding)
            np.savez_compressed(path, **out)

    monkeypatch.setattr(study, 'PrecisionSimulation', FakeCarry)
    monkeypatch.setattr(study, 'context_gate', lambda *args:dict(passed=True))
    monkeypatch.setattr(study.torch.cuda, 'empty_cache', lambda:None)
    pilot = study.advance(tmp_path, job, p, .15)
    assert len(pilot) == 2 and not (tmp_path/job['key']/'result.json').exists()
    fail[0] = True
    with pytest.raises(RuntimeError, match='interrupted'):
        study.advance(tmp_path, job, p, .3)
    final = study.advance(tmp_path, job, p, .3)
    saved = FakeCarry.restore(tmp_path/job['key']/'latest_state.npz')
    np.testing.assert_allclose(saved.activator, host.activator+.001*sum(1.25+.1*k for k in range(80)), rtol=1e-12)
    assert saved.phase_carry.flat[0] == pytest.approx(8.25)
    assert study.advance(tmp_path, job, p, .15) == final[:2]
    result, _ = study.current_evidence(tmp_path, job, p)
    assert result['quality_pass']


def test_paired_pilot_snapshots_are_reproducible_and_immutable(tmp_path):
    host, source = source_fixture(tmp_path)
    seed = host.config.seed; ids = host.ids.tolist(); m = host.volumes()
    selections = {str(s):dict(selected_ids=ids, pair=[0, 1], masses=m.tolist()) for s in (7, 8, 9)}
    jobs = [j for j in study.job_design(tmp_path, selections) if j['seed'] == seed and j['stage'] == 'formation']
    p = dict(interval=.15, pilot_duration=6., late_window=24., contrast_min=.1,
        selections=selections, jobs=jobs, mechanical_limits=study.MECHANICAL_LIMITS, accepted_configs={})
    for job in jobs:
        values = np.array([host.activator, host.inhibitor])
        if job['background'] == 'fresh_exchange':values = exchange(values, m, 0, 1)
        new = study.prepared_checkpoint(source, job['source'], job['dt'], values)
        p['accepted_configs'][f'{seed}_{job["level"]}'] = asdict(new.config)
        graph = new.signaling_graph(); folder = tmp_path/job['key']; folder.mkdir()
        history = [dict(elapsed=float(t), time=390.+float(t), ids=ids, chemistry=values.tolist(),
            volumes=graph.masses.tolist(), delta=graph.delta.tolist(), polarity=new.polarity.tolist(),
            log_activator_sd=float(np.std(np.log(values[0]))), axis_ratio=1.2) for t in np.arange(41)*.15]
        write_json(folder/'history.json', history)
        z = study.payload(job['source']); meta = json.loads(str(z['metadata']))
        meta['time'] = 396.; meta['step_number'] = round(396./job['dt']); z['metadata'] = np.array(json.dumps(meta))
        np.savez_compressed(folder/'latest_state.npz', **z)
    write_json(tmp_path/'protocol.json', p)
    result = study.pair_report(tmp_path, p, jobs, 6.)
    assert result['passed'] and result['same_state_outcome']
    assert study.pair_report(tmp_path, p, jobs, 6.) == result
    snapshot = next(Path(f) for f in result['evidence_sha256'] if f.endswith('.npz'))
    with np.load(snapshot) as z:phi = z['phi'].copy()
    phi.flat[0] += .001; np.savez_compressed(snapshot, phi=phi)
    with pytest.raises(ValueError, match='paired history/fields'):
        study.pair_report(tmp_path, p, jobs, 6.)


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')
@pytest.mark.parametrize('dt', [.00375, .001875])
def test_nonzero_inherited_carry_context_gate_on_gpu(tmp_path, dt):
    host, source = source_fixture(tmp_path)
    expanded = AttributeSimulation(type(host.config)(**{**asdict(host.config), 'grid':24, 'extent':1.2}), mode='direct')
    partition = .5*(1+np.tanh(expanded.xyz[0]/expanded.config.interface_width))
    expanded.phi = np.array([expanded.phi[0]*partition, expanded.phi[0]*(1-partition)], dtype=np.float32)
    expanded.target = expanded.volumes()
    for key in ('ids', 'parents', 'due', 'fate', 'activator', 'inhibitor', 'polarity'):
        setattr(expanded, key, getattr(host, key).copy())
    expanded.time = host.time; expanded.step_number = host.step_number; expanded.checkpoint(source)
    with np.load(source) as z:out={key:z[key].copy() for key in z.files}
    out.update(precision_arm=np.array('phase_carry'), phase_carry=np.full(expanded.phi.shape, 1e-8), rounding=np.array([5, 2], np.int32))
    np.savez_compressed(source, **out)
    edited = tmp_path/'retimed.npz'; new = study.prepared_checkpoint(source, edited, dt)
    job = dict(key='actual-context', seed=host.config.seed, level='coarse', dt=dt, source=str(edited))
    p = dict(accepted_configs={f'{host.config.seed}_coarse':asdict(new.config)},
        amount_error_max=2e-14, prefix_native_threads=2, prefix_duration=.015, interval=.00375)
    result = study.context_gate(tmp_path, job, p)
    assert result['passed'] and result['restart_exact_steps'] == 4
    assert result['cpu_chemical_log_max'] < 1e-11
