from types import SimpleNamespace

import numpy as np
import pytest

from embryo import polarity_robustness_refinement as study
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.resolution import write_json


def physical_source(tmp_path):
    c = Config(grid=16, extent=.8, max_cells=2, dt=.001875, steps=80000,
               save_every=80, differentiation=False, signal_transport='conservative')
    host = AttributeSimulation(c, mode='direct')
    host.step_number = 80000; host.time = 150.
    path = tmp_path/'source.npz'; host.checkpoint(path)
    return host, path


def test_second_halving_preserves_physical_state_and_random_streams(tmp_path):
    host, path = physical_source(tmp_path)
    refined = study.retime(path, tmp_path/'finer.npz', .0009375)
    study.check_halving(host, refined)
    assert refined.step_number == 160000
    assert refined.config.save_every*refined.config.dt == host.config.save_every*host.config.dt
    assert refined.config.steps*refined.config.dt == host.config.steps*host.config.dt
    refined.signal_rng.normal()
    with pytest.raises(ValueError, match='random stream'):
        study.check_halving(host, refined)


@pytest.mark.parametrize('attribute,value,match', [
    ('surface_tension', .9, 'physical parameters'),
    ('save_every', 80, 'physical observation'),
    ('dt', .00046875, 'exactly one')])
def test_halving_rejects_other_parameter_and_clock_changes(tmp_path, attribute, value, match):
    host, path = physical_source(tmp_path)
    refined = study.retime(path, tmp_path/'finer.npz', .0009375)
    setattr(refined.config, attribute, value)
    with pytest.raises(ValueError, match=match): study.check_halving(host, refined)


def traces():
    return [dict(elapsed=t, time=150+t, ids=[3, 4],
        chemistry=[np.exp([-s, s]).tolist(), [1., 1.]], volumes=[1., 1.],
        delta=[[-1., 1.], [1., -1.]], polarity=[[1., 0., 0.], [-1., 0., 0.]],
        axis_ratio=1.2, uniform_growth_max=.01-t/15, log_activator_sd=s)
        for t, s in zip([0., .15, .3], [1e-8, .2, .2])]


def comparison_inputs(tmp_path, monkeypatch):
    p = dict(start=150., interval=.15, pilot_duration=.15, duration=.3,
        pilot_late_window=.15, late_window=.15, criteria=dict(late_log_sd_min=.1),
        refinement_criteria=study.REFINEMENT_CRITERIA)
    reference = tmp_path/'reference'; reference.mkdir()
    fine = traces(); write_json(reference/'history.json', fine)
    chemical = tmp_path/'chemical.npz'
    np.savez(chemical, ids=[3, 4], uniform=np.array(fine[0]['chemistry']))
    job = dict(key='finer', family='uniform', chemical_file=str(chemical), dt=.0009375, level='finer')
    p.update(reference_folder=str(reference), reference_job=dict(job, dt=.001875), jobs=[job])
    folder = tmp_path/'finer'; folder.mkdir(); write_json(folder/'history.json', fine)
    (folder/'latest_state.npz').write_bytes(b'checkpoint')
    write_json(tmp_path/'protocol.json', p)
    monkeypatch.setattr(study, 'restore_checkpoint', lambda *args:
        (SimpleNamespace(time=150.3), {}, study.read(folder/'history.json')))
    return p, folder


def test_transient_error_fails_despite_matching_final_chemistry_and_outcomes(tmp_path, monkeypatch):
    p, folder = comparison_inputs(tmp_path, monkeypatch)
    history = study.read(folder/'history.json'); history[1]['chemistry'][0][0] *= np.exp(.02)
    history[1]['log_activator_sd'] = float(np.std(np.log(history[1]['chemistry'][0])))
    write_json(folder/'history.json', history)
    result = study.compare(tmp_path, p, .3)
    assert not result['passed']
    assert result['errors']['chemical_log_max'] == pytest.approx(.02)
    assert result['persistent_contrast'] == [True, True]
    assert result['compared_timesteps'] == [.001875, .0009375]


def test_pilot_evidence_is_immutable_after_continuation(tmp_path, monkeypatch):
    p, folder = comparison_inputs(tmp_path, monkeypatch)
    result = study.compare(tmp_path, p, .15); assert result['passed']
    pilot_hash = study.digest(tmp_path/'pilot-refinement.json')
    write_json(folder/'history.json', traces()+[dict(elapsed=.45)])
    assert study.compare(tmp_path, p, .15) == result
    assert study.digest(tmp_path/'pilot-refinement.json') == pilot_hash
    write_json(folder/'pilot-history.json', [])
    with pytest.raises(ValueError, match='Changed refinement evidence'):
        study.compare(tmp_path, p, .15)


def test_saved_numerical_decision_is_independently_recomputed(tmp_path, monkeypatch):
    p, _ = comparison_inputs(tmp_path, monkeypatch)
    result = study.compare(tmp_path, p, .3)
    result['passed'] = False; write_json(tmp_path/'long-refinement.json', result)
    with pytest.raises(ValueError, match='Changed refinement decision'):
        study.compare(tmp_path, p, .3)


def stage_protocol():
    return dict(jobs=[dict(key='finer', level='finer', dt=.0009375)],
        accepted_configs=dict(finer={}), pilot_duration=60., duration=240.)


def test_failed_prefix_blocks_all_scientific_steps(tmp_path, monkeypatch):
    monkeypatch.setattr(study, 'prefix', lambda *args: (_ for _ in ()).throw(RuntimeError('backend error')))
    monkeypatch.setattr(study, 'advance', lambda *args: pytest.fail('Backend failure must stop all advances'))
    result = study.run_stages(tmp_path, stage_protocol())
    assert result['state'] == 'stopped_at_backend_gate'


def test_failed_pilot_blocks_long_advance(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(study, 'prefix', lambda *args: dict(passed=True))
    monkeypatch.setattr(study, 'advance', lambda root, job, p, horizon: called.append(horizon))
    monkeypatch.setattr(study, 'compare', lambda *args: dict(passed=False))
    result = study.run_stages(tmp_path, stage_protocol())
    assert called == [60.]
    assert result['state'] == 'stopped_at_numerical_gate' and result['completed'] == 0


@pytest.mark.parametrize('long_pass', [False, True])
def test_only_finer_job_advances_and_full_failure_is_retained(tmp_path, monkeypatch, long_pass):
    called = []
    monkeypatch.setattr(study, 'prefix', lambda *args: dict(passed=True))
    monkeypatch.setattr(study, 'advance', lambda root, job, p, horizon: called.append((job['key'], horizon)))
    monkeypatch.setattr(study, 'compare', lambda *args: dict(passed=True))
    monkeypatch.setattr(study, 'assess', lambda *args: dict(passed=long_pass, prefix_pass=True,
        pilot_refinement_pass=True, long_refinement_pass=long_pass, endpoint_numerical_pass=True))
    study.run_stages(tmp_path, stage_protocol())
    assert called == [('finer', 60.), ('finer', 240.)]
    status = study.read(tmp_path/'status.json')
    assert status['state'] == ('completed' if long_pass else 'completed_with_unresolved_checks')
    assert status['long_refinement_pass'] == long_pass
