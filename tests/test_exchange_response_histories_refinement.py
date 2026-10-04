import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('torch')
from embryo import exchange_response_histories_refinement as study
from embryo.attribute_development import AttributeSimulation
from embryo.cell_response_moving import name
from embryo.feedback_long import digest
from embryo.model import Config
from embryo.resolution import write_json
from test_exchange_response_histories import formation_fixture
from test_cell_exchange_response_refinement import fixture_study as response_fixture


def paired_formation(tmp_path, outcome='donor'):
    base, fine = tmp_path/'base', tmp_path/'fine'
    p, _, _, _ = formation_fixture(base)
    formation_fixture(fine, outcome=outcome)
    study.formation_assessment(base, 8, p)
    for root in (base, fine):
        for family in study.BACKGROUNDS:
            np.savez(root/'seed-8'/'formation'/f'{family}_control'/'latest_state.npz',
                     phi=np.full((2, 4, 4, 4), .5, np.float32))
    p.update(baseline=str(base), formation_criteria=study.FORMATION_CRITERIA)
    return base, fine, p


def test_identical_formation_and_outcomes_pass(tmp_path):
    _, fine, p = paired_formation(tmp_path)
    result = study.formation_comparison(fine, 8, p)
    assert result['passed']
    assert all(r['outcomes_unchanged'] and not any(r['errors'].values()) for r in result['trials'])


def test_changed_collective_outcome_fails_refinement(tmp_path):
    _, fine, p = paired_formation(tmp_path, outcome='destination')
    result = study.formation_comparison(fine, 8, p)
    assert not result['passed']
    assert not result['trials'][1]['outcomes_unchanged']


@pytest.mark.parametrize('kind', ['chemistry', 'polarity', 'volume', 'transport', 'shape', 'field'])
def test_numerical_drift_fails_even_with_passing_worker_flags(tmp_path, kind):
    _, fine, p = paired_formation(tmp_path)
    folder = fine/'seed-8'/'formation'/'fresh_exchange_control'
    if kind == 'field':
        np.savez(folder/'latest_state.npz', phi=np.full((2, 4, 4, 4), .4, np.float32))
    else:
        h = json.loads((folder/'history.json').read_text())
        if kind == 'chemistry': h[1]['chemistry'][0][0] *= 1.03
        elif kind == 'polarity': h[1]['polarity'][0][0] += .02
        elif kind == 'volume': h[1]['volumes'][0] *= 1.01
        elif kind == 'transport': h[1]['delta'][0][0] *= 1.10
        elif kind == 'shape': h[1]['axis_ratio'] *= 1.02
        write_json(folder/'history.json', h)
    result = study.formation_comparison(fine, 8, p)
    assert not result['passed'] and not result['trials'][1]['passed']


@pytest.mark.parametrize('kind', ['initial', 'clock', 'ids', 'incomplete'])
def test_comparison_rejects_incompatible_or_partial_evidence(tmp_path, kind):
    _, fine, p = paired_formation(tmp_path)
    folder = fine/'seed-8'/'formation'/'fresh_exchange_control'
    h = json.loads((folder/'history.json').read_text())
    if kind == 'initial': h[0]['chemistry'][0][0] *= 1.001
    elif kind == 'clock': h[1]['time'] += .01
    elif kind == 'ids': h[1]['ids'].reverse()
    else: h.pop()
    write_json(folder/'history.json', h)
    with pytest.raises(ValueError): study.formation_comparison(fine, 8, p)


def test_retimed_source_preserves_fields_ids_polarity_and_random_streams(tmp_path):
    source = tmp_path/'source'; source.mkdir()
    sim = AttributeSimulation(Config(grid=16, extent=1.4, interface_width=.14,
                                     max_cells=1, dt=.00375), 'direct')
    sim.time = .03; sim.step_number = 8
    sim.activator[:] = 1.2; sim.inhibitor[:] = 1.4; sim.polarity[:] = [.2, .1, 0.]
    sim.checkpoint(source/'source.npz')
    np.savez(source/'initial_states.npz', state=np.array([sim.activator,sim.inhibitor]),
             ids=sim.ids, masses=sim.volumes())
    jobs = [dict(family='state', target=None, factor=1.)]
    cp = dict(checkpoint=str(source/'source.npz'), start=.03, dt=.00375,
              duration=.015, interval=.00375, jobs=jobs, source_sha256={}, input_sha256={})
    write_json(source/'protocol.json', cp)
    before = digest(source/'source.npz')
    output = tmp_path/'fine'
    study.prepare_source(source, output, .001875, jobs, {})
    fine = AttributeSimulation.restore(output/'source.npz')
    assert fine.time == sim.time and fine.step_number == 16 and fine.config.dt == .001875
    for key in ('phi','target','ids','parents','due','activator','inhibitor','polarity'):
        np.testing.assert_array_equal(getattr(fine,key),getattr(sim,key))
    for key in ('rng','signal_rng','fate_rng'):
        assert getattr(fine,key).bit_generator.state == getattr(sim,key).bit_generator.state
    assert digest(source/'source.npz') == before
    assert digest(source/'initial_states.npz') == digest(output/'initial_states.npz')
    study.verify(study.read(output/'protocol.json'))


@pytest.mark.parametrize('state,count,passed', [('prepared',0,False),('running',3,False),
                                              ('failed',0,False),('completed',2,True),
                                              ('completed',4,False)])
def test_partial_or_failed_fine_gpu_gate_stops_scientific_execution(tmp_path, monkeypatch, state, count, passed):
    gate = tmp_path/'gate'; gate.mkdir()
    write_json(gate/'status.json',dict(state=state,completed=count,passed=passed))
    write_json(tmp_path/'protocol.json',dict(validation=str(gate),dt=.001875,contexts=[],
        source_sha256={},input_sha256={}))
    started=[]
    monkeypatch.setattr(study.gpu_validation,'run',lambda _: None)
    monkeypatch.setattr(study,'gpu_run',lambda _: started.append(True))
    with pytest.raises(ValueError): study.run(tmp_path)
    assert not started and study.read(tmp_path/'status.json')['state']=='failed'


def test_existing_direct_response_criteria_detect_transfer_and_recovery_changes(tmp_path):
    _, fine = response_fixture(tmp_path)
    assert study.response_assessment(fine)['passed']
    source = fine/'unexchanged'/'unexchanged_cell-27_factor-0.9'/'history.json'
    target = fine/'fresh_exchange'/'fresh_exchange_cell-27_factor-0.9'/'history.json'
    target.write_text(source.read_text())
    result = study.response_assessment(fine)
    assert not result['passed'] and not result['classifications'][0]['passed']


def test_output_overwrite_rejected_before_accessing_old_evidence(tmp_path):
    with pytest.raises(FileExistsError): study.prepare(tmp_path)
