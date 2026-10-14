"""Protect finer-clock controls, immutable references and failure stopping."""
from dataclasses import asdict
from types import SimpleNamespace

import pytest

from embryo import fresh_zygote_control_completion as f
from embryo.model import Config


def baseline():
    return asdict(Config(seed=1103, grid=120, dt=.000234375, steps=76800,
        save_every=640, feedback=False, fate_tension=0., fate_adhesion=0., polarity_tension=0.))


def protocol(root):
    jobs = [dict(key=key, config=baseline(), reused=True) for key in f.REFERENCES]
    jobs += [dict(j, source=str(root/j['key']/'source.npz')) for j in f.design(baseline())]
    return dict(jobs=jobs, histories=[1103], branches=list(f.BRANCHES),
        time_pairs=[dict(seed=1103, branch=b,
            cases=[f'h1103_{b}_{level}' for level, _ in f.LEVELS]) for b in f.BRANCHES],
        new_scientific_histories=0, new_numerical_paths=6, reused_numerical_paths=2)


def test_three_matched_branches_two_clocks_and_no_physics_drift(tmp_path):
    p = protocol(tmp_path); f.validate_jobs(tmp_path, p)
    assert len(p['jobs']) == 8
    for job in p['jobs'][2:]:
        c = job['config']; f.controls.config_changes(baseline(), c, job['branch'])
        assert c['dt'] in (.000234375, .0001171875)
        assert c['steps']*c['dt'] == 18 and c['save_every']*c['dt'] == .15
        assert {k for k,v in c.items() if v != baseline()[k]} <= {
            'dt','steps','save_every','feedback','fate_tension','fate_adhesion','polarity_tension'}


@pytest.mark.parametrize('change', [dict(seed=1104), dict(grid=104), dict(dt=.00046875),
    dict(signal_dh=.56), dict(graph_contact_cutoff=0.), dict(polarity_tension=.7),
    dict(steps=38400), dict(save_every=1)])
def test_unregistered_history_clock_or_model_drift_is_rejected(tmp_path, change):
    p = protocol(tmp_path); p['jobs'][-1]['config'].update(change)
    with pytest.raises(ValueError): f.validate_jobs(tmp_path, p)


def test_duplicate_or_missing_branch_is_not_complete():
    rows = [dict(seed=1103, branch=b, passed=True) for b in f.BRANCHES]
    assert f.complete_design(rows)
    assert not f.complete_design(rows[:-1])
    assert not f.complete_design(rows[:-1]+[rows[0]])
    assert not f.complete_design([dict(r, passed=False) for r in rows])
    assert not f.complete_design([dict(r, seed=1104) for r in rows])


def test_existing_study_is_never_overwritten(tmp_path):
    marker = tmp_path/'original'; marker.write_text('retain')
    with pytest.raises(FileExistsError): f.prepare(tmp_path)
    assert marker.read_text() == 'retain'


def mock_run(tmp_path, monkeypatch, failure):
    binary = tmp_path/'binary'; binary.write_bytes(b'unchanged')
    p = protocol(tmp_path)
    p.update(runtime_binary_sha256=f.digest(binary), scope='Existing-history controls')
    f.write(tmp_path/'protocol.json', p)
    for name in ('validate_prelaunch','verify','check_references'):
        monkeypatch.setattr(f, name, lambda *args:None)
    monkeypatch.setattr(f.td.gpu, 'library', lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(f.torch.cuda, 'get_device_name', lambda *args:'fixture')
    calls = []
    def prefix(root, job, p):
        calls.append(('prefix', job['key'])); return dict(passed=failure!='prefix')
    def execute(root, job):
        calls.append(('path', job['key'])); return dict(passed=failure!='path', seconds=1.)
    def compare(root, pair, p):
        calls.append(('pair', pair['branch'])); return dict(pair, passed=failure!='comparison')
    monkeypatch.setattr(f.pilot, 'prefix', prefix)
    monkeypatch.setattr(f.controls, 'execute_path', execute)
    monkeypatch.setattr(f.pilot, 'assess_pair', compare)
    return calls


@pytest.mark.parametrize('failure', ['prefix','path','comparison'])
def test_failed_gate_stops_next_branch_and_never_executes_reference(tmp_path, monkeypatch, failure):
    calls = mock_run(tmp_path, monkeypatch, failure)
    with pytest.raises(RuntimeError): f.run(tmp_path)
    status = f.read(tmp_path/'status.json')
    assert status['state'] == 'failed' and not status['qualification_accepted']
    assert not status['post_cleavage_started'] and status['original_pilot_remains_failed']
    assert status['new_scientific_histories'] == 0 and status['reused_paths'] == 2
    assert not any(case in f.REFERENCES for _,case in calls)
    assert all('tension_adhesion' in case for _,case in calls)
    assert not (tmp_path/'summary.json').exists()


def test_passing_controls_require_review_and_preserve_original_failure(tmp_path, monkeypatch):
    calls = mock_run(tmp_path, monkeypatch, None)
    result = f.run(tmp_path)
    assert result['developmental_screens_passed'] and result['independent_review_pending']
    assert not result['qualification_accepted'] and not result['post_cleavage_started']
    assert result['original_pilot_remains_failed'] and result['new_scientific_histories'] == 0
    assert result['new_numerical_paths'] == 6 and result['reused_numerical_paths'] == 2
    assert len([x for x in calls if x[0]=='path']) == 6
    assert not any(case in f.REFERENCES for _,case in calls)


def test_prepared_run_rejects_a_recorded_path(tmp_path, monkeypatch):
    p = protocol(tmp_path); p['gpu_uuid'] = 'fixture'
    monkeypatch.setattr(f, 'verify', lambda *args:None)
    monkeypatch.setattr(f, 'check_references', lambda *args:None)
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', 'fixture')
    f.write(tmp_path/'protocol.json', p)
    f.write(tmp_path/'status.json', dict(state='prepared',protocol_sha256=f.digest(tmp_path/'protocol.json')))
    folder = tmp_path/p['jobs'][2]['key']; folder.mkdir()
    (folder/'gpu-result.json').write_text('original result')
    with pytest.raises(FileExistsError): f.validate_prelaunch(tmp_path, p)
