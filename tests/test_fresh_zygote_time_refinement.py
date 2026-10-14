"""Protect model invariance, failed-pilot preservation and diagnostic scope."""
from dataclasses import asdict
from types import SimpleNamespace

import pytest

from embryo import fresh_zygote_time_refinement as f
from embryo.model import Config


def reference():
    return asdict(Config(seed=1103,grid=120,dt=f.REFERENCE_DT,steps=76800,
        save_every=640,feedback=False,fate_tension=0.,fate_adhesion=0.,polarity_tension=0.))


def test_only_timestep_and_integer_clocks_change():
    base=reference(); candidate=f.candidate_config(base)
    assert {k for k,v in candidate.items() if v!=base[k]}=={'dt','steps','save_every'}
    assert candidate['dt']==.0001171875 and candidate['steps']==153600
    assert candidate['steps']*candidate['dt']==18
    assert candidate['save_every']*candidate['dt']==.15
    f.validate_config(base,candidate)


@pytest.mark.parametrize('change',[dict(seed=1104),dict(grid=104),dict(dt=.00005859375),
    dict(graph_contact_cutoff=0.),dict(signal_dh=.56),dict(feedback=True),
    dict(polarity_tension=.35),dict(steps=76800),dict(save_every=1)])
def test_unregistered_model_clock_or_history_drift_is_rejected(change):
    with pytest.raises(ValueError):
        f.validate_config(reference(),dict(f.candidate_config(reference()),**change))


def test_existing_directory_is_never_overwritten(tmp_path):
    marker=tmp_path/'original';marker.write_text('retain failure')
    with pytest.raises(FileExistsError):f.prepare(tmp_path)
    assert marker.read_text()=='retain failure'


@pytest.mark.parametrize('stage',['prefix','path','comparison'])
def test_failure_stops_without_reference_execution_or_extension(tmp_path,monkeypatch,stage):
    binary=tmp_path/'kernel';binary.write_bytes(b'unchanged')
    jobs=[dict(key=f.REFERENCE,reused=True),dict(key=f.CANDIDATE,reused=False)]
    p=dict(jobs=jobs,time_pair=dict(seed=1103,branch='off',cases=[f.REFERENCE,f.CANDIDATE]),
        runtime_binary_sha256=f.digest(binary))
    f.write(tmp_path/'protocol.json',p)
    monkeypatch.setattr(f,'validate_prelaunch',lambda *args:None)
    monkeypatch.setattr(f.td.gpu,'library',lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(f.torch.cuda,'get_device_name',lambda *args:'fixture')
    calls=[]
    def prefix(root,job,p):
        calls.append(('prefix',job['key']));return dict(passed=stage!='prefix')
    def execute(root,job):
        calls.append(('path',job['key']));return dict(passed=stage!='path',seconds=1.)
    monkeypatch.setattr(f.pilot,'prefix',prefix)
    monkeypatch.setattr(f.controls,'execute_path',execute)
    monkeypatch.setattr(f.pilot,'assess_pair',lambda *args:dict(passed=stage!='comparison'))
    with pytest.raises(RuntimeError):f.run(tmp_path)
    status=f.read(tmp_path/'status.json')
    assert status['state']=='failed' and not status['qualification_accepted']
    assert not status['post_cleavage_started'] and status['new_scientific_histories']==0
    assert status['original_pilot_remains_failed']
    assert all(case==f.CANDIDATE for _,case in calls)
    assert not (tmp_path/'summary.json').exists()


def test_passing_repeat_requires_review_and_does_not_accept_original_cohort(tmp_path,monkeypatch):
    binary=tmp_path/'kernel';binary.write_bytes(b'unchanged')
    p=dict(jobs=[dict(key=f.REFERENCE,reused=True),dict(key=f.CANDIDATE,reused=False)],
        time_pair=dict(seed=1103,branch='off',cases=[f.REFERENCE,f.CANDIDATE]),
        runtime_binary_sha256=f.digest(binary),scope='Existing-history numerical diagnosis')
    f.write(tmp_path/'protocol.json',p)
    for name in ('validate_prelaunch','verify','check_reference'):
        monkeypatch.setattr(f,name,lambda *args:None)
    monkeypatch.setattr(f.td.gpu,'library',lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(f.torch.cuda,'get_device_name',lambda *args:'fixture')
    monkeypatch.setattr(f.pilot,'prefix',lambda *args:dict(passed=True))
    monkeypatch.setattr(f.controls,'execute_path',lambda *args:dict(passed=True,seconds=1.))
    monkeypatch.setattr(f.pilot,'assess_pair',lambda *args:dict(passed=True))
    result=f.run(tmp_path)
    assert result['developmental_screens_passed'] and result['independent_review_pending']
    assert not result['qualification_accepted'] and result['original_pilot_remains_failed']
    assert result['new_numerical_paths']==1 and result['new_scientific_histories']==0
    assert not result['post_cleavage_started']
