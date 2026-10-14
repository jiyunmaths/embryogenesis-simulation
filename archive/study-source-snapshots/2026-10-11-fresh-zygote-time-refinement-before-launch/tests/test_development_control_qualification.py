"""Guard the branch design, read-only baseline and required refinement/replays."""
from dataclasses import asdict
import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from embryo import development_control_qualification as c
from embryo.model import Config
from embryo.resolution import _steps


def baseline():
    return asdict(Config(grid=120, dt=.00046875, steps=38400, save_every=320,
        feedback=False, fate_tension=0., fate_adhesion=0., polarity_tension=0.))


@pytest.mark.parametrize('branch', c.BRANCHES)
def test_only_prespecified_mechanics_and_discretization_can_change(branch):
    base = baseline(); candidate = dict(base)
    candidate.update(zip(('feedback','fate_tension','fate_adhesion','polarity_tension'), c.COUPLINGS[branch]))
    candidate.update(grid=104, dt=.000234375, steps=76800, save_every=640)
    c.config_changes(base, candidate, branch)
    with pytest.raises(ValueError, match='Physical model drift'):
        c.config_changes(base, dict(candidate, signal_dh=.56), branch)
    with pytest.raises(ValueError, match='matched control'):
        c.config_changes(base, dict(candidate, polarity_tension=.7), branch)
    with pytest.raises(ValueError, match='clock drift'):
        c.config_changes(base, dict(candidate, save_every=320), branch)


def test_all_design_clocks_and_candidate_resolution_are_exact():
    assert len(c.LEVELS)==4 and c.BRANCHES==('tension_adhesion','polarity','full')
    for _,grid,dt in c.LEVELS:
        assert grid in (88,104,120)
        for time in (18., .15, .015, 3.):
            assert abs(_steps(time,dt)*dt-time)<1e-12


def test_existing_study_and_baseline_cannot_be_overwritten(tmp_path, monkeypatch):
    marker = tmp_path/'original'; marker.write_text('preserved')
    with pytest.raises(FileExistsError):
        c.prepare(tmp_path)
    assert marker.read_text()=='preserved'
    monkeypatch.setattr(c.td, 'run_path', lambda *a:pytest.fail('Reused evidence would be executed'))
    with pytest.raises(ValueError, match='read-only'):
        c.execute_path(tmp_path, dict(key='off120', reused=True))
    old = tmp_path/'old'; old.mkdir(); (tmp_path/'linked').symlink_to(old)
    with pytest.raises(ValueError, match='read-only'):
        c.execute_path(tmp_path, dict(key='linked', reused=False))


def test_reuse_checks_original_protocol_configuration_and_all_evidence(tmp_path):
    old = tmp_path/'original'; old.mkdir(); (tmp_path/'off120').symlink_to(old)
    evidence = old/'state'; evidence.write_text('original state')
    original_job = dict(key='off120', config={'grid':120})
    report = dict(job=original_job, passed=True, protocol_sha256='parent',
        evidence_sha256={str(evidence):c.digest(evidence)})
    path=old/'gpu-result.json'; path.write_text(json.dumps(report))
    provenance=dict(folder=str(old), job=original_job, protocol_sha256='parent', report_sha256=c.digest(path))
    p=dict(reuse_provenance={'off120':provenance}); job=dict(original_job, reused=True)
    assert c.case_report(tmp_path,job,p)==report
    with pytest.raises(ValueError, match='configuration'):
        c.case_report(tmp_path,dict(job,config={'grid':104}),p)
    evidence.write_text('changed')
    with pytest.raises(ValueError, match='evidence changed'):
        c.case_report(tmp_path,job,p)


@pytest.mark.parametrize('failed', ['physical','space','time','dense','trend'])
def test_each_branch_requires_all_finest_gates(failed):
    paths=[{'passed':True}];spatial=[{'passed':False},{'passed':True}]
    temporal={'passed':True};dense={'peak':{'relative_spectrum':.009}}
    trend=True
    assert c.branch_decision(paths,spatial,temporal,dense,trend,.01)
    if failed=='physical': paths[0]['passed']=False
    if failed=='space': spatial[-1]['passed']=False
    if failed=='time': temporal['passed']=False
    if failed=='dense': dense['peak']['relative_spectrum']=.011
    if failed=='trend': trend=False
    assert not c.branch_decision(paths,spatial,temporal,dense,trend,.01)


def complete_decisions():
    branches=[dict(branch=b, scoped_screens_passed=True) for b in c.BRANCHES]
    replays=[dict(case=key,stage=stage,passed=True)
        for key in ['off120']+[b+'_space120' for b in c.BRANCHES] for stage in (2,4,8)]
    return branches,replays


def test_three_controls_and_all_distinct_live_stages_are_required():
    branches,replays=complete_decisions()
    assert c.final_decision(branches,replays,12)
    assert not c.final_decision(branches[:-1],replays,12)
    assert not c.final_decision(branches,replays[:-1],12)
    replays[-1]=dict(replays[0])
    assert not c.final_decision(branches,replays,12)
    branches,replays=complete_decisions();branches[-1]=dict(branches[0])
    assert not c.final_decision(branches,replays,12)


@pytest.mark.parametrize('failure', ['refinement','restart'])
def test_restart_success_cannot_override_refinement_failure_or_reverse(failure):
    branches,replays=complete_decisions()
    if failure=='refinement': branches[1]['scoped_screens_passed']=False
    else: replays[7]['passed']=False
    assert not c.final_decision(branches,replays,12)


def test_full_restart_refuses_existing_evidence_before_loading_gpu(tmp_path, monkeypatch):
    saved=tmp_path/'replays/full_space120/stage-4';saved.mkdir(parents=True)
    marker=saved/'result.json';marker.write_text('preserved')
    monkeypatch.setattr(c.td.gpu.DevelopmentGpu,'restore',lambda *a:pytest.fail('GPU restore must not happen'))
    with pytest.raises(FileExistsError):
        c.replay_live_state(tmp_path,dict(key='full_space120'),4,{})
    assert marker.read_text()=='preserved'


@pytest.mark.parametrize('fault', ['none', 'observation', 'checkpoint'])
def test_full_replay_records_mismatch_as_failure_without_overwriting_inputs(tmp_path, monkeypatch, fault):
    # A one-step state fixture exercises the replay's failure/recording path
    # without loading CUDA or duplicating the expensive developmental solver.
    cfg=Config(**baseline());folder=tmp_path/'full_space120';folder.mkdir()
    for name in ['restart-stage-2.npz','gpu-result.json']+[f'gpu-t{t:02d}.npz' for t in (0,3,6,9,12,15,18)]:
        (folder/name).write_text('original bytes')
    record=dict(time=18.,parents=[-1]*16,spectrum=[0.]*16)
    (folder/'gpu-history.json').write_text(json.dumps([record]))
    original={p.name:p.read_bytes() for p in folder.iterdir()}

    class State:
        def __init__(self):
            self.config=cfg;self.step_number=38399;self.time=self.step_number*cfg.dt
            self.ids=np.arange(3);self.target=torch.ones(1,dtype=torch.float64)
            self.activator=self.inhibitor=torch.ones(3,dtype=torch.float64)
            self.host=SimpleNamespace(divisions={0:{}},parents=np.full(3,-1));self.event_audits=[]
        def audit(self):
            return dict(max_volume_error=0.,min_radius=9.,max_clipping=0.,dilution_amount_error=0.)
        def step(self):
            self.step_number+=1;self.time=self.step_number*cfg.dt;self.ids=np.arange(16)
            self.host.parents=np.full(16,-1);self.host.divisions={};self.event_audits=[{} for _ in range(15)]
        def checkpoint(self,path):
            path.write_text('separate replay state')
    monkeypatch.setattr(c.td.gpu.DevelopmentGpu,'restore',lambda *a:State())
    monkeypatch.setattr(c.q,'event_passed',lambda *a:True)
    monkeypatch.setattr(c.q,'spectrum',lambda *a:np.zeros(16))
    def observation(sim,*a):
        return dict(time=sim.time+(1 if fault=='observation' else 0),volumes=[1.]*16,delta=np.zeros((16,16)).tolist())
    record.update(volumes=[1.]*16,delta=np.zeros((16,16)).tolist())
    (folder/'gpu-history.json').write_text(json.dumps([record]));original['gpu-history.json']=(folder/'gpu-history.json').read_bytes()
    monkeypatch.setattr(c,'row',observation)
    def compare(*a):
        if fault=='checkpoint':raise RuntimeError('Development restart changed phase_carry')
    monkeypatch.setattr(c,'exact',compare)
    monkeypatch.setattr(torch.cuda,'empty_cache',lambda:None)
    job=dict(key='full_space120',config=asdict(cfg))
    p=dict(physical_limits=c.q.PHYSICAL,checkpoint_interval=3.)
    (tmp_path/'protocol.json').write_text(json.dumps(p))
    result=c.replay_live_state(tmp_path,job,2,p)
    assert result['passed']==(fault=='none')
    assert result['complete_state_bitwise_equal']==(fault=='none')
    if fault!='none':assert result['failure']['type']=='RuntimeError'
    assert {f.name:f.read_bytes() for f in folder.iterdir()}==original
