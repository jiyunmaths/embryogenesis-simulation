from dataclasses import asdict
import os
from pathlib import Path

import numpy as np
import pytest
import torch

from embryo import polarity_conductance_controls as study
from embryo.attribute_development import AttributeSimulation
from embryo.gpu_initiation_controls import InitiationControlSimulation
from embryo.gpu_precision_control import PrecisionSimulation
from embryo.resolution import write_json
from test_gpu_backend import source
from test_polarity_robustness import protocol, trace, small_source


def design():
    return study.job_design({s:dict(coarse=f'{s}-c',fine=f'{s}-f',chemical=f'{s}-chem',transport=f'{s}-g') for s in study.HISTORIES})


def test_balanced_matched_design_and_eight_missing_paths():
    jobs=design()
    assert len(jobs)==24 and len({j['key'] for j in jobs})==24
    assert jobs[0]['seed']==9
    for seed in (7,8,9):
        assert {(j['point']['chi'],j['arm'],j['dt']) for j in jobs if j['seed']==seed}=={
            (chi,arm,dt) for chi in (0.,.35) for arm in study.ARMS for dt in (.00375,.001875)}
    new=[j for j in jobs if j['arm']=='fixed-conductances' and (j['point']['chi']==.35 or j['seed']==9)]
    assert len(new)==8 and all(j['family']=='uniform' and j['point']['ratio']==27.5 for j in jobs)


@pytest.mark.parametrize('values,expected',[
    ((True,False,True,True),'conductance_preservation_removes_observed_polarity_suppression'),
    ((False,False,True,True),'fixed_conductances_enable_both_contrasts_without_original_polarity_specific_loss'),
    ((True,False,True,False),'polarity_suppression_persists_with_fixed_conductances'),
    ((False,False,True,False),'polarity_suppression_persists_with_fixed_conductances'),
    ((True,False,False,False),'neither_fixed_conductance_branch_forms_no_successful_control'),
    ((True,False,False,True),'positive_contrast_only_forms_with_fixed_conductances'),
    ((True,True,True,True),'both_fixed_branches_form_positive_baseline_already_forms')])
def test_interpretation_separates_original_suppression_from_context_rescue(values,expected):
    forms=dict(zip(('chi-0_baseline','chi-0.35_baseline','chi-0_fixed-conductances','chi-0.35_fixed-conductances'),values))
    r=study.interpretation(forms,True)
    assert r['conclusion']==expected
    assert (r['original_polarity_comparison']=='suppression_supported_in_this_history')==(values[0] and not values[1])
    assert r['positive_contrast_conductance_rescue']==(not values[1] and values[3])
    assert study.interpretation(forms,False)==dict(qualified=False,conclusion='unresolved_numerical_or_physical_checks')


def test_reuse_rejects_wrong_polarity_or_contact_intervention():
    job=design()[0];old=dict(job,arm='phase_carry')
    study.check_reuse_context(job,old,'phase')
    bad=dict(old,point=dict(old['point'],chi=.35))
    with pytest.raises(ValueError,match='matched'):study.check_reuse_context(job,bad,'phase')
    with pytest.raises(ValueError,match='matched'):study.check_reuse_context(dict(job,arm='fixed-conductances'),old,'phase')
    fixed=dict(job,arm='fixed-conductances')
    study.check_reuse_context(fixed,fixed,'contact')
    with pytest.raises(ValueError,match='matched'):study.check_reuse_context(fixed,dict(fixed,arm='combined'),'contact')


def test_restart_parameters_allow_history_seed_but_not_wrong_mechanics():
    job=design()[0];job=dict(job,point=dict(job['point'],chi=.35))
    accepted=dict(seed=7,dt=.00375,signal_da=.02,signal_dh=.55,polarity_tension=.35,feedback=True)
    p=dict(accepted_configs=dict(coarse=accepted))
    correct=dict(accepted,seed=9)
    study.check_config(correct,job,p)
    for key,value in (('polarity_tension',0.),('seed',8),('signal_da',.03),('feedback',False)):
        with pytest.raises(ValueError,match='parameters'):study.check_config(dict(correct,**{key:value}),job,p)


def test_failed_pilot_blocks_only_its_own_polarity_contact_pair(tmp_path,monkeypatch):
    jobs=design();reused={j['key']:{} for j in jobs if j['arm']=='baseline' or (j['seed'] in (7,8) and j['point']['chi']==0.)}
    p=dict(jobs=jobs,histories=list(study.HISTORIES),contrasts=list(study.CONTRASTS),arms=list(study.ARMS),reused=reused,pilot_duration=60.,duration=240.,new_moving_jobs=8)
    calls=[]
    def advance(root,job,p,horizon):
        calls.append((job['seed'],job['point']['chi'],job['arm'],horizon));return [dict(elapsed=horizon)]
    monkeypatch.setattr(study,'advance',advance)
    monkeypatch.setattr(study,'pair_report',lambda root,p,job,horizon,h:dict(passed=not(job['seed']==9 and job['point']['chi']==.35 and job['arm']=='fixed-conductances')))
    class Executor:
        def submit(self,*args):return object()
    done,pending,failed=study.run_pairs(tmp_path,p,Executor())
    assert done==6 and len(pending)==6 and failed==[dict(seed=9,chi=.35,arm='fixed-conductances')]
    assert all(h==60 for s,c,a,h in calls if (s,c,a)==(9,.35,'fixed-conductances'))
    assert any(h==240 for s,c,a,h in calls if (s,c,a)==(7,.35,'fixed-conductances'))


def test_pair_evidence_distinguishes_contrasts_and_stays_immutable(tmp_path):
    p=protocol();write_json(tmp_path/'protocol.json',p)
    h={level:[dict(r,geometric_delta=r['delta']) for r in trace()] for level in ('coarse','fine')}
    job=design()[0]
    zero=study.pair_report(tmp_path,p,job,.3,h)
    positive=study.pair_report(tmp_path,p,dict(job,point=dict(job['point'],chi=.35)),.3,h)
    assert zero['passed'] and positive['passed'] and zero['history_sha256']!=positive['history_sha256']
    write_json(Path(next(iter(positive['history_sha256']))),[])
    with pytest.raises(ValueError,match='immutable'):study.pair_report(tmp_path,p,dict(job,point=dict(job['point'],chi=.35)),.3,h)


def test_restart_from_initial_checkpoint_retains_carry_and_matches_continuation(tmp_path,monkeypatch):
    host,cp,chemical=small_source(tmp_path)
    host.config.signal_da=.02;host.config.signal_dh=.55;host.config.polarity_tension=.35
    host.activator=np.ones(2);host.checkpoint(cp)
    graph=host.signaling_graph();g=graph.delta*graph.masses[:,None];np.fill_diagonal(g,0.)
    transport=tmp_path/'transport.npz';np.savez(transport,conductance=g,volumes=graph.masses,ids=host.ids)
    job=dict(key='restart',seed=host.config.seed,level='coarse',dt=host.config.dt,arm='fixed-conductances',family='uniform',point=dict(ratio=27.5,chi=.35),source=str(cp),chemical_file=str(chemical),initial_transport=str(transport))
    p=dict(accepted_configs=dict(coarse=asdict(host.config)),start=150.,duration=.45,pilot_duration=.15,
           interval=.15,checkpoint_interval=.15,late_window=.15,criteria=dict(boundary_max=.01,late_log_sd_min=.1),amount_accounting_max=2e-14,reused={})
    write_json(tmp_path/'protocol.json',p);fail=[40000]
    class Backend:
        # This deliberately simple dynamics checks restart orchestration only.
        def __init__(self,host,arm):
            self.host=host;self.initiation_arm=arm
            self.initial_conductance=torch.from_numpy(g.copy())
            self.phase_carry=np.full(host.phi.shape,.25);self.rounding=np.zeros(2,dtype=np.int32)
        def __getattr__(self,k):return getattr(self.host,k)
        def audit(self):return dict(max_volume_error=0.,min_radius=5.,max_clipping=0.,volume_conversion_amount_error=0.)
        def step(self):
            if self.host.step_number==fail[0]:fail[0]=None;raise RuntimeError('interrupted after initial checkpoint or later checkpoint')
            self.host.activator+=.001*(1+float(self.phase_carry.flat[0]));self.phase_carry+=.1
            self.host.step_number+=1;self.host.time=self.host.step_number*self.host.config.dt
        def observe(self,elapsed):
            graph=self.host.signaling_graph()
            return dict(elapsed=elapsed,time=self.host.time,ids=self.host.ids.tolist(),chemistry=np.array([self.host.activator,self.host.inhibitor]).tolist(),
                        polarity=self.host.polarity.tolist(),volumes=graph.masses.tolist(),delta=graph.delta.tolist(),geometric_delta=graph.delta.tolist(),
                        chemical_conductance=g.tolist(),geometric_conductance=g.tolist(),initiation_arm=self.initiation_arm,
                        cumulative_volume_amount_source=np.zeros((2,2)).tolist(),boundary_occupancy=0.,axis_ratio=1.2,log_activator_sd=0.)
        def precision_diagnostics(self):return dict(marker=float(self.phase_carry.flat[0]))
        def checkpoint(self,path):
            self.host.checkpoint(path)
            with np.load(path) as z:payload={k:z[k].copy() for k in z.files}
            payload.update(precision_arm=np.array('phase_carry'),initiation_arm=np.array(self.initiation_arm),phase_carry=self.phase_carry,
                           rounding=self.rounding,initial_conductance=g,initial_volume=graph.masses,amount_source=np.zeros((2,2)),
                           last_amount_source=np.zeros((2,2)),last_relative_volume_change=np.zeros(2),conversion_error=np.array(0.))
            np.savez_compressed(path,**payload)
        @classmethod
        def restore(cls,path):
            sim=cls(AttributeSimulation.restore(path),'fixed-conductances')
            with np.load(path) as z:sim.phase_carry=z['phase_carry'].copy();sim.rounding=z['rounding'].copy()
            return sim
    monkeypatch.setattr(study,'InitiationControlSimulation',Backend)
    monkeypatch.setattr(study,'initialize',lambda *args:AttributeSimulation.restore(cp))
    monkeypatch.setattr(study.torch.cuda,'empty_cache',lambda:None)
    with pytest.raises(RuntimeError,match='initial checkpoint'):study.advance(tmp_path,job,p,.15)
    initial=Backend.restore(tmp_path/'restart/latest_state.npz')
    assert initial.time==150. and initial.phase_carry.flat[0]==.25
    pilot=study.advance(tmp_path,job,p,.15)
    assert len(pilot)==2 and not (tmp_path/'restart/result.json').exists()
    fail[0]=40080
    with pytest.raises(RuntimeError,match='later checkpoint'):study.advance(tmp_path,job,p,.45)
    saved_hash=study.digest(tmp_path/'restart/latest_state.npz')
    assert study.advance(tmp_path,job,p,.15)==pilot
    assert study.digest(tmp_path/'restart/latest_state.npz')==saved_hash
    final=study.advance(tmp_path,job,p,.45)
    expected=np.ones(2)+.001*sum(1.25+.1*k for k in range(120))
    np.testing.assert_allclose(final[-1]['chemistry'][0],expected,rtol=1e-12)
    assert study.advance(tmp_path,job,p,.45)==final


gpu=pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU opt-in required')


@gpu
@pytest.mark.parametrize('chi',study.CONTRASTS)
@pytest.mark.parametrize('arm',study.ARMS)
def test_polarity_conductance_step_matches_independent_chemistry_and_restarts(tmp_path,chi,arm):
    def make():
        h=source();h.config.polarity_tension=chi;h.config.signal_da=.02;h.config.signal_dh=.55;return h
    sim=InitiationControlSimulation(make(),arm)
    original=PrecisionSimulation(make(),'phase_carry') if arm=='baseline' else None
    for _ in range(8):
        _,d,_,g=sim.matrices();delta,_=sim.chemical_transport(d,g)
        before=np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()]);old=sim.geometry[:,0].cpu().numpy().copy()
        expected=study.contact.cpu_chemical_step(before,delta.cpu().numpy(),sim.config.dt,2.,.02,.55)
        sim.step();expected*=old/sim.geometry[:,0].cpu().numpy()
        np.testing.assert_allclose(np.array([sim.activator.cpu().numpy(),sim.inhibitor.cpu().numpy()]),expected,rtol=1e-13,atol=1e-14)
        if original is not None:original.step();study.contact.exact_state(sim,original)
        assert not torch.any(sim.amount_source) and sim.audit()['volume_conversion_amount_error']<2e-14
    path=tmp_path/'restart.npz';sim.checkpoint(path);restored=InitiationControlSimulation.restore(path)
    for _ in range(3):sim.step();restored.step();study.contact.exact_state(sim,restored,True)
