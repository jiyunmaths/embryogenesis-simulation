from dataclasses import asdict
import os
from pathlib import Path

import numpy as np
import pytest

from embryo import phase_carry_maintenance as study
from embryo.attribute_development import AttributeSimulation
from embryo.resolution import write_json
from test_polarity_robustness import small_source, trace, protocol


def design():
    return study.job_design({s:dict(coarse=f'{s}-c',fine=f'{s}-f',chemical=f'{s}-chem') for s in (7,8,9)})


def test_design_keeps_maintenance_and_initiation_within_three_histories():
    jobs=design()
    assert len(jobs)==24 and len(set(j['key'] for j in jobs))==24
    for seed in (7,8,9):
        assert {(j['point']['chi'],j['dt'],j['family']) for j in jobs if j['seed']==seed}=={
            (chi,dt,family) for chi in (0.,.35) for dt in (.00375,.001875) for family in ('pattern','uniform')}
    assert all(j['arm']=='phase_carry' and j['point']['ratio']==27.5 for j in jobs)


def test_history_checks_actual_family_including_initial_checkpoint(tmp_path):
    _,_,chemical=small_source(tmp_path)
    h=trace();p=protocol();job=dict(family='pattern',arm='phase_carry',chemical_file=str(chemical))
    with np.load(chemical) as z:h[0]['chemistry']=z['pattern'].tolist()
    h[0]['log_activator_sd']=float(np.std(np.log(np.array(h[0]['chemistry'])[0])))
    study.validate_history(h,p,job,.3)
    study.validate_history(h[:1],p,job,0.)
    job['family']='uniform'
    with pytest.raises(ValueError,match='preparation'):study.validate_history(h,p,job,.3)
    job['family']='pattern';h[-1]['log_activator_sd']+=.01
    with pytest.raises(ValueError,match='contrast'):study.validate_history(h,p,job,.3)


def test_late_recovery_is_distinguished_from_continuous_retention():
    p=protocol();p['duration']=.3;p['late_window']=.0
    h=trace();h[0]['log_activator_sd']=.3;h[1]['log_activator_sd']=.05;h[2]['log_activator_sd']=.2
    r=study.maintenance_summary(h,p)
    assert r['persistent_contrast'] and not r['whole_window_retained']
    assert r['first_contrast_loss_elapsed']==pytest.approx(.12)
    assert r['whole_window_min_log_sd']==.05


def test_pair_report_gates_loss_time_and_preserves_immutable_evidence(tmp_path):
    p=protocol();p['loss_time_abs_max']=.01;write_json(tmp_path/'protocol.json',p)
    job=dict(seed=7,point=dict(chi=.35),family='pattern')
    h=trace();h[0]['log_activator_sd']=.3;h[1]['log_activator_sd']=.05
    # Keep chemistry consistent with the large start/weak final contrast.
    for row in h:row['chemistry'][0]=np.exp([-row['log_activator_sd'],row['log_activator_sd']]).tolist()
    r=study.pair_report(tmp_path,p,job,.3,dict(coarse=h,fine=h))
    assert r['passed'] and r['loss_times']==pytest.approx([.12,.12])
    assert 'history 7, pattern' in r['interpretation']
    assert study.pair_report(tmp_path,p,job,.3,dict(coarse=h,fine=h))==r
    write_json(Path(next(iter(r['history_sha256']))),[])
    with pytest.raises(ValueError,match='immutable'):study.pair_report(tmp_path,p,job,.3,dict(coarse=h,fine=h))


def test_unqualified_or_lost_patterns_are_not_accepted_maintenance():
    assert study.maintenance_comparison(True,True,True)=='both_contrasts_maintain'
    assert study.maintenance_comparison(True,False,True)=='positive_contrast_loses_maintenance'
    assert study.maintenance_comparison(True,True,False)=='unresolved_numerical_or_physical_checks'


def test_failed_maintenance_pilot_does_not_launch_its_long_pair(tmp_path,monkeypatch):
    p=dict(jobs=design(),histories=[9,7,8],contrasts=[0.,.35],pilot_duration=60.,duration=240.)
    calls=[]
    def advance(root,job,p,horizon):calls.append((job['seed'],job['point']['chi'],job['family'],horizon));return [dict(elapsed=horizon)]
    monkeypatch.setattr(study,'advance',advance)
    monkeypatch.setattr(study,'pair_report',lambda root,p,job,horizon,h:dict(passed=not(job['seed']==9 and job['point']['chi']==0.)))
    class Executor:
        def submit(self,*args):return object()
    done,pending,failed=study.run_pairs(tmp_path,p,Executor())
    assert done==10 and len(pending)==10 and failed==[dict(seed=9,chi=0.,family='pattern')]
    assert all(family=='pattern' for _,_,family,_ in calls)
    assert all(horizon==60. for seed,chi,_,horizon in calls if seed==9 and chi==0.)


def test_initial_checkpoint_interrupt_and_past_pilot_resume_preserve_residual(tmp_path,monkeypatch):
    host,source,chemical=small_source(tmp_path)
    job=dict(key='maintenance',seed=host.config.seed,level='coarse',dt=host.config.dt,arm='phase_carry',family='pattern',point=dict(key='polarity-0',ratio=27.5,chi=0.),source=str(source),chemical_file=str(chemical))
    p=dict(accepted_configs=dict(coarse=asdict(host.config)),start=150.,duration=.3,pilot_duration=.15,interval=.15,checkpoint_interval=.15,late_window=.15,criteria=dict(boundary_max=.01,dilution_error_max=2e-14,late_log_sd_min=.1),reused={})
    write_json(tmp_path/'protocol.json',p);folder=tmp_path/job['key'];folder.mkdir();fail=[False]
    class CarryBackend:
        def __init__(self,host,arm='phase_carry'):
            self.host=host;self.precision_arm=arm;self.phase_carry=np.full(host.phi.shape,.25,dtype=np.float64);self.rounding=np.zeros(2,dtype=np.int32)
        def __getattr__(self,key):return getattr(self.host,key)
        def audit(self):return dict(max_volume_error=0.,min_radius=5.,max_clipping=0.,dilution_amount_error=0.)
        def step(self):
            if fail[0] and self.host.step_number==40045:fail[0]=False;raise RuntimeError('interrupted')
            self.host.activator+=.001*(1+float(self.phase_carry.flat[0]));self.phase_carry+=.1;self.rounding[0]+=1;self.host.step_number+=1;self.host.time=self.host.step_number*self.host.config.dt
        def observe(self,elapsed):
            g=self.host.signaling_graph()
            return dict(elapsed=elapsed,time=self.host.time,ids=self.host.ids.tolist(),chemistry=np.array([self.host.activator,self.host.inhibitor]).tolist(),polarity=self.host.polarity.tolist(),volumes=g.masses.tolist(),delta=g.delta.tolist(),boundary_occupancy=0.,axis_ratio=1.2,log_activator_sd=float(np.std(np.log(self.host.activator))))
        def precision_diagnostics(self):return dict(marker=float(self.phase_carry.flat[0]))
        def checkpoint(self,path):
            self.host.checkpoint(path)
            with np.load(path) as z:payload={k:z[k].copy() for k in z.files}
            payload.update(precision_arm=np.array(self.precision_arm),phase_carry=self.phase_carry,rounding=self.rounding);np.savez_compressed(path,**payload)
        @classmethod
        def restore(cls,path):
            sim=cls(AttributeSimulation.restore(path))
            with np.load(path) as z:sim.phase_carry=z['phase_carry'].copy();sim.rounding=z['rounding'].copy()
            return sim
    monkeypatch.setattr(study,'PrecisionSimulation',CarryBackend);monkeypatch.setattr(study,'context_gate',lambda *args:dict(passed=True));monkeypatch.setattr(study.torch.cuda,'empty_cache',lambda:None)
    initial=CarryBackend(study.initialize(source,chemical,'pattern',job['point'],p['accepted_configs']['coarse']))
    h=[study.retain(initial.observe(0.),job)]
    audit=dict(max_volume_error=0.,min_radius=5.,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    study.save_checkpoint(initial,folder/'latest_state.npz',audit,h,study.digest(tmp_path/'protocol.json'),job)
    pilot=study.advance(tmp_path,job,p,.15);assert len(pilot)==2 and not (folder/'result.json').exists()
    fail[0]=True
    with pytest.raises(RuntimeError,match='interrupted'):study.advance(tmp_path,job,p,.3)
    final=study.advance(tmp_path,job,p,.3);saved=CarryBackend.restore(folder/'latest_state.npz')
    np.testing.assert_allclose(saved.activator,host.activator+.001*sum(1.25+.1*k for k in range(80)),rtol=1e-12)
    assert saved.phase_carry.flat[0]==pytest.approx(8.25)
    hashes=[study.digest(folder/f) for f in ('latest_state.npz','history.json','result.json')]
    assert study.advance(tmp_path,job,p,.15)==final[:2]
    assert study.advance(tmp_path,job,p,.3)==final
    assert hashes==[study.digest(folder/f) for f in ('latest_state.npz','history.json','result.json')]


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU opt-in required')
@pytest.mark.parametrize('chi,dt',[(0.,.00375),(.35,.00375),(0.,.001875),(.35,.001875)])
def test_developed_state_actual_context_gate_on_gpu(tmp_path,chi,dt):
    host,source,chemical=small_source(tmp_path)
    # Give the synthetic cells empty space while preserving the grid spacing;
    # the smaller CPU bookkeeping fixture fails the physical boundary screen.
    expanded=AttributeSimulation(type(host.config)(**{**asdict(host.config),'grid':24,'extent':1.2}),mode='direct')
    partition=.5*(1+np.tanh(expanded.xyz[0]/expanded.config.interface_width))
    expanded.phi=np.array([expanded.phi[0]*partition,expanded.phi[0]*(1-partition)],dtype=np.float32)
    expanded.target=expanded.volumes()
    for key in ('ids','parents','due','fate','activator','inhibitor','polarity'):
        setattr(expanded,key,getattr(host,key).copy())
    expanded.step_number=host.step_number;expanded.time=host.time
    expanded.checkpoint(source);host=expanded
    if dt!=host.config.dt:
        from embryo.feedback_survival_validation import retime
        source=tmp_path/'fine-source.npz';host=retime(tmp_path/'source.npz',source,dt)
    job=dict(key='context',seed=host.config.seed,level='coarse',dt=dt,arm='phase_carry',family='pattern',point=dict(key=f'polarity-{chi:g}',ratio=27.5,chi=chi),source=str(source),chemical_file=str(chemical))
    p=dict(accepted_configs=dict(coarse=asdict(host.config)),interval=.00375,prefix_duration=.015,prefix_native_threads=2,prefix_criteria=study.GPU_CRITERIA,chemical_reference_log_max=1e-11,criteria=dict(boundary_max=.01,dilution_error_max=2e-14))
    r=study.actual_context_gate(tmp_path/'gate',job,p)
    assert r['passed'] and r['cpu_chemical_log_max']<1e-11 and r['carry_restart_exact_steps']==4
