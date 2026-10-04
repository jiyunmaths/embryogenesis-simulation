from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

from embryo import phase_carry_polarity as study
from embryo.attribute_development import AttributeSimulation
from embryo.neighbor_context import read
from embryo.resolution import write_json
from test_polarity_robustness import small_source, trace, protocol


def design():
    return study.job_design({s:dict(coarse=f'{s}-c',fine=f'{s}-f',chemical=f'{s}-chem') for s in (7,8,9)})


def test_design_balances_histories_and_changes_only_directional_tension_and_dt():
    jobs=design()
    assert len(jobs)==12 and len(set(j['key'] for j in jobs))==12
    for seed in (7,8,9):
        assert {(j['point']['chi'],j['dt']) for j in jobs if j['seed']==seed} == {
            (0.,.00375),(0.,.001875),(.35,.00375),(.35,.001875)}
    assert all(j['family']=='uniform' and j['arm']=='phase_carry' and j['point']['ratio']==27.5 for j in jobs)


def test_uniform_branches_and_numerical_failures_are_not_suppression_evidence():
    assert study.initiation_comparison(True,False,True)=='suppression_supported_in_this_history'
    assert study.initiation_comparison(False,False,True)=='neither_initiates_no_successful_control'
    assert study.initiation_comparison(True,True,True)=='both_initiate'
    assert study.initiation_comparison(False,True,True)=='positive_contrast_only_initiates'
    assert study.initiation_comparison(True,False,False)=='unresolved_numerical_agreement'


def test_failed_pilot_blocks_only_its_own_new_long_continuations(tmp_path,monkeypatch):
    p=dict(jobs=design(),contrasts=[0.,.35],reused={},pilot_duration=60.,duration=240.)
    calls=[]
    def advance(root,job,p,horizon):
        calls.append((job['seed'],job['point']['chi'],horizon));return [dict(elapsed=horizon)]
    monkeypatch.setattr(study,'advance',advance)
    monkeypatch.setattr(study,'pair_report',lambda root,p,seed,chi,horizon,histories:dict(passed=not(seed==9 and chi==0.)))
    class Executor:
        def submit(self,*args):return object()
    done,pending,failed=study.run_pairs(tmp_path,p,Executor())
    assert failed==[dict(seed=9,chi=0.)]
    assert done==10 and len(pending)==10
    assert all(horizon==60. for seed,chi,horizon in calls if seed==9 and chi==0.)
    assert any(horizon==240. for seed,chi,horizon in calls if seed==7)


def test_pair_report_freezes_pilot_evidence_and_rejects_tampering(tmp_path):
    p=protocol();write_json(tmp_path/'protocol.json',p)
    histories=dict(coarse=trace(),fine=trace())
    r=study.pair_report(tmp_path,p,9,.35,p['pilot_duration'],histories)
    assert r['passed']
    saved={f:study.digest(f) for f in r['history_sha256']}
    assert study.pair_report(tmp_path,p,9,.35,p['pilot_duration'],histories)==r
    assert saved=={f:study.digest(f) for f in saved}
    write_json(Path(next(iter(saved))),[])
    with pytest.raises(ValueError,match='immutable'):study.pair_report(tmp_path,p,9,.35,p['pilot_duration'],histories)


def test_pilot_resume_retains_a_numerical_state_that_affects_future_dynamics(tmp_path,monkeypatch):
    host,source,chemical=small_source(tmp_path)
    job=dict(key='pilot',seed=9,level='coarse',dt=host.config.dt,arm='phase_carry',family='uniform',
        point=dict(key='polarity-0',ratio=27.5,chi=0.),source=str(source),chemical_file=str(chemical))
    p=dict(accepted_configs=dict(coarse=asdict(host.config)),start=150.,duration=.3,pilot_duration=.15,
        interval=.15,checkpoint_interval=.15,criteria=dict(boundary_max=.01,dilution_error_max=2e-14),reused={})
    write_json(tmp_path/'protocol.json',p);(tmp_path/'pilot').mkdir();fail=[False]
    class CarryBackend:
        # Deliberately simple dynamics depend on a checkpointed residual marker.
        # This checks resumption, not the numerical accuracy of the real kernel.
        def __init__(self,host,arm='phase_carry'):
            self.host=host;self.precision_arm=arm
            self.phase_carry=np.full(host.phi.shape,.25,dtype=np.float64)
            self.rounding=np.zeros(2,dtype=np.int32)
        def __getattr__(self,key):return getattr(self.host,key)
        def audit(self):return dict(max_volume_error=0.,min_radius=5.,max_clipping=0.,dilution_amount_error=0.)
        def step(self):
            if fail[0] and self.host.step_number==40045:
                fail[0]=False;raise RuntimeError('interrupted')
            self.host.activator+=.001*(1+float(self.phase_carry.flat[0]))
            self.phase_carry+=.1;self.rounding[0]+=1
            self.host.step_number+=1;self.host.time=self.host.step_number*self.host.config.dt
        def observe(self,elapsed):
            graph=self.host.signaling_graph()
            return dict(elapsed=elapsed,time=self.host.time,ids=self.host.ids.tolist(),
                chemistry=np.array([self.host.activator,self.host.inhibitor]).tolist(),polarity=self.host.polarity.tolist(),
                volumes=graph.masses.tolist(),delta=graph.delta.tolist(),boundary_occupancy=0.,axis_ratio=1.2,
                log_activator_sd=float(np.std(np.log(self.host.activator))))
        def precision_diagnostics(self):return dict(residual_marker=float(self.phase_carry.flat[0]))
        def checkpoint(self,path):
            self.host.checkpoint(path)
            with np.load(path) as z:payload={k:z[k].copy() for k in z.files}
            payload.update(precision_arm=np.array(self.precision_arm),phase_carry=self.phase_carry,rounding=self.rounding)
            np.savez_compressed(path,**payload)
        @classmethod
        def restore(cls,path):
            sim=cls(AttributeSimulation.restore(path))
            with np.load(path) as z:sim.phase_carry=z['phase_carry'].copy();sim.rounding=z['rounding'].copy()
            return sim
    monkeypatch.setattr(study,'PrecisionSimulation',CarryBackend)
    monkeypatch.setattr(study,'context_gate',lambda *args:dict(passed=True))
    monkeypatch.setattr(study.torch.cuda,'empty_cache',lambda:None)
    pilot=study.advance(tmp_path,job,p,.15)
    assert len(pilot)==2 and not (tmp_path/'pilot/result.json').exists()
    before=CarryBackend.restore(tmp_path/'pilot/latest_state.npz')
    assert before.phase_carry.flat[0]==pytest.approx(4.25)
    fail[0]=True
    with pytest.raises(RuntimeError,match='interrupted'):study.advance(tmp_path,job,p,.3)
    unchanged=CarryBackend.restore(tmp_path/'pilot/latest_state.npz')
    np.testing.assert_array_equal(unchanged.phase_carry,before.phase_carry)
    final=study.advance(tmp_path,job,p,.3)
    saved=CarryBackend.restore(tmp_path/'pilot/latest_state.npz')
    np.testing.assert_allclose(saved.activator,np.ones(2)+.001*sum(1.25+.1*k for k in range(80)),rtol=1e-12)
    assert saved.phase_carry.flat[0]==pytest.approx(8.25) and len(final)==3
    hashes=[study.digest(tmp_path/'pilot'/f) for f in ('result.json','history.json','latest_state.npz')]
    assert study.advance(tmp_path,job,p,.3)==final
    assert hashes==[study.digest(tmp_path/'pilot'/f) for f in ('result.json','history.json','latest_state.npz')]
