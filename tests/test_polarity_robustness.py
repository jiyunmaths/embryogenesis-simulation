from dataclasses import asdict

import numpy as np
import pytest

from embryo import polarity_robustness as study
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.neighbor_context import read
from embryo.resolution import write_json


def trace(spread=1e-8):
    return [dict(elapsed=t, time=150+t, ids=[3, 4],
        chemistry=[np.exp([-s, s]).tolist(), [1., 1.]],
        volumes=[1., 1.], delta=[[-1., 1.], [1., -1.]],
        polarity=[[1., 0., 0.], [-1., 0., 0.]], axis_ratio=1.2,
        uniform_growth_max=.01-t/15, log_activator_sd=s)
        for t,s in zip([0., .15, .3], [1e-8, spread, spread])]


def protocol():
    return dict(start=150., interval=.15, pilot_duration=.3, duration=.6,
        pilot_late_window=.15, late_window=.3,
        criteria=dict(late_log_sd_min=.1), refinement_criteria=study.REFINEMENT_CRITERIA)


def test_crossings_interpolate_on_physical_clock_and_require_ordered_finite_data():
    assert study.first_crossing([0., 2., 4.], [.3, .1, -.1]) == pytest.approx(3.)
    assert study.first_crossing([0., 2.], [.05, .15], .1, 'up') == pytest.approx(1.)
    assert study.first_crossing([0., 1.], [1., 2.]) is None
    for times,values in [([1., 0.], [0., 1.]), ([0., 1.], [0., np.nan])]:
        with pytest.raises(ValueError): study.first_crossing(times,values)


def test_near_uniform_refinement_uses_raw_chemical_error_without_sd_floor_division():
    result=study.refinement_comparison(trace(),trace(2e-8),protocol(),.3)
    assert result['passed']
    assert result['errors']['chemical_log_max'] < 2e-8
    assert result['persistent_contrast']==[False,False]
    assert result['crossing_times']==pytest.approx([.15,.15])


def test_refinement_rejects_missing_growth_crossing_and_different_late_outcomes():
    fine=trace()
    for row in fine: row['uniform_growth_max']=.01
    r=study.refinement_comparison(trace(),fine,protocol(),.3)
    assert not r['passed'] and r['crossing_time_error'] is None
    fine=trace(.2)
    r=study.refinement_comparison(trace(),fine,protocol(),.3)
    assert not r['passed'] and r['persistent_contrast']==[False,True]
    with pytest.raises(ValueError,match='aligned'):
        study.refinement_comparison(trace()[:-1],fine,protocol(),.3)


def small_source(tmp_path):
    c=Config(grid=16,extent=.8,max_cells=2,dt=.00375,steps=40000,
        save_every=40,differentiation=False,signal_transport='conservative')
    host=AttributeSimulation(c,mode='direct')
    partition=.5*(1+np.tanh(host.xyz[0]/c.interface_width))
    host.phi=np.array([host.phi[0]*partition,host.phi[0]*(1-partition)],dtype=np.float32)
    host.target=host.volumes();host.ids=np.array([3,4]);host.parents=np.array([1,1])
    host.due=np.full(2,np.inf);host.fate=np.zeros(2)
    host.activator=np.array([.2,1.8]);host.inhibitor=np.ones(2)
    host.polarity=np.array([[1.,0.,0.],[-1.,0.,0.]])
    host.step_number=40000;host.time=150.
    cp=tmp_path/'source.npz';host.checkpoint(cp)
    chemical=tmp_path/'chemical.npz'
    np.savez(chemical,ids=host.ids,pattern=np.array([host.activator,host.inhibitor]),uniform=np.ones((2,2)))
    return host,cp,chemical


def test_retiming_changes_only_discrete_clock_and_preserves_physical_state_and_rng(tmp_path):
    host,cp,_=small_source(tmp_path)
    fine=study.retime(cp,tmp_path/'fine.npz',host.config.dt/2)
    study.assert_same_physical_start(host,fine)
    assert fine.step_number==2*host.step_number
    fine.activator[0]+=1e-5
    with pytest.raises(ValueError,match='activator'):study.assert_same_physical_start(host,fine)


def test_pilot_stops_without_final_result_and_resume_never_resets_chemistry(tmp_path,monkeypatch):
    from embryo import gpu_backend
    host,cp,chemical=small_source(tmp_path)
    job=dict(key='pilot',seed=9,level='coarse',dt=host.config.dt,source=str(cp),
        chemical_file=str(chemical),family='pattern',point=dict(ratio=27.5,chi=.7))
    p=dict(accepted_configs=dict(coarse=asdict(host.config)),start=150.,duration=.3,
        pilot_duration=.15,interval=.15,checkpoint_interval=.15,late_window=.15,
        criteria=dict(boundary_max=.01,dilution_error_max=2e-14,late_log_sd_min=.1))
    write_json(tmp_path/'protocol.json',p);fail=[False]
    class CheckpointBackend:
        def __init__(self,host):self.host=host
        def __getattr__(self,key):return getattr(self.host,key)
        def to_cpu(self):return self.host
        def audit(self):return dict(max_volume_error=0.,min_radius=5.,max_clipping=0.,dilution_amount_error=0.)
        def step(self):
            if fail[0] and self.host.step_number==40045:
                fail[0]=False;raise RuntimeError('simulated interruption')
            self.host.activator+=1e-5;self.host.step_number+=1
            self.host.time=self.host.step_number*self.host.config.dt
        def observe(self,elapsed):
            graph=self.host.signaling_graph()
            return dict(elapsed=elapsed,time=self.host.time,ids=self.host.ids.tolist(),
                chemistry=np.array([self.host.activator,self.host.inhibitor]).tolist(),
                volumes=graph.masses.tolist(),delta=graph.delta.tolist(),boundary_occupancy=0.,
                log_activator_sd=float(np.std(np.log(self.host.activator))))
    monkeypatch.setattr(gpu_backend,'GpuSimulation',CheckpointBackend)
    monkeypatch.setattr(study,'prefix',lambda *args:dict(passed=True))
    monkeypatch.setattr(study,'endpoint_assay',lambda *args:dict(numerical_pass=True))
    pilot=study.advance(tmp_path,job,p,.15)
    assert len(pilot)==2 and not (tmp_path/'pilot/result.json').exists()
    assert read(tmp_path/'pilot/status.json')['state']=='awaiting_timestep_gate'
    fail[0]=True
    with pytest.raises(RuntimeError,match='interruption'):study.advance(tmp_path,job,p,.3)
    saved=AttributeSimulation.restore(tmp_path/'pilot/latest_state.npz')
    assert saved.time==pytest.approx(150.15)
    final=study.advance(tmp_path,job,p,.3)
    saved=AttributeSimulation.restore(tmp_path/'pilot/latest_state.npz')
    np.testing.assert_allclose(saved.activator,host.activator+80e-5,rtol=1e-12,atol=1e-14)
    assert saved.time==pytest.approx(150.3) and len(final)==3
    hashes=[study.digest(tmp_path/'pilot'/f) for f in ('latest_state.npz','history.json','result.json')]
    assert study.advance(tmp_path,job,p,.3)==final
    assert hashes==[study.digest(tmp_path/'pilot'/f) for f in ('latest_state.npz','history.json','result.json')]


def test_failed_pilot_gate_blocks_all_long_jobs(tmp_path,monkeypatch):
    p=dict(contrasts=[0.,.35,.7],refinement_seed=9,pilot_duration=60.,duration=240.,
        jobs=[dict(key=f'{chi}-{level}',seed=9,family='uniform',level=level,point=dict(chi=chi))
            for chi in (0.,.35,.7) for level in ('coarse','fine')])
    called=[]
    monkeypatch.setattr(study,'frozen_references',lambda *args:None)
    monkeypatch.setattr(study,'advance',lambda root,job,p,horizon:called.append((job['key'],horizon)))
    monkeypatch.setattr(study,'compare_level',lambda *args:dict(passed=False))
    assert not study.run_stages(tmp_path,p)['passed']
    assert len(called)==6 and all(horizon==60. for _,horizon in called)
    assert read(tmp_path/'status.json')['state']=='stopped_at_numerical_gate'


def test_cached_pilot_hashes_are_immutable_when_long_history_grows(tmp_path,monkeypatch):
    p=protocol();p.update(contrasts=[0.],refinement_seed=9)
    host,_,chemical=small_source(tmp_path)
    # Build matched evidence; checkpoint decoding is mocked so this test focuses on provenance.
    p['jobs']=[dict(key=level,seed=9,family='uniform',level=level,point=dict(chi=0.),chemical_file=str(chemical))
        for level in ('coarse','fine')]
    histories={};write_json(tmp_path/'protocol.json',p)
    with np.load(chemical) as z:initial=z['uniform'].tolist()
    for job in p['jobs']:
        folder=tmp_path/job['key'];folder.mkdir();h=trace()
        h[0]['chemistry']=initial;histories[job['key']]=h
        write_json(folder/'history.json',h)
    host.time=150.3
    monkeypatch.setattr(study,'restore_checkpoint',lambda path,*args:(host,{},histories[path.parent.name]))
    report=study.compare_level(tmp_path,p,.3)
    assert report['passed']
    report_hash=study.digest(tmp_path/'pilot-refinement.json')
    for job in p['jobs']:
        write_json(tmp_path/job['key']/'history.json',histories[job['key']]+[dict(elapsed=.45)])
    assert study.compare_level(tmp_path,p,.3)==report
    assert study.digest(tmp_path/'pilot-refinement.json')==report_hash
    write_json(tmp_path/'fine/pilot-history.json',[])
    with pytest.raises(ValueError,match='Changed refinement evidence'):study.compare_level(tmp_path,p,.3)
