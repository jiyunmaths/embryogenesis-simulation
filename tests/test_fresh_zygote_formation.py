"""Protect inherited starts, sustained-onset semantics and the real long loop."""
from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from embryo import fresh_zygote_formation as s
from embryo.model import Config


def rules():
    return dict(thresholds=[.05, .1, .2], sustained_window=24., late_window=24.,
        physical_limits=s.f.q.PHYSICAL, interval=.15, diagnostic_interval=.015,
        snapshots=list(s.SNAPSHOTS), restart_elapsed=list(s.RESTARTS), restart_steps=8)


def history(values):
    return [dict(elapsed=i * .15, log_activator_sd=float(v)) for i, v in enumerate(values)]


def test_sustained_onset_requires_entire_interval_and_reports_recovery():
    v = np.zeros(401)
    v[10:50] = .3  # Six-unit transient cannot count as sustained.
    v[100:] = .3
    result = s.formation(history(v), rules())
    for r in result:
        assert r['first_crossing_elapsed'] == 1.5
        assert r['first_sustained_onset_elapsed'] == 15.
        assert r['late_persistent'] and r['sustained_formation']
        assert [x['direction'] for x in r['sampled_transitions']] == ['up', 'down', 'up']


def test_late_crossing_and_uniform_history_are_valid_nonformation_outcomes():
    v = np.zeros(401)
    v[300:] = .3
    for r in s.formation(history(v), rules()):
        assert r['first_crossing_elapsed'] == 45.
        assert not r['sustained_formation'] and not r['late_persistent']
    assert all(not r['sustained_formation'] for r in s.formation(history(np.zeros(401)), rules()))


def test_threshold_is_strict_and_windows_include_both_endpoints():
    v = np.ones(401) * .1
    result = s.formation(history(v), rules())
    assert result[0]['first_sustained_onset_elapsed'] == 0
    assert not result[1]['sustained_formation'] and not result[2]['late_persistent']
    v[:] = .3
    v[240] = 0  # Exactly elapsed36, first endpoint of late window.
    assert not s.formation(history(v), rules())[1]['late_persistent']
    assert s.formation(history(v), rules())[1]['sustained_formation']


@pytest.mark.parametrize('bad', [np.zeros(400), np.r_[np.zeros(400), np.nan]])
def test_incomplete_or_nonfinite_formation_history_rejected(bad):
    with pytest.raises(ValueError):
        s.formation(history(bad), rules())


def test_lossless_horizon_transfer_and_model_drift_rejection(tmp_path):
    cfg = asdict(Config(seed=1101, grid=120, dt=.00375, steps=6400, save_every=40))
    rng = np.random.default_rng(1)
    meta = dict(config=cfg, time=18., step_number=4800, divisions={}, rng=rng.bit_generator.state)
    data = dict(metadata=np.array(json.dumps(meta)), phi=rng.random((16,3,3,3)).astype('float32'),
        phase_carry=rng.normal(size=(16,3,3,3))*1e-10, chemical_correction=rng.normal(size=(2,16))*1e-17,
        rounding=np.array([432,23],dtype='int32'), ids=np.arange(16), polarity=rng.random((16,3)),
        event_audits=np.array(json.dumps([dict(event=i) for i in range(15)])),
        chemical_method=np.array(s.f.td.gpu.METHOD), development_schema=np.array(s.f.td.gpu.SCHEMA))
    source=tmp_path/'source.npz';np.savez_compressed(source,**data)
    destination=tmp_path/'transferred.npz';new=dict(cfg,steps=20800)
    assert s.checked_transfer(source,destination,new)==len(data)-1
    with np.load(destination,allow_pickle=False) as z:
        for k in set(data)-{'metadata'}:
            assert z[k].dtype==data[k].dtype and z[k].tobytes()==data[k].tobytes()
        current=json.loads(str(z['metadata']))
        assert current['config']==new and current['time']==18. and current['step_number']==4800
        current['config']=cfg;assert current==meta
    with pytest.raises(FileExistsError):s.checked_transfer(source,destination,new)
    with pytest.raises(ValueError):s.checked_transfer(source,tmp_path/'drift.npz',dict(new,signal_beta=3.))
    with pytest.raises(FileExistsError):s.prepare(tmp_path)


def test_frozen_uniform_growth_excludes_conserved_mode():
    cfg=dict(signal_beta=2.,signal_da=.02,signal_dh=.55)
    assert s.growth([0],cfg)==dict(frozen_uniform_max_growth=None,frozen_uniform_unstable_modes=0)
    r=s.growth([0,3.,12.],cfg)
    assert r['frozen_uniform_max_growth']>0 and r['frozen_uniform_unstable_modes']>0
    json.dumps(r,allow_nan=False)


@pytest.mark.parametrize('dt',[.00375,.001875])
def test_real_full_horizon_loop_initial_checkpoint_both_restarts(tmp_path,monkeypatch,dt):
    cfg=dict(dt=dt,save_every=round(.15/dt),signal_beta=2.,signal_da=.02,signal_dh=.55)
    class Fixture:
        def __init__(self,step):
            self.config=SimpleNamespace(**cfg);self.step_number=step;self.time=step*dt
            self.ids=np.arange(16);self.activator=torch.ones(16);self.inhibitor=torch.ones(16);self.target=torch.ones(16)
            self.event_audits=[dict(event=i) for i in range(15)];rng=np.random.default_rng(1)
            self.host=SimpleNamespace(parents=np.arange(16),due=np.ones(16)*10,next_id=31,lineage=[],divisions={},rng=rng,fate_rng=rng,signal_rng=rng)
        def step(self):self.step_number+=1;self.time=self.step_number*dt
        def audit(self):return dict(max_volume_error=0.,min_radius=8.,max_clipping=0.,dilution_amount_error=0.)
        def checkpoint(self,path):path.write_text(json.dumps(dict(step=self.step_number)))
        @classmethod
        def restore(cls,path):return cls(json.loads(Path(path).read_text())['step'])
    folder=tmp_path/'fixture';folder.mkdir();source=folder/'source.npz';source.write_text(json.dumps(dict(step=round(18/dt))))
    s.write(tmp_path/'protocol.json',{})
    job=dict(key='fixture',source=str(source),config=cfg)
    monkeypatch.setattr(s.f.td.gpu,'DevelopmentGpu',Fixture)
    monkeypatch.setattr(s,'row',lambda sim,gpu:dict(time=sim.time,volumes=[1.]*16,delta=np.zeros((16,16)).tolist(),boundary_occupancy=0.,nondividing_components=[1]*16,log_activator_sd=0.))
    monkeypatch.setattr(s.f.td,'network_trace',lambda sim:dict(time=sim.time))
    replay=[]
    def exact(a,b):assert a.step_number==b.step_number;replay.append(a.step_number)
    monkeypatch.setattr(s,'exact',exact)
    result=s.run_path(tmp_path,job,rules())
    assert result['passed'] and result['final_time']==78 and result['final_cells']==16
    assert result['network_trace_rows']==4001 and len(s.read(folder/'gpu-history.json'))==401
    assert result['exact_restart_steps']=={'12':8,'36':8}
    assert replay==[i for stamp in (30,54) for i in range(round(stamp/dt),round(stamp/dt)+9)]
    assert [json.loads((folder/f'gpu-e{k:02d}.npz').read_text())['step'] for k in s.SNAPSHOTS]==[round((18+k)/dt) for k in s.SNAPSHOTS]


def mock_run(tmp_path,monkeypatch,failure):
    binary=tmp_path/'binary';binary.write_bytes(b'fixed')
    pairs=[dict(seed=h,branch=b,cases=[f'h{h}_{b}_{l}' for l in ('coarse','fine')]) for h in (1101,1102,1103) for b in ('off','tension_adhesion','polarity','full')]
    p=dict(time_pairs=pairs,jobs=[dict(key=k) for pair in pairs for k in pair['cases']],runtime_binary_sha256=s.digest(binary),gpu_uuid='fixture',scope='test')
    s.write(tmp_path/'protocol.json',p);s.write(tmp_path/'status.json',dict(state='prepared',protocol_sha256=s.digest(tmp_path/'protocol.json')))
    monkeypatch.setattr(s,'verify',lambda *args:None);monkeypatch.setenv('CUDA_VISIBLE_DEVICES','fixture')
    monkeypatch.setattr(s.f.td.gpu,'library',lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(s.torch.cuda,'get_device_name',lambda *args:'fixture')
    calls=[]
    def prefix(root,job,p):calls.append(('prefix',job['key']));return dict(passed=failure!='prefix')
    def path(root,job,p):calls.append(('path',job['key']));return dict(passed=failure!='path',seconds=1.)
    def pair(root,pair,p):calls.append(('pair',pair['branch']));return dict(pair,passed=failure!='pair')
    monkeypatch.setattr(s.f.pilot,'prefix',prefix);monkeypatch.setattr(s,'run_path',path);monkeypatch.setattr(s,'compare_pair',pair)
    return calls


@pytest.mark.parametrize('failure',['prefix','path','pair'])
def test_failure_stops_remaining_jobs_without_fallback(tmp_path,monkeypatch,failure):
    calls=mock_run(tmp_path,monkeypatch,failure)
    with pytest.raises(RuntimeError):s.run(tmp_path)
    status=s.read(tmp_path/'status.json')
    assert status['state']=='failed' and not status['qualification_accepted'] and not status['maintenance_extension_started']
    assert all('h1101_off' in key for kind,key in calls if kind!='pair')
    assert not (tmp_path/'summary.json').exists()


def test_completed_coordinator_remains_provisional_and_never_starts_maintenance(tmp_path,monkeypatch):
    calls=mock_run(tmp_path,monkeypatch,None);s.run(tmp_path)
    result=s.read(tmp_path/'summary.json')
    assert result['execution_completed'] and result['independent_review_pending']
    assert not result['scientific_accepted'] and not result['qualification_accepted'] and not result['maintenance_extension_started']
    assert result['independent_histories']==3 and result['new_scientific_histories']==0
    assert len([x for x in calls if x[0]=='path'])==24 and len([x for x in calls if x[0]=='pair'])==12


@pytest.mark.parametrize('contrasts,expected', [((.099,.101),False),((.099,.099),True),((.101,.101),True)])
def test_full_pair_detects_outcome_sensitivity_even_below_chemical_error_gate(tmp_path,monkeypatch,contrasts,expected):
    pair=dict(seed=1101,branch='off',cases=['coarse','fine'])
    parent=s.read(s.PARENT/'protocol.json')
    p=dict(rules(),time_limits=parent['time_limits'],homologous_limits=parent['homologous_limits'],
           dense_time_spectrum_limit=.01,growth_abs_limit=.001,onset_time_limit=.3)
    s.write(tmp_path/'protocol.json',p)
    for key,contrast in zip(pair['cases'],contrasts):
        folder=tmp_path/key;folder.mkdir();s.write(folder/'gpu-result.json',dict(passed=True))
        rows=[]
        for i in range(401):
            u=np.exp(np.r_[np.ones(8),-np.ones(8)]*contrast)
            rows.append(dict(time=18+i*.15,elapsed=i*.15,ids=list(range(16)),parents=[0]*16,
                chemistry=[u.tolist(),[1.]*16],volumes=[1.]*16,polarity=np.zeros((16,3)).tolist(),
                delta=np.zeros((16,16)).tolist(),axis_ratio=1.,log_activator_sd=contrast,
                frozen_uniform_max_growth=None))
        s.write(folder/'gpu-history.json',rows)
        with (folder/'network.jsonl').open('w') as stream:
            for i in range(4001):
                stream.write(json.dumps(dict(time=18+i*.015,ids=list(range(16)),parents=[0]*16,
                    volumes=[1.]*16,delta=np.zeros((16,16)).tolist(),conductance=np.zeros((16,16)).tolist(),edges=0))+'\n')
    monkeypatch.setattr(s.f.q,'checkpoint_density',lambda *args:(np.zeros((6,2,2,2)),[0.,1.]))
    result=s.compare_pair(tmp_path,pair,p)
    assert result['passed']==expected
    assert result['homologous_errors']['chemical_log_max']<.01
    if not expected:assert 'formation_threshold_0.1' in result['failed_gates']
    with pytest.raises(FileExistsError):s.compare_pair(tmp_path,pair,p)
