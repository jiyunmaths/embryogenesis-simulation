"""Protect complete endpoint inheritance and stopping before long science."""
from dataclasses import asdict
import json
from types import SimpleNamespace

import numpy as np
import pytest

from embryo import fresh_zygote_mature_qualification as f
from embryo.model import Config


def baseline():
    return asdict(Config(seed=1101,grid=120,dt=.000234375,steps=76800,save_every=640))


def endpoint(path, **changes):
    # Small arrays deliberately isolate lossless checkpoint transformation from
    # mechanics; the independent prelaunch audit checks actual 120^3 arrays.
    rng=np.random.default_rng(3)
    meta=dict(config=baseline(),time=18.,step_number=76800,divisions={},next_id=32,
        lineage=[dict(parent=1,daughters=[2,3])],rng=rng.bit_generator.state,
        fate_rng=rng.bit_generator.state,signal_rng=rng.bit_generator.state)
    arrays=dict(metadata=np.array(json.dumps(meta)),phi=rng.random((16,3,3,3)).astype('float32'),
        phase_carry=rng.normal(size=(16,3,3,3))*1e-9,rounding=np.arange(16,dtype='int32'),
        chemical_correction=rng.normal(size=(2,16))*1e-16,activator=rng.random(16)+.5,
        inhibitor=rng.random(16)+.5,polarity=rng.random((16,3)),target=rng.random(16)+.1,
        ids=np.arange(16),parents=np.arange(16),due=rng.random(16)+24,
        event_audits=np.array(json.dumps([dict(event=i) for i in range(15)])),
        development_schema=np.array(f.td.gpu.SCHEMA),chemical_method=np.array(f.td.gpu.METHOD))
    arrays.update(changes);np.savez_compressed(path,**arrays)
    return arrays,meta


@pytest.mark.parametrize('dt',[.00375,.001875])
def test_retime_changes_only_clocks_and_preserves_compensation(tmp_path,dt):
    source=tmp_path/'endpoint.npz';destination=tmp_path/'source.npz'
    arrays,meta=endpoint(source);config=f.candidate_config(baseline(),dt)
    f.retime_checkpoint(source,destination,config)
    f.check_start(dict(source_checkpoint=str(source),source=str(destination),config=config))
    with np.load(destination,allow_pickle=False) as z:
        for name in set(arrays)-{'metadata'}:
            assert z[name].dtype==arrays[name].dtype and z[name].tobytes()==arrays[name].tobytes()
        changed=json.loads(str(z['metadata']));assert changed['time']==18.
        assert changed['step_number']*dt==18 and changed['config']['steps']*dt==24
        assert changed['config']['save_every']*dt==.15
        changed.pop('config');changed.pop('step_number');meta.pop('config');meta.pop('step_number')
        assert changed==meta


@pytest.mark.parametrize('change',[dict(seed=1104),dict(grid=104),dict(dt=.0001171875)])
def test_unqualified_source_or_unregistered_clock_rejected(change):
    with pytest.raises(ValueError):f.candidate_config(dict(baseline(),**change),.00375)
    with pytest.raises(ValueError):f.candidate_config(baseline(),.0009375)


@pytest.mark.parametrize('change',[dict(time=17.),dict(step_number=1),dict(divisions={'1':{}})])
def test_incomplete_endpoint_rejected(tmp_path,change):
    source=tmp_path/'endpoint.npz';_,meta=endpoint(source)
    endpoint(source,metadata=np.array(json.dumps(dict(meta,**change))))
    with pytest.raises(ValueError):f.retime_checkpoint(source,tmp_path/'start.npz',f.candidate_config(baseline(),.00375))


def test_model_drift_and_existing_destination_rejected(tmp_path):
    source=tmp_path/'endpoint.npz';endpoint(source)
    config=f.candidate_config(baseline(),.00375);config['signal_dh']*=2
    with pytest.raises(ValueError):f.retime_checkpoint(source,tmp_path/'start.npz',config)
    destination=tmp_path/'recorded.npz';destination.write_bytes(b'original')
    with pytest.raises(FileExistsError):f.retime_checkpoint(source,destination,config)
    assert destination.read_bytes()==b'original'
    with pytest.raises(FileExistsError):f.prepare(tmp_path)


def test_saved_compensation_or_rng_reset_is_detected(tmp_path):
    source=tmp_path/'endpoint.npz';endpoint(source);destination=tmp_path/'start.npz'
    config=f.candidate_config(baseline(),.00375);f.retime_checkpoint(source,destination,config)
    with np.load(destination,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
    arrays['phase_carry'][:]=0;np.savez_compressed(destination,**arrays)
    with pytest.raises(ValueError,match='phase_carry'):f.check_start(dict(source_checkpoint=str(source),source=str(destination),config=config))


def test_invariant_is_an_independent_snapshot():
    rng=np.random.default_rng(4);host=SimpleNamespace(parents=np.array([1]),due=np.array([23.]),
        next_id=32,lineage=[dict(parent=1)],rng=rng,fate_rng=rng,signal_rng=rng)
    sim=SimpleNamespace(ids=np.array([16]),host=host,event_audits=[dict(parent=1)])
    original=f.invariant(sim);host.lineage[0]['parent']=2
    assert original['lineage'][0]['parent']==1 and f.invariant(sim)!=original


def mock_run(root,monkeypatch,failure):
    binary=root/'binary';binary.write_bytes(b'unchanged')
    pairs=[dict(seed=h,branch=b,cases=[f'h{h}_{b}_{l}' for l,_ in f.LEVELS]) for h in f.HISTORIES for b in f.BRANCHES]
    p=dict(time_pairs=pairs,jobs=[dict(key=k) for pair in pairs for k in pair['cases']],
        runtime_binary_sha256=f.digest(binary),gpu_uuid='fixture',scope='Mature clocks only')
    f.write(root/'protocol.json',p);f.write(root/'status.json',dict(state='prepared',protocol_sha256=f.digest(root/'protocol.json')))
    monkeypatch.setattr(f,'verify',lambda *args:None);monkeypatch.setenv('CUDA_VISIBLE_DEVICES','fixture')
    monkeypatch.setattr(f.td.gpu,'library',lambda:SimpleNamespace(binary_path=str(binary)))
    monkeypatch.setattr(f.torch.cuda,'get_device_name',lambda *args:'fixture')
    calls=[]
    def prefix(root,job,p):calls.append(('prefix',job['key']));return dict(passed=failure!='prefix')
    def path(root,job,p):calls.append(('path',job['key']));return dict(passed=failure!='path',seconds=1.)
    def pair(root,pair,p):calls.append(('pair',pair['branch']));return dict(pair,passed=failure!='comparison')
    monkeypatch.setattr(f.pilot,'prefix',prefix);monkeypatch.setattr(f,'run_path',path);monkeypatch.setattr(f,'compare_pair',pair)
    return calls


@pytest.mark.parametrize('failure',['prefix','path','comparison'])
def test_failure_stops_later_histories_without_fallback(tmp_path,monkeypatch,failure):
    calls=mock_run(tmp_path,monkeypatch,failure)
    with pytest.raises(RuntimeError):f.run(tmp_path)
    status=f.read(tmp_path/'status.json')
    assert status['state']=='failed' and not status['qualification_accepted']
    assert not status['long_scientific_continuation_started'] and status['new_scientific_histories']==0
    assert all('h1101_off' in case for kind,case in calls if kind!='pair')
    assert not (tmp_path/'summary.json').exists()


def test_all_paths_passing_still_requires_review(tmp_path,monkeypatch):
    calls=mock_run(tmp_path,monkeypatch,None);f.run(tmp_path)
    summary=f.read(tmp_path/'summary.json');status=f.read(tmp_path/'status.json')
    assert summary['execution_completed'] and summary['mature_screens_passed'] and summary['independent_review_pending']
    assert not summary['qualification_accepted'] and not summary['long_scientific_continuation_started']
    assert summary['new_numerical_paths']==24 and summary['new_scientific_histories']==0
    assert len([x for x in calls if x[0]=='path'])==24 and len([x for x in calls if x[0]=='pair'])==12
    assert status['state']=='completed_provisional'
