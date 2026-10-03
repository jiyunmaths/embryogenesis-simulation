import numpy as np
import pytest
from embryo.cell_exchange_response_moving import local_response,donor_id,assess,BACKGROUNDS,source_ready
from embryo.feedback_long import digest
from embryo.resolution import write_json
from embryo.cell_response_moving import name


def test_local_response_removes_control_drift_and_preserves_sign():
    t=np.arange(4.)[:,None,None];control=np.exp(t*.1)*np.ones((4,2,2));path=control.copy();path[:,0,0]*=.9
    wave=local_response(path,control,0,.9)
    np.testing.assert_allclose(wave[:,0],-1);np.testing.assert_array_equal(wave[:,1],0)
    assert donor_id(27,[27,20])==20 and donor_id(20,[27,20])==27
    with pytest.raises(ValueError):donor_id(3,[27,20])


def test_source_gate_rejects_partial_or_invalid_result(tmp_path):
    source=dict(root=str(tmp_path),folder='arm');folder=tmp_path/'arm';folder.mkdir()
    write_json(tmp_path/'protocol.json',{})
    write_json(tmp_path/'status.json',dict(state='running'))
    assert not source_ready(source)
    write_json(folder/'result.json',dict(quality_pass=False,protocol_sha256=digest(tmp_path/'protocol.json')))
    with pytest.raises(ValueError):source_ready(source)
    write_json(folder/'result.json',dict(quality_pass=True,protocol_sha256=digest(tmp_path/'protocol.json')))
    assert source_ready(source)


def test_full_assessment_detects_transferred_responses(tmp_path):
    p=dict(backgrounds=BACKGROUNDS,targets=[27,20],factors=[.9,1.1],duration=60.,interval=30.,scope='test',interpretation='test',limitations='test')
    times=np.array([0.,30.,60.]);ids=np.array([27,20]);m=np.ones(2)
    for key in BACKGROUNDS:
        child=tmp_path/key;child.mkdir()
        jobs=[dict(family=key,target=None,factor=1.)]+[dict(family=key,target=cell,factor=factor) for cell in p['targets'] for factor in p['factors']]
        write_json(child/'protocol.json',dict(jobs=jobs));np.savez(child/'initial_states.npz',ids=ids,masses=m)
        control=np.exp(.003*times[:,None,None])*np.ones((3,2,2))
        for job in jobs:
            path=control.copy()
            if job['target'] is not None:
                i=list(ids).index(job['target']);tau=[10.,30.][i if key=='unexchanged' else 1-i]
                path[:,0,i]*=np.exp(np.log(job['factor'])*np.exp(-times/tau))
            folder=child/name(job);folder.mkdir()
            write_json(folder/'result.json',dict(quality_pass=True,protocol_sha256=digest(child/'protocol.json')))
            write_json(folder/'history.json',[dict(elapsed=float(t),ids=ids.tolist(),chemistry=x.tolist()) for t,x in zip(times,path)])
    report=assess(tmp_path,p)
    assert len(report['comparisons'])==8
    assert all(x['nearest_reference']=='donor' and x['donor_distance']<1e-12 for x in report['comparisons'])
