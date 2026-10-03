import json
import numpy as np
from embryo.cell_exchange_moving import gate_status,compare,FAMILIES
from embryo.feedback_long import digest
from embryo.resolution import write_json


def test_gate_requires_finished_passing_refinement(tmp_path):
    write_json(tmp_path/'status.json',dict(state='running'))
    assert gate_status(tmp_path)=='waiting'
    write_json(tmp_path/'status.json',dict(state='completed'))
    assert gate_status(tmp_path)=='waiting'
    write_json(tmp_path/'refinement.json',dict(passed=False))
    assert gate_status(tmp_path)=='blocked'
    write_json(tmp_path/'refinement.json',dict(passed=True))
    assert gate_status(tmp_path)=='ready'
    write_json(tmp_path/'status.json',dict(state='failed'))
    assert gate_status(tmp_path)=='blocked'


def test_comparison_distinguishes_destination_from_transferred_state(tmp_path):
    base=np.array([[.2,2.],[.3,1.]])
    swapped=base[:,::-1].copy();m=np.ones(2);ids=[1,2];times=[0.,36.,60.]
    dependency=tmp_path/'dependency';(dependency/'pattern_control').mkdir(parents=True)
    def rows(x):return [dict(elapsed=t,ids=ids,chemistry=x.tolist(),delta=[[-1.,1.],[1.,-1.]],axis_ratio=1.,log_activator_sd=float(np.std(np.log(x[0])))) for t in times]
    write_json(dependency/'pattern_control/history.json',rows(base))
    p=dict(dependency=str(dependency),pair=[0,1],criteria=dict(late_start=36.),jobs=[dict(family=f,target=None,factor=1.) for f in FAMILIES],scope='test',time_refinement='test')
    write_json(tmp_path/'protocol.json',p)
    np.savez(tmp_path/'initial_states.npz',masses=m,ids=ids,baseline=base,fresh_exchange=swapped,relaxed_exchange=swapped)
    np.savez(tmp_path/'frozen_reference.npz',**{f:np.array([swapped]*3) for f in FAMILIES})
    for f,x in zip(FAMILIES,(base,swapped)):
        folder=tmp_path/(f+'_control');folder.mkdir()
        write_json(folder/'history.json',rows(x));write_json(folder/'result.json',dict(quality_pass=True,protocol_sha256=digest(tmp_path/'protocol.json')))
    compare(tmp_path,p);r=json.loads((tmp_path/'comparison.json').read_text())['trials']
    assert r[0]['destination_like'] and not r[0]['transferred_like']
    assert r[1]['transferred_like'] and r[1]['late_pair_order_reversed']
