import json
import numpy as np
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.cell_response import response_metrics
from embryo.cell_response_moving import initialize, worker, observe
from embryo.feedback_long import digest, save_checkpoint


def source(tmp_path):
    sim=AttributeSimulation(Config(grid=24,max_cells=1,interface_width=.12,dt=.0075,signal_transport='conservative'),'direct')
    sim.polarity[:]=[.1,0,.05]
    path=tmp_path/'source.npz';sim.checkpoint(path)
    return sim,path


def test_pulse_preserves_geometry_polarity_clock_and_rng(tmp_path):
    before,path=source(tmp_path);chem=np.array([[1.2],[.9]])
    after=initialize(path,chem,int(before.ids[0]),.9)
    assert after.activator[0]==1.08 and after.inhibitor[0]==.9
    for key in ('phi','polarity','target','ids','parents','due'):
        np.testing.assert_array_equal(getattr(before,key),getattr(after,key))
    assert before.time==after.time
    assert before.signal_rng.bit_generator.state==after.signal_rng.bit_generator.state
    assert before.rng.bit_generator.state==after.rng.bit_generator.state
    assert before.config==after.config


def test_response_subtracts_matched_control_drift():
    times=np.arange(101.)
    control=np.exp(.003*times[:,None,None])*np.ones((101,2,2))
    path=control.copy();path[:,0,0]*=np.exp(np.log(1.1)*np.exp(-times/2))
    moving=response_metrics(path,control,times,np.ones(2),0,1.1)
    normalized=response_metrics(path/control,np.ones_like(control),times,np.ones(2),0,1.1)
    assert moving==normalized
    assert moving['final_log_rms_from_control']<1e-12


def test_worker_resume_matches_uninterrupted_trajectory(tmp_path):
    sim,path=source(tmp_path);root=tmp_path/'run';root.mkdir()
    chem=np.array([sim.activator,sim.inhibitor]);np.savez(root/'initial_states.npz',pattern=chem)
    p=dict(checkpoint=str(path),start=0.,duration=.03,interval=.015,dt=.0075)
    (root/'protocol.json').write_text(json.dumps(p));ph=digest(root/'protocol.json')
    job=dict(family='pattern',target=None,factor=1.)
    folder=root/'pattern_control';folder.mkdir()
    history=[observe(sim,0.)]
    sim.step();sim.step();history.append(observe(sim,sim.time))
    audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    save_checkpoint(sim,folder/'latest_state.npz',audit,history,ph,job)
    sim.step();sim.step()
    result=worker((str(root),job));assert result['quality_pass']
    restored=AttributeSimulation.restore(folder/'latest_state.npz')
    for key in ('phi','activator','inhibitor','polarity'):
        np.testing.assert_array_equal(getattr(sim,key),getattr(restored,key))
    history=json.loads((folder/'history.json').read_text())
    assert [r['time'] for r in history]==[0.,.015,.03]
