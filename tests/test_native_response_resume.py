import numpy as np
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.cell_response_moving import observe
from embryo.feedback_long import save_checkpoint,digest
from embryo.native_response_resume import worker
from embryo.resolution import write_json


def test_original_checkpoint_resumes_with_native_worker(tmp_path):
    sim=AttributeSimulation(Config(grid=24,interface_width=.12,max_cells=1,dt=.0075),'direct');sim.polarity[:]=[.2,0,0]
    source=tmp_path/'source.npz';sim.checkpoint(source)
    root=tmp_path/'experiment';root.mkdir();job=dict(family='pattern',target=None,factor=1.)
    write_json(root/'protocol.json',dict(checkpoint=str(source),start=0.,dt=.0075,duration=.03,interval=.015))
    np.savez(root/'initial_states.npz',pattern=np.array([sim.activator,sim.inhibitor]))
    folder=root/'pattern_control';folder.mkdir();history=[observe(sim,0.)]
    sim.step();sim.step();history.append(observe(sim,sim.time))
    audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    save_checkpoint(sim,folder/'latest_state.npz',audit,history,digest(root/'protocol.json'),job)
    sim.step();sim.step();r=worker((str(root),job,2));assert r['quality_pass']
    restored=AttributeSimulation.restore(folder/'latest_state.npz')
    np.testing.assert_allclose(sim.phi,restored.phi,atol=2e-7,rtol=2e-6)
    np.testing.assert_allclose(sim.activator,restored.activator,atol=1e-7)
    assert sim.time==restored.time
