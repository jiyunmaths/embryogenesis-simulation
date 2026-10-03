import json
import pytest
from embryo.validate_fast_polarity import run as polarity_run
from embryo.validate_native_mechanics import run as native_run
import numpy as np
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.cell_response import pulse
from embryo.cell_response_moving import observe
from embryo.resolution import write_json
from embryo.validate_fast_mechanics import run


@pytest.mark.parametrize('runner', [run, polarity_run, native_run])
def test_full_validation_runner_with_short_matched_references(tmp_path, runner):
    source=tmp_path/'source';source.mkdir();root=tmp_path/'validation';root.mkdir()
    sim=AttributeSimulation(Config(grid=24,interface_width=.12,dt=.0075,max_cells=2),'direct');sim.divide(0,direction=[1,0,0])
    for _ in range(240):
        sim.step()
        if len(sim.ids)==2 and not sim.divisions:break
    assert len(sim.ids)==2
    sim.checkpoint(source/'source.npz');cell=int(sim.ids[0]);start=sim.time
    np.savez(source/'initial_states.npz',ids=sim.ids,masses=sim.volumes())
    jobs=[dict(key='control',cell=None,factor=1.),dict(key='pulse',cell=cell,factor=.9)]
    for job in jobs:
        s=AttributeSimulation.restore(source/'source.npz');folder=source/job['key'];folder.mkdir()
        if job['cell'] is not None:s.activator,s.inhibitor=pulse(np.array([s.activator,s.inhibitor]),0,.9)
        history=[observe(s,0.)]
        for i in range(4):
            s.step()
            if i%2==1:history.append(observe(s,s.time-start))
        write_json(folder/'history.json',history);s.checkpoint(folder/'latest_state.npz')
    criteria=dict(chemical_log_max=1e-5,polarity_abs_max=1e-5,relative_axis_max=1e-4,relative_volume_max=1e-5,final_phi_abs_max=2e-5,normalized_response_max=.001,relative_auc_max=.002,recovery_time_error_max=.15)
    p=dict(source=str(source),checkpoint=str(source/'source.npz'),start=start,duration=.03,interval=.015,dt=.0075,jobs=jobs,criteria=criteria,scope='fixture',source_sha256={},input_sha256={})
    write_json(root/'protocol.json',p);runner(root)
    result=json.loads((root/'comparison.json').read_text());assert result['passed']
    runner(root) # Reusing accepted completed results preserves the comparison.
    assert json.loads((root/'comparison.json').read_text())==result
