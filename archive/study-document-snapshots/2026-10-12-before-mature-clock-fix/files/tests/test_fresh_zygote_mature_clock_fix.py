"""Exercise initial checkpoint scheduling through the real coordinator loop."""
from types import SimpleNamespace
from pathlib import Path
import json

import numpy as np
import pytest
import torch

from embryo import fresh_zygote_mature_clock_fix as c
from embryo import fresh_zygote_mature_qualification as f


@pytest.mark.parametrize('dt',[.00375,.001875])
def test_zero_checkpoint_and_positive_clocks(dt):
    assert c.elapsed_steps(0,dt)==0
    for duration in (.015,.15,3.,6.,18.,24.):
        assert c.elapsed_steps(duration,dt)==c.ORIGINAL_STEPS(duration,dt)


@pytest.mark.parametrize('dt',[0.,-1.,float('nan'),float('inf')])
def test_zero_checkpoint_retains_dt_validation(dt):
    with pytest.raises(ValueError):c.elapsed_steps(0,dt)


def test_real_path_loop_saves_initial_state_and_completes_restart(tmp_path,monkeypatch):
    class FixtureSimulation:
        def __init__(self,step=4800):
            self.config=SimpleNamespace(dt=.00375,save_every=40);self.step_number=step;self.time=step*.00375
            self.ids=np.arange(16);self.activator=torch.ones(16);self.inhibitor=torch.ones(16);self.target=torch.ones(16)
            self.event_audits=[dict(event=i) for i in range(15)]
            rng=np.random.default_rng(1)
            self.host=SimpleNamespace(parents=np.arange(16),due=np.ones(16)*10,next_id=32,lineage=[],
                divisions={},rng=rng,fate_rng=rng,signal_rng=rng)
        def step(self):self.step_number+=1;self.time=self.step_number*self.config.dt
        def audit(self):return dict(max_volume_error=0.,min_radius=8.,max_clipping=0.,dilution_amount_error=0.)
        def checkpoint(self,path):path.write_text(json.dumps(dict(step=self.step_number)))
        @classmethod
        def restore(cls,path):return cls(json.loads(Path(path).read_text())['step'])
    folder=tmp_path/'fixture';folder.mkdir();source=folder/'source.npz';source.write_text('{"step":4800}')
    f.write(tmp_path/'protocol.json',{})
    job=dict(key='fixture',source=str(source),config=dict(dt=.00375))
    p=dict(physical_limits=f.q.PHYSICAL,diagnostic_interval=.015,restart_elapsed=3.,restart_steps=8)
    monkeypatch.setattr(f,'_steps',c.elapsed_steps)
    monkeypatch.setattr(f.td.gpu,'DevelopmentGpu',FixtureSimulation)
    monkeypatch.setattr(f,'row',lambda sim,gpu:dict(time=sim.time,volumes=[1.]*16,delta=np.zeros((16,16)).tolist(),boundary_occupancy=0.,nondividing_components=[1]*16))
    monkeypatch.setattr(f.td,'network_trace',lambda sim:dict(time=sim.time))
    replay=[]
    def exact(a,b):
        assert a.step_number==b.step_number
        replay.append(a.step_number)
    monkeypatch.setattr(f,'exact',exact)
    result=f.run_path(tmp_path,job,p)
    assert result['passed'] and result['final_time']==24 and result['lineage_cycle_rng_unchanged']
    assert result['network_trace_rows']==401 and result['exact_restart_steps']==8
    assert [json.loads((folder/f'gpu-e{k:02d}.npz').read_text())['step'] for k in (0,3,6)]==[4800,5600,6400]
    assert replay==list(range(5600,5609))
    assert len(f.read(folder/'gpu-history.json'))==41


def test_existing_corrected_study_is_never_overwritten(tmp_path):
    with pytest.raises(FileExistsError):c.prepare(tmp_path)
