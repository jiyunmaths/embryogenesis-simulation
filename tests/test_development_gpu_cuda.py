"""Hardware software-contract checks; not developmental convergence evidence."""
import os
import numpy as np
import pytest
import torch
from embryo.attribute_development import AttributeSimulation
from embryo.development_gpu import DevelopmentGpu
from embryo.model import Config
from embryo.qualify_development_gpu import exact

pytestmark=pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='Explicit CUDA opt-in required')


def test_shape_axis_and_division_start_preserve_initial_fields():
    torch.set_num_threads(1)
    h=AttributeSimulation(Config(grid=32,extent=1.6,interface_width=.12,max_cells=2,dt=.00375))
    axis=np.array([1.,2.,-1.]);axis/=np.linalg.norm(axis)
    along=np.einsum('d,dijk->ijk',axis,h.xyz)
    radius=np.sqrt((along/.95)**2+(np.sum(h.xyz**2,axis=0)-along**2)/.45**2)
    h.phi=(.5*(1-np.tanh((radius-1)/.12)))[None].astype(np.float32);h.target=h.volumes();h.initial_volume=float(h.target.sum())
    s=DevelopmentGpu(h);before={n:getattr(s,n).clone() for n in ('phi','phase_carry','activator','inhibitor','chemical_correction','geometry')}
    assert s.divide(0)
    assert abs(np.asarray(s.host.divisions[0]['axis'])@axis)>.999
    for name,value in before.items():assert torch.equal(getattr(s,name),value)
    assert len(s.ids)==1 and not s.divide(0)


def test_naturally_scheduled_four_cells_and_restart_with_another_active_mother(tmp_path):
    torch.set_num_threads(1)
    h=AttributeSimulation(Config(grid=32,extent=1.6,interface_width=.12,max_cells=4,dt=.0075,
        division_interval=.3,signal_dh=.55,signal_partition_noise=.001))
    s=DevelopmentGpu(h);clone=None;observed=set();target=float(s.target.sum())
    for _ in range(800):
        s.step();observed.add(len(s.ids))
        if len(s.ids)==3 and clone is None:
            assert s.host.divisions
            cp=tmp_path/'three-cells-live-furrow.npz';s.checkpoint(cp);clone=DevelopmentGpu.restore(cp);exact(s,clone)
        elif clone is not None:clone.step();exact(s,clone)
        if len(s.ids)==4 and not s.host.divisions:break
    assert observed=={1,2,3,4} and clone is not None and len(s.event_audits)==3
    assert abs(float(s.target.sum())/target-1)<=2e-14 and not s.host.divisions and np.all(s.host.fate==0)
    for ev in s.event_audits:
        assert ev['amount_relative_jump']<=2e-14 and ev['target_relative_jump']<=2e-14 and ev['volume_relative_jump']<=1e-6
    assert s.phase_carry.shape==s.phi.shape==(4,32,32,32) and s.chemical_correction.shape==(2,4)
