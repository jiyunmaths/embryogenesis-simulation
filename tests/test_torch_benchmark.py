"""Opt-in correctness checks for the PyTorch benchmark implementation."""
import os
import numpy as np
import pytest
torch=pytest.importorskip('torch')
from embryo.torch_benchmark import arrays,contacts,attraction,geometry,mechanics
from embryo.native_mechanics import NativeSimulation,kernel
from embryo.cuda_benchmark import arguments
from embryo.model import Config
from embryo.polarity import exposure_cue

pytestmark=pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU benchmarks require explicit opt-in')


@pytest.mark.parametrize('cells',[1,2,4])
@torch.inference_mode()
def test_torch_spatial_boundaries_and_matrix_layout(cells):
    sim=NativeSimulation(Config(grid=12,interface_width=.2),'direct')
    sim.phi=np.random.default_rng(4).random((cells,12,12,12),dtype=np.float32)
    sim.ids=np.arange(cells);sim.activator=np.ones(cells)
    tensor=lambda a:torch.from_numpy(np.ascontiguousarray(a)).cuda()
    phi=tensor(sim.phi);h,shell,shell2,occupied=arrays(phi)
    vol,cent,cue=geometry(phi,h,shell,occupied,tensor(sim.xyz),sim.dx)
    np.testing.assert_allclose(vol.cpu(),sim.volumes(),rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(cent.cpu(),sim.centers(),rtol=0,atol=1e-13)
    np.testing.assert_allclose(cue.cpu(),exposure_cue(sim.phi,sim.dx),rtol=0,atol=2e-7)
    np.testing.assert_allclose(contacts(shell,sim.dx).cpu(),sim.contacts()[0],rtol=2e-6,atol=1e-8)
    adhesion=np.random.default_rng(7).random((cells,cells))
    target=(adhesion@((sim.phi*(1-sim.phi))**2).reshape(cells,-1)).reshape(sim.phi.shape)
    np.testing.assert_allclose(attraction(shell2,tensor(adhesion)).cpu(),target,rtol=1e-13,atol=1e-13)


@pytest.mark.parametrize('precision',[64,32])
@pytest.mark.parametrize('kind',['random','empty','full'])
@torch.inference_mode()
def test_torch_mechanical_forces_match_cpu(precision,kind):
    sim=NativeSimulation(Config(grid=12,interface_width=.2),'direct')
    rng=np.random.default_rng(5);sim.phi=rng.random(sim.phi.shape,dtype=np.float32)
    if kind=='empty':sim.phi.fill(0)
    if kind=='full':sim.phi.fill(1)
    prepared=(sim.phi,rng.random(sim.phi.shape,dtype=np.float32),rng.random(sim.phi.shape),
              np.array([[.2,-.1,.07]]),np.array([[.3,-.2,.1]]),np.array([1.2]),np.array([-.4]))
    expected=np.empty_like(sim.phi);kernel()(*arguments(sim,prepared,expected),1)
    tensors=[torch.from_numpy(a).cuda() for a in prepared];c=sim.config
    actual=mechanics(*tensors,torch.from_numpy(sim.xyz).cuda(),sim.dx,c.interface_width,
                     c.polarity_tension,c.repulsion,c.dt,precision)
    np.testing.assert_allclose(actual.cpu(),expected,rtol=0,atol=2e-7 if precision==32 else 1e-7)
