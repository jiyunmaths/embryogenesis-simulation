import numpy as np
import pytest
from embryo.fast_polarity import compiled_exposure,FastPolaritySimulation
from embryo.polarity import exposure_cue
from embryo.model import Config
from embryo.fast_mechanics import FastAttributeSimulation


def test_cue_matches_surfaces_contacts_and_boundary_gradients():
    rng=np.random.default_rng(8)
    for phi in (rng.random((3,12,12,12),dtype=np.float32),np.ones((1,12,12,12),dtype=np.float32),np.zeros((2,12,12,12),dtype=np.float32)):
        np.testing.assert_allclose(compiled_exposure(phi,.1),exposure_cue(phi,.1),atol=2e-7,rtol=2e-5)
    sim=FastPolaritySimulation(Config(grid=24,interface_width=.12,max_cells=1,dt=.0075))
    np.testing.assert_array_equal(compiled_exposure(sim.phi,sim.dx),np.zeros((1,3)))
    phi=np.concatenate([sim.phi,np.roll(sim.phi,5,axis=1)])
    np.testing.assert_allclose(compiled_exposure(phi,sim.dx),exposure_cue(phi,sim.dx),atol=2e-7,rtol=2e-5)
    with pytest.raises(ValueError):compiled_exposure(phi.astype(float),.1)


def test_trajectory_restart_and_clamp(tmp_path):
    a=FastPolaritySimulation(Config(grid=24,interface_width=.12,max_cells=1,dt=.0075));a.polarity[:]=[.2,0,0]
    path=tmp_path/'start.npz';a.checkpoint(path);b=FastPolaritySimulation.restore(path);old=FastAttributeSimulation.restore(path)
    for _ in range(12):
        for sim in (a,b,old):sim.step(prescribed_signals=(np.array([1.2]),np.array([.9])))
    np.testing.assert_array_equal(a.phi,b.phi)
    np.testing.assert_allclose(a.phi,old.phi,atol=2e-7,rtol=2e-6)
    np.testing.assert_allclose(a.polarity,old.polarity,atol=1e-7)
