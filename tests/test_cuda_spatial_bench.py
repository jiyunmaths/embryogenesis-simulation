"""Opt-in checks for the additional CUDA spatial benchmark components."""
import os
import numpy as np
import pytest
from embryo.model import Config
from embryo.native_mechanics import NativeSimulation
from embryo.cuda_spatial_bench import library, snapshot

pytestmark=pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='Explicit GPU opt-in required')


@pytest.mark.parametrize('cells', [1,2,4])
def test_spatial_reductions_and_matrix_layout(cells):
    sim=NativeSimulation(Config(grid=12,interface_width=.2),'direct')
    rng=np.random.default_rng(19)
    sim.phi=rng.random((cells,12,12,12),dtype=np.float32)
    sim.ids=np.arange(cells);sim.activator=np.ones(cells)
    # Deliberately asymmetric matrix detects an accidental transpose.
    adhesion=rng.random((cells,cells));np.fill_diagonal(adhesion,0)
    sim.material_coefficients=lambda:(np.ones(cells),adhesion)
    result=snapshot(library(),sim,2);errors=result['errors']
    assert errors['volume_relative']<1e-13
    assert errors['center_abs']<1e-13
    assert errors['cue_abs']<2e-7
    assert errors['contact_relative_frobenius']<2e-6
    assert errors['adhesion_abs']<1e-13
    assert all(v>0 for v in result['component_ms'].values())
