"""Opt-in GPU tests: EMBRYO_CUDA_TESTS=1 CUDA_VISIBLE_DEVICES=1 pytest ..."""
import os
import numpy as np
import pytest
from embryo.cuda_benchmark import update, arguments
from embryo.native_mechanics import NativeSimulation, kernel
from embryo.model import Config

pytestmark = pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1',
                               reason='GPU tests explicitly opt in; never occupy GPU by default')


@pytest.mark.parametrize('precision', [64, 32])
@pytest.mark.parametrize('field_kind', ['zero', 'one', 'random'])
def test_cuda_stencils_and_precision(precision, field_kind):
    sim = NativeSimulation(Config(grid=12, interface_width=.2), 'direct')
    rng = np.random.default_rng(12)
    sim.phi = rng.random(sim.phi.shape, dtype=np.float32)
    if field_kind == 'zero': sim.phi.fill(0)
    if field_kind == 'one': sim.phi.fill(1)
    shape = sim.phi.shape
    # Nonuniform inputs exercise all boundaries, attraction, repulsion,
    # volume response and a displaced polar cortical-tension field.
    arrays = (sim.phi, rng.random(shape, dtype=np.float32),
              rng.random(shape), np.array([[.1, -.2, .05]]),
              np.array([[.3, .2, -.1]]), np.array([1.2]), np.array([-.4]))
    expected = np.empty_like(sim.phi)
    kernel()(*arguments(sim, arrays, expected), 1)
    actual, seconds = update(sim, arrays, precision, 2)
    assert seconds > 0
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, atol=2e-7 if precision == 32 else 1e-7, rtol=0)
