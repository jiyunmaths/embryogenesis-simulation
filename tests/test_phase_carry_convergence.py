import numpy as np
import pytest

from embryo.feedback_survival_validation import retime
from embryo.phase_carry_convergence import error_trend
from embryo.polarity_robustness_refinement import check_halving
from test_gpu_backend import source


def test_coarsening_preserves_state_and_physical_clocks(tmp_path):
    fine = source(); fine.time = 150.; fine.config.dt = .001875
    fine.step_number = round(fine.time/fine.config.dt)
    fine.config.steps = 80000; fine.config.save_every = 80
    before = tmp_path/'fine.npz'; after = tmp_path/'coarse.npz'; fine.checkpoint(before)
    coarse = retime(before, after, .00375)
    check_halving(coarse, fine)
    assert coarse.config.save_every*coarse.config.dt == fine.config.save_every*fine.config.dt
    np.testing.assert_array_equal(coarse.phi, fine.phi)
    assert coarse.signal_rng.bit_generator.state == fine.signal_rng.bit_generator.state


def test_declared_trend_does_not_accept_increasing_or_degenerate_errors():
    result = error_trend(.0008, .0002)
    assert result['passed'] and result['observed_pair_order'] == 2.
    assert not error_trend(.0002, .0008)['passed']
    assert not error_trend(0., 0.)['passed']
    assert not error_trend(.00024, .0002)['passed']
    with pytest.raises(ValueError): error_trend(float('nan'), .0002)
    with pytest.raises(ValueError): error_trend(-1., .0002)
