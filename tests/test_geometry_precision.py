import numpy as np
import pytest
from embryo.geometry_precision import errors, checked_history


def rows():
    return [dict(elapsed=t, time=150+t, ids=[1, 2], chemistry=[[1., 2.], [3., 4.]],
                 volumes=[1., 2.], delta=[[-2., 2.], [1., -1.]], polarity=[[1., 0., 0.], [0., 1., 0.]],
                 centers=[[0., 0., 0.], [1., 0., 0.]], uniform_growth_max=.2) for t in (0., .15)]


def test_comparisons_use_absolute_log_error_near_uniform_and_check_clocks():
    a, b = rows(), rows()
    b[-1]['chemistry'][0][0] = np.exp(.02)
    assert errors(a, b)['chemical_log_max'] == pytest.approx(.02)
    assert errors(a, b)['relative_transport_max'] == 0
    b[-1]['elapsed'] += .01
    with pytest.raises(ValueError, match='Unaligned'): errors(a, b)


def test_history_rejects_mismatched_physical_clock_and_nonconservative_graph(tmp_path):
    f = tmp_path/'chemical.npz'; np.savez(f, uniform=rows()[0]['chemistry'], ids=[1, 2])
    p, job = dict(interval=.15, start=150.), dict(chemical_file=str(f))
    assert checked_history(rows(), p, job, .15) == rows()
    a = rows(); a[-1]['time'] += .01
    with pytest.raises(ValueError, match='physical clock'): checked_history(a, p, job, .15)
    a = rows(); a[-1]['delta'][0][0] += .1
    with pytest.raises(ValueError): checked_history(a, p, job, .15)


import os
from dataclasses import asdict
from embryo.feedback_survival_validation import retime
from embryo.geometry_precision import worker
from embryo.feedback_long import digest
from embryo.neighbor_context import read
from embryo.resolution import write_json

gpu = pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')


@gpu
def test_moving_runner_resumes_carry_without_changing_completed_path(tmp_path):
    from test_gpu_backend import source
    from embryo.gpu_precision_control import PrecisionSimulation
    import torch
    host = source(); host.time = 150.; host.step_number = round(150/host.config.dt)
    checkpoint = tmp_path/'source.npz'; host.checkpoint(checkpoint)
    chemical = tmp_path/'chemical.npz'
    np.savez(chemical, uniform=[host.activator, host.inhibitor], ids=host.ids)
    job = dict(key='both_fine', arm='both', level='fine', dt=host.config.dt, source=str(checkpoint),
               chemical_file=str(chemical), family='uniform', point=dict(ratio=27.5, chi=0., key='test'))
    p = dict(start=150., duration=.3, interval=.15, checkpoint_interval=.15,
             accepted_configs=dict(fine=asdict(host.config)), criteria=dict(boundary_max=.01, dilution_error_max=2e-14),
             source_sha256={}, input_sha256={})
    write_json(tmp_path/'protocol.json', p)
    # Simulate an interruption after a complete first-half atomic checkpoint.
    worker(tmp_path, job, dict(p, duration=.15))
    (tmp_path/'both_fine/result.json').unlink()
    resumed = worker(tmp_path, job, p)
    saved = PrecisionSimulation.restore(tmp_path/'both_fine/latest_state.npz')
    other = dict(job, key='uninterrupted')
    full = worker(tmp_path, other, p)
    original = PrecisionSimulation.restore(tmp_path/'uninterrupted/latest_state.npz')
    assert torch.equal(saved.phi, original.phi)
    assert torch.equal(saved.activator, original.activator)
    assert torch.equal(saved.polarity, original.polarity)
    assert torch.equal(saved.phase_carry, original.phase_carry)
    assert read(tmp_path/'both_fine/history.json') == read(tmp_path/'uninterrupted/history.json')
    assert worker(tmp_path, job, p) == resumed
    assert full['quality_pass']
