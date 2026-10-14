"""Scientific diagnostic noninterference, not tests of copied equations."""
import io
from types import SimpleNamespace
import numpy as np
import pytest
import torch
from embryo import development_gpu as gpu
from embryo.development_time_diagnosis import projection_trace


def simulation():
    return SimpleNamespace(ids=np.array([1,2]),host=SimpleNamespace(divisions={2:{}}),
        config=SimpleNamespace(dt=.001875),time=4.,step_number=100)


def extrema():
    return dict(calls=0,correction_magnitude_max=0.,before_relative_volume_error_max=0.,after_relative_volume_error_max=0.)


def test_pressure_logging_does_not_change_inputs_or_returned_phase_pair():
    field=torch.full((4,4,4),.4,dtype=torch.float32);carry=torch.full_like(field,1e-10,dtype=torch.float64)
    before=[field.clone(),carry.clone()];volume=float(gpu.occupancy64(field.double()+carry).sum())*1.001
    expected=gpu.project_phase(field,carry,volume,1.)
    original=gpu.project_phase;stream=io.StringIO();bounds=extrema()
    with projection_trace(simulation(),stream,bounds):actual=gpu.project_phase(field,carry,volume,1.)
    assert gpu.project_phase is original
    for a,b in zip(expected[:2],actual[:2]):assert torch.equal(a,b)
    assert expected[2]==actual[2]
    assert torch.equal(field,before[0]) and torch.equal(carry,before[1])
    assert bounds['calls']==1 and bounds['after_relative_volume_error_max']<=1e-7
    import json
    assert json.loads(stream.getvalue())['parent']==2


def test_failed_projection_restores_original_function_and_lock(monkeypatch):
    def fail(*args):raise FloatingPointError('original failure')
    monkeypatch.setattr(gpu,'project_phase',fail)
    field=torch.ones((2,2,2),dtype=torch.float32);carry=torch.zeros_like(field,dtype=torch.float64)
    with pytest.raises(FloatingPointError,match='original failure'):
        with projection_trace(simulation(),io.StringIO(),extrema()):gpu.project_phase(field,carry,1.,1.)
    assert gpu.project_phase is fail
    # A second context can acquire the guard after the exception.
    with pytest.raises(FloatingPointError):
        with projection_trace(simulation(),io.StringIO(),extrema()):gpu.project_phase(field,carry,1.,1.)


def test_nested_projection_tracing_cannot_mix_parent_logs():
    s=SimpleNamespace(ids=np.array([1]),host=SimpleNamespace(divisions={}),config=SimpleNamespace(dt=.01),time=0.,step_number=0)
    with projection_trace(s,io.StringIO(),extrema()):
        with pytest.raises(RuntimeError,match='serial'):
            with projection_trace(s,io.StringIO(),extrema()):pass
