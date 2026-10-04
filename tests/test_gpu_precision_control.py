import os
import numpy as np
import pytest

torch = pytest.importorskip('torch')
from embryo.gpu_backend import GpuSimulation
from embryo.gpu_precision_control import PrecisionSimulation, ARMS, force
from test_gpu_backend import source

gpu = pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')


def equal_state(a, b):
    for key in ('phi', 'activator', 'inhibitor', 'polarity', 'geometry'):
        assert torch.equal(getattr(a, key), getattr(b, key)), key
    assert a.time == b.time and a.step_number == b.step_number
    assert torch.equal(a.matrices()[1], b.matrices()[1])


@gpu
def test_baseline_multistep_and_zero_carry_match_accepted_kernel():
    original = GpuSimulation(source()); control = PrecisionSimulation(source())
    carry = PrecisionSimulation(source(), 'phase_carry')
    stream = torch.cuda.Stream()
    with torch.cuda.stream(stream):
        stream.wait_stream(torch.cuda.default_stream())
        original.step(); control.step(); carry.step()
        equal_state(original, control); equal_state(original, carry)
        assert torch.any(carry.phase_carry != 0)
        for _ in range(15):
            original.step(); control.step()
            equal_state(original, control)
    stream.synchronize()
    assert control.precision_diagnostics()['phase_carry_abs_max'] == 0
    counts = control.precision_diagnostics()
    assert 0 <= counts['uncarried_stuck_voxels'] <= counts['active_force_voxels']


@gpu
def test_rounding_counter_detects_deliberately_sub_float32_updates():
    host = source(); host.config.dt = 1e-10
    sim = PrecisionSimulation(host, 'phase_carry'); sim.step()
    counts = sim.precision_diagnostics()
    assert counts['uncarried_stuck_voxels'] > 0
    assert counts['phase_carry_abs_max'] > 0


@gpu
@pytest.mark.parametrize('arm', ARMS)
def test_carry_checkpoint_restart_is_exact(tmp_path, arm):
    sim = PrecisionSimulation(source(), arm)
    for _ in range(8): sim.step()
    path = tmp_path/'state.npz'; sim.checkpoint(path)
    restored = PrecisionSimulation.restore(path)
    assert restored.precision_arm == arm
    assert torch.equal(sim.phase_carry, restored.phase_carry)
    for _ in range(4):
        sim.step(); restored.step(); equal_state(sim, restored)
        assert torch.equal(sim.phase_carry, restored.phase_carry)
        assert torch.equal(sim.rounding, restored.rounding)
    assert sim.audit()['dilution_amount_error'] < 2e-14
    # Tampering is checked before uploading the invalid state.
    with np.load(path) as z: payload = {k: z[k].copy() for k in z.files}
    payload['phase_carry'] = payload['phase_carry'].astype(np.float32)
    np.savez_compressed(path, **payload)
    with pytest.raises(ValueError, match='checkpoint'): PrecisionSimulation.restore(path)


@gpu
def test_double_contacts_match_cpu_dot_and_keep_conservation():
    sim = PrecisionSimulation(source(), 'contact64')
    contact, delta, _, _ = sim.matrices()
    shell = sim.shell.cpu().numpy().astype(np.float64).reshape(len(sim.ids), -1)
    expected = shell@shell.T*sim.dx**3; np.fill_diagonal(expected, 0)
    np.testing.assert_allclose(contact.cpu(), expected, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(delta.cpu().numpy().sum(1), 0, atol=1e-13)
    np.testing.assert_allclose(sim.geometry[:, 0].cpu().numpy()@delta.cpu().numpy(), 0, atol=1e-13)


@gpu
def test_precision_pointer_contract():
    sim = PrecisionSimulation(source()); c = sim.config
    args = [sim.phi, sim.excluded, sim.gamma, sim.geometry[:, 1:4].contiguous(),
            sim.polarity, sim.target, sim.target, sim.gamma, sim.updated, sim.quality,
            sim.phase_carry, sim.rounding]
    scalars = (len(sim.ids), c.grid, 0, sim.dx, c.extent, c.interface_width,
               c.polarity_tension, c.repulsion, c.dt)
    for bad, pattern in ((sim.phi.double(), 'dtype'), (sim.phi[:, :-1].contiguous(), 'shape'),
                         (sim.phi.transpose(-1, -2), 'contiguous')):
        with pytest.raises(ValueError, match=pattern): force([bad, *args[1:]], scalars)
