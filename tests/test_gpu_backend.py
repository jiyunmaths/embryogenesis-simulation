import os
import numpy as np
import pytest

torch = pytest.importorskip('torch')
from embryo.gpu_backend import GpuSimulation, call, gm_step, transport_matrices
from embryo.signaling import normalized_graph
from embryo.transport import contact_transport, integrate_gm


def test_conservative_and_polarity_operators_match_cpu_and_conserve_amount():
    contacts = np.array([[0., .3, .001], [.3, 0., .4], [.001, .4, 0.]])
    volumes = np.array([.8, 1.2, .4])
    centers = np.array([[0., 0., 0.], [1., .1, 0.], [1.5, .7, 0.]])
    tensors = [torch.from_numpy(a) for a in (contacts, volumes, centers)]
    delta, normalized, conductance = transport_matrices(*tensors, .12, .02, .1)
    reference = contact_transport(contacts, volumes, centers, .12, .02)
    np.testing.assert_allclose(delta, reference.delta.toarray(), atol=1e-13, rtol=1e-14)
    np.testing.assert_allclose(normalized, normalized_graph(contacts, .1).delta, atol=1e-15)
    np.testing.assert_allclose(conductance, reference.conductance.toarray(), atol=1e-13)
    np.testing.assert_allclose(volumes@delta.numpy(), 0, atol=1e-13)
    np.testing.assert_allclose(delta.numpy().sum(axis=1), 0, atol=1e-13)
    a, h = np.array([.7, 1.4, .9]), np.array([.8, 1.3, 1.1])
    # Large dt forces multiple positivity-preserving substeps.
    expected = integrate_gm(a, h, reference, .03, 2., .02, .4)
    actual = gm_step(torch.from_numpy(a), torch.from_numpy(h), delta, .03, 2., .02, .4)
    for x, y in zip(actual, expected):
        np.testing.assert_allclose(x, y, atol=1e-14, rtol=1e-14)
        assert torch.isfinite(x).all() and (x > 0).all()


def test_disconnected_polarity_graph_has_zero_rows():
    c = torch.zeros((2, 2), dtype=torch.float64)
    v = torch.ones(2, dtype=torch.float64)
    x = torch.tensor([[0., 0., 0.], [1., 0., 0.]], dtype=torch.float64)
    delta, normalized, _ = transport_matrices(c, v, x, .1, .02, .02)
    assert torch.equal(delta, c) and torch.equal(normalized, c)


gpu = pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS') != '1', reason='GPU opt-in required')


def source():
    from embryo.native_mechanics import NativeSimulation
    from embryo.model import Config
    sim = NativeSimulation(Config(grid=24, extent=1.4, interface_width=.12, max_cells=2, dt=.00375), 'direct')
    fields = []
    for x in (-.5, .5):
        radius = np.sqrt(np.sum((sim.xyz-np.array([x, 0., 0.])[:, None, None, None])**2, axis=0))
        fields.append(.5*(1-np.tanh((radius-.55)/(np.sqrt(2)*.12))))
    sim.phi = np.array(fields, dtype=np.float32)
    sim.ids = np.array([0, 1]); sim.parents = np.array([-1, -1])
    sim.due = np.full(2, np.inf); sim.fate = np.zeros(2)
    sim.activator = np.array([.8, 1.6]); sim.inhibitor = np.array([.9, 1.5])
    sim.polarity = np.array([[-.15, 0., 0.], [.2, .1, 0.]])
    sim.target = sim.volumes(); sim.next_id = 2
    sim.lineage = [dict(id=i, parent=-1, birth=0., division=None) for i in (0, 1)]
    return sim


@gpu
def test_resident_trajectory_checkpoint_restart_and_stream_ordering(tmp_path):
    native = source()
    native.native_threads = 2
    device = GpuSimulation(source())
    stream = torch.cuda.Stream()
    with torch.cuda.stream(stream):
        stream.wait_stream(torch.cuda.default_stream())
        for _ in range(8):
            device.step()
            device.audit()
    stream.synchronize()
    for _ in range(8):
        native.step()
    np.testing.assert_allclose(device.phi.cpu(), native.phi, atol=1e-7, rtol=0)
    np.testing.assert_allclose(device.activator.cpu(), native.activator, atol=1e-7, rtol=1e-7)
    np.testing.assert_allclose(device.polarity.cpu(), native.polarity, atol=2e-8, rtol=0)
    device.checkpoint(tmp_path/'gpu.npz')
    restart = GpuSimulation.restore(tmp_path/'gpu.npz')
    device.step(); restart.step()
    assert torch.equal(device.phi, restart.phi)
    assert torch.equal(device.activator, restart.activator)
    assert torch.equal(device.inhibitor, restart.inhibitor)
    assert torch.equal(device.polarity, restart.polarity)
    cpu = restart.to_cpu()
    assert cpu.rng.bit_generator.state == native.rng.bit_generator.state
    assert cpu.signal_rng.bit_generator.state == native.signal_rng.bit_generator.state
    assert cpu.ids.tolist() == native.ids.tolist() and cpu.fate.max() == 0
    assert restart.audit()['dilution_amount_error'] < 2e-14


@gpu
def test_pointer_contract_rejects_wrong_dtype_shape_and_noncontiguous_buffers():
    sim = GpuSimulation(source())
    args = [sim.phi, sim.h, sim.shell, sim.shell2, sim.occupied, sim.excluded]
    dimensions = (len(sim.ids), sim.config.grid)
    for bad, expected in ((sim.phi.double(), 'dtype'), (sim.phi[:, :-1].contiguous(), 'shape'),
                          (sim.phi.transpose(-1, -2), 'contiguous')):
        with pytest.raises(ValueError, match=expected):
            call('gpu_arrays', [bad, *args[1:]], dimensions)


@gpu
def test_unsupported_division_or_nonpolar_states_are_rejected():
    for kind in ('division', 'nonpolar', 'no_feedback'):
        sim = source()
        if kind == 'division':
            sim.divisions[int(sim.ids[0])] = {}
        elif kind == 'nonpolar':
            sim.polarity[0] = 0
        else:
            sim.config.feedback = False
        with pytest.raises(ValueError, match='mature'):
            GpuSimulation(sim)
