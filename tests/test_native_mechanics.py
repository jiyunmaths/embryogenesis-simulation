"""Numerical contract of the opt-in threaded backend."""
import numpy as np
import pytest
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config
from embryo.native_mechanics import NativeSimulation
from embryo.polarity import exposure_cue


@pytest.mark.parametrize('threads', [1, 2, 4])
def test_cue_matches_reference_at_contacts_and_boundaries(threads):
    sim = NativeSimulation(Config(grid=12, interface_width=.2, max_cells=2), 'direct')
    sim.divide(0, direction=[1, 0, 0])
    sim.native_threads = threads
    sim.phi = np.random.default_rng(19).random(sim.phi.shape, dtype=np.float32)
    np.testing.assert_array_equal(sim.native_exposure(), exposure_cue(sim.phi, sim.dx))
    assert sim._fast_cache is None


@pytest.mark.parametrize('threads', [1, 4])
@pytest.mark.parametrize('mode', ['polar', 'nonpolar', 'division', 'prescribed'])
def test_trajectory_and_restart(tmp_path, threads, mode):
    c = Config(grid=24, interface_width=.12, max_cells=2, dt=.0075)
    a = AttributeSimulation(c, 'direct')
    if mode != 'nonpolar':
        a.polarity[:] = [.2, .1, -.1]
    if mode == 'division':
        a.divide(0, direction=[1, 0, 0])
    a.checkpoint(tmp_path/'source.npz')
    b = NativeSimulation.restore(tmp_path/'source.npz')
    b.native_threads = threads
    for _ in range(240 if mode == 'division' else 24):
        kwargs = dict(prescribed_signals=(a.activator.copy(), a.inhibitor.copy())) if mode == 'prescribed' else {}
        a.step(**kwargs)
        b.step(**kwargs)
        if mode == 'division' and len(a.ids) == 2 and not a.divisions:
            break
    if mode == 'division':
        assert len(a.ids) == len(b.ids) == 2 and not b.divisions
    for key in ('phi', 'activator', 'inhibitor', 'polarity', 'ids'):
        np.testing.assert_array_equal(getattr(a, key), getattr(b, key), err_msg=key)
    assert a.lineage == b.lineage
    assert a.rng.bit_generator.state == b.rng.bit_generator.state
    assert b._fast_cache is None
    b.checkpoint(tmp_path/'native.npz')
    restarted = NativeSimulation.restore(tmp_path/'native.npz')
    restarted.native_threads = threads
    b.step(); restarted.step()
    np.testing.assert_array_equal(b.phi, restarted.phi)


def test_invalid_threads_and_layout_fail_before_native_call():
    sim = NativeSimulation(Config(grid=12, interface_width=.2), 'direct')
    sim.native_threads = 0
    with pytest.raises(ValueError, match='threads'):
        sim.step()
    sim.native_threads = 1
    sim.phi = sim.phi[..., ::-1]
    with pytest.raises(ValueError, match='contiguous'):
        sim.native_exposure()


@pytest.mark.parametrize('clamp,contrast', [(False, .35), (True, .35), (False, 0.)])
def test_cue_hooks_and_checkpoint_preserved(tmp_path, clamp, contrast):
    from embryo.fertilization_cue import CueSimulation
    from embryo.native_cue import NativeCueSimulation
    a = CueSimulation(Config(grid=24, interface_width=.12, max_cells=2, dt=.0075,
                             polarity_tension=contrast), 'direct', strength=.4, clamp=clamp)
    a.checkpoint(tmp_path/'cue.npz')
    b = NativeCueSimulation.restore(tmp_path/'cue.npz')
    b.native_threads = 2
    for _ in range(110):
        a.step(); b.step()
    for key in ('phi', 'activator', 'inhibitor', 'polarity'):
        np.testing.assert_array_equal(getattr(a, key), getattr(b, key), err_msg=key)
    assert a.cue_delivered == b.cue_delivered
    assert b.cue_delivered == pytest.approx(.4)
    b.checkpoint(tmp_path/'native-cue.npz')
    restarted = NativeCueSimulation.restore(tmp_path/'native-cue.npz')
    b.step(); restarted.step()
    np.testing.assert_array_equal(b.phi, restarted.phi)
    assert b.cue == restarted.cue
