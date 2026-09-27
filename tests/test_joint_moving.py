import numpy as np
import pytest

from embryo.model import Config
from embryo.resolution import manufactured
from embryo.causal_signaling import reaction
from embryo.joint_fate import rhs, advance, discrepancies
from embryo.moving_causal import initialize, JointMovingSimulation, prepare


def test_deviation_rhs_matches_reactions_and_joint_fate_drift():
    c=Config();delta=np.zeros((3,3));state=np.array([[.1,-.2,.03],[.2,.1,-.02],[.1,0,-.2]])
    for arm in ('full','no_self_activation','no_transport','equal_diffusion','no_signal_to_fate'):
        value=rhs(state,delta,c,arm)
        expected=reaction(1+state[0],1+state[1],c.signal_beta,arm)
        np.testing.assert_allclose(value[:2],expected,atol=1e-15)
        gain=0 if arm=='no_signal_to_fate' else c.signal_fate_gain
        np.testing.assert_allclose(value[2],c.fate_rate*(state[2]-state[2]**3+gain*state[0]))


def test_moving_step_joint_fate_and_single_amount_dilution(tmp_path):
    source=tmp_path/'source.npz';manufactured(Config(grid=32,extent=2.24,dt=.0075)).checkpoint(source)
    sim=initialize(source,'full',7)
    old=sim.volumes();expected=advance(sim.joint_state,sim.signaling_graph().delta,sim.config,'full',sim.config.dt)
    expected_amounts=[old@(1+expected[k]) for k in (0,1)]
    sim.step();new=sim.volumes()
    np.testing.assert_allclose([new@sim.activator,new@sim.inhibitor],expected_amounts,rtol=1e-12,atol=1e-12)
    np.testing.assert_array_equal(sim.fate,expected[2])
    sim.checkpoint(tmp_path/'joint.npz');restored=JointMovingSimulation.restore(tmp_path/'joint.npz')
    np.testing.assert_array_equal(sim.joint_state,restored.joint_state)
    sim.step();restored.step()
    np.testing.assert_array_equal(sim.phi,restored.phi)
    np.testing.assert_array_equal(sim.joint_state,restored.joint_state)


def test_interventions_and_validation_gate(tmp_path):
    source=tmp_path/'source.npz';manufactured(Config(grid=32,extent=2.24)).checkpoint(source)
    off=initialize(source,'no_signal_to_fate',7);off.step()
    np.testing.assert_array_equal(off.fate,0)
    no_self=initialize(source,'no_self_activation',7)
    assert not no_self.graph_snapshot()['unstable_modes']
    assert not initialize(source,'no_feedback',7).config.feedback
    validation=tmp_path/'validation';validation.mkdir();(validation/'comparison.json').write_text('{"passed": false}')
    with pytest.raises(ValueError,match='must pass'):
        prepare(source,validation,tmp_path/'blocked')
