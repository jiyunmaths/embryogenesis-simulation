import numpy as np
import pytest
from embryo.attribute_development import AttributeSimulation
from embryo.feedback_long import ARMS, initialize, save_checkpoint, restore_checkpoint
from embryo.model import Config


def source(tmp_path):
    sim=AttributeSimulation(Config(grid=24,max_cells=1,interface_width=.12))
    sim.polarity[:]=[.1,0,0]
    path=tmp_path/'source.npz';sim.checkpoint(path)
    return sim,path


def test_arms_share_geometry_and_change_only_declared_coefficients(tmp_path):
    original,path=source(tmp_path);signals=np.array([[1.2],[.8]])
    for arm,expected in ARMS.items():
        sim=initialize(path,signals,arm)
        np.testing.assert_array_equal(sim.phi,original.phi)
        np.testing.assert_array_equal(sim.polarity,original.polarity)
        np.testing.assert_array_equal(sim.activator,signals[0])
        assert (sim.config.fate_tension,sim.config.fate_adhesion,sim.config.polarity_tension)==expected
        assert not sim.config.differentiation
        assert np.all(sim.fate==0)


def test_atomic_checkpoint_restores_exact_dynamics_and_rejects_wrong_job(tmp_path):
    original,path=source(tmp_path);sim=initialize(path,np.array([[1.2],[.8]]),'full')
    sim.step();job={'state':'formation','arm':'full'};audit={'steps':1};history=[{'time':sim.time}]
    checkpoint=tmp_path/'latest.npz';save_checkpoint(sim,checkpoint,audit,history,'hash',job)
    restored,a,h=restore_checkpoint(checkpoint,'hash',job)
    assert a==audit and h==history
    for _ in range(2):sim.step();restored.step()
    for field in ['phi','activator','inhibitor','polarity','fate']:
        np.testing.assert_array_equal(getattr(sim,field),getattr(restored,field))
    with pytest.raises(ValueError):restore_checkpoint(checkpoint,'other',job)
    with pytest.raises(ValueError):restore_checkpoint(checkpoint,'hash',{'state':'persistence','arm':'full'})
