from copy import deepcopy
import numpy as np
import pytest
from embryo.attribute_development import AttributeSimulation, attributes
from embryo.model import Config, Simulation


def pair(sim):
    original=sim.phi[0].copy()
    sim.phi=np.stack([np.roll(original,3,axis=0),np.roll(original,-3,axis=0)])
    sim.target=sim.volumes();sim.ids=np.array([0,1]);sim.parents=np.array([-1,-1]);sim.due=np.array([10.,10.])
    sim.activator=np.array([1.2,.8]);sim.inhibitor=np.ones(2);sim.fate=np.zeros(2)
    sim.polarity=np.array([[.1,0,0],[-.1,0,0]])
    return sim


@pytest.mark.parametrize('mode',['direct','no_feedback'])
def test_material_substitution_matches_core_mechanics(mode):
    sim=pair(AttributeSimulation(Config(grid=24,dt=.0075),mode))
    reference=Simulation(deepcopy(sim.config))
    for key in ('phi','target','ids','parents','due','activator','inhibitor','polarity'):
        setattr(reference,key,getattr(sim,key).copy())
    reference.fate=sim.activator-1
    sim.mechanical_step();reference.mechanical_step()
    np.testing.assert_array_equal(sim.phi,reference.phi)


@pytest.mark.parametrize('mode',['direct','no_feedback'])
def test_no_fate_drift_or_labels_and_checkpoint_replay(tmp_path,mode):
    sim=AttributeSimulation(Config(grid=24,dt=.0075),mode)
    sim.fate[:]=.9
    sim.step()
    np.testing.assert_array_equal(sim.fate,0)
    row=attributes(sim)
    assert not {'fate_a','fate_b','uncommitted','fate_separation'}&row['metrics'].keys()
    assert 'fate' not in sim.surfaces(include_mesh=False)[0]
    sim.checkpoint(tmp_path/'state.npz');restored=AttributeSimulation.restore(tmp_path/'state.npz')
    assert restored.attribute_mode==mode
    sim.step();restored.step()
    for name in ('phi','activator','inhibitor','polarity'):np.testing.assert_array_equal(getattr(sim,name),getattr(restored,name))


def test_material_response_is_instantaneous_not_a_stored_identity():
    sim=pair(AttributeSimulation(Config(grid=24,dt=.0075)))
    tension,adhesion=sim.material_coefficients()
    assert tension[0]>sim.config.surface_tension>tension[1]
    sim.fate[:]=[1e6,-1e6]
    np.testing.assert_array_equal(sim.material_coefficients()[0],tension)
    sim.activator[:]=1
    neutral,attraction=sim.material_coefficients()
    np.testing.assert_array_equal(neutral,sim.config.surface_tension)
    assert attraction[0,1]==sim.config.adhesion


def test_active_cleavage_restores_and_daughters_have_no_fate(tmp_path):
    sim=AttributeSimulation(Config(grid=24,dt=.0075,max_cells=2))
    sim.divide(0,np.array([1.,0.,0.]))
    for _ in range(8):sim.step()
    sim.checkpoint(tmp_path/'dividing.npz');restored=AttributeSimulation.restore(tmp_path/'dividing.npz')
    sim.step();restored.step()
    np.testing.assert_array_equal(sim.phi,restored.phi)
    for _ in range(500):
        sim.step()
        if len(sim.phi)==2:break
    assert len(sim.phi)==2 and not sim.divisions
    np.testing.assert_array_equal(sim.fate,0)
