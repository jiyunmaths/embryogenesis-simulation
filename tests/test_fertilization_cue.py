import numpy as np
import pytest
from embryo.fertilization_cue import CueSimulation,cue_increment,observe
from embryo.attribute_development import AttributeSimulation
from embryo.model import Config


def config(**kw):return Config(grid=24,interface_width=.12,dt=.0075,max_cells=2,**kw)


def test_cue_integral_cutoff_and_duration():
    for dt in (.0075,.00375,.0021):
        assert sum(cue_increment(t,dt,.4,.75) for t in np.arange(0,1.5,dt))==pytest.approx(.4,abs=1e-13)
    assert cue_increment(.75,.1,.4,.75)==0
    assert cue_increment(10,.1,.4,.75)==0
    with pytest.raises(ValueError):CueSimulation(config(),duration=1.9)


def test_zero_cue_exactly_preserves_existing_dynamics():
    a=AttributeSimulation(config(),'direct');b=CueSimulation(config(),'direct')
    for _ in range(4):a.step();b.step()
    for key in ('phi','activator','inhibitor','polarity','due'):np.testing.assert_array_equal(getattr(a,key),getattr(b,key))
    assert a.rng.bit_generator.state==b.rng.bit_generator.state


def test_no_polarity_mechanics_control_is_identical():
    a=CueSimulation(config(polarity_tension=0.),strength=0.,duration=.03)
    b=CueSimulation(config(polarity_tension=0.),strength=.4,duration=.03)
    for _ in range(8):a.step();b.step()
    for key in ('phi','activator','inhibitor'):np.testing.assert_array_equal(getattr(a,key),getattr(b,key))
    assert b.cue_delivered==pytest.approx(.4)
    assert np.linalg.norm(b.polarity)>0


def test_checkpoint_preserves_cue_and_does_not_reapply_past_pulse(tmp_path):
    a=CueSimulation(config(),strength=.1,duration=.03,direction=(1,2,3))
    a.step();a.step();path=tmp_path/'state.npz';a.checkpoint(path);b=CueSimulation.restore(path)
    for _ in range(6):a.step();b.step()
    for key in ('phi','activator','inhibitor','polarity'):np.testing.assert_array_equal(getattr(a,key),getattr(b,key))
    assert a.cue==b.cue and a.cue_delivered==b.cue_delivered
    assert b.cue_delivered==pytest.approx(.1)


def test_spherical_zygote_does_not_report_an_axis():
    r=observe(CueSimulation(config()))
    assert r['shape_axis_alignment'] is None and r['chemical_dipole_alignment'] is None
