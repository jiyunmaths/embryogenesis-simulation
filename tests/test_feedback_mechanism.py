import numpy as np
import pytest
from embryo.feedback_response import ResponseSimulation
from embryo.model import Config
from embryo.signaling import mode_growth
from embryo.transport import conservative_transport,transport_graph


def test_uniform_signals_annihilate_geometry_operator_changes():
    v=np.array([1.,2.,3.]);g=np.array([[0.,1.,.3],[1.,0.,2.],[.3,2.,0.]])
    first=transport_graph(conservative_transport(g,v))
    second=transport_graph(conservative_transport(g*2,v))
    np.testing.assert_allclose((second.delta-first.delta)@np.ones(3),0,atol=1e-15)
    np.testing.assert_allclose(second.eigenvalues,2*first.eigenvalues,atol=1e-14)
    # The constitutive coefficients do not enter this fixed-geometry dispersion relation.
    assert mode_growth([5.41],2,.02,.4)[0]<0
    assert mode_growth([6.636],2,.02,.4)[0]>0


def test_extended_law_is_positive_and_linear_high_strength_is_rejected():
    sim=ResponseSimulation(Config(grid=16, interface_width=.15))
    sim.config.fate_tension=2.;sim.config.fate_adhesion=2.8
    sim.activator[:]=.1
    with pytest.raises(ValueError):sim.material_coefficients()
    sim.law='exponential'
    tension,adhesion=sim.material_coefficients()
    assert np.all(tension>0)
    sim.activator[:]=1.
    np.testing.assert_allclose(sim.material_coefficients()[0],sim.config.surface_tension)
