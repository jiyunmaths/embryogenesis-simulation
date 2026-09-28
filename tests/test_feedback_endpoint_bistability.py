import numpy as np
from embryo.feedback_endpoint_bistability import chemical_jacobian
from embryo.attribute_persistence import rhs
from embryo.signaling import mode_growth


def test_pattern_jacobian_matches_derivative_and_uniform_dispersion():
    delta=np.array([[-1.,1.],[1.,-1.]])
    y=np.array([[.3,1.4],[.6,1.2]])
    fun=rhs(delta,2.,.02,.4);eps=1e-6
    numerical=np.column_stack([(fun(0,y.ravel()+eps*e)-fun(0,y.ravel()-eps*e))/(2*eps) for e in np.eye(4)])
    np.testing.assert_allclose(chemical_jacobian(y,delta,2.,.02,.4),numerical,rtol=1e-8,atol=1e-9)
    spectrum=np.linalg.eigvals(chemical_jacobian(np.ones((2,2)),delta,2.,.02,.4))
    assert abs(spectrum.real.max()-mode_growth([0,2],2,.02,.4).max())<1e-12
