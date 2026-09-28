import numpy as np
from embryo.attribute_persistence import perturb, rhs, trajectory, distance


def test_perturbation_preserves_species_amounts():
    values=np.array([[.1,2.,3.],[2.,.5,4.]])
    masses=np.array([1.,2.,5.])
    altered=perturb(values,masses,.05,np.random.default_rng(7))
    assert np.all(altered>0)
    np.testing.assert_allclose(altered@masses,values@masses,rtol=1e-14)
    assert not np.allclose(altered,values)


def test_diffusion_conserves_amount_and_homogeneous_state_is_stationary():
    masses=np.array([1.,2.])
    delta=np.array([[-1.,1.],[.5,-.5]])
    y=np.array([[.8,1.2],[1.3,.7]])
    fun=rhs(delta,2.,.02,.4)
    reaction=rhs(np.zeros((2,2)),2.,.02,.4)
    np.testing.assert_allclose((fun(0,y.ravel())-reaction(0,y.ravel())).reshape(2,2)@masses,0,atol=1e-15)
    path=trajectory(fun,np.ones((2,2)),np.linspace(0,10,21))
    np.testing.assert_allclose(path,1,atol=1e-13)
    np.testing.assert_allclose(distance(path,np.ones((2,2)),masses),0,atol=1e-13)
