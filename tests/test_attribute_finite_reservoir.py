import numpy as np
from scipy.linalg import expm
from embryo.attribute_finite_reservoir import system, integrate


def test_exchange_conserves_amount_and_matches_exact_linear_evolution():
    v=np.array([1.,2.,3.]);rv=4.;initial=np.array([[.3,2.,1.,.8],[2.,.2,.7,1.3]])
    fun,jac=system(v,rv,2.,.1,2.,reactions=False)
    mass=np.append(v,rv)
    np.testing.assert_allclose(fun(0,initial.ravel()).reshape(2,4)@mass,0,atol=1e-14)
    times=np.linspace(0,20,21);path=integrate(fun,jac,initial,times,True)
    exact=np.array([(expm(t*jac(0,initial.ravel()))@initial.ravel()).reshape(2,4) for t in times])
    np.testing.assert_allclose(path,exact,rtol=1e-9,atol=1e-11)
    np.testing.assert_allclose(path@mass,np.broadcast_to(initial@mass,(len(times),2)),atol=1e-10)


def test_full_jacobian_and_reaction_amount_balance():
    v=np.array([1.,2.]);rv=.3;state=np.array([[.2,2.,1.],[.7,3.,.8]])
    fun,jac=system(v,rv,2.,.05,.9)
    y=state.ravel();eps=1e-6
    numeric=np.column_stack([(fun(0,y+eps*e)-fun(0,y-eps*e))/(2*eps) for e in np.eye(len(y))])
    np.testing.assert_allclose(jac(0,y),numeric,atol=2e-9,rtol=1e-8)
    a,b=state[:,:2];reaction=np.array([a*a/b-a,2*(a*a-b)])
    np.testing.assert_allclose(fun(0,y).reshape(2,3)@np.append(v,rv),reaction@v,atol=1e-14)


def test_uniform_equilibrium_and_fixed_reservoir_limit():
    fun,jac=system(np.array([1.,2.]),3.,2.,.1,2.)
    np.testing.assert_array_equal(fun(0,np.ones(6)),0)
    fixed,_=system(np.array([1.,2.]),3.,2.,.1,2.,clamped=True)
    slow,_=system(np.array([1.,2.]),3e12,2.,.1,2.)
    state=np.array([.5,2.,1.,.3,3.,1.])
    np.testing.assert_allclose(slow(0,state),fixed(0,state),atol=1e-11)


def test_weighted_zero_mean_contrast_has_reservoir_independent_growth():
    v=np.array([1.,2.]);contrast=np.array([2.,-1.])
    ka=.2621471462181868;kb=5.242942924363736
    local=np.array([[1-ka,-1],[4,-2-kb]])
    assert np.linalg.eigvals(local).real.max()>0
    perturbation=np.array([*contrast,0.,*(.3*contrast),0.])
    expected=local@np.array([1.,.3])
    for rv in [.1,1.,100.]:
        _,jac=system(v,rv,2.,ka,kb)
        out=(jac(0,np.ones(6))@perturbation).reshape(2,3)
        np.testing.assert_allclose(out[:,:2],expected[:,None]*contrast,atol=1e-14)
        np.testing.assert_allclose(out[:,-1],0,atol=1e-14)
