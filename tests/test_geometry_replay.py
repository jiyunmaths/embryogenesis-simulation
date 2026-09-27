import numpy as np
from embryo.geometry_replay import geometry, trajectory
from embryo.joint_fate import advance
from embryo.model import Config


def test_interpolation_retains_symmetric_flux_and_positive_capacity():
    times=np.array([0.,1.]);g=np.array([[[0,1],[1,0]],[[0,3],[3,0]]],float);v=np.array([[1,2],[2,4]],float)
    delta,mass=geometry(times,g,v,.5)
    np.testing.assert_allclose(delta@np.ones(2),0,atol=1e-15)
    np.testing.assert_allclose(mass@delta,0,atol=1e-15)
    np.testing.assert_allclose(mass,[1.5,3.])


def test_constant_history_matches_verified_joint_solver():
    times=np.array([0.,.06]);g=np.array([[[0,1],[1,0]]]*2,float);v=np.ones((2,2));state=np.array([[.001,-.001],[0.,0.],[0.,0.]])
    c=Config(dt=.0075);delta,_=geometry(times,g,v,0);expected=state.copy()
    for _ in range(8):expected=advance(expected,delta,c,'full',c.dt)
    for arm in ('frozen','replay','no_dilution','frozen_transport'):
        result=trajectory(times,g,v,state,c,arm,c.dt,times)
        np.testing.assert_array_equal(result[-1],expected)


def test_single_dilution_matches_amount_balance_after_joint_reaction():
    times=np.array([0.,.0075]);g=np.zeros((2,2,2));v=np.array([[1.,2.],[1.1,2.4]])
    state=np.array([[.01,-.02],[.005,.003],[0.,0.]])
    c=Config(dt=.0075);reacted=advance(state,np.zeros((2,2)),c,'full',c.dt)
    result=trajectory(times,g,v,state,c,'replay',c.dt,times)
    np.testing.assert_allclose((1+result[-1,:2])*v[-1],(1+reacted[:2])*v[0],rtol=1e-14)
    np.testing.assert_array_equal(result[-1,2],reacted[2])
