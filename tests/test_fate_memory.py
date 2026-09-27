import numpy as np
from scipy.integrate import solve_ivp
from embryo.fate_memory import neutral_fate, withdrawal_history, joint_rhs
from embryo.joint_fate import rhs
from embryo.model import Config


def test_exact_neutral_flow_matches_independent_ode_and_preserves_zero():
    initial=np.array([-1.4,-.4,-1e-5,0.,1e-5,.4,1.4]);times=np.linspace(0,30,61)
    for bistable in [True,False]:
        solution=solve_ivp(lambda t,z:.8*(z-z**3 if bistable else -z),(0,30),initial,t_eval=times,method='DOP853',rtol=1e-11,atol=1e-13)
        exact=neutral_fate(initial,times[:,None],.8,bistable)
        np.testing.assert_allclose(exact,solution.y.T,atol=3e-8,rtol=3e-8)
        np.testing.assert_array_equal(exact[:,3],0)
        np.testing.assert_array_equal(exact[0],initial)


def test_withdrawal_has_no_state_jump_and_keeps_preceding_history():
    times=np.array([0.,.6,1.2]);driven=np.arange(3*4*2*2,dtype=float).reshape(3,4,2,2)/100
    for bistable,index in [(True,2),(False,3)]:
        result=withdrawal_history(driven,times,.6,.8,bistable)
        np.testing.assert_array_equal(result[:2],driven[:2,index])
        np.testing.assert_allclose(result[-1],neutral_fate(driven[1,index],.6,.8,bistable))


def test_joint_control_does_not_change_chemistry_or_bistable_rhs():
    state=np.array([[.01,-.01],[.02,-.02],[.1,-.1],[.3,-.3]])
    delta=np.array([[-1.,1.],[1.,-1.]])
    c=Config();value=joint_rhs(state,delta,c)
    np.testing.assert_array_equal(value[:3],rhs(state[:3],delta,c,'full'))
    np.testing.assert_allclose(value[3],c.fate_rate*(-state[3]+c.signal_fate_gain*state[0]))
