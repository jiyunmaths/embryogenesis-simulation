"""Independent accuracy, amount balance and accumulation-state checks."""
from decimal import Decimal,localcontext
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import torch
from embryo.compensated_chemistry import NumpyStepper,TorchStepper,numpy_add


def system():
    v=np.array([1.,2.,.4,1.6]);g=np.zeros((4,4))
    for i,j,c in [(0,1,.3),(1,2,.7),(2,3,.2),(3,0,.4),(1,3,.15)]:g[i,j]=g[j,i]=c
    d=g.copy();np.fill_diagonal(d,-g.sum(axis=1));d/=v[:,None]
    u=np.array([[.8,1.3,.7,1.1],[1.1,.7,1.2,.8]])
    return d,v,u


def integrate(stepper,u,steps):
    e=np.zeros_like(u)
    for _ in range(steps):u,e=stepper(u,e)
    return u,e


def test_subulp_increments_are_retained_against_decimal():
    u=np.array([1.]);e=np.zeros_like(u);increment=np.array([2.**-56])
    naive=u.copy()
    for _ in range(128):u,e=numpy_add(u,e,increment);naive+=increment
    with localcontext() as ctx:
        ctx.prec=80;exact=Decimal(1)+128*Decimal.from_float(float(increment[0]))
    assert naive[0]==1. and u[0]==float(exact) and e[0]==0.


def test_ssp_order_against_independent_nonlinear_solution():
    d,v,u=system()
    def rhs(t,y):
        a,b=y.reshape(2,4);return np.stack([a*a/b-a+.02*(d@a),2*(a*a-b)+.55*(d@b)]).ravel()
    reference=solve_ivp(rhs,(0,2),u.ravel(),method='DOP853',rtol=1e-13,atol=1e-15,max_step=.015)
    assert reference.success
    exact=reference.y[:,-1].reshape(u.shape);errors=[]
    for dt in (.02,.01,.005):
        result,_=integrate(NumpyStepper(d,dt),u.copy(),round(2/dt))
        errors.append(float(abs(np.log(result/exact)).max()))
    orders=np.log2(np.array(errors[:-1])/errors[1:])
    assert np.all((orders>1.8)&(orders<2.2))


def test_conservative_exchange_matches_exponential_with_unequal_volumes():
    d,v,u=system();step=NumpyStepper(d,.005,reaction=False)
    result,e=integrate(step,u.copy(),2000)
    np.testing.assert_allclose(result@v,u@v,rtol=2e-13,atol=0)
    exact=np.stack([expm(10*D*d)@a for D,a in zip((.02,.55),u)])
    assert abs(np.log(result/exact)).max()<1e-7
    assert np.all(result>0) and np.isfinite(e).all()


def test_restart_requires_and_preserves_correction(tmp_path):
    d,v,u=system();step=NumpyStepper(d,.001875);e=np.zeros_like(u)
    for _ in range(71):u,e=step(u,e)
    assert np.any(e!=0)
    np.savez(tmp_path/'checkpoint.npz',state=u,correction=e)
    continued=u.copy();ce=e.copy()
    for _ in range(151):continued,ce=step(continued,ce)
    with np.load(tmp_path/'checkpoint.npz') as z:resumed=z['state'].copy();re=z['correction'].copy()
    for _ in range(151):resumed,re=step(resumed,re)
    np.testing.assert_array_equal(continued,resumed);np.testing.assert_array_equal(ce,re)


def test_torch_cpu_and_independent_numpy_trajectories_agree():
    d,v,u=system();numpy=NumpyStepper(d,.00375);torch_step=TorchStepper(torch.from_numpy(d),.00375)
    x=u.copy();e=np.zeros_like(u);y=torch.from_numpy(u.copy());c=torch.zeros_like(y)
    for _ in range(800):x,e=numpy(x,e);y,c=torch_step(y,c)
    assert abs(np.log(x/y.numpy())).max()<1e-12
    assert np.isfinite(e).all() and torch.isfinite(c).all()


def test_large_timestep_is_positivity_substepped():
    d,v,u=system();step=NumpyStepper(d,2.)
    assert step.count>1
    result,e=step(u,np.zeros_like(u));assert np.all(result>0)


@pytest.mark.parametrize('bad',['integer','negative','nan_correction','wrong_shape','negative_dt'])
def test_invalid_or_nonpositive_states_are_not_clipped(bad):
    d,v,u=system();e=np.zeros_like(u)
    if bad=='negative_dt':
        with pytest.raises(ValueError):NumpyStepper(d,-.01)
        return
    if bad=='integer':u=u.astype(np.int64)
    if bad=='negative':u[0,0]=-1
    if bad=='nan_correction':e[0,0]=np.nan
    if bad=='wrong_shape':e=e[:,:-1]
    with pytest.raises((ValueError,FloatingPointError)):NumpyStepper(d,.01)(u,e)
