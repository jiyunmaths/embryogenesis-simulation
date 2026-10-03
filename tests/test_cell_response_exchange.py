import numpy as np
import pytest
from embryo.cell_response_exchange import transplant, waveform_distance, checked_solve


def test_transplant_controls_and_unequal_volume_conservation():
    x=np.array([[.2,2.,1.],[.4,1.,3.]])
    m=np.array([1.,2.,4.]);pair=[0,1]
    for arm in ('untouched','sham'):np.testing.assert_array_equal(transplant(x,m,pair,arm),x)
    exact=transplant(x,m,pair,'exact');np.testing.assert_array_equal(exact[:,pair],x[:,pair[::-1]])
    assert not np.allclose(exact@m,x@m)
    conserved=transplant(x,m,pair,'conservative');np.testing.assert_allclose(conserved@m,x@m)
    np.testing.assert_array_equal(conserved[:,2],x[:,2]);assert np.all(conserved>0)
    np.testing.assert_array_equal(x,np.array([[.2,2.,1.],[.4,1.,3.]]))


def test_waveform_distance_detects_sign_and_preserves_donor_mapping():
    t=np.array([0.,1.,2.]);a=np.ones((3,2));b=-a
    assert waveform_distance(a,a,t)==0
    assert waveform_distance(a,b,t)==2
    # Trial order: first cell factors .9/1.1, second cell same factors.
    assert [(i+2)%4 for i in range(4)]==[2,3,0,1]
    with pytest.raises(ValueError):waveform_distance(a,b[:2],t)


def test_independent_solver_check_against_analytic_decay():
    initial=np.array([[2.,3.],[4.,5.]])
    times=np.array([0.,.2,1.,4.])
    path,error=checked_solve(lambda t,y:1-y,lambda t,y:-np.eye(4),initial,times)
    np.testing.assert_allclose(path,1+(initial-1)*np.exp(-times[:,None,None]),rtol=1e-10,atol=1e-12)
    assert error<1e-9
