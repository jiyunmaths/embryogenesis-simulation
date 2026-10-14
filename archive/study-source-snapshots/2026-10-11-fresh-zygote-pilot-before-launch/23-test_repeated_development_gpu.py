"""Tests of label-invariant numerical comparisons and zero-polarity handling."""
import numpy as np
import pytest
from embryo.qualify_repeated_development import (spectrum, compare_observations,
    densities, density_errors, homologous_errors)


def example():
    v=np.array([.2,.8]); k=np.array([[1.,-1.],[-1.,1.]])
    return dict(time=0.,ids=[5,8],volumes=v.tolist(),chemistry=[[.2,1.4],[.5,1.1]],
        polarity=[[.1,0.,0.],[-.025,0.,0.]],delta=(-k/v[:,None]).tolist(),axis_ratio=1.2)


def test_graph_spectrum_uses_volume_weighted_symmetry():
    r=example(); eigen=spectrum(r['volumes'],r['delta'])
    np.testing.assert_allclose(eigen,[0.,6.25],atol=1e-14)
    assert not np.allclose(r['delta'],np.array(r['delta']).T)


def test_signaling_and_geometry_comparison_does_not_depend_on_numeric_cell_ids():
    a=example(); order=np.array([1,0]); b=dict(a,ids=[100,200],
        chemistry=np.array(a['chemistry'])[:,order].tolist(),
        volumes=np.array(a['volumes'])[order].tolist(),
        polarity=np.array(a['polarity'])[order].tolist(),
        delta=np.array(a['delta'])[np.ix_(order,order)].tolist())
    error,unmatched=compare_observations([a],[b])
    assert unmatched==[] and max(error.values())<1e-14


def test_count_mismatch_keeps_whole_amount_and_distribution_checks():
    a=dict(time=0.,ids=[0],volumes=[1.],chemistry=[[1.],[2.]],polarity=[[0.,0.,0.]],delta=[[0.]],axis_ratio=1.)
    b=dict(time=0.,ids=[1,2],volumes=[.5,.5],chemistry=[[1.,1.],[2.,2.]],
        polarity=[[0.,0.,0.],[0.,0.,0.]],delta=[[-1.,1.],[1.,-1.]],axis_ratio=1.)
    error,unmatched=compare_observations([a],[b])
    assert unmatched==[0.] and error['cell_count_difference_max']==1
    assert error['whole_amount_relative_max']==error['concentration_wasserstein_max']==0
    b['chemistry'][0]=[2.,2.]
    error,_=compare_observations([a],[b]);assert error['whole_amount_relative_max']==.5


def test_field_comparison_is_permutation_invariant_and_keeps_spatial_information():
    phi=np.zeros((2,8,8,8),dtype=np.float32);phi[0,:4]=.8;phi[1,4:]=.8
    u=np.array([[.2,1.4],[.5,1.1]]);p=np.array([[.1,0,0],[-.025,0,0]])
    a=densities(phi,u,p);b=densities(phi[::-1],u[:,::-1],p[::-1])
    assert max(density_errors(a,b).values())==0
    c=densities(phi,u[:,::-1],p)
    assert density_errors(a,c)['activator_density']>.1
    # Equal distributions alone cannot detect this spatial exchange.
    r=example();s=dict(r,chemistry=np.array(r['chemistry'])[:,::-1].tolist(),volumes=[.5,.5])
    r=dict(r,volumes=[.5,.5],delta=[[-1.,1.],[1.,-1.]])
    s=dict(s,delta=r['delta']);error,_=compare_observations([r],[s])
    assert error['concentration_wasserstein_max']==0


def test_zero_polarity_field_has_finite_error_without_normalizing_by_zero():
    phi=np.ones((1,4,4,4),dtype=np.float32);u=np.ones((2,1))
    a=densities(phi,u,np.zeros((1,3)));b=densities(phi,u,np.array([[.1,0,0]]))
    error=density_errors(a,b)
    assert error['polarity_density']==pytest.approx(.1)
    assert all(np.isfinite(list(error.values())))


def test_corresponding_cells_reveal_a_spatial_exchange_hidden_by_equal_distributions():
    a=dict(example(),volumes=[.5,.5],delta=[[-1.,1.],[1.,-1.]],parents=[0,0])
    b=dict(a,chemistry=np.array(a['chemistry'])[:,::-1].tolist())
    aggregate,_=compare_observations([a],[b]);raw,unmatched=homologous_errors([a],[b])
    assert aggregate['concentration_wasserstein_max']==0 and raw['chemical_log_max']>1 and unmatched==[]
    b['parents']=[1,1];raw,unmatched=homologous_errors([a],[b])
    assert unmatched==[0.] and raw['chemical_log_max']==0


@pytest.mark.parametrize('delta',[[[-1.,1.],[.1,-.1]],[[-1.,.9],[.25,-.25]],[[1.,-1.],[-.25,.25]]])
def test_nonconservative_or_negative_exchange_cannot_supply_an_accepted_spectrum(delta):
    with pytest.raises(FloatingPointError):spectrum([.2,.8],delta)
