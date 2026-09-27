import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from embryo.skew_boundary import closed_transport, manufactured, boundary_reference, evolve_manufactured, prepare, run


def test_closed_operator_has_no_periodic_wrap_and_preserves_mass():
    n=8;t=closed_transport(n,.5)
    i=np.arange(n**3).reshape((n,n,n))
    assert t.conductance[i[0,3,3],i[-1,3,3]]==0
    assert t.conductance[i[3,0,3],i[3,-1,3]]==0
    assert t.conductance.data.min()>=0
    np.testing.assert_allclose(t.delta@np.ones(n**3),0,atol=1e-12)
    np.testing.assert_allclose(t.volumes@t.delta,0,atol=1e-12)


@pytest.mark.parametrize('shear',[0.,.5,1.])
def test_exact_boundary_reference_exercises_tangential_gradients(shear):
    reference=boundary_reference(shear)
    assert reference['max_conormal_residual']<1e-12
    assert reference['max_logical_tangential_derivative_on_xi_wall']>.001
    f,laplace=manufactured(16,shear)
    assert abs(f.mean())<1e-12
    assert abs(laplace.mean())<1e-12


def test_wall_pulse_positive_and_constant_equilibrium_stationary():
    n=8;t=closed_transport(n,1.)
    pulse=np.zeros(n**3);pulse[0]=1.
    final=expm_multiply(t.delta*.001,pulse)
    assert final.min()>=-1e-12
    assert final.sum()==pytest.approx(1.,abs=1e-12)
    np.testing.assert_allclose(expm_multiply(t.delta*.01,np.ones(n**3)),1.,atol=1e-12)


def test_orthogonal_manufactured_error_decreases():
    errors=[]
    for n in (8,16):
        t=closed_transport(n,0.)
        _,final,exact,_=evolve_manufactured(t,n,0.,.02)
        errors.append(np.linalg.norm(final-exact)/np.linalg.norm(exact))
    assert errors[1]<errors[0]/3


def test_protocol_refuses_overwrite_and_modified_config(tmp_path):
    path=tmp_path/'study';prepare(path)
    with pytest.raises(FileExistsError):prepare(path)
    (path/'protocol.json').write_text((path/'protocol.json').read_text()+' ')
    with pytest.raises(ValueError,match='unchanged'):run(path)
