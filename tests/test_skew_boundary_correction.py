import numpy as np
import pytest

from embryo.skew_boundary import closed_transport
from embryo.skew_boundary_correction import corrected_transport, reference, boundary_residual, evolve, prepare, run


@pytest.mark.parametrize('shear',[-1.,-.5,0.,.25,1.])
def test_correction_is_nonnegative_symmetric_and_preserves_invariants(shear):
    base=closed_transport(8,shear);fixed=corrected_transport(8,shear)
    extra=fixed.conductance-base.conductance
    assert not extra.nnz or extra.data.min()>=0
    assert (fixed.conductance-fixed.conductance.T).nnz==0
    np.testing.assert_allclose(fixed.delta@np.ones(8**3),0,atol=1e-12)
    np.testing.assert_allclose(fixed.volumes@fixed.delta,0,atol=1e-12)
    if shear==0:assert extra.nnz==0


@pytest.mark.parametrize('family',['xi','eta'])
def test_both_manufactured_families_have_exact_conormal_walls(family):
    for shear in (-1.,-.5,0.,.25,1.):
        assert boundary_residual(shear,family)<1e-12
        f,laplace=reference(16,shear,family)
        assert abs(f.mean())<1e-12
        assert abs(laplace.mean())<1e-12


def test_correction_improves_original_boundary_case_without_fitting():
    errors=[]
    for builder in (closed_transport,corrected_transport):
        t=builder(16,.5)
        _,final,exact,_=evolve(t,16,.5,'xi',.02)
        errors.append(np.linalg.norm(final-exact)/np.linalg.norm(exact))
    assert errors[1]<errors[0]


def test_protocol_immutable_and_no_overwrite(tmp_path):
    path=tmp_path/'study';prepare(path)
    with pytest.raises(FileExistsError):prepare(path)
    (path/'protocol.json').write_text((path/'protocol.json').read_text()+' ')
    with pytest.raises(ValueError,match='unchanged'):run(path)
