import numpy as np
import pytest

from embryo.skew_flux import linear_patch, divergence, symbol, evolve, manufactured, prepare, run


def test_face_flux_exact_for_affine_fields_on_all_skewed_faces():
    for shear in (0.,.25,.5,1.):
        assert max(r['absolute_error_per_area'] for r in linear_patch(shear))<1e-12


def test_flux_pair_conservation_and_operator_spectrum_agree():
    rng=np.random.default_rng(123)
    field=rng.normal(size=(8,8,8))
    derivative=divergence(field,1.)
    expected=np.fft.ifftn(np.fft.fftn(field)*symbol(8,1.)).real
    np.testing.assert_allclose(derivative,expected,atol=1e-11)
    assert abs(derivative.sum())<1e-10
    np.testing.assert_allclose(divergence(np.ones_like(field),1.),0,atol=1e-12)
    assert symbol(8,1.).max()<=1e-12


def test_exact_semigroup_exposes_spatial_positivity_failure():
    pulse=np.zeros((16,16,16));pulse[8,8,8]=1
    time=.01/16**2
    assert evolve(pulse,0.,time).min()>-1e-12
    final=evolve(pulse,1.,time)
    assert final.min() < -1e-3
    assert final.sum()==pytest.approx(1.,abs=1e-12)


def test_manufactured_solution_errors_decrease_with_refinement():
    errors=[]
    for n in (8,16,32):
        start,exact=manufactured(n,.5,.02)
        errors.append(np.linalg.norm((evolve(start,.5,.02)-exact).ravel())/np.linalg.norm(exact.ravel()))
    assert errors[0]>errors[1]>errors[2]
    assert np.log2(errors[1]/errors[2])>1.8


def test_full_report_retains_failed_positivity_and_refuses_overwrite(tmp_path):
    output=tmp_path/'study';prepare(output)
    frozen=(output/'protocol.json').read_bytes()
    result=run(output)
    assert sum(result['checks'].values())==7
    assert not result['checks']['all_positive_pulses_remain_nonnegative']
    assert not result['passed'] and not result['live_integration_ready']
    assert (output/'protocol.json').read_bytes()==frozen
    with pytest.raises(FileExistsError):prepare(output)
    with pytest.raises(ValueError):run(output)
