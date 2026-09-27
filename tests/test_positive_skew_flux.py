import numpy as np
import pytest

from embryo.skew_flux import metric
from embryo.positive_skew_flux import directions, divergence, face_fluxes, linear_patch, positivity_step, ssprk2, prepare, run


@pytest.mark.parametrize('shear',[-1.,-.5,0.,.25,1.])
def test_directional_tensor_and_face_fluxes_are_consistent(shear):
    terms=directions(shear)
    reconstructed=sum(w*np.outer(d,d) for d,w in terms)
    np.testing.assert_allclose(reconstructed,metric(shear)[1],atol=1e-14)
    assert all(w>=0 for _,w in terms)
    assert max(r['error_per_area'] for r in linear_patch(shear))<1e-12
    u=np.random.default_rng(7).random((8,8,8))
    q=face_fluxes(u,shear)
    np.testing.assert_allclose(8**3*sum(f-np.roll(f,1,a) for a,f in enumerate(q)),divergence(u,shear),atol=1e-11)


def test_positive_generator_and_step_bound_preserve_mass_and_bounds():
    pulse=np.zeros((8,8,8));pulse[4,4,4]=1
    rates=divergence(pulse,1.);rates[4,4,4]=0
    assert rates.min()>=0
    bound=positivity_step(8,1.)
    value=ssprk2(pulse,1.,10*.9*bound,10)
    assert value.min()>=0 and value.max()<=1
    assert value.sum()==pytest.approx(1.,abs=1e-12)
    with pytest.raises(ValueError,match='positivity bound'):
        ssprk2(pulse,1.,1.1*bound,1)
    with pytest.raises(ValueError,match='shear'):
        directions(1.1)


def test_complete_benchmark_passes_without_claiming_live_readiness(tmp_path):
    out=tmp_path/'study';prepare(out)
    frozen=(out/'protocol.json').read_bytes()
    result=run(out)
    assert result['passed'] and all(result['checks'].values())
    assert not result['live_integration_ready']
    assert (out/'protocol.json').read_bytes()==frozen
    with pytest.raises(ValueError):run(out)
    with pytest.raises(FileExistsError):prepare(out)
