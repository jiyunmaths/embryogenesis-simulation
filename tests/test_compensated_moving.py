"""Independent multiplication precision and evolving-volume amount checks."""
from decimal import Decimal,localcontext
import numpy as np
import pytest
import torch
from embryo.compensated_moving import numpy_scale,torch_scale,validate_correction,MovingSimulation


def test_scaled_represented_value_matches_decimal():
    u=np.array([[.013,1.17,12.3],[.41,1.001,5.73]]);e=np.spacing(u)*np.array([[.23,-.41,.38],[-.2,.31,.44]])
    r=np.array([1.000031,.999963,1.0027]);v,c=numpy_scale(u,e,r)
    with localcontext() as context:
        context.prec=90
        for i in range(2):
            for j in range(3):
                exact=(Decimal.from_float(u[i,j])-Decimal.from_float(e[i,j]))*Decimal.from_float(r[j])
                reconstructed=Decimal.from_float(v[i,j])-Decimal.from_float(c[i,j])
                assert abs(exact-reconstructed)<=Decimal.from_float(np.spacing(v[i,j]))*Decimal('1e-15')
    validate_correction(v,c)


def test_torch_and_numpy_scale_agree_exactly():
    rng=np.random.default_rng(15);u=np.exp(rng.normal(size=(2,16)));e=.3*np.spacing(u);r=1+.001*rng.normal(size=16)
    a,b=numpy_scale(u,e,r);x,y=torch_scale(torch.from_numpy(u),torch.from_numpy(e),torch.from_numpy(r))
    np.testing.assert_array_equal(a,x.numpy());np.testing.assert_array_equal(b,y.numpy())


def test_pure_dilution_amounts_on_evolving_unequal_volumes():
    volumes=np.array([.8,1.7,.43]);u=np.array([[.27,1.83,4.1],[.7,1.9,.93]]);e=np.zeros_like(u);initial=u*volumes
    for k in range(1500):
        new=volumes*np.exp(.0001*np.sin(k+np.arange(3)));u,e=numpy_scale(u,e,volumes/new);volumes=new
    np.testing.assert_allclose(u*volumes,initial,rtol=2e-14,atol=0);validate_correction(u,e)


def test_checkpoint_without_chemical_correction_rejected_before_gpu(tmp_path):
    f=tmp_path/'bad.npz';np.savez_compressed(f,chemical_method=np.array('bad'))
    with pytest.raises(ValueError,match='mandatory'):MovingSimulation.restore(f)


@pytest.mark.parametrize('bad',['nan','wrong_dtype','wrong_shape','too_large'])
def test_invalid_correction_is_rejected(bad):
    u=np.ones((2,3));e=np.zeros_like(u)
    if bad=='nan':e[0,0]=np.nan
    if bad=='wrong_dtype':e=e.astype(np.float32)
    if bad=='wrong_shape':e=e[:,:2]
    if bad=='too_large':e[0,0]=1e-3
    with pytest.raises(ValueError):validate_correction(u,e)
