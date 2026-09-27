import numpy as np
from embryo.fate_robustness import pulse, noise_paths, wilson
from embryo.fate_memory import neutral_fate


def test_pulses_are_sign_symmetric_and_subcritical_equilibria_persist():
    initial=np.array([-1.,1.])
    sub=pulse(initial,.25,24.,.0075)
    np.testing.assert_allclose(sub[0],-sub[1],atol=1e-14)
    assert sub[1]>1/np.sqrt(3)
    over=pulse(initial,.75,12.,.0075)
    recovered=neutral_fate(over,40.,.8)
    np.testing.assert_allclose(recovered,-initial,atol=1e-10)


def test_shared_brownian_increments_and_zero_noise():
    signs,coarse,fine,_,_=noise_paths([0.,.2],per_sign=4,until=.06,dt=.0075,interval=.015,rate=0.)
    np.testing.assert_allclose(coarse,fine,atol=1e-15)
    np.testing.assert_array_equal(coarse[:,0],np.broadcast_to(signs,coarse[:,0].shape))
    low,high=wilson(0,512)
    assert abs(low)<1e-15 and 0<high<.01
