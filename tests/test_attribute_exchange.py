import numpy as np
from embryo.attribute_exchange import exchange, classify


def test_exchange_conserves_pair_amount_and_is_reversible():
    x=np.array([[.1,2.,3.],[2.,.5,4.]])
    m=np.array([1.,2.,5.])
    swapped=exchange(x,m,0,2)
    assert np.all(swapped>0)
    np.testing.assert_allclose(swapped@m,x@m,rtol=1e-14)
    np.testing.assert_array_equal(swapped[:,1],x[:,1])
    np.testing.assert_allclose(exchange(swapped,m,0,2),x,rtol=1e-14)
    np.testing.assert_allclose(exchange(x,np.ones(3),0,2),x[:,[2,1,0]])


def test_classification_requires_contrast_and_stationarity():
    assert classify(0.,1.,True,False)=='uninformative'
    assert classify(.01,1.,False,True)=='unsettled'
    assert classify(.01,1.,True,True)=='destination'
    assert classify(1.,.01,True,True)=='transferred'
    assert classify(.5,.5,True,True)=='reorganized'
