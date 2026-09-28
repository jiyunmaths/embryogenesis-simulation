import numpy as np
from embryo.attribute_common_environment import bath_rhs, equilibria


def test_isolation_has_unit_equilibrium_and_no_cross_cell_coupling():
    f=bath_rhs(2,0,0,[1,1])
    assert equilibria(2,0,0,[1,1])[0]['locally_stable']
    np.testing.assert_array_equal(f(0,np.ones(4)),0)
    y=np.array([.5,2.,1.,3.]);changed=y.copy();changed[0]=.8
    np.testing.assert_array_equal(f(0,y)[[1,3]],f(0,changed)[[1,3]])


def test_analytic_reservoir_equilibria_satisfy_odes_and_include_bistability():
    f=bath_rhs(2,.4,8,[1,1])
    roots=equilibria(2,.4,8,[1,1])
    assert len(roots)==3
    assert sum(x['locally_stable'] for x in roots)==2
    for x in roots:
        np.testing.assert_allclose(f(0,np.array([x['activator'],x['inhibitor']])),0,atol=1e-12)
