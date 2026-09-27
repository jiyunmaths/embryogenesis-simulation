import numpy as np
import pytest

from embryo.model import Config
from embryo.transport import conservative_transport, transport_graph, integrate_gm
from embryo.causal_signaling import ARMS, reaction, jacobian, initial_signals, simulate, preflight


def graph():
    return transport_graph(conservative_transport(np.array([[0.,2.,.5],[2.,0.,1.],[.5,1.,0.]]),np.array([1.,2.,3.])))


@pytest.mark.parametrize('arm',ARMS)
def test_ablations_preserve_equilibrium_and_have_correct_jacobian(arm):
    a,h=reaction(np.ones(1),np.ones(1),2.,arm)
    assert a[0]==h[0]==0
    eps=1e-6;columns=[]
    for perturb in ([eps,0],[0,eps]):
        plus=np.array(reaction(1+perturb[0],1+perturb[1],2.,arm))
        minus=np.array(reaction(1-perturb[0],1-perturb[1],2.,arm))
        columns.append((plus-minus)/(2*eps))
    np.testing.assert_allclose(np.stack(columns,axis=1),jacobian(2.,arm),rtol=1e-8,atol=1e-8)


def test_perturbations_are_paired_and_volume_centered():
    volumes=np.array([1.,2.,3.]);a,h=initial_signals(volumes,[0,1,2],.001)
    weights=volumes/volumes.sum()
    np.testing.assert_allclose(a@weights,1,atol=1e-15)
    np.testing.assert_allclose(np.sqrt(((a-1)**2)@weights),.001,rtol=1e-12)
    np.testing.assert_array_equal(a,initial_signals(volumes,[0,1,2],.001)[0])


def test_full_signaling_matches_existing_positive_integrator():
    g=graph();config=Config()
    result=simulate(g,config,[7],'full',until=.02,interval=.01)
    a,h=initial_signals(g.masses,[7],.001)
    transport=conservative_transport(g.weights,g.masses)
    for _ in range(2):a[0],h[0]=integrate_gm(a[0],h[0],transport,.01)
    np.testing.assert_allclose(result['history'][-1]['activator'],a,atol=1e-14)
    np.testing.assert_allclose(result['history'][-1]['inhibitor'],h,atol=1e-14)


def test_downstream_ablation_changes_fate_not_signals_and_runaway_is_reported():
    g=graph();config=Config()
    full=simulate(g,config,[7],'full',until=.1,interval=.05)
    off=simulate(g,config,[7],'no_signal_to_fate',until=.1,interval=.05)
    np.testing.assert_array_equal(full['history'][-1]['activator'],off['history'][-1]['activator'])
    assert np.max(np.abs(off['history'][-1]['fate']))==0
    assert np.max(np.abs(full['history'][-1]['fate']))>0
    stopped=simulate(g,config,[7],'no_inhibitor_action',until=.1,interval=.05,ceiling=1.)
    assert not stopped['outcomes'][0]['completed']
    assert stopped['outcomes'][0]['stopped'] is not None
    assert not stopped['outcomes'][0]['persistent_contrast_above_0_1']
    assert not preflight(g.eigenvalues,2.,.02,.4,'no_inhibitor_action')['local_stable']
