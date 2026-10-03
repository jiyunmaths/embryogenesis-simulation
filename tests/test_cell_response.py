import numpy as np
import pytest
from embryo.cell_response import pulse, sustained_recovery, response_metrics


def test_pulse_changes_only_target_and_external_amount():
    initial=np.array([[1,2,3],[4,5,6]])
    masses=np.array([.1,.3,.2]);result=pulse(initial,1,.5)
    assert np.array_equal(initial,[[1,2,3],[4,5,6]])
    assert np.array_equal(result,[[1,1,3],[4,5,6]])
    assert (result[0]-initial[0])@masses==pytest.approx(-.3)
    for factor in (0,-1,np.nan):
        with pytest.raises(ValueError):pulse(initial,1,factor)


def test_recovery_requires_no_late_rebound_and_observation_window():
    t=np.arange(101.)
    y=np.exp(-t/5);y[60]=.5
    assert sustained_recovery(t,y,.1)==61
    y[90]=.5
    assert sustained_recovery(t,y,.1) is None


def test_response_metrics_against_exponential_recovery():
    t=np.arange(0,101.,.1);control=np.ones((len(t),2,2));path=control.copy()
    path[:,0,0]=np.exp(np.log(1.1)*np.exp(-t/2))
    r=response_metrics(path,control,t,np.ones(2),0,1.1)
    assert r['target_activator_peak_gain']==pytest.approx(1)
    assert r['target_inhibitor_peak_gain']==0
    assert r['target_activator_log_auc_per_log_pulse']==pytest.approx(2,rel=.001)
    assert r['target_recovery_time']==pytest.approx(4.7)
    assert r['other_cells_peak_log_rms_per_log_pulse']==0
