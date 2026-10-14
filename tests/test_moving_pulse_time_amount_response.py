"""Preparation/replay safeguards for immediate moving matched-start responses."""
import numpy as np
import pytest
from embryo import moving_pulse_time_amount_response as study
from test_moving_recipient_amount_response import response_fixture


def state():
    return dict(initial=np.array([[1.,2.,3.],[2.,3.,4.]]),volumes=np.array([1.,2.,3.]),ids=np.array([19,20,21]),delta=np.array([[-1.,1.,0.],[.5,-1.,.5],[0.,1/3,-1/3]]))


def test_alignment_rejects_endpoint_transport_or_capacity_changes():
    a=state();row=dict(time=570.,elapsed=60.,**{k:v.tolist() for k,v in a.items() if k!='initial'})
    study.validate_alignment(a,a['ids'],row)
    changed={**row,'volumes':[1.,2.,3.001]}
    with pytest.raises(ValueError,match='capacity/operator'):study.validate_alignment(a,a['ids'],changed)
    with pytest.raises(ValueError,match='clock'):study.validate_alignment(a,a['ids'],dict(row,time=570.1))


def test_matched_amounts_allow_allocation_difference_but_not_total_drift():
    a=state();b={k:v.copy() for k,v in a.items()};b['initial'][:,1]+=.15;b['initial'][:,2]-=.1
    assert max(study.matched_closure(a,b,19).values())<1e-15
    b['initial'][0,1]+=.01
    with pytest.raises(ValueError,match='budgets'):study.matched_closure(a,b,19)


def test_equal_fractional_pulse_preserves_matched_budgets_between_backgrounds():
    a=state();b={k:v.copy() for k,v in a.items()};b['initial'][:,1]+=.15;b['initial'][:,2]-=.1
    for x in (a,b):x['initial'][0,0]*=.9
    assert max(study.matched_closure(a,b,19).values())<1e-15
    assert a['initial'][1,0]==b['initial'][1,0]==2.


def test_zero_correction_only_at_explicit_first_intervention(tmp_path,monkeypatch):
    a=state();source=tmp_path/'initial.npz';study.save_npz(source,**a,correction=np.full_like(a['initial'],1e-17))
    p=dict(contexts=[dict(key='a',physical_source='unused')]);job=dict(context='a',source=str(source))
    def forbidden(*args,**kwargs):raise AssertionError('GPU should not be allocated for invalid reset')
    monkeypatch.setattr(study.backend.MovingSimulation,'from_physical',forbidden)
    with pytest.raises(ValueError,match='zero initial chemical correction'):study.initialize(tmp_path,job,p)


def test_moving_response_uses_own_controls_and_signed_species_vectors(tmp_path):
    p,pair=response_fixture(tmp_path);p['jobs']=[j for j in p['jobs'] if j['dose']!='amount']
    report=study.response_contrast(tmp_path,p,pair,60.)
    assert report['passed'] and len(report['rows'])==1 and report['rows'][0]['fine_above_screen']
    assert report['rows'][0]['direct_signed_waveform_error']==0.
    with np.load(tmp_path/'contrasts/history-8_cell-19_fractional_full.npz') as z:
        # R responds more slowly, so R-S must be negative after the pulse.
        assert z['fine'][1,0]<0 and np.all(z['fine'][:,1]==0)
