"""Protect reset specificity, exact carry reuse and own-background responses."""
from dataclasses import asdict
import json
import os

import numpy as np
import pytest

from embryo import phase_carry_network_context as study
from embryo import phase_carry_exchange_response as moving
from embryo.attribute_development import AttributeSimulation
from embryo.resolution import write_json
from test_phase_carry_exchange_response import source_fixture


def test_context_reset_preserves_recipient_and_records_actual_external_amounts():
    x=np.array([[.1,2.,.3],[.5,3.,.8]])
    ref=np.array([[2.,.2,.4],[2.,.4,.9]])
    masses=np.array([1.,3.,2.]);ids=np.array([11,17,23])
    reset,record=study.reset_context(x,ref,masses,ids,17)
    np.testing.assert_array_equal(reset[:,1],x[:,1])
    np.testing.assert_array_equal(reset[:,[0,2]],ref[:,[0,2]])
    np.testing.assert_array_equal(x,[[.1,2.,.3],[.5,3.,.8]])
    np.testing.assert_allclose(record['added_amounts'],(reset-x)@masses,rtol=0,atol=1e-15)
    np.testing.assert_allclose(record['reset_total_amounts'],np.asarray(record['initial_total_amounts'])+record['added_amounts'])
    assert record['reset_ids']==[11,23] and record['recipient_exact']
    assert record['surrounding_log_rms']>.05 and not record['sham_exact']


def test_sham_is_exact_and_duplicate_ids_or_invalid_reference_are_rejected():
    x=np.array([[.1,2.],[.5,3.]]);m=np.array([1.,3.]);ids=[11,17]
    sham,record=study.reset_context(x,x*2,m,ids)
    np.testing.assert_array_equal(sham,x)
    assert record['sham_exact'] and record['surrounding_log_rms']==0.
    assert record['added_amounts']==[0.,0.] and record['reset_ids']==[]
    with pytest.raises(ValueError):study.reset_context(x,x,m,[11,11],11)
    with pytest.raises(ValueError):study.reset_context(x,x,m,ids,12)
    with pytest.raises(ValueError):study.reset_context(x,-x,m,ids,11)


def test_nonchemical_arrays_carry_and_rngs_survive_context_reset_and_retiming(tmp_path):
    host,source=source_fixture(tmp_path,start=450.)
    original=moving.payload(source);x=np.array([host.activator,host.inhibitor]);ref=x[:,::-1].copy()
    values,record=study.reset_context(x,ref,host.volumes(),host.ids,int(host.ids[0]))
    edited=tmp_path/'edited.npz';new=moving.prepared_checkpoint(source,edited,host.config.dt/2,values)
    after=moving.payload(edited)
    for key in original:
        if key not in ('metadata','activator','inhibitor','long_experiment'):
            np.testing.assert_array_equal(after[key],original[key])
    old_meta,new_meta=[json.loads(str(z['metadata'])) for z in (original,after)]
    for key in ('rng','fate_rng','signal_rng','time','lineage','divisions'):
        assert old_meta[key]==new_meta[key]
    assert new.step_number==host.step_number*2 and record['recipient_exact']


def test_reuse_requires_exact_full_initial_state_not_approximate_agreement(tmp_path):
    host,source=source_fixture(tmp_path,start=450.)
    dest=tmp_path/'copy.npz';moving.prepared_checkpoint(source,dest,host.config.dt)
    study.source_equivalence(source,dest)
    z=moving.payload(dest);z['phase_carry'].flat[0]+=1e-16
    np.savez_compressed(dest,**z)
    with pytest.raises(ValueError,match='phase_carry'):study.source_equivalence(source,dest)


def histories():
    times=np.arange(401)*.15;ids=[11,17,23];m=[1.,2.,3.]
    init=np.array([[.1,2.,.4],[.5,3.,.9]])
    result={};refs={}
    for name in ('sham','reset-cell-11','reset-cell-17','reference'):
        initial=init.copy()
        if name.startswith('reset'):initial[:,2]*=2
        control=np.repeat(initial[None],len(times),axis=0)*np.exp(.002*times[:,None,None])
        for target in (None,11,17):
            if name.startswith('reset') and target not in (None,int(name.split('-')[-1])):continue
            path=control.copy()
            if target is not None:
                i=ids.index(target)
                # Reference has 1/4 time constants; sham represents a transferred 4/1 response.
                tau=(1. if target==11 else 4.) if name=='reference' else (4. if target==11 else 1.)
                if name.startswith('reset'):tau*=2
                path[:,0,i]*=np.exp(np.log(.9)*np.exp(-times/tau))
            h=[dict(elapsed=float(t),time=450.+float(t),ids=ids,chemistry=x.tolist(),
                volumes=m,log_activator_sd=float(np.std(np.log(x[0])))) for t,x in zip(times,path)]
            if name=='reference':refs[target]=h
            else:result[name,target]=h
    p=dict(interval=.15,late_window=24.,contrast_min=.1,effect_min=.01,response_limits=moving.RESPONSE_CRITERIA)
    return result,refs,p


def test_context_effect_uses_own_drifting_control_and_changes_response_not_self_initial_state():
    h,ref,p=histories();rows=study.context_summary(h,ref,[11,17],p,60.)
    assert all(r['response_shift_detected'] for r in rows)
    assert all(not r['late_target_state_shift_detected'] for r in rows)
    assert all(r['contexts']['sham']['nearest_reference']=='donor' for r in rows)
    assert study.response_comparison(h,h,ref,[11,17],p,60.)['passed']


def test_uninformative_references_stay_unresolved():
    h,ref,p=histories()
    control=np.asarray([r['chemistry'] for r in ref[None]])
    for target in (11,17):
        i=ref[None][0]['ids'].index(target)
        t=np.arange(401)*.15;x=control.copy();x[:,0,i]*=np.exp(np.log(.9)*np.exp(-t/2))
        for row,chem in zip(ref[target],x):row['chemistry']=chem.tolist()
    rows=study.context_summary(h,ref,[11,17],p,60.)
    assert all(not c['informative_reference'] and c['nearest_reference']=='unresolved' for r in rows for c in r['contexts'].values())


def test_refinement_rejects_normalized_response_error_even_with_small_raw_difference():
    h,ref,p=histories();fine=json.loads(json.dumps({str(k):v for k,v in h.items()}))
    fine={k:fine[str(k)] for k in h};fine['reset-cell-11',11][5]['chemistry'][0][0]*=np.exp(.005)
    report=study.response_comparison(h,fine,ref,[11,17],p,60.)
    assert not report['passed'] and max(r['normalized_response_error'] for r in report['rows'])>.01


def test_shared_sham_control_and_context_specific_reset_controls(tmp_path):
    selections={str(seed):dict(selected_ids=[11,17]) for seed in study.HISTORIES}
    contexts,jobs=study.job_design(tmp_path,selections)
    assert len(contexts)==9 and len(jobs)==42
    assert len({j['key'] for j in jobs})==42
    assert sum(j['recipient'] is None and j['target'] is None for j in jobs)==6
    assert sum(j['recipient'] is not None for j in jobs)==24
    assert all(j['target'] in (None,j['recipient']) for j in jobs if j['recipient'] is not None)
    assert sum(j['recipient'] is None and j['level']=='fine' for j in jobs)==9


def test_failed_pilot_cannot_launch_long_runs(tmp_path,monkeypatch):
    selections={str(seed):dict(selected_ids=[11,17]) for seed in study.HISTORIES}
    _,jobs=study.job_design(tmp_path,selections);p=dict(jobs=jobs,pilot_duration=6.,duration=60.,total_jobs=42,new_jobs=33,reused_jobs=9)
    calls=[];monkeypatch.setattr(study,'advance',lambda root,job,p,h:calls.append(h))
    monkeypatch.setattr(study,'pair_report',lambda *args:dict(passed=False))
    done,passed=study.schedule_history(tmp_path,p,7,0)
    assert not passed and done==0 and calls==[6.]*14


def test_reused_pilot_uses_its_pilot_field_not_the_completed_field(tmp_path):
    pilot,full=tmp_path/'pilot.npz',tmp_path/'full.npz'
    np.savez_compressed(pilot,phi=np.zeros(3,dtype=np.float32))
    np.savez_compressed(full,phi=np.ones(3,dtype=np.float32))
    p=dict(pilot_duration=6.,reused=dict(reference=dict(pilot_phi=str(pilot),full_phi=str(full))))
    np.testing.assert_array_equal(study.field_snapshot(tmp_path,dict(key='reference'),p,6.),np.zeros(3))
    np.testing.assert_array_equal(study.field_snapshot(tmp_path,dict(key='reference'),p,60.),np.ones(3))


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU opt-in required')
@pytest.mark.parametrize('dt',[.00375,.001875])
def test_reset_pulse_context_and_nonzero_carry_on_gpu(tmp_path,dt):
    host,source=source_fixture(tmp_path,start=450.)
    config=type(host.config)(**{**asdict(host.config),'grid':24,'extent':1.2})
    expanded=AttributeSimulation(config,mode='direct')
    partition=.5*(1+np.tanh(expanded.xyz[0]/config.interface_width))
    expanded.phi=np.array([expanded.phi[0]*partition,expanded.phi[0]*(1-partition)],dtype=np.float32)
    expanded.target=expanded.volumes()
    for key in ('ids','parents','due','fate','activator','inhibitor','polarity'):
        setattr(expanded,key,getattr(host,key).copy())
    expanded.time=host.time;expanded.step_number=host.step_number;expanded.checkpoint(source)
    with np.load(source) as z:out={k:z[k].copy() for k in z.files}
    out.update(precision_arm=np.array('phase_carry'),phase_carry=np.full(expanded.phi.shape,1e-8),rounding=np.array([5,2],np.int32))
    np.savez_compressed(source,**out)
    x=np.array([expanded.activator,expanded.inhibitor])
    changed,_=study.reset_context(x,x[:,::-1],expanded.volumes(),expanded.ids,int(expanded.ids[0]))
    changed[0,0]*=.9;dest=tmp_path/'reset-pulse.npz'
    new=moving.prepared_checkpoint(source,dest,dt,changed)
    job=dict(key='reset-context',seed=host.config.seed,level='coarse',dt=dt,source=str(dest))
    p=dict(accepted_configs={f'{host.config.seed}_coarse':asdict(new.config)},amount_error_max=2e-14,
        prefix_native_threads=2,prefix_duration=.015,interval=.00375)
    gate=moving.context_gate(tmp_path,job,p)
    assert gate['passed'] and gate['restart_exact_steps']==4 and gate['cpu_chemical_log_max']<=1e-11
