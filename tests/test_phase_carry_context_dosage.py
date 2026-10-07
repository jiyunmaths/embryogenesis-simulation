"""Test distribution/dosage separation and signed factorial comparisons."""
from dataclasses import asdict
import json
import os

import numpy as np
import pytest

from embryo import phase_carry_context_dosage as study
from embryo import phase_carry_exchange_response as moving
from embryo.attribute_development import AttributeSimulation
from embryo.feedback_long import digest
from embryo.resolution import write_json
from test_phase_carry_exchange_response import source_fixture


def example():
    return (np.array([[.1, 2., .3, 4.], [.5, 3., .8, .2]]),
        np.array([[2., .2, .4, .6], [2., .4, .9, 2.]]), np.array([1., 3., 2., 4.]), [11, 17, 23, 29])


@pytest.mark.parametrize('arm', study.ARMS)
def test_factorial_controls_exact_recipient_and_species_amounts(arm):
    x, ref, m, ids = example(); old = x.copy(); i = 1; mask = np.arange(4) != i
    changed, record = study.surrounding_state(x, ref, m, ids, 17, arm)
    np.testing.assert_array_equal(changed[:, i], x[:, i]); np.testing.assert_array_equal(x, old)
    distribution, dosage = study.FACTORS[arm]
    base = ref if distribution == 'reference' else x
    desired = (ref if dosage == 'reference' else x)[:, mask]@m[mask]
    np.testing.assert_allclose(changed[:, mask]@m[mask], desired, rtol=2e-14)
    # Every within-species concentration ratio preserves its assigned distribution.
    np.testing.assert_allclose(changed[:, mask]/changed[:, mask][:, :1],
        base[:, mask]/base[:, mask][:, :1], rtol=3e-15)
    np.testing.assert_allclose(record['added_amounts'], (changed-x)@m, atol=1e-14)
    if arm == 'sham':
        np.testing.assert_array_equal(changed, x)
    if arm == 'reset':
        np.testing.assert_array_equal(changed[:, mask], ref[:, mask])
    if arm == 'redistributed':
        np.testing.assert_allclose(changed@m, x@m, rtol=2e-14)


@pytest.mark.parametrize('bad', ['recipient', 'ids', 'masses', 'positive', 'finite', 'arm'])
def test_invalid_factorial_preparations_are_rejected(bad):
    x, ref, m, ids = example(); target=17; arm='redistributed'
    if bad == 'recipient': target=99
    if bad == 'ids': ids=[11, 17, 17, 29]
    if bad == 'masses': m[0]=0
    if bad == 'positive': ref[0, 0]=-1
    if bad == 'finite': x[0, 0]=np.nan
    if bad == 'arm': arm='unknown'
    with pytest.raises(ValueError):
        study.surrounding_state(x, ref, m, ids, target, arm)


def test_species_scaling_is_separate_and_single_surrounding_compartment_is_valid():
    x, ref, m, ids = example()
    changed, record = study.surrounding_state(x, ref, m, ids, 17, 'redistributed')
    assert record['scale_factors'][0] != record['scale_factors'][1]
    # With only one other cell, amount preservation leaves no spatial redistribution.
    changed, record = study.surrounding_state(x[:, :2], ref[:, :2], m[:2], ids[:2], 17, 'redistributed')
    np.testing.assert_allclose(changed, x[:, :2], rtol=2e-15)


def test_unrepresentable_amount_is_rejected():
    x=np.full((2, 2), 1e308)
    with np.errstate(over='ignore'):
        with pytest.raises(ValueError, match='representable'):
            study.surrounding_state(x, x, [2., 2.], [1, 2], 1, 'bulk')


def test_design_has_24_new_and_18_reusable_paths_nested_in_three_histories(tmp_path):
    choices={str(seed):dict(selected_ids=ids) for seed, ids in zip(study.HISTORIES, ([22,23], [19,20], [29,30]))}
    jobs=study.job_design(tmp_path, choices)
    assert len(jobs)==42 and len({j['key'] for j in jobs})==42
    assert sum(j['arm'] in ('bulk', 'redistributed') for j in jobs)==24
    assert sum(j['arm']=='sham' for j in jobs)==6
    assert all(j['target'] is None and j['factor']==1. for j in jobs)


def factorial_histories(spatial=.02, dosage=.04, interaction=0.):
    times=np.arange(401)*.15; targets=[11,17]; histories={}
    for arm in study.ARMS:
        s,d=study.FACTORS[arm]
        delta=(spatial if s=='reference' else 0)+(dosage if d=='reference' else 0)
        if arm=='reset': delta+=interaction
        for target in ([None] if arm=='sham' else targets):
            path=np.ones((len(times),2,2))*np.exp(delta)
            histories[arm,target]=[dict(ids=targets, chemistry=x.tolist()) for x in path]
    p=dict(interval=.15,late_window=24.,effect_min=.01,state_metric_error_max=.01)
    return histories,targets,p


def test_signed_factorial_interaction_is_zero_for_additive_effects():
    histories,targets,p=factorial_histories()
    row=study.state_summary(histories,targets,p,60.)[0]['contrasts']
    assert row['spatial_original_dose']['late_log_rms']==pytest.approx(.02)
    assert row['dosage_original_distribution']['late_log_rms']==pytest.approx(.04)
    assert row['total_reset']['late_log_rms']==pytest.approx(.06)
    assert row['interaction']['late_log_rms']<1e-15 and not row['interaction']['effect_detected']
    assert study.state_comparison(histories,histories,targets,p,60.)['passed']


def test_interaction_and_effect_decision_disagreement_fail_refinement():
    a,targets,p=factorial_histories(interaction=.03)
    b,_,_=factorial_histories(spatial=.001, interaction=.03)
    assert study.state_summary(a,targets,p,60.)[0]['contrasts']['interaction']['late_log_rms']==pytest.approx(.03)
    report=study.state_comparison(a,b,targets,p,60.)
    assert not report['passed'] and not report['same_effect_decisions']


def test_late_screen_excludes_early_transient_but_full_window_retains_it():
    histories,targets,p=factorial_histories(spatial=0,dosage=0)
    for row in histories['redistributed',11][:40]:
        row['chemistry'][0][0]=np.exp(.3)
    row=study.state_summary(histories,targets,p,60.)[0]['contrasts']['spatial_original_dose']
    assert not row['effect_detected'] and row['late_log_rms']==0 and row['full_log_rms']>.01


def test_retiming_keeps_physical_arrays_carry_and_random_streams(tmp_path):
    host,source=source_fixture(tmp_path,start=450.)
    x=np.array([host.activator,host.inhibitor]); ref=x[:,::-1].copy()
    changed,_=study.surrounding_state(x,ref,host.volumes(),host.ids,int(host.ids[0]),'bulk')
    before=moving.payload(source); dest=tmp_path/'new.npz'
    new=moving.prepared_checkpoint(source,dest,host.config.dt/2,changed); after=moving.payload(dest)
    for key in before:
        if key not in ('metadata','activator','inhibitor','long_experiment'):
            assert before[key].dtype==after[key].dtype
            np.testing.assert_array_equal(before[key],after[key])
    a,b=[json.loads(str(z['metadata'])) for z in (before,after)]
    for key in ('rng','fate_rng','signal_rng','time','lineage','divisions'):
        assert a[key]==b[key]
    assert new.step_number==host.step_number*2


@pytest.mark.parametrize('volume_error', [.01,.06])
def test_bad_physical_pilot_blocks_long_continuation(tmp_path,volume_error):
    host,source=source_fixture(tmp_path,start=450.); seed=host.config.seed
    job=dict(key='pilot',seed=seed,start=450.)
    p=dict(jobs=[job],reused={},pilot_duration=6.,amount_error_max=2e-14)
    write_json(tmp_path/'protocol.json',p)
    state=moving.payload(source); meta=json.loads(str(state['metadata']))
    meta['time']=456.; meta['step_number']=round(456./host.config.dt)
    state['metadata']=np.array(json.dumps(meta))
    audit=dict(max_volume_error=volume_error,min_radius=5.,max_clipping=0.,boundary_max=0.,dilution_error_max=0.)
    state['long_experiment']=np.array(json.dumps(dict(protocol_hash=digest(tmp_path/'protocol.json'),
        job=job,audit=audit,history=[dict(elapsed=6.)])))
    folder=tmp_path/'pilot';folder.mkdir();np.savez_compressed(folder/'latest_state.npz',**state)
    if volume_error<.05:
        assert study.pilot_quality(tmp_path,p,seed)['passed']
    else:
        with pytest.raises(RuntimeError,match='long continuation blocked'):
            study.pilot_quality(tmp_path,p,seed)
        record=json.loads((tmp_path/'pairs'/f'seed-{seed}'/'pilot-quality.json').read_text())
        assert not record['passed']


@pytest.mark.skipif(os.environ.get('EMBRYO_CUDA_TESTS')!='1',reason='GPU opt-in required')
@pytest.mark.parametrize('arm', ['redistributed','bulk'])
@pytest.mark.parametrize('dt', [.00375,.001875])
def test_new_surrounding_arm_passes_native_carry_chemical_and_restart_gate(tmp_path,arm,dt):
    host,source=source_fixture(tmp_path,start=450.)
    config=type(host.config)(**{**asdict(host.config),'grid':24,'extent':1.2,'max_cells':3})
    expanded=AttributeSimulation(config,mode='direct')
    split=.5*(1+np.tanh(expanded.xyz[0]/config.interface_width))
    second=.5*(1+np.tanh(expanded.xyz[1]/config.interface_width))
    expanded.phi=np.array([expanded.phi[0]*split,expanded.phi[0]*(1-split)*second,
        expanded.phi[0]*(1-split)*(1-second)],dtype=np.float32)
    expanded.target=expanded.volumes()
    for key in ('ids','parents','due','fate','activator','inhibitor','polarity'):
        a=getattr(host,key)
        setattr(expanded,key,np.concatenate([a,a[-1:]],axis=0))
    expanded.ids[-1]=expanded.ids[-2]+1
    expanded.time,expanded.step_number=host.time,host.step_number;expanded.checkpoint(source)
    with np.load(source) as z: out={k:z[k].copy() for k in z.files}
    out.update(precision_arm=np.array('phase_carry'),phase_carry=np.full(expanded.phi.shape,1e-8),rounding=np.array([5,2],np.int32))
    np.savez_compressed(source,**out)
    x=np.array([expanded.activator,expanded.inhibitor]); ref=x.copy(); ref[:,1]*=[1.4,.7];ref[:,2]*=[.8,1.4]
    changed,_=study.surrounding_state(x,ref,expanded.volumes(),expanded.ids,int(expanded.ids[0]),arm)
    assert not np.array_equal(changed,x)
    dest=tmp_path/'prepared.npz';new=moving.prepared_checkpoint(source,dest,dt,changed)
    job=dict(key='factorial-arm',seed=host.config.seed,level='fine',dt=dt,source=str(dest))
    p=dict(accepted_configs={f'{host.config.seed}_fine':asdict(new.config)},amount_error_max=2e-14,
        prefix_native_threads=2,prefix_duration=.015,interval=.00375)
    record=moving.context_gate(tmp_path,job,p)
    assert record['passed'] and record['restart_exact_steps']==4 and record['cpu_chemical_log_max']<=1e-11
