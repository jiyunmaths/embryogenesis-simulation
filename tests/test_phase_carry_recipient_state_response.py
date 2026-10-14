"""Qualify local matching, independent backgrounds and conditional responses."""
from pathlib import Path
import numpy as np
import pytest
import torch
from embryo import phase_carry_recipient_state_response as study
from test_phase_carry_shared_context_response import fixture as parent_fixture


def fixture(tmp_path, duration=.6):
    parent=parent_fixture(tmp_path/'parent-inputs',duration)
    saved={}
    for c in parent['contexts']:
        with np.load(c['source'],allow_pickle=False) as z:saved[c['key']]={k:z[k].copy() for k in z.files}
    contexts,trials,factorials,arrays=study.design(parent,saved,tmp_path)
    for c in contexts:study.save_npz(c['source'],**arrays[c['key']])
    for t in trials:
        z=arrays[t['context']];u=z['initial'].copy();factor=1.
        if t['target'] is not None:
            i=z['ids'].tolist().index(t['target'])
            u,record=study.shared.shared_pulse(u,z['volumes'],z['ids'],t['target'],'fractional',[u[0,i],u[0,i]])
            t.update(record);factor=record['factor']
        t.update(source=str(tmp_path/'starts'/(t['key']+'.npz')),factor=factor,
            initial_log_amplitude=float(abs(np.log(factor))))
        study.save_npz(t['source'],initial=u)
    p={**parent,'parent':parent['parent'],'contexts':contexts,'trials':trials,'factorials':factorials,
       'primary_paths':186,'diagnostic_paths':93,'total_production_paths':279,
       'cpu_verification_paths':279,'logical_reference_cases':93,'reference_trajectories':186,
       'source_sha256':{},'input_sha256':study.frozen.pin([c['source'] for c in contexts]+[t['source'] for t in trials])}
    p['unique_initial_value_problems']=study.frozen.unique_problem_count(p)
    study.write_json(tmp_path/'protocol.json',p)
    return p,parent,saved,arrays


def test_realistic_factorial_counts_and_exact_source_isolation(tmp_path):
    p,parent,saved,arrays=fixture(tmp_path)
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(45,93,12)
    assert sum(t['target'] is None for t in p['trials'])==45
    assert sum(t['target'] is not None for t in p['trials'])==48
    assert p['unique_initial_value_problems']==93
    for f in p['factorials']:
        origin=next(x for x in parent['factorials'] if (x['seed'],x['cell'])==(f['seed'],f['cell']))
        ref=saved[origin['contexts'][f['matching_reference']+'_'+f['matching_reference']]]
        index=ref['ids'].tolist().index(f['cell']);pair=ref['initial'][:,index]
        for condition,key in f['contexts'].items():
            source=saved[origin['contexts'][condition]];result=arrays[key]
            np.testing.assert_array_equal(result['initial'][:,index],pair)
            others=np.arange(len(ref['ids']))!=index
            np.testing.assert_array_equal(result['initial'][:,others],source['initial'][:,others])
            for name in ('delta','volumes','ids'):np.testing.assert_array_equal(result[name],source[name])
    for seed in study.HISTORIES:
        families=[f for f in p['factorials'] if f['seed']==seed and f['matching_reference']=='S']
        assert families[0]['contexts']['S_S']==families[1]['contexts']['S_S']


def test_matching_changes_only_two_target_species_and_accounts_amounts():
    data=dict(initial=np.array([[.4,2.],[.8,3.]]),volumes=np.array([1.,2.]),
              delta=np.array([[-.2,.2],[.1,-.1]]),ids=np.array([11,17]))
    result,record=study.match_recipient(data,17,np.array([.3,1.1]))
    np.testing.assert_array_equal(result['initial'][:,0],data['initial'][:,0])
    np.testing.assert_array_equal(record['matching_per_cell_species_amount_change'],[[0.,-3.4],[0.,-3.8]])
    np.testing.assert_array_equal(record['post_matching_global_species_amounts'],[1.,3.])
    np.testing.assert_array_equal(data['initial'],[[.4,2.],[.8,3.]])


@pytest.mark.parametrize('pair',[np.array([0.,1.]),np.array([-1.,1.]),np.array([np.nan,1.]),np.array([1.]),np.array([1.,1.],dtype=np.float32)])
def test_invalid_measured_matching_pair_rejected(pair):
    data=dict(initial=np.array([[.4,2.],[.8,3.]]),volumes=np.array([1.,2.]),
              delta=np.array([[-.2,.2],[.1,-.1]]),ids=np.array([11,17]))
    with pytest.raises(ValueError):study.match_recipient(data,17,pair)


def test_equal_local_state_means_equal_pulse_amount_within_physics(tmp_path):
    p,_,_,arrays=fixture(tmp_path)
    for f in p['factorials']:
        for P in ('S','R'):
            pair=[]
            for C in ('S','R'):
                t=next(t for t in p['trials'] if t['context']==f['contexts'][C+'_'+P] and t['target']==f['cell'])
                z=arrays[t['context']];i=z['ids'].tolist().index(f['cell'])
                pair.append(t)
                assert t['factor']==.9
                expected=z['initial'].copy();expected[0,i]*=.9
                with np.load(t['source']) as q:np.testing.assert_array_equal(q['initial'],expected)
                equivalent,_=study.shared.shared_pulse(z['initial'],z['volumes'],z['ids'],f['cell'],'amount',[z['initial'][0,i]]*2)
                np.testing.assert_array_equal(equivalent,expected)
            assert pair[0]['requested_removed_amount']==pair[1]['requested_removed_amount']
            assert pair[0]['removed_activator_amount']==pair[1]['removed_activator_amount']


def test_block_layout_is_uncoupled_and_preserves_actual_cases(tmp_path):
    p,_,_,arrays=fixture(tmp_path)
    initial,matrix,mapping,count=study.frozen.batch_layout(p)
    assert initial.shape==(2,90,3) and count==2
    assert [m['trial'] for m in mapping]==[t['key'] for t in p['trials']]
    for index,c in enumerate(p['contexts']):
        np.testing.assert_array_equal(matrix[2*index:2*index+2,2*index:2*index+2],arrays[c['key']]['delta'])
        mask=np.ones(matrix.shape[1],bool);mask[2*index:2*index+2]=False
        assert not np.any(matrix[2*index:2*index+2,mask])
    for level in p['levels']:
        global_count=study.frozen.substeps(matrix,level['dt'],p['beta'],p['da'],p['db'])
        assert all(study.frozen.substeps(x['delta'],level['dt'],p['beta'],p['da'],p['db'])==global_count for x in arrays.values())


def synthetic(p,times,label):
    paths={}
    for ci,c in enumerate(p['contexts']):
        ctrl=next(t for t in p['trials'] if t['context']==c['key'] and t['target'] is None)
        _,_,_,u,ids=study.frozen.load_context(p,ctrl)
        baseline=u[None]*np.exp(.0002*(ci+1)*times[:,None,None])
        paths[label,ctrl['key']]=baseline
        for t in [t for t in p['trials'] if t['context']==c['key'] and t['target'] is not None]:
            i=ids.tolist().index(t['target']);v=baseline.copy()
            wave=np.column_stack((-np.exp(-times/2),.1*times*np.exp(-times/2)))
            wave+=(.004*ci)*times[:,None]*np.exp(-times[:,None]/4)
            v[:,:,i]*=np.exp(wave*t['initial_log_amplitude']);paths[label,t['key']]=v
    return paths


def test_reference_labels_and_signed_interaction_are_not_collapsed(tmp_path):
    p,_,_,_=fixture(tmp_path,60.)
    times=np.arange(401)*.15;paths=synthetic(p,times,'fine')
    summary,waves=study.summarize(paths,'fine',p,times)
    assert len(summary['metrics'])==48 and len(summary['effects'])==60 and len(waves)==60
    assert len(summary['native_diagonal'])==12
    assert sum('response_shift_detected' in row for row in summary['effects'])==24
    assert {row['matching_reference'] for row in summary['effects']}=={'S','R'}
    for f in p['factorials']:
        prefix=f'{f["seed"]}_{f["cell"]}_{f["matching_reference"]}_fractional_'
        np.testing.assert_allclose(waves[prefix+'interaction_difference_of_differences'],
            waves[prefix+'chemical_redistributed_physics']-waves[prefix+'chemical_sham_physics'],atol=1e-15)
    desc=study.baseline_descriptions(paths,p,times)
    assert len(desc['paired_chemical_states'])==24
    assert all(row['initial_species_ratios']==[1.,1.] for row in desc['paired_chemical_states'])


def test_direct_four_way_reference_gate_is_retained(tmp_path):
    p,_,_,_=fixture(tmp_path,60.);times=np.arange(401)*.15
    paths=synthetic(p,times,'dop853');paths.update({('radau',k):v.copy() for (l,k),v in list(paths.items())})
    f=p['factorials'][0]
    for name,sign in [('R_R',1.),('S_R',-1.),('R_S',-1.),('S_S',1.)]:
        t=next(t for t in p['trials'] if t['context']==f['contexts'][name] and t['target']==f['cell'])
        _,_,_,_,ids=study.frozen.load_context(p,t);i=ids.tolist().index(f['cell'])
        paths['radau',t['key']][:,:,i]*=np.exp(sign*6e-9*(1-np.exp(-times[:,None]))*t['initial_log_amplitude'])
    unpacked={l:study.summarize(paths,l,p,times) for l in ('dop853','radau')}
    result=study.compare(paths,'dop853','radau',p,times,{l:r[0] for l,r in unpacked.items()},{l:r[1] for l,r in unpacked.items()},True)
    assert all(r['passed'] for r in result['pulse_rows'])
    assert not result['passed']
    assert any(r['kind']=='interaction' and r['direct_waveform_error']>1e-8 for r in result['factorial_contrasts'])


def test_modified_prepared_arrays_rejected_before_stepping(tmp_path):
    p,_,_,_=fixture(tmp_path);f=Path(p['contexts'][0]['source'])
    with np.load(f) as z:q={k:z[k].copy() for k in z.files}
    q['initial'][0,0]*=1.01;study.save_npz(f,**q)
    with pytest.raises(ValueError,match='preparation'):study.frozen.load_context(p,p['trials'][0])


def test_full_concurrent_cpu_fixture_qualifies_all_cases(tmp_path):
    p,_,_,_=fixture(tmp_path)
    result=study.run(tmp_path,workers=2)
    assert result['passed']
    assert len(result['production_accuracy'])==279 and len(result['orders'])==93
    assert len(result['raw_timestep_comparisons'])==186
    assert len(result['response_comparisons'])==9
    assert all(len(r['pulse_rows'])==48 and len(r['factorial_contrasts'])==60 for r in result['response_comparisons'])
    assert len(result['history_outcomes'])==12 and len(result['reciprocal_history_outcomes'])==3
    assert study.read(tmp_path/'timing.json')['cpu_verification_overlaps_gpu']
