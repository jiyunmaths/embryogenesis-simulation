"""Qualify delayed source, dose controls, grouping and numerical accounting."""
from pathlib import Path
import numpy as np
import pytest

from embryo import phase_carry_recipient_delayed_response as study
from test_phase_carry_recipient_state_response import fixture as matched_fixture, synthetic


def fixture(tmp_path,duration=.3):
    old,_,_,_=matched_fixture(tmp_path/'accepted',duration)
    arrays={}
    for i,c in enumerate(old['contexts']):
        with np.load(c['source']) as z:arrays[c['key']]={k:z[k].copy() for k in z.files}
        # A delayed target may have diverged again: never enforce old matching.
        arrays[c['key']]['initial'] *= np.array([1+.01*(i+1),1+.005*(i+1)])[:,None]
    contexts,trials,factorials,starts=study.design(old,arrays,tmp_path)
    for c in contexts:study.save_npz(c['source'],**arrays[c['parent_context']])
    for t in trials:study.save_npz(t['source'],initial=starts[t['key']])
    times=np.arange(round(duration/.15)+1)*.15
    immediate=study.matched.summarize(synthetic(old,times,'fine'),'fine',old,times)[0]['effects']
    p={**old,'parent':str(tmp_path/'accepted'),'start':570.,'wait':60.,'contexts':contexts,'trials':trials,'factorials':factorials,
        'primary_paths':282,'diagnostic_paths':141,'total_production_paths':423,'cpu_verification_paths':423,
        'logical_reference_cases':141,'reference_trajectories':282,'source_checks':{},
        'accepted_immediate_fractional_effects':immediate,
        'source_sha256':{},'input_sha256':study.frozen.pin([c['source'] for c in contexts]+[t['source'] for t in trials])}
    p['unique_initial_value_problems']=study.frozen.unique_problem_count(p)
    study.write_json(tmp_path/'protocol.json',p)
    return p,old,arrays


def test_delayed_design_preserves_whole_evolved_states_without_rematching(tmp_path):
    p,old,arrays=fixture(tmp_path)
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(45,141,12)
    assert sum(t['target'] is None for t in p['trials'])==45
    assert sum(t['target'] is not None for t in p['trials'])==96
    assert 93<=p['unique_initial_value_problems']<141
    for c in p['contexts']:
        with np.load(c['source']) as z:
            for name in ('initial','delta','volumes','ids'):
                np.testing.assert_array_equal(z[name],arrays[c['parent_context']][name])
    for f in p['factorials']:
        for P in ('S','R'):
            pair=[next(c for c in p['contexts'] if c['key']==f['contexts'][C+'_'+P]) for C in ('S','R')]
            a,b=[arrays[c['parent_context']] for c in pair]
            for name in ('delta','volumes','ids'):np.testing.assert_array_equal(a[name],b[name])
            i=a['ids'].tolist().index(f['cell'])
            assert not np.array_equal(a['initial'][:,i],b['initial'][:,i])


def test_delayed_amount_dose_matches_removal_not_fraction(tmp_path):
    p,_,_=fixture(tmp_path)
    for f in p['factorials']:
        for P in ('S','R'):
            pair={dose:[next(t for t in p['trials'] if (t['context'],t['target'],t['dose'])==(f['contexts'][C+'_'+P],f['cell'],dose)) for C in ('S','R')] for dose in study.DOSES}
            a,b=pair['amount']
            assert a['requested_removed_amount']==b['requested_removed_amount']
            assert a['factor']!=b['factor']
            assert pair['fractional'][0]['requested_removed_amount']!=pair['fractional'][1]['requested_removed_amount']
            for dose in study.DOSES:
                for t in pair[dose]:
                    c=next(c for c in p['contexts'] if c['key']==t['context'])
                    with np.load(c['source']) as z:initial=z['initial'].copy();ids=z['ids'];volumes=z['volumes']
                    i=ids.tolist().index(t['target']);expected=initial.copy();expected[0,i]*=t['factor']
                    with np.load(t['source']) as z:np.testing.assert_array_equal(z['initial'],expected)
                    removal=(initial[0,i]-expected[0,i])*volumes[i]
                    assert abs(removal/t['requested_removed_amount']-1)<2e-13


def test_source_drift_is_reported_without_dropping_unsettled_background():
    delta=np.array([[-.2,.2],[.1,-.1]])
    initial=np.ones((2,2));times=np.array([0.,36.,48.,60.]);path=np.stack([initial*1.2,initial*1.1,initial*1.05,initial])
    p=dict(beta=2.,da=.02,db=.55)
    r=study.source_check(path,initial,delta,times,p,{'reference':initial.copy()})
    assert r['source_gate_passed'] and not r['descriptive_settling_screen']
    assert r['preceding_max_raw_log_drift']>0.09 and r['challenge_max_absolute_log_rate']==0.
    r=study.source_check(np.repeat(initial[None],4,axis=0),initial,delta,times,p,{'reference':initial.copy()})
    assert r['descriptive_settling_screen']


def test_source_gate_rejects_wrong_endpoint_or_inaccurate_donor():
    delta=np.array([[-.2,.2],[.1,-.1]]);initial=np.ones((2,2));p=dict(beta=2.,da=.02,db=.55)
    times=np.array([0.,60.]);path=np.repeat(initial[None],2,axis=0)
    with pytest.raises(ValueError,match='endpoint'):study.source_check(path,initial*1.01,delta,times,p,{'reference':initial})
    with pytest.raises(ValueError,match='source gate'):study.source_check(path,initial,delta,times,p,{'reference':initial*1.001})


def test_grouping_preserves_all_141_delayed_labels_and_operator_isolation(tmp_path):
    p,_,_=fixture(tmp_path);u,d,m,n=study.grouped.batch_layout(p)
    assert len(m)==141 and d.shape==(18,18) and n==2
    assert [x['trial'] for x in m]==[t['key'] for t in p['trials']]
    np.testing.assert_array_equal(study.frozen.unpack(u,m,n),np.stack([study.frozen.load_context(p,t)[3] for t in p['trials']]))
    for item,t in zip(m,p['trials']):
        _,delta,_,_,_=study.frozen.load_context(p,t);start=item['context_index']*n
        np.testing.assert_array_equal(d[start:start+n,start:start+n],delta)
        mask=np.ones(d.shape[1],bool);mask[start:start+n]=False
        assert not np.any(d[start:start+n,mask])


def test_analysis_keeps_both_doses_and_original_matching_references(tmp_path):
    p,_,_=fixture(tmp_path,60.);times=np.arange(401)*.15
    paths=synthetic(p,times,'fine')
    summary,waves=study.summarize(paths,'fine',p,times)
    assert len(summary['metrics'])==96 and len(summary['effects'])==120 and len(waves)==120
    assert {x['dose'] for x in summary['effects']}=={'fractional','amount'}
    assert {x['matching_reference'] for x in summary['effects']}=={'S','R'}
    assert sum('response_shift_detected' in x for x in summary['effects'])==48
    assert len({x['key'] for x in summary['effects']})==120
    for f in p['factorials']:
        for dose in study.DOSES:
            prefix=f'{f["seed"]}_{f["cell"]}_{f["matching_reference"]}_{dose}_'
            np.testing.assert_allclose(waves[prefix+'interaction_difference_of_differences'],
                waves[prefix+'chemical_redistributed_physics']-waves[prefix+'chemical_sham_physics'],atol=1e-15)


def test_tampered_delayed_start_rejected(tmp_path):
    p,_,_=fixture(tmp_path);t=p['trials'][1]
    with np.load(t['source']) as z:changed=z['initial'].copy()
    changed[0,0]*=1.01;study.save_npz(t['source'],initial=changed)
    with pytest.raises(ValueError,match='preparation'):study.frozen.load_context(p,t)


def test_full_small_delayed_pipeline_keeps_every_numerical_gate(tmp_path):
    p,_,_=fixture(tmp_path)
    result=study.run(tmp_path,workers=2)
    assert result['passed']
    assert len(result['production_accuracy'])==423 and len(result['orders'])==141
    assert len(result['raw_timestep_comparisons'])==282 and len(result['response_comparisons'])==9
    assert all(len(x['pulse_rows'])==96 and len(x['factorial_contrasts'])==120 for x in result['response_comparisons'])
    assert len(result['history_outcomes'])==24 and len(result['immediate_vs_delayed'])==24
    assert len(result['matching_reference_response_contrasts'])==48
    assert study.read(tmp_path/'timing.json')['cpu_verification_overlaps_gpu']
