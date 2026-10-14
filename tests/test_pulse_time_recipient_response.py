import numpy as np
import pytest

from embryo import pulse_time_recipient_response as study


def endpoints():
    ids=np.array([19,20,21]);v=np.array([1.,2.,.5]);w=np.array([1.3,1.8,.7])
    g=np.array([[0.,.2,.1],[.2,0.,.3],[.1,.3,0.]])
    k=g.copy();np.fill_diagonal(k,-g.sum(axis=1))
    a=np.array([[1.1,.05,.8],[1.,.2,.9]]);b=np.array([[.9,.004,1.2],[.95,.07,1.1]])
    return [dict(seed=8,cell=19,contexts=dict(S='S',R='R'))],dict(
        S=dict(initial=a,volumes=v,delta=k/v[:,None],ids=ids),
        R=dict(initial=b,volumes=w,delta=1.4*k/w[:,None],ids=ids))


def test_crossed_geometry_and_local_reset_preserve_surrounding_concentrations(tmp_path):
    pairs,saved=endpoints();contexts,trials,factorials,states=study.design(pairs,saved,tmp_path)
    assert (len(contexts),len(trials),len(factorials))==(12,24,3)
    for c in contexts:
        z=states[c['key']];chemical=saved[c['chemical_origin']];physical=saved[c['physical_origin']]
        np.testing.assert_array_equal(z['delta'],physical['delta'])
        np.testing.assert_array_equal(z['volumes'],physical['volumes'])
        np.testing.assert_array_equal(z['initial'][:,1:],chemical['initial'][:,1:])
        reference=c['matching_reference']
        expected=chemical['initial'][:,0] if reference=='native' else saved[reference]['initial'][:,0]
        np.testing.assert_array_equal(z['initial'][:,0],expected)
        original=chemical['initial']*chemical['volumes'];reassigned=chemical['initial']*physical['volumes']
        np.testing.assert_allclose(c['reassignment_global_amount_change'],(reassigned-original).sum(axis=1))
        np.testing.assert_allclose(c['post_reset_global_species_amounts'],z['initial']@z['volumes'])
    for f in factorials:
        for physical in ('S','R'):
            a,b=[states[f['contexts'][C+'_'+physical]] for C in ('S','R')]
            np.testing.assert_array_equal(a['delta'],b['delta'])
            np.testing.assert_array_equal(a['volumes'],b['volumes'])
            if f['matching_reference']!='native':np.testing.assert_array_equal(a['initial'][:,0],b['initial'][:,0])


def test_equal_local_state_and_capacity_make_matched_removed_amounts_identical(tmp_path):
    pairs,saved=endpoints();contexts,trials,factorials,states=study.design(pairs,saved,tmp_path)
    for f in factorials:
        for P in ('S','R'):
            amounts=[]
            for C in ('S','R'):
                key=f['contexts'][C+'_'+P];t=next(t for t in trials if t['context']==key and t['target']==19)
                before=states[key]['initial'];after=states[t['key']]['initial'];v=states[key]['volumes']
                expected=before.copy();expected[0,0]*=.9
                np.testing.assert_array_equal(after,expected)
                measured=(before[0,0]-after[0,0])*v[0]
                assert measured==t['removed_activator_amount']
                assert t['pulse_amount_relative_error']<=2e-13
                amounts.append(measured)
            if f['matching_reference']!='native':assert amounts[0]==amounts[1]
            else:assert amounts[0]!=amounts[1]


def test_design_rejects_misaligned_cell_ids(tmp_path):
    pairs,saved=endpoints();saved['R']['ids']=np.array([20,19,21])
    with pytest.raises(AssertionError):study.design(pairs,saved,tmp_path)


def test_summary_measures_matching_effect_against_the_same_geometry(tmp_path):
    pairs,saved=endpoints();contexts,trials,factorials,states=study.design(pairs,saved,tmp_path)
    for c in contexts:study.save_npz(c['source'],**states[c['key']])
    for t in trials:study.save_npz(t['source'],**states[t['key']])
    p=dict(contexts=contexts,trials=trials,factorials=factorials,effect_min=.01)
    times=np.linspace(0,60,401);paths={}
    for c in contexts:
        control=next(t for t in trials if t['context']==c['key'] and t['target'] is None)
        pulse=next(t for t in trials if t['context']==c['key'] and t['target'] is not None)
        # Different control relaxation must not be mistaken for pulse response.
        base=states[c['key']]['initial'][None]*np.exp(.01*times[:,None,None])
        perturbed=base.copy();tau=1. if c['chemical_origin']=='S' else 4.
        if c['matching_reference']!='native':tau=2.
        perturbed[:,0,0]*=np.exp(np.log(.9)*np.exp(-times/tau))
        paths['fine',control['key']]=base;paths['fine',pulse['key']]=perturbed
    summary,waves=study.summarize(paths,'fine',p,times)
    chemical=[r for r in summary['effects'] if r['kind']=='chemical']
    assert all(r['response_shift_detected'] for r in chemical if r['matching_reference']=='native')
    assert all(not r['response_shift_detected'] for r in chemical if r['matching_reference']!='native')
    changes=[r for r in summary['effects'] if r['kind']=='recipient_matching']
    assert len(changes)==4 and all(r['normalized_response_shift_rms']>0 for r in changes)
    for r in changes:
        P=r['condition'];M=r['matching_reference']
        np.testing.assert_allclose(waves[r['key']],waves[f'8_19_{M}_fractional_chemical_{P}']-waves[f'8_19_native_fractional_chemical_{P}'])
