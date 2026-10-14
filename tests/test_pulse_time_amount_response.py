import numpy as np
import pytest

from embryo import pulse_time_amount_response as study
from embryo import pulse_time_recipient_response as local


def endpoints():
    ids=np.array([19,20,21]);v=np.array([1.,2.,.5]);w=np.array([1.3,1.8,.7])
    g=np.array([[0.,.2,.1],[.2,0.,.3],[.1,.3,0.]])
    k=g.copy();np.fill_diagonal(k,-g.sum(axis=1))
    a=np.array([[1.1,.05,.8],[1.,.2,.9]]);b=np.array([[.9,.004,1.2],[.95,.07,1.1]])
    return [dict(seed=8,cell=19,contexts=dict(S='S',R='R'))],dict(
        S=dict(initial=a,volumes=v,delta=k/v[:,None],ids=ids),
        R=dict(initial=b,volumes=w,delta=1.4*k/w[:,None],ids=ids))


def setup(tmp_path):
    pairs,saved=endpoints();contexts,trials,factorials,arrays=local.design(pairs,saved,tmp_path/'parent')
    old=dict(contexts=contexts,trials=trials,factorials=factorials)
    result=study.design(old,{c['key']:arrays[c['key']] for c in contexts},tmp_path/'new')
    return old,arrays,result


def test_both_measured_budgets_preserve_target_and_surrounding_allocation(tmp_path):
    old,arrays,(contexts,trials,factorials,states)=setup(tmp_path)
    assert (len(contexts),len(trials),len(factorials))==(28,56,7)
    for c in contexts:
        z=states[c['key']];before=arrays[c['parent_context']]
        for k in ('delta','volumes','ids'):np.testing.assert_array_equal(z[k],before[k])
        np.testing.assert_array_equal(z['initial'][:,0],before['initial'][:,0])
        if c['budget_reference']=='native':np.testing.assert_array_equal(z['initial'],before['initial']);continue
        f=next(f for f in old['factorials'] if f['matching_reference']==c['matching_reference'])
        ref=arrays[f['contexts'][c['budget_reference']+'_'+c['physical_origin']]]
        common=ref['initial'][:,1:]@z['volumes'][1:]
        np.testing.assert_allclose(z['initial'][:,1:]@z['volumes'][1:],common,rtol=2e-13,atol=0)
        # A within-species concentration ratio is unchanged by budget matching.
        np.testing.assert_allclose(z['initial'][:,1]/z['initial'][:,2],before['initial'][:,1]/before['initial'][:,2],rtol=2e-13,atol=0)
        np.testing.assert_allclose(c['amount_matching']['per_cell_species_amount_change'],(z['initial']-before['initial'])*z['volumes'])
    for f in factorials:
        for P in ('S','R'):
            a,b=[states[f['contexts'][C+'_'+P]] for C in ('S','R')]
            if f['budget_reference']!='native':
                np.testing.assert_array_equal(a['initial'][:,0],b['initial'][:,0])
                np.testing.assert_allclose(a['initial']@a['volumes'],b['initial']@b['volumes'],rtol=2e-13,atol=0)
                pulse=[next(t for t in trials if t['context']==f['contexts'][C+'_'+P] and t['target']==19) for C in ('S','R')]
                assert pulse[0]['removed_activator_amount']==pulse[1]['removed_activator_amount']
                np.testing.assert_allclose(states[pulse[0]['key']]['initial']@a['volumes'],states[pulse[1]['key']]['initial']@b['volumes'],rtol=2e-13,atol=0)


@pytest.mark.parametrize('bad',[[0.,1.],[1.,np.inf],[1.,np.nan],[1.],np.array([1,2])])
def test_rescaling_rejects_invalid_budget(bad):
    _,saved=endpoints()
    with pytest.raises(ValueError):study.rescale_surrounding(saved['S'],19,np.asarray(bad))


def test_rescaling_rejects_missing_target_and_misaligned_ids(tmp_path):
    _,saved=endpoints()
    with pytest.raises(ValueError):study.rescale_surrounding(saved['S'],99,np.array([1.,1.]))
    old,arrays,_=setup(tmp_path)
    key=old['factorials'][0]['contexts']['R_S'];arrays[key]['ids']=np.array([20,19,21])
    with pytest.raises(AssertionError):study.design(old,arrays,tmp_path/'bad')


def test_response_uses_own_reset_control_and_signed_budget_change(tmp_path):
    _,_,(contexts,trials,factorials,states)=setup(tmp_path)
    for c in contexts:study.save_npz(c['source'],**states[c['key']])
    for t in trials:study.save_npz(t['source'],**states[t['key']])
    p=dict(contexts=contexts,trials=trials,factorials=factorials,effect_min=.01)
    times=np.linspace(0,60,401);paths={}
    for c in contexts:
        control=next(t for t in trials if t['context']==c['key'] and t['target'] is None)
        pulse=next(t for t in trials if t['context']==c['key'] and t['target'] is not None)
        base=states[c['key']]['initial'][None]*np.exp(.01*times[:,None,None]);perturbed=base.copy()
        tau=(1. if c['chemical_origin']=='S' else 4.) if c['budget_reference']=='native' else 2.
        perturbed[:,0,0]*=np.exp(np.log(.9)*np.exp(-times/tau))
        paths['fine',control['key']]=base;paths['fine',pulse['key']]=perturbed
    summary,waves=study.summarize(paths,'fine',p,times)
    chemical=[r for r in summary['effects'] if r['kind']=='chemical']
    assert all(r['response_shift_detected'] for r in chemical if r['budget_reference']=='native')
    assert all(not r['response_shift_detected'] for r in chemical if r['budget_reference']!='native')
    changes=[r for r in summary['effects'] if r['kind']=='budget_change'];assert len(changes)==8
    for r in changes:
        assert 'response_shift_detected' not in r
        f=dict(seed=8,cell=19,matching_reference=r['matching_reference'],budget_reference=r['budget_reference'])
        np.testing.assert_array_equal(waves[r['key']],waves[study.effect_key(f,'chemical',r['condition'])]-waves[study.effect_key({**f,'budget_reference':'native'},'chemical',r['condition'])])
