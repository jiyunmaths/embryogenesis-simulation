"""Qualify symmetric capacity-weighted preparation and both phase pipelines."""
from pathlib import Path
import numpy as np
import pytest
import torch

from embryo import phase_carry_recipient_amount_response as study


def matched_input(root,duration=.3):
    ids=np.array([11,17,23,29]);selections={str(seed):dict(selected_ids=[11,17]) for seed in (7,8,9)}
    original={}
    for seed in (7,8,9):
        for cell in (None,11,17):
            name='sham' if cell is None else f'redistributed-cell-{cell}'
            u=np.array([[.8,1.1,.9,.6],[.9,1.2,.7,1.]])*(1.+.02*(seed-7))
            v=np.array([1.,2.,1.3,.8])
            if cell is not None:
                i=ids.tolist().index(cell);u[:,i]*=np.array([.7,1.1]) if cell==11 else np.array([1.2,.85])
                u[:,2]*=np.array([1.2,.9]);u[:,3]*=np.array([.8,1.15])
                v*=np.array([1.1,.9,1.05,.95]) if cell==11 else np.array([.9,1.1,.95,1.05])
            g=np.zeros((4,4))
            for i in range(3):g[i,i+1]=g[i+1,i]=.1+.015*(seed-7)+(0. if cell is None else .01*(ids.tolist().index(cell)+1))
            d=g.copy();np.fill_diagonal(d,-g.sum(axis=1));d/=v[:,None]
            original[f'seed-{seed}_{name}']=dict(delta=d,volumes=v,ids=ids,initial=u)
    contexts,_,factorials=study.shared.design(selections,root/'raw');saved={};assignments={}
    for c in contexts:
        saved[c['key']],assignments[c['key']]=study.shared.transferred_context(original[c['chemical_context']],original[c['physical_context']])
    previous=dict(contexts=contexts,factorials=factorials,assignments=assignments)
    contexts,trials,factorials,arrays=study.matched.design(previous,saved,root/'matched')
    old=dict(contexts=contexts,trials=trials,factorials=factorials,histories=[7,8,9],selections=selections,
        levels=[dict(level=n,dt=dt) for n,dt in study.frozen.LEVELS],beta=2.,da=.02,db=.55,
        start=510.,duration=duration,interval=.15,checkpoint_interval=.3,late_window=24.,effect_min=.01,
        criteria=dict(study.frozen.CRITERIA),reference_settings=dict(rtol=1e-13,atol=1e-15,max_step=.15),
        runtime=dict(torch=str(torch.__version__),torch_cuda=torch.version.cuda,expected_gpu_name='unused'),device='cpu',
        source_sha256={},scope='Small four-compartment test fixture')
    return old,arrays


def fixture(root,phase='immediate',duration=.3):
    old,saved=matched_input(root,duration)
    contexts,trials,factorials,arrays,starts=study.design(old,saved,root)
    p={**old,'phase':'immediate','contexts':contexts,'trials':trials,'factorials':factorials,
        'primary_paths':192,'diagnostic_paths':96,'total_production_paths':288,'cpu_verification_paths':288,
        'logical_reference_cases':96,'reference_trajectories':192,'parent':str(root/'legacy')}
    if phase=='delayed':
        immediate=[]
        for f in factorials:
            for kind,condition in [('chemical','sham_physics'),('chemical','redistributed_physics')]:
                immediate.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],dose='fractional',kind=kind,condition=condition,normalized_response_shift_rms=0.,response_shift_detected=False))
        delayed_arrays={}
        for index,c in enumerate(contexts):
            delayed_arrays[c['key']]={k:v.copy() for k,v in arrays[c['key']].items()}
            delayed_arrays[c['key']]['initial']*=np.array([1+.001*(index+1),1+.0005*(index+1)])[:,None]
        contexts,trials,factorials,starts=study.delayed.design(p,delayed_arrays,root/'delayed')
        arrays={c['key']:delayed_arrays[c['parent_context']] for c in contexts}
        p={**p,'phase':'delayed','contexts':contexts,'trials':trials,'factorials':factorials,
            'start':570.,'wait':60.,'primary_paths':288,'diagnostic_paths':144,'total_production_paths':432,
            'cpu_verification_paths':432,'logical_reference_cases':144,'reference_trajectories':288,
            'accepted_immediate_fractional_effects':immediate,'source_checks':{}}
    for c in contexts:study.save_npz(c['source'],**arrays[c['key']])
    for t in trials:study.save_npz(t['source'],initial=starts[t['key']])
    times=np.arange(round(duration/.15)+1)*.15
    for level in ('coarse','fine'):study.save_npz(Path(p['parent'])/'production'/level/'paths.npz',times=times,trial_keys=np.array([],dtype=str),trajectory=np.empty((len(times),0,2,4)))
    p['input_sha256']=study.frozen.pin([c['source'] for c in contexts]+[t['source'] for t in trials])
    p['unique_initial_value_problems']=study.frozen.unique_problem_count(p)
    study.write_json(root/'protocol.json',p)
    return p,old,saved,arrays


def test_equalization_preserves_target_physics_and_species_allocation():
    ids=np.array([11,17,23,29]);volumes=np.array([.5,2.,1.2,.8]);g=np.ones((4,4))*.1;np.fill_diagonal(g,0.)
    delta=g.copy();np.fill_diagonal(delta,-g.sum(axis=1));delta/=volumes[:,None]
    first=dict(delta=delta,volumes=volumes,ids=ids,initial=np.array([[2.,4.,8.,1.],[3.,5.,2.,4.]]))
    second={**first,'initial':np.array([[2.,3.,7.,5.],[3.,8.,1.,2.]])}
    output,records=study.equalize_surrounding(first,second,11);mask=ids!=11
    common=np.sqrt(first['initial'][:,mask]@volumes[mask])*np.sqrt(second['initial'][:,mask]@volumes[mask])
    for source,result,r in zip((first,second),output,records):
        np.testing.assert_array_equal(result['initial'][:,0],source['initial'][:,0])
        for k in ('delta','volumes','ids'):np.testing.assert_array_equal(result[k],source[k])
        np.testing.assert_allclose(result['initial'][:,mask]@volumes[mask],common,rtol=2e-13,atol=0)
        np.testing.assert_allclose(r['post_surrounding_allocation_fractions'],r['pre_surrounding_allocation_fractions'],rtol=2e-13,atol=0)
        assert r['maximum_relative_amount_error']<=2e-13
    np.testing.assert_allclose(output[0]['initial']@volumes,output[1]['initial']@volumes,rtol=2e-13,atol=0)
    assert not np.array_equal(output[0]['initial'][:,mask],output[1]['initial'][:,mask])
    swapped,_=study.equalize_surrounding(second,first,11)
    for a,b in zip(output,swapped[::-1]):np.testing.assert_array_equal(a['initial'],b['initial'])
    np.testing.assert_array_equal(first['initial'],[[2.,4.,8.,1.],[3.,5.,2.,4.]])


@pytest.mark.parametrize('problem',['local','physics','integer','negative','missing'])
def test_invalid_matching_sources_rejected(problem):
    first=dict(delta=np.array([[-.2,.2],[.1,-.1]]),volumes=np.array([1.,2.]),ids=np.array([11,17]),initial=np.array([[1.,2.],[1.,3.]]))
    second={k:v.copy() for k,v in first.items()};target=11
    if problem=='local':second['initial'][0,0]*=1.1
    if problem=='physics':second['volumes'][0]*=1.1
    if problem=='integer':second['initial']=second['initial'].astype(int)
    if problem=='negative':second['initial'][0,1]=-1.
    if problem=='missing':target=29
    with pytest.raises(ValueError):study.equalize_surrounding(first,second,target)


def test_all_48_backgrounds_control_both_budgets_and_keep_distinct_allocations(tmp_path):
    p,old,saved,arrays=fixture(tmp_path)
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(48,96,12)
    assert p['unique_initial_value_problems']==96
    for f in p['factorials']:
        for P in ('S','R'):
            a,b=[arrays[f['contexts'][C+'_'+P]] for C in ('S','R')];i=a['ids'].tolist().index(f['cell']);others=a['ids']!=f['cell']
            np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i])
            np.testing.assert_allclose(a['initial'][:,others]@a['volumes'][others],b['initial'][:,others]@b['volumes'][others],rtol=2e-13,atol=0)
            np.testing.assert_allclose(a['initial']@a['volumes'],b['initial']@b['volumes'],rtol=2e-13,atol=0)
            assert not np.array_equal(a['initial'][:,others],b['initial'][:,others])
            ts=[next(t for t in p['trials'] if (t['context'],t['target'])==(f['contexts'][C+'_'+P],f['cell'])) for C in ('S','R')]
            assert ts[0]['requested_removed_amount']==ts[1]['requested_removed_amount']
    u,d,m,n=study.grouped.batch_layout(p)
    assert u.shape==(2,36,16) and d.shape==(36,36) and len(m)==96 and n==4
    np.testing.assert_array_equal(study.frozen.unpack(u,m,n),np.stack([study.frozen.load_context(p,t)[3] for t in p['trials']]))


def test_delayed_design_does_not_repeat_local_or_amount_matching(tmp_path):
    p,_,_,arrays=fixture(tmp_path,'delayed')
    assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(48,144,12)
    assert p['unique_initial_value_problems']==120
    u,d,m,n=study.grouped.batch_layout(p)
    assert u.shape==(2,36,24) and len(m)==144
    for f in p['factorials']:
        for P in ('S','R'):
            a,b=[arrays[f['contexts'][C+'_'+P]] for C in ('S','R')];i=a['ids'].tolist().index(f['cell'])
            assert not np.array_equal(a['initial'][:,i],b['initial'][:,i])
            assert not np.allclose(a['initial']@a['volumes'],b['initial']@b['volumes'],rtol=2e-13,atol=0)
            ts=[next(t for t in p['trials'] if (t['context'],t['target'],t['dose'])==(f['contexts'][C+'_'+P],f['cell'],'amount')) for C in ('S','R')]
            assert ts[0]['requested_removed_amount']==ts[1]['requested_removed_amount']


def test_modified_amount_matched_input_rejected_before_stepping(tmp_path):
    p,_,_,_=fixture(tmp_path);c=p['contexts'][0]
    with np.load(c['source']) as z:data={k:z[k].copy() for k in z.files}
    data['initial'][0,2]*=1.01;study.save_npz(c['source'],**data)
    with pytest.raises(ValueError,match='preparation'):study.frozen.load_context(p,p['trials'][0])


@pytest.mark.parametrize('phase',['immediate','delayed'])
def test_complete_small_phase_pipeline_preserves_all_gates(tmp_path,phase):
    p,_,_,_=fixture(tmp_path,phase);result=study.run(tmp_path,workers=2)
    assert result['passed'] and result['phase']==phase
    count=96 if phase=='immediate' else 144
    assert len(result['production_accuracy'])==3*count and len(result['orders'])==count
    assert len(result['raw_timestep_comparisons'])==2*count
    assert len(result['response_comparisons'])==9
    pulses=48 if phase=='immediate' else 96;contrasts=60 if phase=='immediate' else 120
    assert all(len(x['pulse_rows'])==pulses and len(x['factorial_contrasts'])==contrasts for x in result['response_comparisons'])
    assert len(result['amount_evolution'])==24
    if phase=='immediate':assert all(x['initial_paired_amount_relative_error']<=2e-13 for x in result['amount_evolution'])
    else:assert any(x['initial_paired_amount_relative_error']>1e-4 for x in result['amount_evolution'])
    assert study.read(tmp_path/'timing.json')['cpu_verification_overlaps_gpu']
