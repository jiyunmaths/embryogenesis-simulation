"""Match initial local chemistry and surrounding species amounts, then challenge."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import fcntl
from multiprocessing import get_context
import os
from pathlib import Path
import time

import numpy as np
import torch

from . import phase_carry_frozen_response as frozen
from . import phase_carry_shared_context_response as shared
from . import phase_carry_recipient_state_response as matched
from . import phase_carry_recipient_state_response_grouped as grouped
from . import phase_carry_recipient_delayed_response as delayed
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

AMOUNT_RELATIVE_MAX=2e-13
write_json=frozen.write_json
save_npz=frozen.save_npz


def equalize_surrounding(first,second,target):
    """Preserve target state and each species' surrounding allocation fractions."""
    for z in (first,second):frozen.validate_context(z['delta'],z['volumes'],z['initial'],z['ids'])
    for key in ('delta','volumes','ids'):
        if not np.array_equal(first[key],second[key]):raise ValueError('Shared physics/capacities required')
    if target not in first['ids']:raise ValueError('Recipient absent')
    i=first['ids'].tolist().index(target);others=np.arange(len(first['ids']))!=i
    if not np.any(others):raise ValueError('At least one surrounding compartment required')
    if not np.array_equal(first['initial'][:,i],second['initial'][:,i]):raise ValueError('Initially matched recipient pair required')
    volumes=first['volumes'];budgets=[z['initial'][:,others]@volumes[others] for z in (first,second)]
    common=np.sqrt(budgets[0])*np.sqrt(budgets[1])
    if not np.isfinite(common).all() or np.any(common<=0):raise ValueError('Positive finite surrounding budgets required')
    output=[];records=[]
    for z,before in zip((first,second),budgets):
        initial=z['initial'];changed=initial.copy();scale=common/before
        changed[:,others]*=scale[:,None];frozen.valid_state(changed)
        after=changed[:,others]@volumes[others]
        fraction_before=initial[:,others]*volumes[others]/before[:,None]
        fraction_after=changed[:,others]*volumes[others]/after[:,None]
        error=float(np.max(abs(after/common-1.)))
        allocation_error=float(np.max(abs(fraction_after/fraction_before-1.)))
        if max(error,allocation_error)>AMOUNT_RELATIVE_MAX:raise ValueError('Amount/allocation preparation closure failed')
        if not np.array_equal(changed[:,i],initial[:,i]):raise ValueError('Target changed during amount matching')
        output.append(dict(delta=z['delta'].copy(),volumes=volumes.copy(),ids=z['ids'].copy(),initial=changed))
        records.append(dict(target=target,target_pair=initial[:,i].tolist(),surrounding_indices=np.flatnonzero(others).tolist(),
            pre_surrounding_species_amounts=before.tolist(),common_surrounding_species_amounts=common.tolist(),
            scale_factors=scale.tolist(),post_surrounding_species_amounts=after.tolist(),
            pre_per_cell_species_amounts=(initial*volumes).tolist(),post_per_cell_species_amounts=(changed*volumes).tolist(),
            per_cell_species_amount_change=((changed-initial)*volumes).tolist(),
            pre_global_species_amounts=(initial@volumes).tolist(),post_global_species_amounts=(changed@volumes).tolist(),
            pre_surrounding_allocation_fractions=fraction_before.tolist(),post_surrounding_allocation_fractions=fraction_after.tolist(),
            maximum_relative_amount_error=error,maximum_relative_allocation_error=allocation_error))
    if max(float(np.max(abs(output[1]['initial'][:,others]@volumes[others]/(output[0]['initial'][:,others]@volumes[others])-1.))),
           float(np.max(abs((output[1]['initial']@volumes)/(output[0]['initial']@volumes)-1.))))>AMOUNT_RELATIVE_MAX:
        raise ValueError('Paired surrounding/global amounts differ')
    return output,records


def design(old,saved,root):
    """Recipient-specific budgets prevent sharing formerly native controls."""
    root=Path(root);contexts=[];factorials=[];arrays={}
    for f in old['factorials']:
        combination={}
        for P in ('S','R'):
            names=[f['contexts'][C+'_'+P] for C in ('S','R')]
            output,records=equalize_surrounding(*(saved[k] for k in names),f['cell'])
            for C,previous,data,record in zip(('S','R'),names,output,records):
                origin=next(c for c in old['contexts'] if c['key']==previous)
                key=f'seed-{f["seed"]}_cell-{f["cell"]}_matched-{f["matching_reference"]}_physics-{P}_chemistry-{C}_amountmatched'
                contexts.append(dict(key=key,seed=f['seed'],background=key,targets=[f['cell']],matching_reference=f['matching_reference'],
                    chemical_origin=C,physics_label=P,parent_context=previous,
                    source=str(root/'contexts'/(key+'.npz')),amount_matching=record,
                    original_matching_accounting=origin['matching_accounting'],target_pair=record['target_pair']))
                arrays[key]=data;combination[C+'_'+P]=key
        factorials.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],contexts=combination))
    trials=[];starts={}
    for c in contexts:
        for target in (None,c['targets'][0]):
            key=c['key']+('_control' if target is None else f'_cell-{target}_fractional')
            z=arrays[c['key']];changed=z['initial'].copy()
            record=dict(factor=1.,initial_log_amplitude=0.,requested_removed_amount=0.,removed_activator_amount=0.,
                initial_target_activator_amount=None,initial_target_volume=None,pulse_amount_relative_error=0.)
            if target is not None:
                i=z['ids'].tolist().index(target)
                changed,record=shared.shared_pulse(z['initial'],z['volumes'],z['ids'],target,'fractional',[z['initial'][0,i]]*2)
            t=dict(key=key,context=c['key'],seed=c['seed'],background=c['background'],target=target,
                dose=None if target is None else 'fractional',matching_reference=c['matching_reference'],
                source=str(root/'starts'/(key+'.npz')),**record)
            trials.append(t);starts[key]=changed
    return contexts,trials,factorials,arrays,starts


def finalize(root,p,inputs,protocol_doc):
    p['unique_initial_value_problems']=frozen.unique_problem_count(p)
    signatures=set()
    for trial in p['trials']:
        _,delta,volumes,initial,ids=frozen.load_context(p,trial)
        signatures.add(tuple((a.dtype.str,a.shape,a.tobytes()) for a in (delta,volumes,ids,initial)))
    p['unique_physical_initial_value_problems']=len(signatures)
    initial,matrix,mapping,count=grouped.batch_layout(p)
    p['batching']=dict(mode='accepted identical-physics grouping with independent chemical columns',
        physical_blocks=matrix.shape[0]//count,width=initial.shape[-1],logical_paths=len(mapping),
        padding_paths=initial.shape[1]*initial.shape[-1]//count-len(mapping))
    if p['batching']['physical_blocks']!=9:raise ValueError('Unexpected measured physical contexts')
    for level in p['levels']:
        n=frozen.substeps(matrix,level['dt'],p['beta'],p['da'],p['db'])
        for c in p['contexts']:
            with np.load(c['source']) as z:
                if frozen.substeps(z['delta'],level['dt'],p['beta'],p['da'],p['db'])!=n:raise ValueError('Batch changes positivity substeps')
    p['source_sha256'].update(frozen.pin([Path(__file__),Path('tests/test_phase_carry_recipient_amount_response.py')]))
    p['input_sha256']=frozen.pin([*inputs,Path(protocol_doc)])
    write_json(Path(root)/'protocol.json',p);write_json(Path(root)/'status.json',dict(state='prepared',phase=p['phase'],logical_starts=len(p['trials']),
        unique_initial_value_problems=p['unique_initial_value_problems'],independent_histories=3,new_histories=0))
    return p


def prepare(root,parent=Path('output/studies/phase-carry-recipient-state-response-grouped'),
            review_file=Path('docs/phase_carry_recipient_state_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_recipient_amount_response.md')):
    root,parent=Path(root).resolve(),Path(parent).absolute()
    if root.exists():raise FileExistsError(root)
    old,audit,inputs=shared.accepted_parent(parent,review_file)
    latest,_,latest_inputs=shared.accepted_parent(Path('output/studies/phase-carry-recipient-delayed-response'),
        Path('docs/phase_carry_recipient_delayed_response_assessment.json'))
    if (len(old['contexts']),len(old['trials']),len(old['factorials']))!=(45,93,12) or old['start']!=510. or latest['parent_protocol_sha256']!=digest(parent/'protocol.json'):
        raise ValueError('Accepted immediate and delayed recipient-matching lineage required')
    inputs += latest_inputs;saved={}
    for c in old['contexts']:
        file=Path(c['source'])
        if audit['inputs_sha256'].get(str(file))!=digest(file):raise ValueError('Unaccepted t510 matched source')
        with np.load(file) as z:saved[c['key']]={k:z[k].copy() for k in ('delta','volumes','ids','initial')}
        inputs.append(file)
    contexts,trials,factorials,arrays,starts=design(old,saved,root)
    if (len(contexts),len(trials),len(factorials))!=(48,96,12):raise ValueError('Unexpected amount-controlled design')
    root.mkdir(parents=True)
    try:
        for c in contexts:save_npz(c['source'],**arrays[c['key']]);inputs.append(Path(c['source']))
        for t in trials:save_npz(t['source'],initial=starts[t['key']]);inputs.append(Path(t['source']))
        p=dict(parent=str(parent),parent_protocol_sha256=digest(parent/'protocol.json'),phase='immediate',
            preceding_delayed_parent='output/studies/phase-carry-recipient-delayed-response',
            contexts=contexts,trials=trials,factorials=factorials,histories=old['histories'],selections=old['selections'],
            levels=old['levels'],beta=old['beta'],da=old['da'],db=old['db'],start=510.,duration=60.,interval=.15,
            late_window=24.,checkpoint_interval=3.,effect_min=.01,criteria=old['criteria'],reference_settings=old['reference_settings'],
            runtime=old['runtime'],device=old['device'],amount_relative_max=AMOUNT_RELATIVE_MAX,
            primary_paths=192,diagnostic_paths=96,total_production_paths=288,cpu_verification_paths=288,
            logical_reference_cases=96,reference_trajectories=192,independent_histories=3,new_histories=0,parent_paths_reused=0,
            source_sha256={**old['source_sha256'],**latest['source_sha256']},
            scope='Initial recipient pairs and surrounding/whole species amounts equal within each shared physics/reference pair. Symmetric geometric-mean budgets preserve species allocation fractions. Immediate pulse and own controls, frozen physics; later amounts free. Three mature histories, no new development or autonomous identity.')
        # The accepted analysis reads these old curves only for optional native checks.
        inputs += [parent/'production'/level/f for level in ('coarse','fine') for f in ('paths.npz','result.json')]
        return finalize(root,p,inputs,protocol_doc)
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def prepare_delayed(root,parent=Path('output/studies/phase-carry-recipient-amount-response'),
                    review_file=Path('docs/phase_carry_recipient_amount_response_assessment.json'),
                    protocol_doc=Path('docs/phase_carry_recipient_amount_response.md')):
    root,parent=Path(root).resolve(),Path(parent).absolute()
    if root.exists():raise FileExistsError(root)
    old,audit,inputs=shared.accepted_parent(parent,review_file)
    if old.get('phase')!='immediate' or (len(old['contexts']),len(old['trials']),len(old['factorials']))!=(48,96,12) or old['start']!=510. or old['duration']!=60.:
        raise ValueError('Accepted complete amount-controlled immediate phase required')
    production={}
    for level in ('coarse','fine','diagnostic'):
        file=parent/'production'/level/'paths.npz';result=parent/'production'/level/'result.json'
        for f in (file,result):
            if audit['inputs_sha256'].get(str(f))!=digest(f):raise ValueError('Changed qualified waiting paths')
        with np.load(file) as z:
            production[level]={str(k):v.copy() for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
            times=z['times'].copy()
            if not np.array_equal(times,np.arange(401)*.15):raise ValueError('Changed waiting clock')
        inputs += [file,result]
    arrays={};checks={}
    for c in old['contexts']:
        ctrl=next(t for t in old['trials'] if t['context']==c['key'] and t['target'] is None)
        files=[Path(c['source']),parent/'references'/ctrl['key']/'paths.npz',parent/'references'/ctrl['key']/'result.json']
        for file in files:
            if audit['inputs_sha256'].get(str(file))!=digest(file):raise ValueError('Unaccepted waiting source/reference')
        with np.load(files[0]) as z:data={k:z[k].copy() for k in ('delta','volumes','ids')}
        path=production['fine'][ctrl['key']];data['initial']=path[-1].copy()
        endpoints={level:production[level][ctrl['key']][-1] for level in ('coarse','diagnostic')}
        with np.load(files[1]) as z:
            if not np.array_equal(z['times'],times):raise ValueError('Waiting reference clock differs')
            endpoints.update({label:z[label][-1].copy() for label in ('dop853','radau')})
        checks[c['key']]=delayed.source_check(path,data['initial'],data['delta'],times,old,endpoints)
        checks[c['key']]['parent_control']=ctrl['key'];arrays[c['key']]=data;inputs += files
    contexts,trials,factorials,starts=delayed.design(old,arrays,root)
    if (len(contexts),len(trials),len(factorials))!=(48,144,12):raise ValueError('Unexpected delayed amount-control design')
    root.mkdir(parents=True)
    try:
        for c in contexts:
            z=arrays[c['parent_context']];save_npz(c['source'],**z);inputs.append(Path(c['source']))
            c['settling']=checks[c['parent_context']]
            c['challenge_global_species_amounts']=(z['initial']@z['volumes']).tolist()
            c['challenge_per_cell_species_amounts']=(z['initial']*z['volumes']).tolist()
        for t in trials:save_npz(t['source'],initial=starts[t['key']]);inputs.append(Path(t['source']))
        p={**{k:v for k,v in old.items() if k not in ('input_sha256','batching')},
            'parent':str(parent),'parent_protocol_sha256':digest(parent/'protocol.json'),'phase':'delayed',
            'contexts':contexts,'trials':trials,'factorials':factorials,'start':570.,'wait':60.,'original_matching_time':510.,
            'source_log_max':delayed.SOURCE_LOG_MAX,'settling_criteria':delayed.SETTLING,'source_checks':checks,
            'primary_paths':288,'diagnostic_paths':144,'total_production_paths':432,'cpu_verification_paths':432,
            'logical_reference_cases':144,'reference_trajectories':288,
            'accepted_immediate_fractional_effects':read(parent/'summary.json')['response_summaries']['fine']['effects'],
            'scope':'Fixed60-unit wait after initial local and surrounding-amount matching; actual qualified fine unpulsed endpoint, no rematching or continued clamp. Fractional/equal-amount pulses; frozen shared physics. Later local/global chemistry may differ. Three mature histories, no new development or autonomous identity.'}
        return finalize(root,p,inputs,protocol_doc)
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def amount_evolution(root,p):
    """Describe endogenous amount changes; impose no continuing budget clamp."""
    with np.load(Path(root)/'production/fine/paths.npz') as z:
        times=z['times'].copy();paths={str(k):v for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
    rows=[]
    for f in p['factorials']:
        for P in ('S','R'):
            data=[]
            for C in ('S','R'):
                c=next(c for c in p['contexts'] if c['key']==f['contexts'][C+'_'+P])
                t=next(t for t in p['trials'] if t['context']==c['key'] and t['target'] is None)
                with np.load(c['source']) as z:v=z['volumes'].copy();ids=z['ids'].copy()
                others=ids!=f['cell'];u=paths[t['key']]
                data.append((u@v,u[:,:,others]@v[others]))
            global_ratio=data[1][0]/data[0][0];surrounding_ratio=data[1][1]/data[0][1]
            closure=max(float(np.max(abs(global_ratio[0]-1))),float(np.max(abs(surrounding_ratio[0]-1))))
            if p['phase']=='immediate' and closure>AMOUNT_RELATIVE_MAX:raise ValueError('Initial amount closure lost in production')
            rows.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],physics=P,
                times=times.tolist(),global_species_amounts_S=data[0][0].tolist(),global_species_amounts_R=data[1][0].tolist(),
                surrounding_species_amounts_S=data[0][1].tolist(),surrounding_species_amounts_R=data[1][1].tolist(),
                global_species_amount_ratios=global_ratio.tolist(),surrounding_species_amount_ratios=surrounding_ratio.tolist(),
                initial_global_species_amount_ratios=global_ratio[0].tolist(),final_global_species_amount_ratios=global_ratio[-1].tolist(),
                maximum_global_species_log_amount_separation=np.max(abs(np.log(global_ratio)),axis=0).tolist(),
                initial_paired_amount_relative_error=closure))
    return rows


def assess(root):
    root=Path(root);p=read(root/'protocol.json')
    result=(matched.assess if p['phase']=='immediate' else delayed.assess)(root)
    result['phase']=p['phase'];result['amount_evolution']=amount_evolution(root,p)
    write_json(root/'summary.json',result);return result


def run(root,workers=8):
    root=Path(root).resolve()
    if not 1<=workers<=8:raise ValueError('Use one to eight single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
        start=time.perf_counter()
        qualify=matched.reference_qualification if p['phase']=='immediate' else delayed.reference_qualification
        try:
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                pending=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,f in enumerate(as_completed(pending),1):
                    if not f.result()['passed']:raise RuntimeError('Raw reference gate failed')
                    write_json(root/'status.json',dict(state='running',phase=p['phase'],stage='independent_references',completed=i,total=len(pending)))
                if not qualify(root,p)['passed']:raise RuntimeError('Normalized reference gate failed; production blocked')
                pending=[pool.submit(frozen.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for level in p['levels']:
                    for f in pending:
                        if f.done():f.result()
                    write_json(root/'status.json',dict(state='running',phase=p['phase'],stage='production_and_cpu_verification',level=level['level'],cpu_total=len(pending)))
                    grouped.production(root,p,level)
                for i,f in enumerate(as_completed(pending),1):
                    f.result();write_json(root/'status.json',dict(state='running',phase=p['phase'],stage='cpu_verification',completed=i,total=len(pending)))
            result=assess(root)
            write_json(root/'timing.json',dict(coordinator_wall_seconds=time.perf_counter()-start,cpu_verification_overlaps_gpu=True,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed'],phase=p['phase'],
                total_production_paths=p['total_production_paths'],independent_histories=3,new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',phase=p['phase'],error=str(error)));raise


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('action',choices=('prepare','prepare-delayed','run','assess'))
    a.add_argument('--output',type=Path);a.add_argument('--parent',type=Path);a.add_argument('--review',type=Path)
    a.add_argument('--protocol-doc',type=Path,default=Path('docs/phase_carry_recipient_amount_response.md'));a.add_argument('--workers',type=int,default=8)
    args=a.parse_args();root=args.output or Path('output/studies/phase-carry-recipient-amount-'+('delayed-response' if args.action=='prepare-delayed' else 'response'))
    options={k:v for k,v in dict(parent=args.parent,review_file=args.review,protocol_doc=args.protocol_doc).items() if v is not None}
    result=prepare(root,**options) if args.action=='prepare' else prepare_delayed(root,**options) if args.action=='prepare-delayed' else run(root,args.workers) if args.action=='run' else assess(root)
    print({k:v for k,v in result.items() if k in ('passed','phase','primary_paths','diagnostic_paths','total_production_paths','unique_initial_value_problems')})
