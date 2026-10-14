"""Challenge freely evolved matched-recipient networks after a fixed wait."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

DOSES = ('fractional', 'amount')
WAIT = 60.
SOURCE_LOG_MAX = 1e-8
SETTLING = dict(window=24., drift_log_max=.001, log_rate_max=.0001)
write_json = frozen.write_json
save_npz = frozen.save_npz


def source_check(path, initial, delta, times, p, other_endpoints):
    """Qualify the measured source and describe drift without choosing timing."""
    frozen.valid_state(path)
    if path.shape != (len(times), *initial.shape) or not np.array_equal(path[-1], initial):
        raise ValueError('Delayed source must be the actual accepted fine control endpoint')
    discrepancies = {label: float(np.max(np.abs(np.log(initial/u))))
                     for label, u in other_endpoints.items()}
    if any(error > SOURCE_LOG_MAX for error in discrepancies.values()):
        raise ValueError('Delayed source endpoint fails the prespecified numerical source gate')
    mask = times >= times[-1]-SETTLING['window']
    drift = float(np.max(np.abs(np.log(path[mask]/initial[None]))))
    rate = frozen.rhs(delta, p['beta'], p['da'], p['db'])(0., initial.ravel()).reshape(initial.shape)/initial
    maximum_rate = float(np.max(np.abs(rate)))
    return dict(endpoint_log_discrepancies=discrepancies, source_gate_passed=True,
        preceding_window=SETTLING['window'], preceding_max_raw_log_drift=drift,
        challenge_max_absolute_log_rate=maximum_rate,
        descriptive_settling_screen=bool(drift <= SETTLING['drift_log_max'] and maximum_rate <= SETTLING['log_rate_max']))


def design(old, arrays, root):
    """Retain all original references; never rematch the freely evolved target."""
    root = Path(root); contexts = []
    for previous in old['contexts']:
        key = previous['key']+'_wait-60'
        contexts.append(dict(**{k:v for k,v in previous.items() if k not in ('key','background','source','parent_context')},
            key=key, background=key, source=str(root/'contexts'/(key+'.npz')),
            parent_context=previous['key'], original_shared_context=previous['parent_context'], original_matching_time=510., wait=WAIT,
            challenge_time=570.))
    by_parent = {c['parent_context']:c for c in contexts}
    factorials = [dict(**{k:v for k,v in f.items() if k != 'contexts'},
        contexts={name:by_parent[key]['key'] for name,key in f['contexts'].items()}) for f in old['factorials']]
    paired = {}
    for f in factorials:
        for physics in ('S','R'):
            c0,c1 = [next(c for c in contexts if c['key']==f['contexts'][chem+'_'+physics]) for chem in ('S','R')]
            z0,z1 = arrays[c0['parent_context']],arrays[c1['parent_context']]
            for key in ('delta','volumes','ids'):
                if not np.array_equal(z0[key],z1[key]):
                    raise ValueError('Shared delayed physics or capacity differs between chemical origins')
            i = z0['ids'].tolist().index(f['cell'])
            values = [float(z0['initial'][0,i]),float(z1['initial'][0,i])]
            for c in (c0,c1):
                key = c['key'],f['cell']
                if key in paired and paired[key] != values:
                    raise ValueError('Shared control has ambiguous paired amount')
                paired[key] = values
    trials = []
    for c in contexts:
        for target,dose in [(None,None),*[(target,dose) for target in c['targets'] for dose in DOSES]]:
            key = c['key']+('_control' if target is None else f'_cell-{target}_{dose}')
            z = arrays[c['parent_context']];initial=z['initial'];changed=initial.copy()
            record = dict(factor=1.,initial_log_amplitude=0.,requested_removed_amount=0.,
                removed_activator_amount=0.,initial_target_activator_amount=None,initial_target_volume=None,
                pulse_amount_relative_error=0.)
            if target is not None:
                changed,record=shared.shared_pulse(initial,z['volumes'],z['ids'],target,dose,paired[c['key'],target])
            trial = dict(key=key,context=c['key'],seed=c['seed'],background=c['background'],target=target,dose=dose,
                matching_reference=c['matching_reference'],source=str(root/'starts'/(key+'.npz')),**record)
            if target is not None:trial['paired_activators']=paired[c['key'],target]
            trials.append((trial,changed))
    return contexts,[t for t,u in trials],factorials,{t['key']:u for t,u in trials}


def prepare(root,parent=Path('output/studies/phase-carry-recipient-state-response-grouped'),
            review_file=Path('docs/phase_carry_recipient_state_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_recipient_delayed_response.md')):
    root,parent=Path(root).resolve(),Path(parent).absolute()
    if root.exists():raise FileExistsError(root)
    old,audit,inputs=shared.accepted_parent(parent,review_file)
    if (len(old['contexts']),len(old['trials']),len(old['factorials'])) != (45,93,12) or old['start'] != 510. or old['duration'] != WAIT:
        raise ValueError('Accepted actual matched-recipient 60-unit parent required')
    production={}
    for level in ('coarse','fine','diagnostic'):
        file=parent/'production'/level/'paths.npz';result=parent/'production'/level/'result.json'
        for f in (file,result):
            if audit['inputs_sha256'].get(str(f)) != digest(f):raise ValueError('Unaccepted parent production evidence')
        with np.load(file,allow_pickle=False) as z:
            production[level]={str(k):v.copy() for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
            times=z['times'].copy()
            if not np.array_equal(times,np.arange(401)*.15):raise ValueError('Changed accepted wait clock')
        inputs += [file,result]
    arrays={};checks={}
    for context in old['contexts']:
        ctrl=next(t for t in old['trials'] if t['context']==context['key'] and t['target'] is None)
        consumed=[Path(context['source']),parent/'references'/ctrl['key']/'paths.npz',parent/'references'/ctrl['key']/'result.json']
        for file in consumed:
            if audit['inputs_sha256'].get(str(file)) != digest(file):raise ValueError('Unaccepted parent source or control reference')
        with np.load(consumed[0],allow_pickle=False) as z:
            data={k:z[k].copy() for k in ('delta','volumes','ids')}
        path=production['fine'][ctrl['key']];data['initial']=path[-1].copy()
        endpoints={level:production[level][ctrl['key']][-1] for level in ('coarse','diagnostic')}
        with np.load(consumed[1],allow_pickle=False) as z:
            if not np.array_equal(z['times'],times):raise ValueError('Changed accepted control-reference clock')
            endpoints.update({method:z[method][-1].copy() for method in ('dop853','radau')})
        checks[context['key']]=source_check(path,data['initial'],data['delta'],times,old,endpoints)
        checks[context['key']]['parent_control']=ctrl['key']
        frozen.validate_context(data['delta'],data['volumes'],data['initial'],data['ids'])
        arrays[context['key']]=data;inputs += consumed
    contexts,trials,factorials,starts=design(old,arrays,root)
    if (len(contexts),len(trials),len(factorials)) != (45,141,12):raise ValueError('Unexpected delayed design counts')
    root.mkdir(parents=True)
    try:
        for c in contexts:
            z=arrays[c['parent_context']];save_npz(c['source'],**z);inputs.append(Path(c['source']))
            c['settling']=checks[c['parent_context']]
            c['challenge_global_species_amounts']=(z['initial']@z['volumes']).tolist()
            c['challenge_per_cell_species_amounts']=(z['initial']*z['volumes']).tolist()
        for t in trials:save_npz(t['source'],initial=starts[t['key']]);inputs.append(Path(t['source']))
        immediate=read(parent/'summary.json')['response_summaries']['fine']['effects']
        p=dict(parent=str(parent),parent_protocol_sha256=digest(parent/'protocol.json'),
            contexts=contexts,trials=trials,factorials=factorials,histories=old['histories'],selections=old['selections'],
            levels=old['levels'],beta=old['beta'],da=old['da'],db=old['db'],start=570.,duration=60.,interval=.15,
            wait=WAIT,original_matching_time=510.,late_window=24.,checkpoint_interval=3.,effect_min=.01,
            criteria=old['criteria'],reference_settings=old['reference_settings'],runtime=old['runtime'],device=old['device'],
            source_log_max=SOURCE_LOG_MAX,settling_criteria=SETTLING,source_checks=checks,
            primary_paths=282,diagnostic_paths=141,total_production_paths=423,cpu_verification_paths=423,
            logical_reference_cases=141,reference_trajectories=282,independent_histories=3,new_histories=0,parent_paths_reused=0,
            accepted_immediate_fractional_effects=immediate,
            scope='Fixed60-unit wait after original matching; challenge on accepted fine t570 controls, no rematching. Fractional and common-amount pulses. Shared frozen physical contexts; target chemistry and surrounding amounts can differ. Three existing mature histories, no new development or autonomous identity.')
        p['unique_initial_value_problems']=frozen.unique_problem_count(p)
        initial,matrix,mapping,count=grouped.batch_layout(p)
        p['batching']=dict(mode='accepted identical-physics grouping with independent chemical columns',
            physical_blocks=matrix.shape[0]//count,width=initial.shape[-1],logical_paths=len(mapping),
            padding_paths=initial.shape[1]*initial.shape[-1]//count-len(mapping))
        if p['batching']['physical_blocks'] != 9:raise ValueError('Unexpected delayed physical grouping')
        for level in p['levels']:
            n=frozen.substeps(matrix,level['dt'],p['beta'],p['da'],p['db'])
            if any(frozen.substeps(z['delta'],level['dt'],p['beta'],p['da'],p['db']) != n for z in arrays.values()):
                raise ValueError('Grouping changes constituent positivity substeps')
        p['source_sha256']={**old['source_sha256'],**frozen.pin([Path(__file__),Path('tests/test_phase_carry_recipient_delayed_response.py')])}
        p['input_sha256']=frozen.pin([*inputs,Path(protocol_doc)])
        write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',logical_starts=141,
            unique_initial_value_problems=p['unique_initial_value_problems'],independent_histories=3,new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def summarize(paths,label,p,times):
    # Reuse the accepted one-dose analysis separately, then retain dose labels.
    metrics={};effects=[];native=[];contrast={}
    for dose in DOSES:
        one={**p,'trials':[t if t['dose'] is None else {**t,'dose':'fractional'} for t in p['trials'] if t['dose'] in (None,dose)]}
        summary,waves=matched.summarize(paths,label,one,times)
        metrics.update(summary['metrics'])
        for row in summary['effects']:
            row={**row,'dose':dose,'key':row['key'].replace('_fractional_',f'_{dose}_')}
            effects.append(row)
        native += [{**row,'dose':dose} for row in summary['native_diagonal']]
        contrast.update({key.replace('_fractional_',f'_{dose}_'):wave for key,wave in waves.items()})
    return dict(metrics=metrics,effects=effects,native_diagonal=native),contrast


def reference_qualification(root,p):
    paths,times,errors,evidence=frozen.reference_path_sets(root,p)
    unpacked={label:summarize(paths,label,p,times) for label in ('dop853','radau')}
    summaries={k:v[0] for k,v in unpacked.items()};waves={k:v[1] for k,v in unpacked.items()}
    result=shared.compare(paths,'dop853','radau',p,times,summaries,waves,reference=True)
    r=dict(protocol_sha256=digest(Path(root)/'protocol.json'),
        passed=bool(all(e<=p['criteria']['reference_log_max'] for e in errors.values()) and result['passed']),
        raw_errors=errors,responses=result,evidence_sha256=evidence)
    write_json(Path(root)/'reference-qualification.json',r);return r


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p)
    paths,times,ref_errors,evidence=frozen.path_sets(root,p)
    labels=['dop853','radau',*[x['level'] for x in p['levels']]]
    unpacked={label:summarize(paths,label,p,times) for label in labels}
    summaries={k:v[0] for k,v in unpacked.items()};waves={k:v[1] for k,v in unpacked.items()}
    comparisons=[shared.compare(paths,'dop853','radau',p,times,summaries,waves,True)]
    comparisons += [shared.compare(paths,a,b,p,times,summaries,waves) for a,b in [('coarse','fine'),('fine','diagnostic')]]
    comparisons += [shared.compare(paths,level['level'],ref,p,times,summaries,waves) for level in p['levels'] for ref in ('dop853','radau')]
    accuracy=[];raw_pairs=[];orders=[]
    for t in p['trials']:
        errors=[]
        for level in p['levels']:
            name=level['level'];error=max(float(np.max(abs(np.log(paths[name,t['key']]/paths[ref,t['key']])))) for ref in ('dop853','radau'))
            cpu=float(np.max(abs(np.log(paths[name,t['key']]/paths['cpu-'+name,t['key']]))))
            accuracy.append(dict(trial=t['key'],level=name,raw_reference_log_error=error,cpu_torch_log_error=cpu,
                passed=bool(error<=p['criteria']['raw_log_max'] and cpu<=p['criteria']['cpu_torch_log_max'])))
            errors.append(error)
        for a,b in [('coarse','fine'),('fine','diagnostic')]:
            error=float(np.max(abs(np.log(paths[a,t['key']]/paths[b,t['key']]))))
            raw_pairs.append(dict(trial=t['key'],first=a,second=b,raw_log_error=error,passed=bool(error<=p['criteria']['raw_log_max'])))
        noise=ref_errors[t['key']];resolved=errors[-1]>p['criteria']['reference_noise_factor']*noise
        values=np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved and all(x>0 for x in errors) else None
        orders.append(dict(trial=t['key'],errors=errors,reference_noise=noise,order_resolved=bool(resolved),orders=values,
            passed=bool(not resolved or (values is not None and all(p['criteria']['order_min']<=x<=p['criteria']['order_max'] for x in values)))))
    outcomes=[]
    for seed in p['histories']:
        for reference in ('S','R'):
            for dose in DOSES:
                for physics in ('sham_physics','redistributed_physics'):
                    rows=[x for x in summaries['fine']['effects'] if (x['seed'],x['matching_reference'],x['dose'],x['kind'],x['condition'])==(seed,reference,dose,'chemical',physics)]
                    if len(rows)!=2:raise ValueError('Missing delayed chemical comparison')
                    outcomes.append(dict(seed=seed,matching_reference=reference,dose=dose,physics=physics,recipients=2,
                        above_screen=sum(x['response_shift_detected'] for x in rows),cells=[dict(cell=x['cell'],above_screen=x['response_shift_detected']) for x in rows]))
    immediate=[]
    for row in summaries['fine']['effects']:
        if row['kind']=='chemical' and row['dose']=='fractional':
            old=next(x for x in p['accepted_immediate_fractional_effects'] if (x['seed'],x['cell'],x['matching_reference'],x['kind'],x['condition'])==
                (row['seed'],row['cell'],row['matching_reference'],row['kind'],row['condition']))
            immediate.append(dict(seed=row['seed'],cell=row['cell'],matching_reference=row['matching_reference'],physics=row['condition'],
                immediate_rms=old['normalized_response_shift_rms'],delayed_rms=row['normalized_response_shift_rms'],
                immediate_above_screen=old['response_shift_detected'],delayed_above_screen=row['response_shift_detected']))
    reference_contrasts=[]
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            fs={f['matching_reference']:f for f in p['factorials'] if (f['seed'],f['cell'])==(seed,cell)}
            for name in ('S_S','R_S','S_R','R_R'):
                for dose in DOSES:
                    pair=[next(t for t in p['trials'] if (t['context'],t['target'],t['dose'])==(fs[L]['contexts'][name],cell,dose)) for L in ('S','R')]
                    i=frozen.load_context(p,pair[0])[4].tolist().index(cell)
                    target_waves=[]
                    for t in pair:
                        c=next(u for u in p['trials'] if u['context']==t['context'] and u['target'] is None)
                        target_waves.append(np.log(paths['fine',t['key']][:,:,i]/paths['fine',c['key']][:,:,i])/t['initial_log_amplitude'])
                    reference_contrasts.append(dict(seed=seed,cell=cell,chemical_physics=name,dose=dose,
                        rms=frozen.waveform_distance(*target_waves,times),maximum_normalized_difference=float(np.max(abs(target_waves[0]-target_waves[1])))))
    result=dict(protocol_sha256=digest(root/'protocol.json'),passed=bool(all(e<=p['criteria']['reference_log_max'] for e in ref_errors.values()) and all(x['passed'] for x in comparisons+accuracy+raw_pairs+orders)),
        primary_paths=p['primary_paths'],diagnostic_paths=p['diagnostic_paths'],total_production_paths=p['total_production_paths'],
        cpu_verification_paths=p['cpu_verification_paths'],logical_reference_cases=p['logical_reference_cases'],
        unique_initial_value_problems=p['unique_initial_value_problems'],reference_trajectories=p['reference_trajectories'],independent_histories=3,new_histories=0,
        references_raw_errors=ref_errors,production_accuracy=accuracy,raw_timestep_comparisons=raw_pairs,orders=orders,
        response_comparisons=comparisons,response_summaries=summaries,baseline_descriptions=matched.baseline_descriptions(paths,p,times),
        history_outcomes=outcomes,immediate_vs_delayed=immediate,matching_reference_response_contrasts=reference_contrasts,
        source_checks=p['source_checks'],evidence_sha256=evidence,scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=4):
    root=Path(root).resolve()
    if not 1<=workers<=4:raise ValueError('Use one to four single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
        start=time.perf_counter()
        try:
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                pending=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,f in enumerate(as_completed(pending),1):
                    if not f.result()['passed']:raise RuntimeError('Delayed raw reference gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=len(pending)))
                if not reference_qualification(root,p)['passed']:raise RuntimeError('Delayed normalized reference gate failed')
                pending=[pool.submit(frozen.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for level in p['levels']:
                    for f in pending:
                        if f.done():f.result()
                    write_json(root/'status.json',dict(state='running',stage='production_and_cpu_verification',level=level['level'],cpu_total=len(pending)))
                    grouped.production(root,p,level)
                for i,f in enumerate(as_completed(pending),1):
                    f.result();write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=len(pending)))
            result=assess(root)
            write_json(root/'timing.json',dict(coordinator_wall_seconds=time.perf_counter()-start,cpu_verification_overlaps_gpu=True,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed'],
                total_production_paths=p['total_production_paths'],independent_histories=3,new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/phase-carry-recipient-delayed-response'))
    parser.add_argument('--parent',type=Path,default=Path('output/studies/phase-carry-recipient-state-response-grouped'))
    parser.add_argument('--review',type=Path,default=Path('docs/phase_carry_recipient_state_response_assessment.json'))
    parser.add_argument('--protocol-doc',type=Path,default=Path('docs/phase_carry_recipient_delayed_response.md'))
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    result=prepare(args.output,args.parent,args.review,args.protocol_doc) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('passed','primary_paths','diagnostic_paths','total_production_paths','unique_initial_value_problems')})
