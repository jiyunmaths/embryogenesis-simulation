"""Match challenge-time surrounding budgets on accepted endpoint graphs.

No shared solver is changed. Measured S/R budgets are both crossed with S/R
local matching references. Native and unmatched preparations remain controls.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from multiprocessing import get_context
import os
from pathlib import Path
import time

import numpy as np
import torch

from . import compensated_chemistry as chemistry
from . import compensated_chemistry_validation as qualification
from . import compensated_recipient_response as branch
from . import phase_carry_frozen_response as frozen
from . import phase_carry_recipient_state_response as matched
from . import phase_carry_recipient_state_response_grouped as grouped
from . import phase_carry_shared_context_response as shared
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

write_json=frozen.write_json
save_npz=frozen.save_npz
AMOUNT_MAX=2e-13


def accepted_parent(parent, review_file):
    parent=Path(parent).resolve();review_file=Path(review_file)
    p=read(parent/'protocol.json');review=read(review_file);audit_file=Path(review['independent_review'])
    if not review['full_study_accepted'] or review['protocol_sha256']!=digest(parent/'protocol.json'):
        raise ValueError('Independently accepted challenge-time parent required')
    if digest(audit_file)!=review['independent_review_sha256']:
        raise ValueError('Changed independent parent review')
    audit=read(audit_file)
    if not audit['passed'] or not audit['accepted_at_tested_scope'] or audit['protocol_sha256']!=review['protocol_sha256']:
        raise ValueError('Parent saved arrays not independently accepted')
    verify(p);inputs=[review_file,audit_file]
    for name in ('protocol.json','summary.json','status.json'):
        file=parent/name
        if audit['evidence_sha256'].get(str(file))!=digest(file):raise ValueError('Changed accepted parent report')
        inputs.append(file)
    if read(parent/'status.json')['state']!='completed' or not read(parent/'summary.json')['passed']:
        raise ValueError('Completed passing parent required')
    saved={}
    for c in p['contexts']:
        file=Path(c['source'])
        if audit['evidence_sha256'].get(str(file.resolve()))!=digest(file):raise ValueError('Changed accepted challenge context')
        with np.load(file,allow_pickle=False) as z:
            d,v,u,ids=frozen.validate_context(*(z[k] for k in ('delta','volumes','initial','ids')))
            saved[c['key']]=dict(delta=d.copy(),volumes=v.copy(),initial=u.copy(),ids=ids.copy())
        inputs.append(file)
    if (len(p['contexts']),len(p['trials']),len(p['factorials']))!=(72,144,18) or p['histories']!=[7,8,9] or p['start']!=570.:
        raise ValueError('Wrong challenge-time parent design')
    return p,saved,inputs


def rescale_surrounding(data, target, common):
    """Uniform species-specific scale preserves spatial amount fractions."""
    d,v,u,ids=frozen.validate_context(*(data[k] for k in ('delta','volumes','initial','ids')))
    if target not in ids:raise ValueError('Missing recipient')
    common=np.asarray(common)
    if common.dtype!=np.float64 or common.shape!=(2,) or not np.isfinite(common).all() or np.any(common<=0):
        raise ValueError('Positive finite float64 two-species budget required')
    i=ids.tolist().index(target);others=ids!=target;before=u[:,others]@v[others]
    if not others.any() or np.any(before<=0):raise ValueError('Positive surrounding capacity and amounts required')
    changed=u.copy();scale=common/before;changed[:,others]*=scale[:,None]
    frozen.valid_state(changed)
    np.testing.assert_array_equal(changed[:,i],u[:,i])
    after=changed[:,others]@v[others]
    fractions_before=u[:,others]*v[others]/before[:,None]
    fractions_after=changed[:,others]*v[others]/after[:,None]
    amount_error=float(abs(after/common-1).max());allocation_error=float(abs(fractions_after/fractions_before-1).max())
    if max(amount_error,allocation_error)>AMOUNT_MAX:raise ValueError('Budget/allocation closure failed')
    record=dict(target=target,target_pair=u[:,i].tolist(),surrounding_indices=np.flatnonzero(others).tolist(),
        pre_surrounding_species_amounts=before.tolist(),common_surrounding_species_amounts=common.tolist(),
        scale_factors=scale.tolist(),post_surrounding_species_amounts=after.tolist(),
        pre_global_species_amounts=(u@v).tolist(),post_global_species_amounts=(changed@v).tolist(),
        pre_per_cell_species_amounts=(u*v).tolist(),post_per_cell_species_amounts=(changed*v).tolist(),
        per_cell_species_amount_change=((changed-u)*v).tolist(),
        pre_surrounding_allocation_fractions=fractions_before.tolist(),post_surrounding_allocation_fractions=fractions_after.tolist(),
        maximum_relative_amount_error=amount_error,maximum_relative_allocation_error=allocation_error)
    return dict(delta=d.copy(),volumes=v.copy(),ids=ids.copy(),initial=changed),record


def design(old, saved, root):
    root=Path(root);contexts=[];trials=[];factorials=[];states={}
    for f in old['factorials']:
        seed,cell,M=f['seed'],f['cell'],f['matching_reference']
        # Native recipient chemistry is a no-reset comparator, not a matched
        # local-state test; do not add budget-only native arms unnecessarily.
        for B in (('native',) if M=='native' else ('native','S','R')):
            combinations={}
            for P in ('S','R'):
                pair=[saved[f['contexts'][C+'_'+P]] for C in ('S','R')]
                for k in ('delta','volumes','ids'):np.testing.assert_array_equal(pair[0][k],pair[1][k])
                ids=pair[0]['ids'];others=ids!=cell;i=ids.tolist().index(cell)
                if M!='native':np.testing.assert_array_equal(pair[0]['initial'][:,i],pair[1]['initial'][:,i])
                common=None if B=='native' else pair[0 if B=='S' else 1]['initial'][:,others]@pair[0]['volumes'][others]
                outputs=[]
                for C,z in zip(('S','R'),pair):
                    key=f'seed-{seed}_cell-{cell}_at570_match-{M}_budget-{B}_chemistry-{C}_physics-{P}'
                    parent_key=f['contexts'][C+'_'+P]
                    if common is None:data={k:v.copy() for k,v in z.items()};account=None
                    else:data,account=rescale_surrounding(z,cell,common)
                    outputs.append(data)
                    contexts.append(dict(key=key,seed=seed,background=key,targets=[cell],matching_reference=M,budget_reference=B,
                        chemical_origin=C,physical_origin=P,parent_context=parent_key,
                        amount_matching=account,post_reset_global_species_amounts=(data['initial']@data['volumes']).tolist(),
                        source=str(root/'contexts'/(key+'.npz'))))
                    states[key]=data;combinations[C+'_'+P]=key
                    for target in (None,cell):
                        suffix='control' if target is None else f'cell-{cell}_fractional';initial=data['initial'].copy()
                        record=dict(factor=1.,initial_log_amplitude=0.,requested_removed_amount=0.,removed_activator_amount=0.,
                            initial_target_activator_amount=None,initial_target_volume=None,pulse_amount_relative_error=0.)
                        if target is not None:
                            initial,record=shared.shared_pulse(initial,data['volumes'],data['ids'],cell,'fractional',[initial[0,i]]*2)
                        trial=dict(key=key+'_'+suffix,context=key,seed=seed,background=key,target=target,
                            dose=None if target is None else 'fractional',matching_reference=M,budget_reference=B,
                            source=str(root/'starts'/(key+'_'+suffix+'.npz')),**record)
                        trials.append(trial);states[trial['key']]=dict(initial=initial)
                if common is not None:
                    if max(float(abs((outputs[1]['initial']@outputs[1]['volumes'])/(outputs[0]['initial']@outputs[0]['volumes'])-1).max()),
                           float(abs((outputs[1]['initial'][:,others]@outputs[1]['volumes'][others])/(outputs[0]['initial'][:,others]@outputs[0]['volumes'][others])-1).max()))>AMOUNT_MAX:
                        raise ValueError('Paired initial surrounding/global budgets differ')
            factorials.append(dict(seed=seed,cell=cell,matching_reference=M,budget_reference=B,contexts=combinations))
    return contexts,trials,factorials,states


def prepare(root,parent=Path('output/studies/pulse-time-recipient-response'),
            review_file=Path('docs/pulse_time_recipient_response_assessment.json')):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    old,saved,inputs=accepted_parent(parent,review_file)
    contexts,trials,factorials,states=design(old,saved,root)
    if (len(contexts),len(trials),len(factorials))!=(168,336,42):raise ValueError('Unexpected amount factorial size')
    root.mkdir(parents=True)
    try:
        for c in contexts:save_npz(c['source'],**states[c['key']]);inputs.append(Path(c['source']))
        for t in trials:save_npz(t['source'],**states[t['key']]);inputs.append(Path(t['source']))
        n=len(trials)
        p=dict(parent=str(Path(parent).resolve()),parent_protocol_sha256=digest(Path(parent)/'protocol.json'),
            contexts=contexts,trials=trials,factorials=factorials,histories=old['histories'],selections=old['selections'],
            levels=old['levels'],beta=old['beta'],da=old['da'],db=old['db'],start=570.,duration=60.,interval=.15,
            late_window=24.,checkpoint_interval=3.,effect_min=.01,criteria=old['criteria'],reference_settings=old['reference_settings'],
            runtime=old['runtime'],device=old['device'],primary_paths=2*n,diagnostic_paths=n,total_production_paths=3*n,
            cpu_verification_paths=3*n,logical_reference_cases=n,reference_trajectories=2*n,independent_histories=3,new_histories=0,parent_paths_reused=0,
            amount_relative_max=AMOUNT_MAX,physical_context_definitions=old['physical_context_definitions'],
            parent_historical_export_note=old['parent_historical_export_note'],
            solver=dict(method=chemistry.METHOD,correction_initialization='zero at explicit rounded-state t570 matching/rescaling/pulse',mandatory_checkpoint_correction=True),
            intervention_design='Native/unmatched comparators; two local references crossed with both measured surrounding-total references; uniform per-species surrounding rescaling only.',
            scope='Frozen actual t570 co-evolved graphs, immediate recipient/pulse reset with own controls. Matched arms have equal initial local chemistry, surrounding species amounts and global amounts under shared physical context. Relative per-species surrounding allocations preserved. Later amounts evolve freely. Three existing mature histories; no moving reset, settled challenge, nearest-neighbor specificity, autonomy/inheritance, new development or spatial convergence.')
        p['unique_initial_value_problems']=frozen.unique_problem_count(p)
        signatures={tuple((a.dtype.str,a.shape,a.tobytes()) for a in (*[states[t['context']][k] for k in ('delta','volumes','ids')],states[t['key']]['initial'])) for t in trials}
        p['unique_physical_initial_value_problems']=len(signatures)
        initial,d,m,count=grouped.batch_layout(p)
        p['batching']=dict(physical_blocks=d.shape[0]//count,width=initial.shape[-1],logical_paths=len(m),padding_paths=initial.shape[1]*initial.shape[-1]//count-len(m))
        if p['batching']!=dict(physical_blocks=12,width=28,logical_paths=336,padding_paths=0):raise ValueError('Wrong physical grouping')
        for level in p['levels']:
            ns=frozen.substeps(d,level['dt'],p['beta'],p['da'],p['db'])
            for c in contexts:
                if ns!=frozen.substeps(states[c['key']]['delta'],level['dt'],p['beta'],p['da'],p['db']):raise ValueError('Batch changes positivity substeps')
        p['source_sha256']={**old['source_sha256'],**frozen.pin([Path(__file__),Path('tests/test_pulse_time_amount_response.py')])}
        p['input_sha256']={**old['input_sha256'],**frozen.pin([*inputs,Path('docs/pulse_time_amount_response.md')])}
        write_json(root/'protocol.json',p);verify(p)
        write_json(root/'prefix/protocol.json',{**p,'duration':.6,'checkpoint_interval':.3,'scope':'Every new start at three timesteps; GPU/NumPy plus exact compensated restart qualification only.'})
        write_json(root/'status.json',dict(state='prepared',contexts=len(contexts),logical_starts=n,production_paths=3*n))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def effect_key(f, kind, condition):
    return f"{f['seed']}_{f['cell']}_{f['matching_reference']}_budget-{f['budget_reference']}_fractional_{kind}_{condition}"


def summarize(paths,label,p,times):
    metrics={};responses={};effects=[];waves={};diagonals=[]
    for t in p['trials']:
        if t['target'] is None:continue
        ctrl=next(c for c in p['trials'] if c['context']==t['context'] and c['target'] is None)
        _,_,v,_,ids=frozen.load_context(p,t);i=ids.tolist().index(t['target'])
        values,baseline=paths[label,t['key']],paths[label,ctrl['key']]
        responses[t['key']]=np.log(values[:,:,i]/baseline[:,:,i])/t['initial_log_amplitude']
        metrics[t['key']]=frozen.response_metrics(values,baseline,times,v,i,t['factor'])
    for f in p['factorials']:
        four={name:responses[next(t['key'] for t in p['trials'] if t['context']==context and t['target']==f['cell'])] for name,context in f['contexts'].items()}
        common=dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],budget_reference=f['budget_reference'],dose='fractional')
        diagonals.append(dict(**common,normalized_response_shift_rms=frozen.waveform_distance(four['R_R'],four['S_S'],times)))
        for kind,condition,wave in [('chemical','sham_physics',four['R_S']-four['S_S']),('chemical','redistributed_physics',four['R_R']-four['S_R']),
            ('physical','sham_chemistry',four['S_R']-four['S_S']),('physical','redistributed_chemistry',four['R_R']-four['R_S']),
            ('interaction','difference_of_differences',four['R_R']-four['S_R']-four['R_S']+four['S_S'])]:
            key=effect_key(f,kind,condition);value=frozen.waveform_distance(wave,np.zeros_like(wave),times)
            row=dict(key=key,**common,kind=kind,condition=condition,normalized_response_shift_rms=value)
            if kind=='chemical':row['response_shift_detected']=bool(value>p['effect_min'])
            effects.append(row);waves[key]=wave
    # Descriptive changes are signed vector differences before taking RMS.
    for f in p['factorials']:
        if f['matching_reference']=='native' or f['budget_reference']=='native':continue
        original={**f,'budget_reference':'native'}
        for condition in ('sham_physics','redistributed_physics'):
            wave=waves[effect_key(f,'chemical',condition)]-waves[effect_key(original,'chemical',condition)]
            key=effect_key(f,'budget_change',condition)
            effects.append(dict(key=key,seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],budget_reference=f['budget_reference'],dose='fractional',kind='budget_change',condition=condition,normalized_response_shift_rms=frozen.waveform_distance(wave,np.zeros_like(wave),times)))
            waves[key]=wave
    return dict(metrics=metrics,effects=effects,native_diagonal=diagonals),waves


def reference_gate(root,p):
    paths,times,errors,evidence=frozen.reference_path_sets(root,p)
    data={label:summarize(paths,label,p,times) for label in ('dop853','radau')}
    report=shared.compare(paths,'dop853','radau',p,times,{l:x[0] for l,x in data.items()},{l:x[1] for l,x in data.items()},reference=True)
    result=dict(passed=bool(report['passed'] and all(e<=p['criteria']['reference_log_max'] for e in errors.values())),protocol_sha256=digest(Path(root)/'protocol.json'),responses=report,raw_errors=errors,evidence_sha256=evidence)
    write_json(Path(root)/'reference-qualification.json',result);return result


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p)
    paths,times,noise,evidence=frozen.path_sets(root,p);labels=['dop853','radau',*[l['level'] for l in p['levels']]]
    data={label:summarize(paths,label,p,times) for label in labels};summaries={l:x[0] for l,x in data.items()};waves={l:x[1] for l,x in data.items()}
    specs=[('dop853','radau'),('coarse','fine'),('fine','diagnostic')]+[(l['level'],ref) for l in p['levels'] for ref in ('dop853','radau')]
    comparisons=[shared.compare(paths,a,b,p,times,summaries,waves,reference=a=='dop853') for a,b in specs]
    accuracy=[];raw=[];orders=[]
    for t in p['trials']:
        errors=[]
        for l in p['levels']:
            label=l['level'];error=max(float(abs(np.log(paths[label,t['key']]/paths[ref,t['key']])).max()) for ref in ('dop853','radau'))
            cpu=float(abs(np.log(paths[label,t['key']]/paths['cpu-'+label,t['key']])).max());errors.append(error)
            accuracy.append(dict(trial=t['key'],level=label,raw_reference_log_error=error,cpu_torch_log_error=cpu,passed=bool(error<=p['criteria']['raw_log_max'] and cpu<=p['criteria']['cpu_torch_log_max'])))
        for a,b in (('coarse','fine'),('fine','diagnostic')):
            error=float(abs(np.log(paths[a,t['key']]/paths[b,t['key']])).max());raw.append(dict(trial=t['key'],first=a,second=b,raw_log_error=error,passed=bool(error<=p['criteria']['raw_log_max'])))
        resolved=errors[-1]>p['criteria']['reference_noise_factor']*noise[t['key']]
        values=np.log2(np.asarray(errors[:-1])/errors[1:]).tolist() if resolved and all(e>0 for e in errors) else None
        orders.append(dict(trial=t['key'],errors=errors,reference_noise=noise[t['key']],order_resolved=bool(resolved),orders=values,passed=bool(not resolved or values is not None and all(p['criteria']['order_min']<=x<=p['criteria']['order_max'] for x in values))))
    outcomes=[]
    for seed in p['histories']:
        for M,B in [('native','native'),*[(M,B) for M in ('S','R') for B in ('native','S','R')]]:
            for P in ('sham_physics','redistributed_physics'):
                rows=[r for r in summaries['fine']['effects'] if (r['seed'],r['matching_reference'],r['budget_reference'],r['kind'],r['condition'])==(seed,M,B,'chemical',P)]
                if len(rows)!=2:raise ValueError('Missing nested recipient comparison')
                outcomes.append(dict(history=seed,matching_reference=M,budget_reference=B,physical_context=P,recipients=2,above_screen=sum(r['response_shift_detected'] for r in rows),cells=[dict(cell=r['cell'],above_screen=r['response_shift_detected']) for r in rows]))
    baseline=matched.baseline_descriptions(paths,p,times)
    for row,(f,P) in zip(baseline['paired_chemical_states'],[(f,P) for f in p['factorials'] for P in ('S','R')]):row['budget_reference']=f['budget_reference']
    amount_evolution=[]
    for f in p['factorials']:
        for P in ('S','R'):
            ts=[next(t for t in p['trials'] if t['context']==f['contexts'][C+'_'+P] and t['target'] is None) for C in ('S','R')]
            _,_,v,_,ids=frozen.load_context(p,ts[0]);others=ids!=f['cell'];a,b=[paths['fine',t['key']] for t in ts]
            ga,gb=a@v,b@v;sa,sb=a[:,:,others]@v[others],b[:,:,others]@v[others]
            amount_evolution.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],budget_reference=f['budget_reference'],physics=P,times=times.tolist(),global_species_amounts_S=ga.tolist(),global_species_amounts_R=gb.tolist(),surrounding_species_amounts_S=sa.tolist(),surrounding_species_amounts_R=sb.tolist(),initial_global_species_amount_ratios=(gb[0]/ga[0]).tolist(),final_global_species_amount_ratios=(gb[-1]/ga[-1]).tolist(),initial_surrounding_species_amount_ratios=(sb[0]/sa[0]).tolist(),final_surrounding_species_amount_ratios=(sb[-1]/sa[-1]).tolist()))
    prefix=read(root/'prefix/summary.json')
    result=dict(passed=bool(prefix['passed'] and all(e<=p['criteria']['reference_log_max'] for e in noise.values()) and all(r['passed'] for r in comparisons+accuracy+raw+orders)),independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),method=chemistry.METHOD,
        independent_histories=3,new_histories=0,total_production_paths=p['total_production_paths'],cpu_verification_paths=p['cpu_verification_paths'],reference_trajectories=p['reference_trajectories'],unique_initial_value_problems=p['unique_initial_value_problems'],unique_physical_initial_value_problems=p['unique_physical_initial_value_problems'],
        references_raw_errors=noise,production_accuracy=accuracy,raw_timestep_comparisons=raw,orders=orders,response_comparisons=comparisons,response_summaries=summaries,history_outcomes=outcomes,baseline_descriptions=baseline,amount_evolution=amount_evolution,prefix=prefix,evidence_sha256=evidence,scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=8,qualify_only=False):
    root=Path(root).resolve()
    if not 1<=workers<=8:raise ValueError('Use one to eight single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
        clock=time.perf_counter();n=len(p['trials'])
        try:
            prefix_file=root/'prefix/summary.json'
            if prefix_file.exists():prefix=read(prefix_file)
            else:
                prefix=qualification.prefix(root)
                # Correct the inherited helper's fixed descriptive case count
                # on this new report; its stepping already covers every trial.
                for row in prefix['levels']:row['logical_cases']=n
                prefix['actual_start_cases']=n;prefix['actual_context_timestep_cases']=3*n
                write_json(prefix_file,prefix)
            if not prefix['passed'] or prefix['method']!=chemistry.METHOD or prefix['protocol_sha256']!=digest(root/'prefix/protocol.json'):raise ValueError('Actual-start qualification failed')
            if prefix['actual_start_cases']!=n or prefix['actual_context_timestep_cases']!=3*n or any(row['logical_cases']!=n for row in prefix['levels']):raise ValueError('Wrong actual-start counts')
            if qualify_only:
                write_json(root/'status.json',dict(state='qualified',actual_start_cases=n,timestep_prefix_cases=3*n,scientific_results_pending=True));return prefix
            audit=read(root/'prefix/independent-audit.json')
            if not audit['passed'] or audit['protocol_sha256']!=digest(root/'prefix/protocol.json') or audit['summary_sha256']!=digest(root/'prefix/summary.json') or audit['actual_context_timestep_cases']!=3*n:raise ValueError('Separate actual-array prefix audit required')
            for file,h in audit['evidence_sha256'].items():
                if digest(file)!=h:raise ValueError('Changed separately audited prefix')
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                tasks=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,f in enumerate(as_completed(tasks),1):
                    if not f.result()['passed']:raise RuntimeError('Dual-reference raw gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=n))
                if not reference_gate(root,p)['passed']:raise RuntimeError('Normalized reference gate failed')
                tasks=[pool.submit(branch.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for l in p['levels']:
                    for f in tasks:
                        if f.done():f.result()
                    write_json(root/'status.json',dict(state='running',stage='gpu_and_cpu_production',level=l['level'],cpu_total=3*n));branch.production(root,p,l)
                for i,f in enumerate(as_completed(tasks),1):
                    f.result();write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=3*n))
            result=assess(root);write_json(root/'timing.json',dict(seconds=time.perf_counter()-clock,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed'],independent_review_pending=True,total_production_paths=3*n));return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','qualify','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/pulse-time-amount-response'));parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args();r=prepare(args.output) if args.action=='prepare' else assess(args.output) if args.action=='assess' else run(args.output,args.workers,args.action=='qualify')
    print({k:v for k,v in r.items() if k in ('passed','total_production_paths','unique_physical_initial_value_problems')})
