"""Challenge-time recipient matching on actual co-evolved endpoint graphs.

Uses unchanged qualified frozen chemistry/assessment helpers. Shared endpoint
physics, two measured reset references and a no-reset comparison are explicit.
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

write_json = frozen.write_json
save_npz = frozen.save_npz
REFERENCES = ('native', 'S', 'R')


def accepted_endpoints(parent, review_file):
    """Consume the separately reviewed moving parent without altering its flag."""
    parent = Path(parent).resolve();review_file = Path(review_file)
    p,review = read(parent/'protocol.json'),read(review_file)
    if not review['full_study_accepted'] or read(parent/'status.json')['state']!='completed' or not read(parent/'summary.json')['passed']:
        raise ValueError('Independently accepted completed moving parent required')
    inputs = [review_file]
    for key,name in (('protocol_sha256','protocol.json'),('original_summary_sha256','summary.json'),('original_status_sha256','status.json')):
        file = parent/name
        if digest(file)!=review[key]:raise ValueError('Accepted moving parent changed')
        inputs.append(file)
    verify(p)
    audit_file = Path(review['independent_review'])
    if digest(audit_file)!=review['independent_review_sha256']:raise ValueError('Changed independent moving review')
    audit = read(audit_file)
    if not audit['passed'] or not audit['accepted_at_tested_scope'] or audit['protocol_sha256']!=review['protocol_sha256']:
        raise ValueError('Moving saved-array review did not qualify the parent')
    inputs.append(audit_file);saved={}
    for c in p['contexts']:
        job = next(j for j in p['jobs'] if (j['context'],j['stage'],j['level'])==(c['key'],'wait','fine'))
        folder = parent/job['key'];files = [folder/name for name in ('history.json','result.json','latest_state.npz')]
        for file in files:
            if audit['evidence_sha256'].get(str(file.resolve()))!=digest(file):raise ValueError('Changed actual fine endpoint')
        inputs.extend(files);h=read(files[0]);r=read(files[1])
        if r['job']!=job or not r['quality_pass'] or h[-1]['time']!=570. or h[-1]['elapsed']!=60.:
            raise ValueError('Wrong accepted moving endpoint')
        with np.load(files[2],allow_pickle=False) as z:
            u=np.array([z['activator'],z['inhibitor']]);ids=z['ids'].copy()
            if not np.array_equal(u,h[-1]['chemistry']) or not np.array_equal(ids,h[-1]['ids']):raise ValueError('Moving checkpoint/observation mismatch')
        d,v,u,ids=frozen.validate_context(np.asarray(h[-1]['delta']),np.asarray(h[-1]['volumes']),u,ids)
        saved[c['key']]=dict(delta=d,volumes=v,initial=u,ids=ids)
    if len(saved)!=12 or p['histories']!=[7,8,9] or p['new_histories']!=0:raise ValueError('Wrong history/recipient design')
    return p,saved,inputs,review['historical_guard_audit']


def design(pairs, saved, root):
    root=Path(root);contexts=[];trials=[];factorials=[];states={}
    for pair in pairs:
        seed,cell=pair['seed'],pair['cell']
        first,second=[saved[pair['contexts'][C]] for C in ('S','R')]
        np.testing.assert_array_equal(first['ids'],second['ids'])
        i=first['ids'].tolist().index(cell)
        for reference in REFERENCES:
            combinations={}
            for C in ('S','R'):
                chemical=saved[pair['contexts'][C]]
                for P in ('S','R'):
                    physical=saved[pair['contexts'][P]];values=chemical['initial'].copy()
                    reassigned=values*physical['volumes'];original=values*chemical['volumes']
                    data=dict(delta=physical['delta'].copy(),volumes=physical['volumes'].copy(),initial=values,ids=physical['ids'].copy())
                    accounting=None
                    if reference!='native':
                        target_pair=saved[pair['contexts'][reference]]['initial'][:,i]
                        data,accounting=matched.match_recipient(data,cell,target_pair)
                    key=f'seed-{seed}_cell-{cell}_at570_match-{reference}_chemistry-{C}_physics-{P}'
                    context=dict(key=key,seed=seed,background=key,targets=[cell],matching_reference=reference,
                        chemical_origin=C,physical_origin=P,chemical_parent_context=pair['contexts'][C],physical_parent_context=pair['contexts'][P],
                        matching_accounting=accounting,reassignment_per_cell_amount_change=(reassigned-original).tolist(),
                        reassignment_global_amount_change=(reassigned-original).sum(axis=1).tolist(),
                        post_reset_global_species_amounts=(data['initial']@data['volumes']).tolist(),source=str(root/'contexts'/(key+'.npz')))
                    contexts.append(context);states[key]=data;combinations[C+'_'+P]=key
                    for target in (None,cell):
                        suffix='control' if target is None else f'cell-{cell}_fractional'
                        initial=data['initial'].copy()
                        record=dict(factor=1.,initial_log_amplitude=0.,requested_removed_amount=0.,removed_activator_amount=0.,
                            initial_target_activator_amount=None,initial_target_volume=None,pulse_amount_relative_error=0.)
                        if target is not None:
                            initial,record=shared.shared_pulse(initial,data['volumes'],data['ids'],target,'fractional',[initial[0,i],initial[0,i]])
                        trial=dict(key=key+'_'+suffix,context=key,seed=seed,background=key,target=target,dose=None if target is None else 'fractional',
                            matching_reference=reference,source=str(root/'starts'/(key+'_'+suffix+'.npz')),**record)
                        trials.append(trial);states[trial['key']]=dict(initial=initial)
            factorials.append(dict(seed=seed,cell=cell,matching_reference=reference,contexts=combinations))
    return contexts,trials,factorials,states


def prepare(root,parent=Path('output/studies/moving-recipient-amount-response'),
            review_file=Path('docs/moving_recipient_amount_response_assessment.json')):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    old,saved,inputs,guard_note=accepted_endpoints(parent,review_file)
    contexts,trials,factorials,states=design(old['pairs'],saved,root)
    if (len(contexts),len(trials),len(factorials))!=(72,144,18):raise ValueError('Wrong factorial size')
    root.mkdir(parents=True)
    try:
        for c in contexts:save_npz(c['source'],**states[c['key']]);inputs.append(Path(c['source']))
        for t in trials:save_npz(t['source'],**states[t['key']]);inputs.append(Path(t['source']))
        frozen_parent=read(Path(old['accepted_delayed_parent'])/'protocol.json')
        p=dict(parent=str(Path(parent).resolve()),parent_protocol_sha256=digest(Path(parent)/'protocol.json'),
            contexts=contexts,trials=trials,factorials=factorials,histories=[7,8,9],
            selections={str(seed):dict(selected_ids=[pair['cell'] for pair in old['pairs'] if pair['seed']==seed]) for seed in (7,8,9)},
            levels=[dict(level=l,dt=d) for l,d in frozen.LEVELS],beta=2.,da=.02,db=.55,start=570.,duration=60.,interval=.15,late_window=24.,checkpoint_interval=3.,effect_min=.01,
            criteria=dict(frozen.CRITERIA),reference_settings=frozen_parent['reference_settings'],runtime=old['runtime'],device=old['device'],
            primary_paths=288,diagnostic_paths=144,total_production_paths=432,cpu_verification_paths=432,logical_reference_cases=144,reference_trajectories=288,
            independent_histories=3,new_histories=0,parent_paths_reused=0,
            solver=dict(method=chemistry.METHOD,correction_initialization='zero at explicit t570 rounded-state reassignment/reset',mandatory_checkpoint_correction=True),
            physical_context_definitions=dict(S='Actual fine t570 endpoint from S allocation',R='Actual fine t570 endpoint from R allocation; both had initial sham physics'),
            parent_historical_export_note=guard_note,
            scope='Frozen actual co-evolved t570 endpoint graphs. Recipient-local matching at challenge time using S/R measured pairs, plus no-reset arms on the identical physical contexts. Both whole-network chemical allocations; no total-budget matching, no ongoing clamp. Immediate negative pulses against own resetting controls, not a settled delayed challenge. Three existing mature histories; no moving post-reset continuation, nearest-neighbor specificity, autonomous/inherited identity, spatial convergence or new development.')
        p['unique_initial_value_problems']=frozen.unique_problem_count(p)
        physical_signatures=set()
        for t in trials:
            _,d,v,u,ids=frozen.load_context(p,t)
            physical_signatures.add(tuple((a.dtype.str,a.shape,a.tobytes()) for a in (d,v,ids,u)))
        p['unique_physical_initial_value_problems']=len(physical_signatures)
        initial,d,m,n=grouped.batch_layout(p)
        p['batching']=dict(physical_blocks=d.shape[0]//n,width=initial.shape[-1],logical_paths=len(m),padding_paths=initial.shape[1]*initial.shape[-1]//n-len(m))
        for level in p['levels']:
            count=frozen.substeps(d,level['dt'],p['beta'],p['da'],p['db'])
            for c in contexts:
                if count!=frozen.substeps(states[c['key']]['delta'],level['dt'],p['beta'],p['da'],p['db']):raise ValueError('Batch changes positivity substeps')
        p['source_sha256']={**old['source_sha256'],**frozen.pin([Path(__file__),Path(qualification.__file__),Path('tests/test_pulse_time_recipient_response.py')])}
        p['input_sha256']={**old['input_sha256'],**frozen.pin([*inputs,Path('docs/pulse_time_recipient_response.md')])}
        write_json(root/'protocol.json',p);verify(p)
        prefix={**p,'duration':.6,'checkpoint_interval':.3,'scope':'Every new logical start at three timesteps, GPU/NumPy and bitwise compensated restart qualification only.'}
        write_json(root/'prefix/protocol.json',prefix)
        write_json(root/'status.json',dict(state='prepared',logical_starts=144,contexts=72,production_paths=432))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def summarize(paths,label,p,times):
    summary,waves=matched.summarize(paths,label,p,times)
    for pair in [f for f in p['factorials'] if f['matching_reference']=='native']:
        for reference in ('S','R'):
            for physical in ('sham_physics','redistributed_physics'):
                native_key=f"{pair['seed']}_{pair['cell']}_native_fractional_chemical_{physical}"
                matched_key=f"{pair['seed']}_{pair['cell']}_{reference}_fractional_chemical_{physical}"
                wave=waves[matched_key]-waves[native_key];key=f"{pair['seed']}_{pair['cell']}_{reference}_matching_change_{physical}"
                summary['effects'].append(dict(key=key,seed=pair['seed'],cell=pair['cell'],dose='fractional',matching_reference=reference,
                    kind='recipient_matching',condition=physical,normalized_response_shift_rms=frozen.waveform_distance(wave,np.zeros_like(wave),times)))
                waves[key]=wave
    return summary,waves


def reference_gate(root,p):
    paths,times,errors,evidence=frozen.reference_path_sets(root,p)
    data={label:summarize(paths,label,p,times) for label in ('dop853','radau')}
    summaries={label:x[0] for label,x in data.items()};waves={label:x[1] for label,x in data.items()}
    comparison=shared.compare(paths,'dop853','radau',p,times,summaries,waves,reference=True)
    result=dict(passed=bool(comparison['passed'] and all(e<=p['criteria']['reference_log_max'] for e in errors.values())),
        protocol_sha256=digest(Path(root)/'protocol.json'),responses=comparison,raw_errors=errors,evidence_sha256=evidence)
    write_json(Path(root)/'reference-qualification.json',result);return result


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p)
    paths,times,noise,evidence=frozen.path_sets(root,p);labels=['dop853','radau',*[l['level'] for l in p['levels']]]
    data={label:summarize(paths,label,p,times) for label in labels}
    summaries={label:x[0] for label,x in data.items()};waves={label:x[1] for label,x in data.items()}
    specs=[('dop853','radau'),('coarse','fine'),('fine','diagnostic')]+[(l['level'],ref) for l in p['levels'] for ref in ('dop853','radau')]
    comparisons=[shared.compare(paths,a,b,p,times,summaries,waves,reference=a=='dop853') for a,b in specs]
    accuracy=[];raw=[];orders=[]
    for t in p['trials']:
        errors=[]
        for l in p['levels']:
            label=l['level'];error=max(float(abs(np.log(paths[label,t['key']]/paths[ref,t['key']])).max()) for ref in ('dop853','radau'))
            cpu=float(abs(np.log(paths[label,t['key']]/paths['cpu-'+label,t['key']])).max());errors.append(error)
            accuracy.append(dict(trial=t['key'],level=label,raw_reference_log_error=error,cpu_torch_log_error=cpu,
                passed=bool(error<=p['criteria']['raw_log_max'] and cpu<=p['criteria']['cpu_torch_log_max'])))
        for a,b in (('coarse','fine'),('fine','diagnostic')):
            error=float(abs(np.log(paths[a,t['key']]/paths[b,t['key']])).max())
            raw.append(dict(trial=t['key'],first=a,second=b,raw_log_error=error,passed=bool(error<=p['criteria']['raw_log_max'])))
        resolved=errors[-1]>p['criteria']['reference_noise_factor']*noise[t['key']]
        values=np.log2(np.asarray(errors[:-1])/errors[1:]).tolist() if resolved and all(e>0 for e in errors) else None
        orders.append(dict(trial=t['key'],errors=errors,reference_noise=noise[t['key']],order_resolved=bool(resolved),orders=values,
            passed=bool(not resolved or values is not None and all(p['criteria']['order_min']<=x<=p['criteria']['order_max'] for x in values))))
    outcomes=[]
    for seed in p['histories']:
        for reference in REFERENCES:
            for physical in ('sham_physics','redistributed_physics'):
                effects=[r for r in summaries['fine']['effects'] if (r['seed'],r['matching_reference'],r['kind'],r['condition'])==(seed,reference,'chemical',physical)]
                if len(effects)!=2:raise ValueError('Missing history-level recipient comparison')
                outcomes.append(dict(history=seed,matching_reference=reference,physical_context=physical,recipients=2,
                    above_screen=sum(r['response_shift_detected'] for r in effects),cells=[dict(cell=r['cell'],above_screen=r['response_shift_detected']) for r in effects]))
    prefix=read(root/'prefix/summary.json')
    result=dict(passed=bool(prefix['passed'] and all(e<=p['criteria']['reference_log_max'] for e in noise.values()) and all(r['passed'] for r in comparisons+accuracy+raw+orders)),
        independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),method=chemistry.METHOD,
        independent_histories=3,new_histories=0,total_production_paths=432,cpu_verification_paths=432,reference_trajectories=288,
        unique_initial_value_problems=p['unique_initial_value_problems'],unique_physical_initial_value_problems=p['unique_physical_initial_value_problems'],
        references_raw_errors=noise,production_accuracy=accuracy,raw_timestep_comparisons=raw,orders=orders,response_comparisons=comparisons,
        response_summaries=summaries,history_outcomes=outcomes,baseline_descriptions=matched.baseline_descriptions(paths,p,times),
        prefix=prefix,evidence_sha256=evidence,scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=8,qualify_only=False):
    root=Path(root).resolve()
    if not 1<=workers<=8:raise ValueError('Use one to eight single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
        clock=time.perf_counter()
        try:
            prefix_file=root/'prefix/summary.json'
            prefix=read(prefix_file) if prefix_file.exists() else qualification.prefix(root)
            if (not prefix['passed'] or prefix['method']!=chemistry.METHOD
                    or prefix['protocol_sha256']!=digest(root/'prefix/protocol.json')):
                raise RuntimeError('Actual-start prefix qualification failed')
            if qualify_only:
                write_json(root/'status.json',dict(state='qualified',actual_start_cases=144,timestep_prefix_cases=432,scientific_results_pending=True));return prefix
            audit=read(root/'prefix/independent-audit.json')
            if not audit['passed'] or audit['protocol_sha256']!=digest(root/'prefix/protocol.json') or audit['summary_sha256']!=digest(root/'prefix/summary.json'):
                raise ValueError('Separate prefix saved-array review required before full production')
            for file,h in audit['evidence_sha256'].items():
                if digest(file)!=h:raise ValueError('Changed separately reviewed prefix')
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                tasks=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,f in enumerate(as_completed(tasks),1):
                    if not f.result()['passed']:raise RuntimeError('Dual-reference raw gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=144))
                if not reference_gate(root,p)['passed']:raise RuntimeError('Normalized reference gate failed')
                tasks=[pool.submit(branch.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for l in p['levels']:
                    write_json(root/'status.json',dict(state='running',stage='gpu_and_cpu_production',level=l['level'],cpu_total=432))
                    branch.production(root,p,l)
                for i,f in enumerate(as_completed(tasks),1):
                    f.result();write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=432))
            result=assess(root);write_json(root/'timing.json',dict(seconds=time.perf_counter()-clock,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed'],independent_review_pending=True,total_production_paths=432))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','qualify','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/pulse-time-recipient-response'));parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    result=prepare(args.output) if args.action=='prepare' else assess(args.output) if args.action=='assess' else run(args.output,args.workers,args.action=='qualify')
    print({k:v for k,v in result.items() if k in ('passed','total_production_paths','unique_physical_initial_value_problems')})
