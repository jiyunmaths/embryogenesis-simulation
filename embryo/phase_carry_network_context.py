"""One-time surrounding-chemistry reset on qualified moving carry states.

Reuse exact fine untouched-neighbor controls; repeat them from the identical
fine physical start at the coarse timestep. Keep scientific kernels unchanged.
"""
import argparse
from dataclasses import asdict
import fcntl
from pathlib import Path
import time

import numpy as np
import torch

from . import phase_carry_exchange_response as moving
from .attribute_persistence import distance
from .cell_exchange_response_moving import local_response
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .cell_response_moving_refinement import response_error
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify
from .resolution import write_json

HISTORIES = (7, 8, 9)
LEVELS = (('coarse', .00375), ('fine', .001875))
EFFECT_MIN = .01
CHALLENGE_MIN = .05


def reset_context(initial, reference, masses, ids, recipient=None):
    """Preserve one recipient, reset all other cells; record external amounts."""
    x, ref, m, ids = [np.asarray(a) for a in (initial, reference, masses, ids)]
    if (x.shape != (2, len(m)) or ref.shape != x.shape or m.ndim != 1 or
            ids.shape != m.shape or len(set(ids.tolist())) != len(m) or
            any(not np.isfinite(a).all() for a in (x, ref, m)) or
            any(np.any(a <= 0) for a in (x, ref, m)) or
            (recipient is not None and recipient not in ids)):
        raise ValueError('Invalid reset chemistry, masses, IDs or recipient')
    result = np.array(x, dtype=np.float64, copy=True)
    selected = np.zeros(len(m), bool)
    if recipient is not None:
        i = ids.tolist().index(recipient); selected[:] = True; selected[i] = False
        result[:, selected] = ref[:, selected]
        if not np.array_equal(result[:, i], x[:, i]):
            raise ValueError('Recipient chemistry changed')
    amounts = (result-x)@m
    rms = float(distance(result[:, selected], x[:, selected], m[selected])) if selected.any() else 0.
    return result, dict(recipient=recipient, reset_ids=ids[selected].astype(int).tolist(),
        surrounding_log_rms=rms, added_amounts=amounts.tolist(),
        initial_total_amounts=(x@m).tolist(), reset_total_amounts=(result@m).tolist(),
        recipient_exact=True, sham_exact=bool(np.array_equal(result, x)))


def source_equivalence(first, second):
    """No reuse based only on matching parameters or approximately equal fields."""
    a, b = moving.payload(first), moving.payload(second)
    if a.keys() != b.keys():
        raise ValueError('Reused source schema differs')
    for key in a:
        if a[key].dtype != b[key].dtype or not np.array_equal(a[key], b[key]):
            raise ValueError('Reused source differs: '+key)


def job_design(root, selections):
    root = Path(root); contexts, jobs = [], []
    for seed in HISTORIES:
        targets = selections[str(seed)]['selected_ids']
        for recipient in (None, *targets):
            name = 'sham' if recipient is None else f'reset-cell-{recipient}'
            context = dict(key=f'seed-{seed}_{name}', seed=seed, recipient=recipient,
                targets=targets if recipient is None else [recipient])
            contexts.append(context)
            for level, dt in LEVELS:
                for target in (None, *context['targets']):
                    suffix = 'control' if target is None else f'cell-{target}_neg10'
                    key = f'seed-{seed}_{level}_{name}_{suffix}'
                    jobs.append(dict(key=key, seed=seed, level=level, dt=dt,
                        stage='network_context', background=context['key'], recipient=recipient,
                        target=target, factor=1. if target is None else .9,
                        source=str(root/key/'source.npz'), start=450., duration=60.,
                        point=dict(ratio=27.5, chi=.35)))
    return contexts, jobs


def current_evidence(root, job, p):
    if job['key'] not in p['reused']:
        return moving.current_evidence(root, job, p)
    item = p['reused'][job['key']]; parent = Path(p['parent'])
    if digest(parent/'protocol.json') != item['protocol_sha256']:
        raise ValueError('Changed reused protocol')
    old = read(parent/'protocol.json')
    source_equivalence(job['source'], item['job']['source'])
    result, history = moving.current_evidence(parent, item['job'], old)
    moving.validate_history(history, job, p, job['duration'])
    for path, sha in item['evidence_sha256'].items():
        if digest(path) != sha:
            raise ValueError('Changed reused context evidence')
    return result, history


def advance(root, job, p, horizon):
    if job['key'] in p['reused']:
        _, h = current_evidence(root, job, p)
        return [r for r in h if r['elapsed'] <= horizon+1e-9]
    return moving.advance(root, job, p, horizon)


def selected_history(root, job, p, horizon):
    if job['key'] in p['reused']:
        _, history = current_evidence(root, job, p)
    else:
        history = read(Path(root)/job['key']/'history.json')
    result = [r for r in history if r['elapsed'] <= horizon+1e-9]
    moving.validate_history(result, job, p, horizon)
    return result


def reference_histories(p, seed, horizon):
    refs = p['references'][str(seed)]; result = {}
    for item in refs:
        h = read(item['history']); history = [r for r in h if r['elapsed'] <= horizon+1e-9]
        if len(history) != round(horizon/p['interval'])+1:
            raise ValueError('Incomplete fixed fine reference waveform')
        result[item['target']] = history
    return result


def context_summary(histories, reference, targets, p, horizon):
    """Each target reset has its own control; sham control is shared."""
    t = np.arange(round(horizon/p['interval'])+1)*p['interval']
    sham = histories[('sham', None)]; ids = sham[0]['ids']
    ref_control = np.array([r['chemistry'] for r in reference[None]])
    rows = []
    for target in targets:
        i = ids.index(target); donor = next(cell for cell in targets if cell != target)
        refwaves = {cell:local_response(np.array([r['chemistry'] for r in reference[cell]]),
            ref_control, ids.index(cell), .9) for cell in (target, donor)}
        reference_sep = waveform_distance(refwaves[target], refwaves[donor], t)
        contexts = {}; waves = {}
        sham_state = np.array([r['chemistry'] for r in sham])
        late = t >= horizon-min(horizon,p['late_window'])-1e-9
        for name in ('sham', f'reset-cell-{target}'):
            h = histories[name, None]; control = np.array([r['chemistry'] for r in h])
            path = np.array([r['chemistry'] for r in histories[name, target]])
            weights = np.array(h[0]['volumes']); wave = local_response(path,control,i,.9)
            waves[name] = wave
            dd, ds = [waveform_distance(wave, refwaves[cell], t) for cell in (target,donor)]
            state_shift = np.sqrt(np.mean((np.log(control[:,:,i])-np.log(sham_state[:,:,i]))**2,axis=1))
            contexts[name] = dict(**moving.contrast_summary(h,p,horizon),
                target_state_shift_late_max=float(state_shift[late].max()),
                target_state_shift_final=float(state_shift[-1]),
                metrics=dict(injected_activator_amount=float((.9-1)*control[0,0,i]*weights[i]),
                    **response_metrics(path,control,t,weights,i,.9)),
                destination_distance=dd, donor_distance=ds,
                informative_reference=bool(reference_sep>p['effect_min']),
                nearest_reference=('donor' if ds<dd else 'destination') if reference_sep>p['effect_min'] else 'unresolved')
        shifted = contexts[f'reset-cell-{target}']; shift = waveform_distance(waves['sham'],waves[f'reset-cell-{target}'],t)
        rows.append(dict(cell=target,donor=donor,reference_separation=reference_sep,contexts=contexts,
            response_shift_rms=shift,response_shift_detected=bool(shift>p['effect_min']),
            late_target_state_shift_detected=bool(shifted['target_state_shift_late_max']>p['effect_min'])))
    return rows


def response_comparison(coarse, fine, reference, targets, p, horizon):
    summaries=[context_summary(h,reference,targets,p,horizon) for h in (coarse,fine)]
    rows=[]
    for name,target in [(name,t) for t in targets for name in ('sham',f'reset-cell-{t}')]:
        controls=[np.array([r['chemistry'] for r in h[name,None]]) for h in (coarse,fine)]
        paths=[np.array([r['chemistry'] for r in h[name,target]]) for h in (coarse,fine)]
        error=response_error(paths[0],controls[0],paths[1],controls[1],.9)
        a,b=[next(r for r in s if r['cell']==target)['contexts'][name] for s in summaries]
        metric='target_activator_log_auc_per_log_pulse'
        auc_error=abs(b['metrics'][metric]/a['metrics'][metric]-1)
        recovery={}
        for key in ('target_recovery_time','network_recovery_time'):
            x,y=a['metrics'][key],b['metrics'][key]
            recovery[key]=0. if x is None and y is None else None if None in (x,y) else abs(x-y)
        limits=p['response_limits']; passed=error<=limits['max_normalized_response_error'] and auc_error<=limits['max_relative_auc_error'] and all(v is not None and v<=limits['max_recovery_time_error']+1e-9 for v in recovery.values())
        rows.append(dict(context=name,cell=target,normalized_response_error=error,relative_auc_error=auc_error,
            recovery_time_errors=recovery,passed=bool(passed)))
    same=True; state_error=shift_error=0.
    for a,b in zip(*summaries):
        same &= all(a[k]==b[k] for k in ('response_shift_detected','late_target_state_shift_detected'))
        shift_error=max(shift_error,abs(a['response_shift_rms']-b['response_shift_rms']))
        for name in a['contexts']:
            x,y=a['contexts'][name],b['contexts'][name]
            same &= x['nearest_reference']==y['nearest_reference'] and x['informative_reference']==y['informative_reference']
            state_error=max(state_error,abs(x['target_state_shift_late_max']-y['target_state_shift_late_max']))
    passed=all(r['passed'] for r in rows) and same and state_error<=.01 and shift_error<=.01
    return dict(passed=bool(passed),rows=rows,same_effect_and_reference_decisions=bool(same),
        target_state_shift_error_max=state_error,response_shift_error_max=shift_error,
        coarse=summaries[0],fine=summaries[1])


def field_snapshot(root, job, p, horizon):
    if job['key'] in p['reused']:
        item=p['reused'][job['key']]
        key='pilot_phi' if horizon==p['pilot_duration'] else 'full_phi'
        with np.load(item[key]) as z:return z['phi'].copy()
    out=moving.payload(Path(root)/job['key']/'latest_state.npz')
    import json
    if abs(json.loads(str(out['metadata']))['time']-job['start']-horizon)>1e-9:
        raise ValueError('Field snapshot clock mismatch')
    return out['phi']


def pair_report(root, p, seed, horizon):
    root=Path(root);label='pilot' if horizon==p['pilot_duration'] else 'full'
    folder=root/'pairs'/f'seed-{seed}';folder.mkdir(parents=True,exist_ok=True)
    file=folder/f'{label}-refinement.json';old=read(file) if file.exists() else None
    ph=digest(root/'protocol.json')
    if old:
        if old['protocol_sha256']!=ph:raise ValueError('Changed paired protocol')
        for path,sha in old['evidence_sha256'].items():
            if digest(path)!=sha:raise ValueError('Changed paired evidence')
    jobs=[j for j in p['jobs'] if j['seed']==seed];histories={};mechanical=[];files=[]
    for level in ('coarse','fine'):
        histories[level]={}
        for job in [j for j in jobs if j['level']==level]:
            name='sham' if job['recipient'] is None else f'reset-cell-{job["recipient"]}'
            history=selected_history(root,job,p,horizon);histories[level][name,job['target']]=history
            stem=folder/f'{label}_{job["key"]}';hf=stem.with_suffix('.json');pf=stem.with_suffix('.npz')
            if old:
                if read(hf)!=history:raise ValueError('Changed paired history')
            else:
                write_json(hf,history);np.savez_compressed(pf,phi=field_snapshot(root,job,p,horizon))
            files.extend((hf,pf))
    for name,target in histories['coarse']:
        phis={}
        for level in ('coarse','fine'):
            job=next(j for j in jobs if j['level']==level and j['target']==target and ('sham' if j['recipient'] is None else f'reset-cell-{j["recipient"]}')==name)
            with np.load(folder/f'{label}_{job["key"]}.npz') as z:phis[level]=z['phi'].copy()
        mechanical.append(dict(context=name,target=target,**moving.mechanical_comparison(
            histories['coarse'][name,target],histories['fine'][name,target],phis['coarse'],phis['fine'],p,horizon)))
    refs=reference_histories(p,seed,horizon)
    responses=response_comparison(histories['coarse'],histories['fine'],refs,p['selections'][str(seed)]['selected_ids'],p,horizon)
    result=dict(seed=seed,horizon=horizon,protocol_sha256=ph,
        passed=bool(all(r['passed'] for r in mechanical) and responses['passed']),mechanical=mechanical,
        response=responses,evidence_sha256={str(f.resolve()):digest(f) for f in files})
    if old:
        if old!=result:raise ValueError('Paired decisions do not reproduce')
    else:write_json(file,result)
    return result


def prepare(root,parent=Path('outputs/phase-carry-exchange-response')):
    root,parent=Path(root).resolve(),Path(parent).resolve()
    if root.exists():raise FileExistsError(root)
    old=read(parent/'protocol.json');verify(old)
    assessment_file=Path('docs/phase_carry_exchange_response_assessment.json').resolve();review=read(assessment_file)
    if (read(parent/'status.json')['state']!='completed' or not read(parent/'summary.json')['passed'] or
            not review['full_study_accepted'] or review['protocol_sha256']!=digest(parent/'protocol.json') or
            review['original_summary_sha256']!=digest(parent/'summary.json') or
            review['full_review_sha256']!=digest(review['full_review'])):
        raise ValueError('Requires the completed independently reviewed carry exchange study')
    root.mkdir(parents=True);inputs=[assessment_file,Path(review['full_review']),parent/'protocol.json',parent/'status.json',parent/'summary.json',Path('docs/phase_carry_network_context.md').resolve()]
    contexts,jobs=job_design(root,old['selections']);accepted={};reused={};references={};interventions={};estimate=0.
    p=dict(parent=str(parent),device=old['device'],selections=old['selections'],contexts=contexts,jobs=jobs,
        accepted_configs=accepted,reused=reused,references=references,histories=list(HISTORIES),independent_histories=3,new_histories=0,
        start=450.,duration=60.,pilot_duration=6.,interval=.15,late_window=24.,checkpoint_interval=3.,
        prefix_duration=.6,prefix_native_threads=4,amount_error_max=2e-14,contrast_min=.1,
        mechanical_limits=moving.MECHANICAL_LIMITS,response_limits=moving.RESPONSE_CRITERIA,
        effect_min=EFFECT_MIN,challenge_min=CHALLENGE_MIN,total_jobs=42,new_jobs=33,reused_jobs=9,
        design='Same fine t=450 developed exchanged physical state at both timesteps. Preserve one recipient and reset all other cells to same-age untouched chemistry once. Immediate -10% activator pulses with own moving controls; exact fine sham control/pulses reused. Live chemistry, mechanics, polarity, conservative transport and dilution.',
        scope='Three existing mature histories; 33 new and nine exact reused paths. Reset changes amounts, all other cells, chemistry-mediated mechanics and geometry; not neighbor-only topology, autonomous/committed/inherited identity, fresh carry zygote development, spatial convergence or general ledger promotion.')
    try:
        for seed in HISTORIES:
            original={}
            for bg in ('fresh_exchange','unexchanged'):
                j=next(j for j in old['jobs'] if (j['seed'],j['level'],j['stage'],j['background'])==(seed,'fine','formation',bg))
                result,h=moving.current_evidence(parent,j,old);original[bg]=(j,result,h[-1])
                inputs.extend(parent/j['key']/f for f in ('result.json','history.json','latest_state.npz'))
            source=parent/original['fresh_exchange'][0]['key']/'latest_state.npz'
            last=original['fresh_exchange'][2];initial=np.array(last['chemistry']);reference=np.array(original['unexchanged'][2]['chemistry']);m=np.array(last['volumes']);ids=last['ids']
            references[str(seed)]=[]
            for target in (None,*old['selections'][str(seed)]['selected_ids']):
                j=next(j for j in old['jobs'] if (j['seed'],j['level'],j['stage'],j['background'],j['target'])==(seed,'fine','response','unexchanged',target))
                moving.current_evidence(parent,j,old)
                refs=[parent/j['key']/f for f in ('result.json','history.json','latest_state.npz')];inputs.extend(refs)
                references[str(seed)].append(dict(target=target,history=str(refs[1]),result=str(refs[0]),job=j))
            for ctx in [c for c in contexts if c['seed']==seed]:
                changed,record=reset_context(initial,reference,m,ids,ctx['recipient'])
                if ctx['recipient'] is not None and record['surrounding_log_rms']<=CHALLENGE_MIN:raise ValueError('Unresolved surrounding challenge')
                record.update(ids=ids,masses=m.tolist(),initial_chemistry=initial.tolist(),reset_chemistry=changed.tolist(),reference_chemistry=reference.tolist())
                interventions[ctx['key']]=record
                for j in [j for j in jobs if j['background']==ctx['key']]:
                    values=changed.copy();pulse_amount=0.
                    if j['target'] is not None:
                        i=ids.index(j['target']);pulse_amount=float((j['factor']-1)*values[0,i]*m[i]);values[0,i]*=j['factor']
                    host=moving.prepared_checkpoint(source,j['source'],j['dt'],values)
                    accepted[f'{seed}_{j["level"]}']=asdict(host.config)
                    j['reset_added_amounts']=record['added_amounts'];j['pulse_added_activator_amount']=pulse_amount
                    inputs.append(Path(j['source']))
                    if ctx['recipient'] is None and j['level']=='fine':
                        prior=next(a for a in old['jobs'] if (a['seed'],a['level'],a['stage'],a['background'],a['target'])==(seed,'fine','response','fresh_exchange',j['target']))
                        source_equivalence(j['source'],prior['source']);moving.current_evidence(parent,prior,old)
                        evidence={};files=[parent/prior['key']/f for f in ('result.json','history.json','latest_state.npz','prefix/result.json')]
                        gate=read(files[-1]);files.extend(Path(f) for f in gate['evidence_sha256'])
                        suffix='control' if j['target'] is None else f'cell-{j["target"]}'
                        snapshots={label:str(parent/'pairs'/f'seed-{seed}_response'/f'{label}_fresh_exchange_{suffix}_fine.npz') for label in ('pilot','full')}
                        files.extend(Path(f) for f in snapshots.values());inputs.extend(files)
                        evidence={str(f.resolve()):digest(f) for f in files}
                        reused[j['key']]=dict(job=prior,protocol_sha256=digest(parent/'protocol.json'),evidence_sha256=evidence,
                            pilot_phi=snapshots['pilot'],full_phi=snapshots['full'])
                    else:
                        prior=next(a for a in old['jobs'] if (a['seed'],a['level'],a['stage'],a['background'],a['target'])==(seed,j['level'],'response','fresh_exchange',j['target']))
                        estimate+=read(parent/prior['key']/'result.json')['audit']['wall_seconds']
            print('Prepared exact carry states and chemical-context accounting:',seed,flush=True)
        p['interventions']=interventions;p['estimated_moving_seconds']=estimate
        files=[Path(__file__),Path('tests/test_phase_carry_network_context.py')]
        p['source_sha256']={**old['source_sha256'],**{str(f.resolve()):digest(f) for f in files}}
        p['input_sha256']={**old['input_sha256'],**{str(f.resolve()):digest(f) for f in inputs}}
        write_json(root/'protocol.json',p)
        write_json(root/'status.json',dict(state='prepared',completed=0,total=42,new_jobs=33,reused_jobs=9,
            prefix_contexts_completed=0,estimated_moving_seconds=estimate,independent_histories=3,new_histories=0))
        return p
    except Exception as error:
        write_json(root/'status.json',dict(state='preparation_failed',error=str(error)));raise


def schedule_history(root,p,seed,done):
    jobs=[j for j in p['jobs'] if j['seed']==seed]
    for horizon in (p['pilot_duration'],p['duration']):
        for j in jobs:
            write_json(Path(root)/'status.json',dict(state='running',stage='pilot' if horizon==p['pilot_duration'] else 'full',
                current_job=j['key'],target_horizon=horizon,completed=done,total=p['total_jobs'],
                new_jobs=p['new_jobs'],reused_jobs=p['reused_jobs'],independent_histories=3,new_histories=0))
            advance(root,j,p,horizon)
            if horizon==p['duration']:done+=1
        report=pair_report(root,p,seed,horizon)
        if not report['passed']:return done,False
    return done,True


def preflight(root):
    """Full-size reset/pulse contexts in both directions and both timesteps."""
    root=Path(root).resolve();p=read(root/'protocol.json');verify(p)
    moving.require_device(p['device']);torch.set_num_threads(1)
    jobs=[j for j in p['jobs'] if j['seed']==HISTORIES[0] and j['recipient'] is not None and j['target'] is not None]
    records=[]
    for job in jobs:
        write_json(root/'status.json',dict(state='preflighting',current_job=job['key'],
            prefix_contexts_completed=len(records),prefix_contexts_total=len(jobs),completed=0,total=42))
        records.append(moving.context_gate(root,job,p))
    result=dict(passed=bool(all(r['passed'] for r in records)),protocol_sha256=digest(root/'protocol.json'),
        actual_reset_pulse_contexts=len(records),contexts=records,
        evidence_sha256={str((root/j['key']/'prefix/result.json').resolve()):digest(root/j['key']/'prefix/result.json') for j in jobs},
        scope='Four full 72^3 reset/pulse contexts, both recipients and both timesteps in history 7. All remaining new contexts still require their own gate before pilots; all paired pilots gate long runs.')
    write_json(root/'preflight.json',result)
    write_json(root/'status.json',dict(state='prepared_gpu_preflight_passed',completed=0,total=42,
        new_jobs=33,reused_jobs=9,actual_preflight_contexts=len(records),independent_histories=3,new_histories=0,
        estimated_moving_seconds=p['estimated_moving_seconds']))
    return result


def assess(root,p=None):
    root=Path(root);p=read(root/'protocol.json') if p is None else p;verify(p)
    results=[];missing=[];by_history=[]
    for j in p['jobs']:
        if j['key'] not in p['reused'] and not (root/j['key']/'result.json').exists():missing.append(j['key']);continue
        result,_=current_evidence(root,j,p)
        results.append(dict(job=j,reused=j['key'] in p['reused'],evidence=result))
    for seed in p['histories']:
        reports={}
        if all(j['key'] not in missing for j in p['jobs'] if j['seed']==seed):
            reports={label:pair_report(root,p,seed,horizon) for label,horizon in (('pilot',p['pilot_duration']),('full',p['duration']))}
        by_history.append(dict(seed=seed,qualified=bool(len(reports)==2 and all(r['passed'] for r in reports.values())),reports=reports))
    result=dict(protocol_sha256=digest(root/'protocol.json'),completed=len(results),total=p['total_jobs'],
        passed=bool(not missing and all(r['qualified'] for r in by_history)),missing=missing,
        independent_histories=3,new_histories=0,new_jobs=33,reused_jobs=9,results=results,by_history=by_history,scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);moving.require_device(p['device']);torch.set_num_threads(1)
        done=0;failed=[]
        try:
            for seed in p['histories']:
                done,passed=schedule_history(root,p,seed,done)
                if not passed:failed.append(dict(seed=seed,reason='prespecified_timestep_gate_failed'));break
            result=assess(root,p)
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'],completed=done,total=42,new_jobs=33,reused_jobs=9,failed_histories=failed))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',completed=done,total=42,error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','preflight','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('outputs/phase-carry-network-context'))
    args=parser.parse_args()
    actions=dict(prepare=prepare,preflight=preflight,run=run,assess=assess)
    result=actions[args.action](args.output)
    print({k:v for k,v in result.items() if k in ('passed','completed','total','estimated_moving_seconds')})
