"""Moving confirmation of challenge-time local chemistry and amount matching."""
import argparse,fcntl,json,os,time
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from multiprocessing import get_context
from pathlib import Path
import shutil
import numpy as np
import torch
from . import compensated_chemistry as chemical
from . import compensated_moving as backend
from . import phase_carry_frozen_response as frozen
from . import phase_carry_shared_context_response as shared
from . import phase_carry_exchange_response as moving
from . import pulse_time_recipient_response as endpoint
import hashlib
import subprocess
import sys
from .attribute_development import AttributeSimulation
from .native_mechanics import NativeSimulation
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .cell_response_moving import observe
from .cell_response_moving_refinement import response_error
from .feedback_long import digest,save_checkpoint,restore_checkpoint
from .neighbor_context import read,validate_graph
from .parameter_robustness import verify
from .resolution import _steps
from .validate_gpu_backend import CRITERIA as GPU_LIMITS,discrepancies

write_json=frozen.write_json
save_npz=frozen.save_npz
LEVELS=moving.LEVELS
METHOD=backend.METHOD



def accepted_challenge(parent,review_file):
    parent=Path(parent).resolve();review_file=Path(review_file)
    p=read(parent/'protocol.json');r=read(review_file);audit_file=Path(r['independent_review'])
    if not r['full_study_accepted'] or r['protocol_sha256']!=digest(parent/'protocol.json') or digest(audit_file)!=r['independent_review_sha256']:
        raise ValueError('Independently accepted amount-matched challenge parent required')
    audit=read(audit_file)
    if not audit['passed'] or not audit['accepted_at_tested_scope'] or not audit['inputs_unchanged'] or audit['protocol_sha256']!=r['protocol_sha256']:
        raise ValueError('Independent challenge-time review not accepted')
    verify(p);inputs=[review_file,audit_file]
    for name in ('protocol.json','summary.json','status.json'):
        f=parent/name
        if audit['evidence_sha256'].get(str(f))!=digest(f):raise ValueError('Changed accepted challenge report')
        inputs.append(f)
    if read(parent/'status.json')['state']!='completed' or not read(parent/'summary.json')['passed']:
        raise ValueError('Completed accepted challenge parent required')
    return p,audit,inputs


def validate_alignment(data,ids,row):
    frozen.validate_context(*(data[k] for k in ('delta','volumes','initial','ids')))
    if row['time']!=570. or row['elapsed']!=60. or not np.array_equal(data['ids'],ids):
        raise ValueError('Wrong accepted physical endpoint clock/IDs')
    for k in ('volumes','delta','ids'):
        if not np.array_equal(data[k],row[k]):raise ValueError('Physical endpoint capacity/operator/IDs mismatch')


def matched_closure(a,b,cell):
    for k in ('delta','volumes','ids'):np.testing.assert_array_equal(a[k],b[k])
    i=a['ids'].tolist().index(cell);others=a['ids']!=cell;v=a['volumes']
    np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i])
    errors={name:float(abs(x/y-1).max()) for name,x,y in (
        ('surrounding',a['initial'][:,others]@v[others],b['initial'][:,others]@v[others]),
        ('global',a['initial']@v,b['initial']@v))}
    if max(errors.values())>2e-13:raise ValueError('Unequal matched initial chemical budgets')
    return errors


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    parent=Path('output/studies/pulse-time-amount-response').resolve()
    old,audit,inputs=accepted_challenge(parent,'docs/pulse_time_amount_response_assessment.json')
    moving_parent=Path('output/studies/moving-recipient-amount-response').resolve()
    mp,physical,more,note=endpoint.accepted_endpoints(moving_parent,'docs/moving_recipient_amount_response_assessment.json');inputs+=more
    local_parent=Path('output/studies/pulse-time-recipient-response').resolve();local=read(local_parent/'protocol.json');verify(local)
    inputs.append(local_parent/'protocol.json')
    if old['parent']!=str(local_parent) or local['parent']!=str(moving_parent):raise ValueError('Wrong accepted endpoint lineage')
    selected=[dict(c) for c in old['contexts'] if (c['matching_reference'],c['budget_reference'],c['physical_origin'])==('S','S','S')]
    if len(selected)!=12:raise ValueError('Expected all six recipients and two chemical backgrounds')
    root.mkdir(parents=True);contexts=[];jobs=[];configs={};pairs=[];states={};estimate=0.;physical_fingerprints={};signatures=set()
    for c in selected:
        original=Path(c['source']);c['accepted_context_source']=str(original)
        if audit['evidence_sha256'].get(str(original.resolve()))!=digest(original):raise ValueError('Changed accepted matched context')
        with np.load(original,allow_pickle=False) as z:data={k:z[k].copy() for k in ('delta','volumes','initial','ids')}
        pc=next(x for x in local['contexts'] if x['key']==c['parent_context']);key=pc['physical_parent_context']
        j=next(j for j in mp['jobs'] if (j['context'],j['stage'],j['level'])==(key,'wait','fine'));folder=moving_parent/j['key']
        row=read(folder/'history.json')[-1];cp=folder/'latest_state.npz';host=AttributeSimulation.restore(cp)
        validate_alignment(data,host.ids,row)
        if host.time!=570. or host.divisions or host.attribute_mode!='direct' or host.config.max_cells!=len(host.ids):raise ValueError('Wrong actual mature physical source')
        # Source is referenced, never copied/projected. Preserve physical carry,
        # RNG, polarity, target volumes and material parameters across reset.
        c.update(recipient=c['targets'][0],physical_source=str(cp),physical_parent_context=key,source=str(root/'contexts'/(c['key']+'.npz')))
        save_npz(c['source'],**data);states[c['key']]=data;contexts.append(c)
        inputs.extend([original,cp,folder/'history.json',folder/'result.json',Path(c['source'])])
        if str(cp) not in physical_fingerprints:
            h=hashlib.sha256()
            with np.load(cp,allow_pickle=False) as z:
                for name in ('phi','phase_carry','rounding','polarity','target','ids','parents','due','fate','metadata','attribute_mode','precision_arm'):
                    x=z[name];h.update(name.encode());h.update(x.dtype.str.encode());h.update(str(x.shape).encode());h.update(x.tobytes())
            physical_fingerprints[str(cp)]=h.hexdigest()
        for level,dt in LEVELS:
            copied=AttributeSimulation.restore(cp);backend.retime(copied,dt);configs[c['key']+'_'+level]=asdict(copied.config)
            for dose in (None,'fractional'):
                target=None if dose is None else c['recipient']
                trial=next(t for t in old['trials'] if t['context']==c['key'] and t['target']==target)
                initial=data['initial'].copy();record={k:trial[k] for k in ('factor','initial_log_amplitude','requested_removed_amount','removed_activator_amount','initial_target_activator_amount','initial_target_volume','pulse_amount_relative_error')}
                if dose:initial[0,data['ids'].tolist().index(target)]*=.9
                with np.load(trial['source'],allow_pickle=False) as z:
                    if not np.array_equal(initial,z['initial']):raise ValueError('Changed exact accepted challenge/pulse')
                if audit['evidence_sha256'].get(str(Path(trial['source']).resolve()))!=digest(trial['source']):raise ValueError('Changed accepted pulse')
                source=root/'starts'/(c['key']+'_'+(dose or 'control')+'.npz')
                if not source.exists():
                    save_npz(source,**{**data,'initial':initial,'correction':np.zeros_like(initial)})
                    write_json(source.with_suffix('.json'),dict(record=record,parent_trial=trial['key'],parent_source=trial['source'],parent_source_sha256=digest(trial['source']),physical_source=str(cp),physical_source_sha256=digest(cp),paths_sha256=digest(source),correction_initialization='Zero only at explicit rounded-state t570 intervention'))
                    inputs.extend([Path(trial['source']),source,source.with_suffix('.json')]);signatures.add((physical_fingerprints[str(cp)],initial.tobytes()))
                jobs.append(dict(key=c['key']+'_'+level+'_'+(dose or 'control'),context=c['key'],seed=c['seed'],recipient=c['recipient'],chemical_origin=c['chemical_origin'],level=level,dt=dt,stage='response',dose=dose,target=target,source=str(source),start=570.,duration=60.))
                # Matched legacy response uses same dt and horizon: timing only.
                timing=next(x for x in mp['jobs'] if (x['context'],x['stage'],x['dose'],x['level'])==(key,'response',dose,level))
                estimate+=read(moving_parent/timing['key']/'result.json')['audit']['wall_seconds']
    for seed,cell in sorted({(c['seed'],c['recipient']) for c in contexts}):
        pair=dict(seed=seed,cell=cell,contexts={C:next(c['key'] for c in contexts if (c['seed'],c['recipient'],c['chemical_origin'])==(seed,cell,C)) for C in ('S','R')})
        matched_closure(*(states[pair['contexts'][C]] for C in ('S','R')),cell)
        if len({next(c['physical_source'] for c in contexts if c['key']==pair['contexts'][C]) for C in ('S','R')})!=1:raise ValueError('Paired physical sources differ')
        pairs.append(pair)
    launch=Path('archive/study-launches/2026-10-09-moving-pulse-time-amount-response')
    sources=[Path(__file__),Path('tests/test_moving_pulse_time_amount_response.py'),launch/'audit_inputs.py',launch/'audit_qualification.py']
    inputs+=[Path('docs/moving_pulse_time_amount_response.md'),launch/'tests.log']
    p=dict(contexts=contexts,pairs=pairs,jobs=jobs,histories=[7,8,9],independent_histories=3,new_histories=0,
        levels=[dict(level=l,dt=d) for l,d in LEVELS],start=570.,wait=0.,response_duration=60.,pilot_duration=6.,interval=.15,late_window=24.,checkpoint_interval=3.,prefix_duration=.6,
        method=METHOD,device='cuda:0',cuda_validation=mp['cuda_validation'],runtime=old['runtime'],accepted_configs=configs,
        amount_error_max=2e-14,preparation_amount_max=2e-13,mechanical_limits=moving.MECHANICAL_LIMITS,response_limits=moving.RESPONSE_CRITERIA,prefix_limits=GPU_LIMITS,
        effect_min=.01,contrast_min=.1,primary_reference='S',budget_reference='S',initial_physics='S',total_jobs=48,response_jobs=48,wait_jobs=0,reused_jobs=0,
        logical_starts=24,distinct_initial_value_problems=len(signatures),estimated_moving_seconds=estimate,storage_estimate_bytes=12*2**30,storage_reserve_bytes=2*2**30,
        parent=str(parent),physical_parent=str(moving_parent),local_parent=str(local_parent),parent_historical_export_note=note,
        preparation_audit=str(launch/'preparation-audit.json'),qualification_audit=str(launch/'qualification-audit.json'),qualification_auditor=str(launch/'audit_qualification.py'),
        source_sha256={**mp['source_sha256'],**old['source_sha256'],**frozen.pin(sources)},input_sha256={**old['input_sha256'],**frozen.pin(inputs)},
        scope='Primary S local reference, S surrounding-amount reference and S actual endpoint physics; all six recipients nested in three existing mature histories. Immediate t570 reset/pulse with same physical fields/carry and exact accepted initial concentrations. Initial matched local chemistry and both surrounding/global species amounts only; subsequent chemistry, mechanics, polarity, contacts, volumes, conservative transport and dilution co-evolve freely to t630. Own moving controls and two timesteps, all-prefix/pilot gates before full continuation. No settling, new development, nearest-neighbor specificity, inherited/autonomous identity, spatial convergence or extension to reciprocal references. Chemical correction zero at explicit reset, mandatory thereafter. Moving splitting remains first order; no second-order moving claim.')
    if (len(jobs),len(pairs))!=(48,6):raise ValueError('Wrong complete recipient design')
    if shutil.disk_usage(root).free<p['storage_estimate_bytes']+p['storage_reserve_bytes']:raise OSError('Insufficient local storage')
    write_json(root/'protocol.json',p);verify(p);write_json(root/'status.json',dict(state='prepared',total=48,estimated_moving_seconds=estimate));return p


def initialize(root,job,p):
    c=context(p,job)
    with np.load(job['source'],allow_pickle=False) as z:initial=z['initial'].copy();e=z['correction'].copy()
    if np.any(e!=0):raise ValueError('Explicit first reset requires zero initial chemical correction')
    record=read(Path(job['source']).with_suffix('.json'))
    if digest(c['physical_source'])!=record['physical_source_sha256'] or digest(job['source'])!=record['paths_sha256']:raise ValueError('Changed source commitment')
    sim=backend.MovingSimulation.from_physical(c['physical_source'],job['dt'],initial,p['device'])
    if asdict(sim.config)!=p['accepted_configs'][c['key']+'_'+job['level']] or sim.time!=570. or sim.host.divisions:raise ValueError('Changed configuration/clock')
    return sim


def context(p,job):return next(c for c in p['contexts'] if c['key']==job['context'])

def check_history(history,job,p,horizon):
    times=np.arange(_steps(horizon,p['interval'])+1)*p['interval']
    if len(history)!=len(times) or not np.allclose([x['elapsed'] for x in history],times,rtol=0,atol=1e-9):raise ValueError('Incomplete observation clock')
    with np.load(job['source']) as z:initial=z['initial'];ids=z['ids']
    if not np.array_equal(history[0]['chemistry'],initial):raise ValueError('Changed exact intervention/pulse')
    for row in history:
        u=np.asarray(row['chemistry'])
        if not np.array_equal(row['ids'],ids) or abs(row['time']-job['start']-row['elapsed'])>1e-9 or not np.isfinite(u).all() or np.any(u<=0):raise ValueError('Invalid alignment/chemistry')
        validate_graph(row['delta'],row['volumes'])
    return history

def current_evidence(root,job,p):
    folder=Path(root)/job['key'];r=read(folder/'result.json');h=read(folder/'history.json');ph=digest(Path(root)/'protocol.json')
    if r['job']!=job or r['protocol_sha256']!=ph or not r['quality_pass'] or r['source_sha256']!=digest(job['source']) or r['history_sha256']!=digest(folder/'history.json') or r['checkpoint_sha256']!=digest(folder/'latest_state.npz'):raise ValueError('Changed completed moving evidence')
    check_history(h,job,p,job['duration'])
    with np.load(folder/'latest_state.npz') as z:
        meta=json.loads(str(z['long_experiment']));backend.validate_correction(np.array([z['activator'],z['inhibitor']]),z['chemical_correction'])
        if str(z['chemical_method'])!=METHOD or meta['protocol_hash']!=ph or meta['job']!=job or meta['history']!=h or meta['audit']!=r['audit'] or not np.array_equal(h[-1]['chemistry'],[z['activator'],z['inhibitor']]):raise ValueError('Checkpoint/history mismatch')
    return r,h

def cpu_prefix(args):
    root,job,threads=args;root=Path(root);p=read(root/'protocol.json');folder=root/job['key']/'prefix';folder.mkdir(parents=True,exist_ok=True)
    c=context(p,job)
    physical=c['physical_source']
    cpu=AttributeSimulation.restore(physical);backend.retime(cpu,job['dt'])
    with np.load(job['source']) as z:cpu.activator,cpu.inhibitor=z['initial'].copy()
    cpu.__class__=NativeSimulation;cpu.native_threads=threads;history=[]
    steps=_steps(p['prefix_duration'],job['dt']);every=_steps(p['interval'],job['dt'])
    for k in range(steps+1):
        if k%every==0:history.append(observe(cpu,k*job['dt']))
        if k<steps:cpu.step()
    write_json(folder/'cpu-history.json',history);save_npz(folder/'cpu-phi.npz',phi=cpu.phi)
    return dict(job=job,history_sha256=digest(folder/'cpu-history.json'),phi_sha256=digest(folder/'cpu-phi.npz'))

def gpu_prefix(root,job,p):
    folder=Path(root)/job['key']/'prefix';folder.mkdir(parents=True,exist_ok=True);sim=initialize(root,job,p);error=0.;amount_error=0.
    x=torch.stack([sim.activator,sim.inhibitor]).cpu().numpy();e=sim.chemical_correction.cpu().numpy()
    for _ in range(8):
        d=sim.matrices()[1].cpu().numpy();old=sim.geometry[:,0].cpu().numpy().copy()
        step=chemical.NumpyStepper(d,job['dt'],sim.config.signal_beta,sim.config.signal_da,sim.config.signal_dh);x,e=step(x,e)
        sim.step();x,e=backend.numpy_scale(x,e,old/sim.geometry[:,0].cpu().numpy())
        actual=torch.stack([sim.activator,sim.inhibitor]).cpu().numpy();error=max(error,float(abs(np.log(actual/x)).max()));amount_error=max(amount_error,sim.audit()['dilution_amount_error'])
    if error>1e-11 or amount_error>p['amount_error_max']:raise RuntimeError('Independent evolving-geometry chemical check failed')
    checkpoint=folder/'restart.npz';sim.checkpoint(checkpoint);restart=backend.MovingSimulation.restore(checkpoint,p['device'])
    for _ in range(4):
        sim.step();restart.step()
        for name in ('phi','activator','inhibitor','polarity','geometry','phase_carry','rounding','chemical_correction'):
            if not torch.equal(getattr(sim,name),getattr(restart,name)):raise RuntimeError('Restart differs: '+name)
    del sim,restart;torch.cuda.empty_cache();sim=initialize(root,job,p);history=[];every=_steps(p['interval'],job['dt'])
    for k in range(_steps(p['prefix_duration'],job['dt'])+1):
        q=sim.audit();amount_error=max(amount_error,q['dilution_amount_error'])
        if amount_error>p['amount_error_max']:raise RuntimeError('Prefix dilution closure failed')
        if k%every==0:history.append(sim.observe(k*job['dt']))
        if k<_steps(p['prefix_duration'],job['dt']):sim.step()
    write_json(folder/'gpu-history.json',history);save_npz(folder/'gpu-phi.npz',phi=sim.phi.cpu().numpy())
    del sim;torch.cuda.empty_cache();return dict(independent_chemical_steps=8,cpu_chemical_log_max=error,amount_accounting_max=amount_error,restart_exact_steps=4)

def qualify_jobs(root,p,jobs,workers=4):
    pending=[]
    for job in jobs:
        file=Path(root)/job['key']/'prefix/result.json'
        if file.exists():
            r=read(file)
            if not r['passed'] or r['protocol_sha256']!=digest(Path(root)/'protocol.json') or r['job']!=job:raise ValueError('Changed prefix acceptance')
            for f,h in r['evidence_sha256'].items():
                if digest(f)!=h:raise ValueError('Changed prefix evidence')
        else:pending.append(job)
    with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        futures={j['key']:pool.submit(cpu_prefix,(root,j,2)) for j in pending}
        for number,job in enumerate(pending,1):
            write_json(Path(root)/'status.json',dict(state='qualifying',stage=job['stage'],completed=number-1,total=len(pending),active_job=job['key']))
            chemical_check=gpu_prefix(root,job,p);futures[job['key']].result();folder=Path(root)/job['key']/'prefix';a=read(folder/'gpu-history.json');b=read(folder/'cpu-history.json')
            errors={key:max(discrepancies(x,y)[key] for x,y in zip(a,b)) for key in ('chemical_log_max','polarity_abs_max','relative_axis_max','relative_volume_max','relative_transport_max')}
            with np.load(folder/'gpu-phi.npz') as x,np.load(folder/'cpu-phi.npz') as y:phi_error=float(abs(x['phi']-y['phi']).max())
            if any(max(x['boundary_occupancy'],y['boundary_occupancy'])>=.01 for x,y in zip(a,b)):raise RuntimeError('Prefix boundary failed')
            passed=all(value<=GPU_LIMITS[key] for key,value in errors.items()) and phi_error<=GPU_LIMITS['final_phi_abs_max']
            files=[folder/name for name in ('gpu-history.json','cpu-history.json','gpu-phi.npz','cpu-phi.npz','restart.npz')]
            r=dict(passed=bool(passed),job=job,protocol_sha256=digest(Path(root)/'protocol.json'),source_sha256=digest(job['source']),method=METHOD,errors=errors,final_phi_abs_max=phi_error,**chemical_check,evidence_sha256=frozen.pin(files),scope='Actual-start native CPU/GPU prefix; independent live-geometry chemical stepping and exact checkpoint correction restart. Not full-horizon native equivalence.')
            write_json(folder/'result.json',r)
            if not passed:raise RuntimeError('Unchanged strict native prefix gate failed: '+job['key'])
            print('Moving context gate passed',number,len(pending),job['key'],flush=True)
    return dict(passed=True,jobs=len(jobs),method=METHOD)

def advance(root,job,p,horizon):
    root=Path(root);folder=root/job['key'];folder.mkdir(parents=True,exist_ok=True);cp=folder/'latest_state.npz';ph=digest(root/'protocol.json')
    gate=read(folder/'prefix/result.json')
    if not gate['passed'] or gate['source_sha256']!=digest(job['source']):raise ValueError('Actual start not qualified')
    if (folder/'result.json').exists():return [x for x in current_evidence(root,job,p)[1] if x['elapsed']<=horizon+1e-9]
    if cp.exists():
        host,audit,history=restore_checkpoint(cp,ph,job);check_history(history,job,p,host.time-job['start']);sim=backend.MovingSimulation.restore(cp,p['device'])
        if not np.array_equal(history[-1]['chemistry'],[host.activator,host.inhibitor]):raise ValueError('Restart differs from history')
    else:
        sim=initialize(root,job,p);history=[];audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,boundary_max=0.,dilution_error_max=0.,wall_seconds=0.)
    start=_steps(job['start'],job['dt']);stop=_steps(job['start']+horizon,job['dt']);every=_steps(p['interval'],job['dt']);save_every=_steps(p['checkpoint_interval'],job['dt'])
    if sim.step_number>stop:return [x for x in history if x['elapsed']<=horizon+1e-9]
    if not start<=sim.step_number<=stop:raise ValueError('Invalid continuation clock')
    clock=time.perf_counter();base=audit['wall_seconds']
    while sim.step_number<=stop:
        q=sim.audit();audit['max_volume_error']=max(audit['max_volume_error'],q['max_volume_error']);audit['max_clipping']=max(audit['max_clipping'],q['max_clipping']);audit['min_radius']=min(audit['min_radius'],q['min_radius']);audit['dilution_error_max']=max(audit['dilution_error_max'],q['dilution_amount_error'])
        if audit['dilution_error_max']>p['amount_error_max']:raise RuntimeError('Dilution gate failed')
        k=sim.step_number-start
        if k%every==0:
            elapsed=k*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed)>1e-10:
                row=sim.observe(elapsed);row['precision']=sim.precision_diagnostics();history.append(row)
            audit['boundary_max']=max(audit['boundary_max'],history[-1]['boundary_occupancy'])
            if audit['boundary_max']>=.01:raise RuntimeError('Boundary gate failed')
            audit['wall_seconds']=base+time.perf_counter()-clock
            write_json(folder/'status.json',dict(state='running',elapsed=elapsed,target_horizon=horizon,audit=audit))
            if k%save_every==0 or sim.step_number==stop:
                if shutil.disk_usage(root).free<p['storage_reserve_bytes']:raise OSError('Storage reserve reached; retain committed checkpoints')
                save_checkpoint(sim,cp,audit,history,ph,job);write_json(folder/'history.json',history)
                print(job['key'],elapsed,horizon,round(audit['wall_seconds'],1),flush=True)
        if sim.step_number==stop:break
        sim.step()
    check_history(history,job,p,horizon)
    if horizon==p['pilot_duration']:
        save_npz(folder/'pilot-phi.npz',phi=sim.phi.cpu().numpy());write_json(folder/'pilot-history.json',history)
    if horizon==job['duration']:
        if not moving.quality_pass(audit,p):raise RuntimeError('Final quality failed')
        write_json(folder/'result.json',dict(job=job,protocol_sha256=ph,source_sha256=digest(job['source']),history_sha256=digest(folder/'history.json'),checkpoint_sha256=digest(cp),method=METHOD,quality_pass=True,audit=audit,**moving.contrast_summary(history,p,horizon)))
        write_json(folder/'status.json',dict(state='completed',elapsed=horizon))
    del sim;torch.cuda.empty_cache();return history

def numerical_pair(root,p,c,stage,horizon):
    root=Path(root);label='pilot' if horizon==p['pilot_duration'] else 'full';folder=root/'pairs'/c['key'];folder.mkdir(parents=True,exist_ok=True);file=folder/(stage+'-'+label+'.json')
    if file.exists():
        r=read(file)
        if r['protocol_sha256']!=digest(root/'protocol.json') or not r['passed']:raise ValueError('Changed/failed paired report')
        for f,h in r['evidence_sha256'].items():
            if digest(f)!=h:raise ValueError('Changed paired evidence')
        return r
    mechanical=[];histories={};files=[]
    for dose in (None,'fractional'):
        pair=[next(j for j in p['jobs'] if (j['context'],j['stage'],j['dose'],j['level'])==(c['key'],stage,dose,level)) for level,_ in LEVELS];hs=[];phis=[]
        for j in pair:
            d=root/j['key'];hf=d/('pilot-history.json' if label=='pilot' else 'history.json');pf=d/('pilot-phi.npz' if label=='pilot' else 'latest_state.npz')
            h=read(hf);check_history(h,j,p,horizon)
            with np.load(pf) as z:phi=z['phi'].copy()
            hs.append(h);phis.append(phi);histories[j['level'],dose]=h;files.extend([hf,pf])
        r=moving.mechanical_comparison(*hs,*phis,p,horizon);mechanical.append(dict(dose=dose,**r))
    response=[]
    if stage=='response':
        for dose in ('fractional',):
            job=next(j for j in p['jobs'] if (j['context'],j['stage'],j['dose'])==(c['key'],stage,dose));record=read(Path(job['source']).with_suffix('.json'))['record'];factor=record['factor'];metrics=[];u=[];controls=[]
            for level,_ in LEVELS:
                h=histories[level,dose];ctrl=histories[level,None];path=np.array([x['chemistry'] for x in h]);base=np.array([x['chemistry'] for x in ctrl]);times=np.array([x['elapsed'] for x in h]);v=np.asarray(h[0]['volumes']);i=h[0]['ids'].index(c['recipient'])
                metrics.append(response_metrics(path,base,times,v,i,factor));u.append(path);controls.append(base)
            err=response_error(u[0],controls[0],u[1],controls[1],factor);auc=abs(metrics[1]['target_activator_log_auc_per_log_pulse']/metrics[0]['target_activator_log_auc_per_log_pulse']-1);recover={}
            for name in ('target_recovery_time','network_recovery_time'):
                a,b=[x[name] for x in metrics];recover[name]=0. if a is None and b is None else None if a is None or b is None else abs(a-b)
            response.append(dict(dose=dose,normalized_response_error=err,relative_auc_error=auc,recovery_time_errors=recover,passed=bool(err<=.01 and auc<=.02 and all(x is not None and x<=.3+1e-9 for x in recover.values()))))
    result=dict(passed=bool(all(x['passed'] for x in mechanical+response)),stage=stage,label=label,context=c['key'],protocol_sha256=digest(root/'protocol.json'),mechanical=mechanical,response=response,evidence_sha256=frozen.pin(files))
    write_json(file,result)
    if not result['passed']:raise RuntimeError('Unchanged timestep gate failed: '+str(file))
    return result

def response_contrast(root,p,pair,horizon):
    """Signed differences normalized by each background's own pulse."""
    root=Path(root);label='pilot' if horizon==p['pilot_duration'] else 'full';rows=[];evidence=[]
    for dose in ('fractional',):
        waves={};metrics={}
        for level,_ in LEVELS:
            ws=[];ms=[]
            for C in ('S','R'):
                key=pair['contexts'][C];job=next(j for j in p['jobs'] if (j['context'],j['stage'],j['dose'],j['level'])==(key,'response',dose,level));control=next(j for j in p['jobs'] if (j['context'],j['stage'],j['dose'],j['level'])==(key,'response',None,level))
                files=[root/j['key']/('pilot-history.json' if label=='pilot' else 'history.json') for j in (job,control)];h,b=map(read,files);evidence+=files
                u=np.array([x['chemistry'] for x in h]);base=np.array([x['chemistry'] for x in b]);times=np.array([x['elapsed'] for x in h]);i=h[0]['ids'].index(pair['cell']);record=read(Path(job['source']).with_suffix('.json'))['record'];amp=abs(np.log(record['factor']))
                ws.append(np.log(u[:,:,i]/base[:,:,i])/amp);ms.append(response_metrics(u,base,times,np.array(h[0]['volumes']),i,record['factor']))
            waves[level]=ws[1]-ws[0];metrics[level]=ms
        effects={l:waveform_distance(w,np.zeros_like(w),times) for l,w in waves.items()};error=float(abs(waves['coarse']-waves['fine']).max());same=(effects['coarse']>.01)==(effects['fine']>.01)
        rows.append(dict(seed=pair['seed'],cell=pair['cell'],dose=dose,effects=effects,direct_signed_waveform_error=error,same_decision=bool(same),fine_above_screen=bool(effects['fine']>.01),metrics=metrics,passed=bool(same and error<=.01 and abs(effects['fine']-effects['coarse'])<=.01)))
        save_npz(root/'contrasts'/f'history-{pair["seed"]}_cell-{pair["cell"]}_{dose}_{label}.npz',times=times,coarse=waves['coarse'],fine=waves['fine'])
    result=dict(passed=all(x['passed'] for x in rows),rows=rows,protocol_sha256=digest(root/'protocol.json'),evidence_sha256=frozen.pin(evidence))
    write_json(root/'contrasts'/f'history-{pair["seed"]}_cell-{pair["cell"]}_{label}.json',result)
    if not result['passed']:raise RuntimeError('Signed response contrast fails refinement')
    return result


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);files=[]
    for job in p['jobs']:
        current_evidence(root,job,p);files+=[root/job['key']/f for f in ('result.json','history.json','latest_state.npz','prefix/result.json')]
    for c in p['contexts']:
        for horizon in (p['pilot_duration'],60.):numerical_pair(root,p,c,'response',horizon)
    rows=[]
    for pair in p['pairs']:rows+=response_contrast(root,p,pair,60.)['rows']
    outcomes=[dict(history=seed,recipients=2,above_screen=sum(x['fine_above_screen'] for x in rows if x['seed']==seed)) for seed in p['histories']]
    amounts=[]
    for pair in p['pairs']:
        curves={};surrounding={};transport={}
        for C in ('S','R'):
            job=next(j for j in p['jobs'] if (j['context'],j['level'],j['dose'])==(pair['contexts'][C],'fine',None));h=read(root/job['key']/'history.json');i=h[0]['ids'].index(pair['cell']);others=np.arange(len(h[0]['ids']))!=i
            curves[C]=[(np.asarray(x['chemistry'])@np.asarray(x['volumes'])).tolist() for x in h]
            surrounding[C]=[(np.asarray(x['chemistry'])[:,others]@np.asarray(x['volumes'])[others]).tolist() for x in h]
            transport[C]=[float(np.linalg.norm(x['delta'])) for x in h]
        amounts.append(dict(history=pair['seed'],cell=pair['cell'],elapsed=[x['elapsed'] for x in h],whole_species_amounts=curves,surrounding_species_amounts=surrounding,operator_frobenius_norm=transport))
    result=dict(passed=True,independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),method=METHOD,total_paths=48,logical_starts=24,distinct_initial_value_problems=p['distinct_initial_value_problems'],independent_histories=3,new_histories=0,response_rows=rows,history_outcomes=outcomes,amount_evolution=amounts,evidence_sha256=frozen.pin(files),scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=4,qualify_only=False):
    root=Path(root).resolve()
    if not 1<=workers<=4:raise ValueError('Use one to four two-thread native-prefix workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);moving.require_device(p['cuda_validation']);torch.set_num_threads(1)
        for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
        pre=read(p['preparation_audit'])
        if not pre['passed'] or pre['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Independent preparation audit required')
        clock=time.perf_counter()
        try:
            q=qualify_jobs(root,p,p['jobs'],workers)
            write_json(root/'qualification.json',dict(**q,protocol_sha256=digest(root/'protocol.json'),actual_gpu=torch.cuda.get_device_name(0),seconds=time.perf_counter()-clock))
            subprocess.run([sys.executable,p['qualification_auditor']],check=True)
            qa=read(p['qualification_audit'])
            if not qa['passed'] or qa['protocol_sha256']!=digest(root/'protocol.json') or qa['actual_start_gates']!=48:raise ValueError('Independent complete prefix audit required')
            if qualify_only:
                write_json(root/'status.json',dict(state='qualified',qualified_jobs=48,total=48));return q
            for horizon in (p['pilot_duration'],60.):
                for c in p['contexts']:
                    for j in [j for j in p['jobs'] if j['context']==c['key']]:
                        write_json(root/'status.json',dict(state='running',stage='pilot' if horizon==p['pilot_duration'] else 'full',horizon=horizon,active_job=j['key'],completed=sum((root/x['key']/'result.json').exists() for x in p['jobs']),total=48));advance(root,j,p,horizon)
                    numerical_pair(root,p,c,'response',horizon)
                for pair in p['pairs']:response_contrast(root,p,pair,horizon)
                if horizon==p['pilot_duration']:
                    write_json(root/'pilot-gate.json',dict(passed=True,protocol_sha256=digest(root/'protocol.json'),contexts=12,contrast_pairs=6,evidence_sha256=frozen.pin([root/'pairs'/c['key']/'response-pilot.json' for c in p['contexts']]+[root/'contrasts'/f'history-{pair["seed"]}_cell-{pair["cell"]}_pilot.json' for pair in p['pairs']])))
            result=assess(root);write_json(root/'timing.json',dict(seconds=time.perf_counter()-clock,prefix_cpu_workers=workers));write_json(root/'status.json',dict(state='completed',passed=True,independent_review_pending=True,total=48));return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('action',choices=('prepare','qualify','run','assess'));a.add_argument('--output',type=Path,default=Path('output/studies/moving-pulse-time-amount-response'));a.add_argument('--workers',type=int,default=4);args=a.parse_args()
    result=prepare(args.output) if args.action=='prepare' else assess(args.output) if args.action=='assess' else run(args.output,args.workers,args.action=='qualify')
    print({k:v for k,v in result.items() if k in ('passed','total_jobs','estimated_moving_seconds','total_paths','history_outcomes')})
