"""Moving confirmation of a one-time matched local state and chemical budget."""
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
from . import phase_carry_delayed_response as delayed
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


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    immediate=Path('output/studies/phase-carry-recipient-amount-response')
    old,_,inputs=shared.accepted_parent(immediate,Path('docs/phase_carry_recipient_amount_response_assessment.json'))
    accepted=Path('output/studies/compensated-recipient-amount-delayed-response')
    recent,_,more=shared.accepted_parent(accepted,Path('docs/compensated_recipient_response_assessment.json'));inputs+=more
    if recent['parent']!=str(immediate.resolve()):raise ValueError('Wrong accepted delayed lineage')
    fp=read('output/studies/phase-carry-frozen-response/protocol.json');verify(fp)
    physical={c['seed']:c for c in fp['contexts'] if c['background']=='sham'}
    contexts=[dict(c) for c in old['contexts'] if c['matching_reference']=='S' and c['physics_label']=='S']
    if len(contexts)!=12:raise ValueError('Expected twelve primary moving backgrounds')
    pairs=[dict(seed=f['seed'],cell=f['cell'],contexts={C:f['contexts'][C+'_S'] for C in ('S','R')}) for f in old['factorials'] if f['matching_reference']=='S']
    root.mkdir(parents=True);jobs=[];configs={};estimate=0.
    for c in contexts:
        c['recipient']=c['targets'][0];c['physical_source']=physical[c['seed']]['original_source']
        if fp['input_sha256'].get(str(Path(c['physical_source']).resolve()))!=digest(c['physical_source']):raise ValueError('Unaccepted actual physical source')
        inputs.extend([Path(c['source']),Path(c['physical_source']),Path(physical[c['seed']]['physical_history'])])
        with np.load(c['source']) as z:initial=z['initial'].copy();ids=z['ids'].copy()
        host=AttributeSimulation.restore(c['physical_source'])
        if host.time!=510. or host.divisions or not np.array_equal(host.ids,ids):raise ValueError('Changed physical source clock/IDs')
        for level,dt in LEVELS:
            copied=AttributeSimulation.restore(c['physical_source']);backend.retime(copied,dt);configs[f'{c["seed"]}_{level}']=asdict(copied.config)
            for stage in ('wait','response'):
                for dose in ([None] if stage=='wait' else [None,*delayed.DOSES]):
                    key=c['key']+'_'+level+'_'+stage+'_'+(dose or 'control')
                    source=c['source'] if stage=='wait' else str(root/'starts'/(c['key']+'_'+(dose or 'control')+'.npz'))
                    jobs.append(dict(key=key,context=c['key'],seed=c['seed'],recipient=c['recipient'],chemical_origin=c['chemical_origin'],level=level,dt=dt,stage=stage,dose=dose,
                        target=c['recipient'] if dose else None,source=source,start=510. if stage=='wait' else 570.,duration=60.))
        # Actual older moving path gives a prospective runtime estimate only.
        hist=read(physical[c['seed']]['physical_history']);result=read(Path(physical[c['seed']]['physical_history']).with_name('result.json'))
        if hist[-1]['elapsed']!=60. or not result['quality_pass']:raise ValueError('Invalid timing-only moving reference')
        estimate+=result['audit']['wall_seconds']*6
    if len(jobs)!=96 or len(pairs)!=6:raise ValueError('Wrong nested design')
    physical_protocol=read(Path(fp['parent'])/'protocol.json')
    inputs+=[Path(fp['parent'])/'protocol.json',Path('output/studies/phase-carry-frozen-response/protocol.json'),Path('docs/moving_recipient_amount_response.md'),Path('docs/compensated_chemistry_validation.json'),Path('archive/study-launches/2026-10-08-moving-recipient-amount-response/tests.log')]
    sources=[Path(__file__),Path(backend.__file__),Path('tests/test_compensated_moving.py'),Path('tests/test_moving_recipient_amount_response.py')]
    p=dict(contexts=contexts,pairs=pairs,jobs=jobs,histories=[7,8,9],independent_histories=3,new_histories=0,levels=[dict(level=l,dt=d) for l,d in LEVELS],
        start=510.,wait=60.,response_duration=60.,pilot_duration=6.,interval=.15,late_window=24.,checkpoint_interval=3.,prefix_duration=.6,
        method=METHOD,device='cuda:0',cuda_validation=physical_protocol['device'],runtime=old['runtime'],accepted_configs=configs,amount_error_max=2e-14,preparation_amount_max=2e-13,
        mechanical_limits=moving.MECHANICAL_LIMITS,response_limits=moving.RESPONSE_CRITERIA,prefix_limits=GPU_LIMITS,effect_min=.01,contrast_min=.1,
        primary_reference='S',initial_physics='S',total_jobs=96,wait_jobs=24,response_jobs=72,reused_jobs=0,
        estimated_moving_seconds=estimate,storage_estimate_bytes=18*2**30,storage_reserve_bytes=2*2**30,
        parent=str(immediate.resolve()),accepted_delayed_parent=str(accepted.resolve()),
        source_sha256={**old['source_sha256'],**recent['source_sha256'],**frozen.pin(sources)},input_sha256=frozen.pin(inputs),
        scope='Sham matching reference and initial sham physical context, six recipients nested in three existing mature histories. Twelve exact accepted t510 chemical preparations, matched recipient pair and species budgets initially only. Live mechanics/polarity/contact transport/dilution. Two timesteps from identical physical starts; response starts from each actual new fine waited endpoint at both timesteps. Both pulse doses, own unpulsed controls. No rematching, no new development, divisions, spatial convergence, autonomy/inheritance or full frozen factorial generalization.')
    if shutil.disk_usage(root).free<p['storage_estimate_bytes']+p['storage_reserve_bytes']:raise OSError('Insufficient local storage for preserved evidence')
    write_json(root/'protocol.json',p);verify(p);write_json(root/'status.json',dict(state='prepared',total=96,estimated_moving_seconds=estimate));return p


def context(p,job):return next(c for c in p['contexts'] if c['key']==job['context'])
def wait_job(p,c,level):return next(j for j in p['jobs'] if (j['context'],j['level'],j['stage'])==(c['key'],level,'wait'))


def materialize(root,p,c):
    root=Path(root);fine=wait_job(p,c,'fine');coarse=wait_job(p,c,'coarse')
    for j in (fine,coarse):current_evidence(root,j,p)
    pair=read(root/'pairs'/c['key']/'wait-full.json')
    if not pair['passed'] or pair['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Waited pair not qualified')
    sim=backend.MovingSimulation.restore(root/fine['key']/'latest_state.npz');ids=sim.ids;values=torch.stack([sim.activator,sim.inhibitor]).cpu().numpy();e=sim.chemical_correction.cpu().numpy();v=sim.geometry[:,0].cpu().numpy()
    opposite=next(x for x in p['contexts'] if (x['seed'],x['recipient'])==(c['seed'],c['recipient']) and x['chemical_origin']!=c['chemical_origin'])
    other=wait_job(p,opposite,'fine');current_evidence(root,other,p)
    with np.load(root/other['key']/'latest_state.npz') as z:other_values=np.array([z['activator'],z['inhibitor']])
    other_history=read(root/other['key']/'history.json');ov=np.asarray(other_history[-1]['volumes']);i=ids.tolist().index(c['recipient'])
    common=delayed.matched_removal(values[0,i]*v[i],other_values[0,i]*ov[i])
    for dose in (None,*delayed.DOSES):
        file=root/'starts'/(c['key']+'_'+(dose or 'control')+'.npz');commit=file.with_suffix('.json')
        changed=values.copy();correction=e.copy();record=dict(factor=1.,initial_log_amplitude=0.,removed_activator_amount=0.,requested_removed_amount=0.)
        if dose:
            changed,record=delayed.pulse(values,v,ids,c['recipient'],dose,common);correction[0,i]=0.
        expected=dict(method=METHOD,context=c['key'],dose=dose,parent_checkpoint=str(root/fine['key']/'latest_state.npz'),parent_checkpoint_sha256=digest(root/fine['key']/'latest_state.npz'),opposite_checkpoint_sha256=digest(root/other['key']/'latest_state.npz'),record=record,protocol_sha256=digest(root/'protocol.json'))
        if file.exists():
            saved=read(commit)
            if {k:v for k,v in saved.items() if k!='paths_sha256'}!=expected or saved['paths_sha256']!=digest(file):raise ValueError('Changed delayed preparation')
            with np.load(file) as z:
                if not np.array_equal(z['initial'],changed) or not np.array_equal(z['correction'],correction):raise ValueError('Changed actual waited concentrations/correction')
        else:
            save_npz(file,initial=changed,correction=correction,ids=ids,volumes=v);write_json(commit,dict(**expected,paths_sha256=digest(file)))
    del sim;torch.cuda.empty_cache()


def initialize(root,job,p):
    c=context(p,job)
    with np.load(job['source']) as z:initial=z['initial'].copy()
    if job['stage']=='wait':sim=backend.MovingSimulation.from_physical(c['physical_source'],job['dt'],initial,p['device'])
    else:
        record=read(Path(job['source']).with_suffix('.json'));cp=Path(record['parent_checkpoint'])
        if digest(cp)!=record['parent_checkpoint_sha256'] or record['paths_sha256']!=digest(job['source']) or record['protocol_sha256']!=digest(Path(root)/'protocol.json'):raise ValueError('Changed response source commitment')
        sim=backend.MovingSimulation.restore(cp,p['device']);backend.retime(sim.host,job['dt']);sim.step_number=sim.host.step_number;sim.config=sim.host.config
        with np.load(job['source']) as z:e=z['correction'].copy()
        backend.validate_correction(initial,e);sim.activator,sim.inhibitor=torch.from_numpy(initial).to(sim.device);sim.chemical_correction=torch.from_numpy(e).to(sim.device)
    if asdict(sim.config)!=p['accepted_configs'][f'{job["seed"]}_{job["level"]}'] or sim.time!=job['start']:raise ValueError('Changed model configuration/clock')
    return sim


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
    if job['stage']=='wait':physical=c['physical_source']
    else:physical=read(Path(job['source']).with_suffix('.json'))['parent_checkpoint']
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
    for dose in ([None] if stage=='wait' else [None,*delayed.DOSES]):
        pair=[next(j for j in p['jobs'] if (j['context'],j['stage'],j['dose'],j['level'])==(c['key'],stage,dose,level)) for level,_ in LEVELS];hs=[];phis=[]
        for j in pair:
            d=root/j['key'];hf=d/('pilot-history.json' if label=='pilot' else 'history.json');pf=d/('pilot-phi.npz' if label=='pilot' else 'latest_state.npz')
            h=read(hf);check_history(h,j,p,horizon)
            with np.load(pf) as z:phi=z['phi'].copy()
            hs.append(h);phis.append(phi);histories[j['level'],dose]=h;files.extend([hf,pf])
        r=moving.mechanical_comparison(*hs,*phis,p,horizon);mechanical.append(dict(dose=dose,**r))
    response=[]
    if stage=='response':
        for dose in delayed.DOSES:
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
    for dose in delayed.DOSES:
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
    root=Path(root);p=read(root/'protocol.json');verify(p);files=[];results=[]
    for job in p['jobs']:
        results.append(current_evidence(root,job,p)[0]);files+=[root/job['key']/f for f in ('result.json','history.json','latest_state.npz','prefix/result.json')]
    for c in p['contexts']:
        for stage in ('wait','response'):
            numerical_pair(root,p,c,stage,p['pilot_duration']);numerical_pair(root,p,c,stage,60.)
    rows=[]
    for pair in p['pairs']:rows+=response_contrast(root,p,pair,60.)['rows']
    outcomes=[dict(history=seed,recipients=2,above_screen_at_both_doses=sum(all(x['fine_above_screen'] for x in rows if (x['seed'],x['cell'])==(seed,pair['cell'])) for pair in p['pairs'] if pair['seed']==seed)) for seed in p['histories']]
    amounts=[]
    for pair in p['pairs']:
        curves={}
        for C in ('S','R'):
            c=next(x for x in p['contexts'] if x['key']==pair['contexts'][C]);job=wait_job(p,c,'fine');h=read(root/job['key']/'history.json');curves[C]=[(np.asarray(x['chemistry'])@np.asarray(x['volumes'])).tolist() for x in h]
        amounts.append(dict(history=pair['seed'],cell=pair['cell'],elapsed=[x['elapsed'] for x in h],whole_species_amounts=curves))
    result=dict(passed=True,independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),method=METHOD,total_paths=96,independent_histories=3,new_histories=0,response_rows=rows,history_outcomes=outcomes,wait_amount_evolution=amounts,evidence_sha256=frozen.pin(files),scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=4,qualify_only=False):
    root=Path(root).resolve()
    if not 1<=workers<=4:raise ValueError('Use one to four two-thread native-prefix workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);moving.require_device(p['cuda_validation']);torch.set_num_threads(1)
        for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
        clock=time.perf_counter()
        try:
            wait=[j for j in p['jobs'] if j['stage']=='wait'];q=qualify_jobs(root,p,wait,workers)
            if not (root/'qualification.json').exists():write_json(root/'qualification.json',dict(**q,protocol_sha256=digest(root/'protocol.json'),actual_gpu=torch.cuda.get_device_name(0),seconds=time.perf_counter()-clock))
            elif read(root/'qualification.json')['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Changed qualification protocol')
            if qualify_only:
                write_json(root/'status.json',dict(state='qualified_wait_starts',qualified_jobs=24,total=96));return q
            for stage in ('wait','response'):
                if stage=='response':
                    for c in p['contexts']:materialize(root,p,c)
                    qualify_jobs(root,p,[j for j in p['jobs'] if j['stage']==stage],workers)
                for horizon in (p['pilot_duration'],60.):
                    for c in p['contexts']:
                        for j in [j for j in p['jobs'] if j['stage']==stage and j['context']==c['key']]:
                            write_json(root/'status.json',dict(state='running',stage=stage,horizon=horizon,active_job=j['key'],completed=sum((root/x['key']/'result.json').exists() for x in p['jobs']),total=96));advance(root,j,p,horizon)
                        numerical_pair(root,p,c,stage,horizon)
                    if stage=='response':
                        for pair in p['pairs']:response_contrast(root,p,pair,horizon)
            result=assess(root);write_json(root/'timing.json',dict(seconds=time.perf_counter()-clock,prefix_cpu_workers=workers));write_json(root/'status.json',dict(state='completed',passed=True,independent_review_pending=True,total=96));return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('action',choices=('prepare','qualify','run','assess'));a.add_argument('--output',type=Path,default=Path('output/studies/moving-recipient-amount-response'));a.add_argument('--workers',type=int,default=4);args=a.parse_args()
    result=prepare(args.output) if args.action=='prepare' else assess(args.output) if args.action=='assess' else run(args.output,args.workers,args.action=='qualify')
    print({k:v for k,v in result.items() if k in ('passed','total_jobs','estimated_moving_seconds','total_paths','history_outcomes')})
