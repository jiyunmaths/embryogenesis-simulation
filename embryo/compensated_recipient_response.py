"""Isolated compensated chemistry repeat of the delayed amount-matching assay."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import fcntl,json,os,time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
import torch
from . import compensated_chemistry as chemistry
from . import phase_carry_frozen_response as frozen
from . import phase_carry_recipient_state_response_grouped as grouped
from . import phase_carry_recipient_delayed_response as delayed
from . import phase_carry_recipient_amount_response as amount
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

write_json=frozen.write_json
save_npz=frozen.save_npz


def prepare(root,qualification=Path('docs/compensated_chemistry_validation.json'),
            protocol_doc=Path('docs/compensated_recipient_response.md')):
    q=read(qualification)
    if not q['passed'] or q['method']!=chemistry.METHOD:raise ValueError('Accepted compensated backend qualification required')
    for file,h in q['files_sha256'].items():
        if digest(file)!=h:raise ValueError('Changed backend qualification asset: '+file)
    # Derive sources again from the ACCEPTED immediate phase, not from the
    # unaccepted delayed trajectories. This retains the exact t570 starts.
    p=amount.prepare_delayed(root,protocol_doc=protocol_doc)
    previous=Path('output/studies/phase-carry-recipient-amount-delayed-response')
    old=read(previous/'protocol.json');status=read(previous/'status.json')
    if status['state']!='completed_with_unresolved_checks' or status['passed']:
        raise ValueError('Original unresolved study must remain preserved')
    for field in ('criteria','levels','reference_settings','beta','da','db','start','duration','interval'):
        if p[field]!=old[field]:raise ValueError('Repeat changes equations/clock/acceptance: '+field)
    originals=[]
    for now,before in zip(p['contexts'],old['contexts']):
        if now['key']!=before['key']:raise ValueError('Changed context ordering')
        with np.load(now['source']) as a,np.load(before['source']) as b:
            for k in ('delta','volumes','ids','initial'):
                if not np.array_equal(a[k],b[k]):raise ValueError('Changed actual delayed context')
        originals.append(Path(before['source']))
    for now,before in zip(p['trials'],old['trials']):
        if {k:v for k,v in now.items() if k!='source'}!={k:v for k,v in before.items() if k!='source'}:
            raise ValueError('Changed delayed pulse/control preparation')
        with np.load(now['source']) as a,np.load(before['source']) as b:
            if not np.array_equal(a['initial'],b['initial']):raise ValueError('Changed exact delayed initial values')
        originals.append(Path(before['source']))
    p['solver']=dict(method=chemistry.METHOD,dtype='float64',correction_initialization='zero at each identical accepted t570 start',
        correction_retention='every SSP substep and timestep; mandatory checkpoint field',
        observation='rounded positive high concentration arrays; correction is numerical bookkeeping',
        scope='Frozen chemical branch only; historical shared kernels and custom CUDA mechanics/polarity unchanged')
    p['historical_failed_study']=str(previous);p['historical_failed_protocol_sha256']=digest(previous/'protocol.json')
    p['backend_qualification']=str(qualification);p['backend_qualification_sha256']=digest(qualification)
    p['source_sha256'].update(frozen.pin([Path(__file__),Path(chemistry.__file__),Path('tests/test_compensated_chemistry.py'),Path('tests/test_compensated_recipient_response.py')]))
    extra=[qualification,*map(Path,q['files_sha256']),*originals,previous/'protocol.json',previous/'summary.json',previous/'status.json',
        Path('docs/phase_carry_recipient_amount_delayed_response_assessment.json'),Path('docs/phase_carry_recipient_amount_delayed_response_assessment.md')]
    p['input_sha256'].update(frozen.pin(extra))
    p['scope']='Compensated float64 SSP-RK2 repeat of identical delayed amount-matched starts and unchanged scientific/numerical gates. Independent accepted immediate sources, no rematching, no source projection, freely evolving amounts. Three existing mature histories; zero new histories. Preserve original failed assay, shared kernels, ledger and manuscript.'
    write_json(Path(root)/'protocol.json',p);verify(p)
    return p


def _validate_correction(state,correction):
    if correction.dtype!=np.float64 or correction.shape!=state.shape or not np.isfinite(correction).all():
        raise ValueError('Missing/invalid float64 checkpoint correction')
    if np.any(abs(correction)>2*np.spacing(state)):
        raise ValueError('Checkpoint correction exceeds its representational bound')


def cpu_worker(args):
    root,t,level=args;root=Path(root);p=read(root/'protocol.json');ph=digest(root/'protocol.json')
    folder=root/'cpu-verification'/level['level']/t['key'];file=folder/'paths.npz';report=folder/'result.json'
    if report.exists():
        r=read(report)
        if r['protocol_sha256']!=ph or r['trial']!=t or r['level']!=level or r['paths_sha256']!=digest(file) or r['method']!=chemistry.METHOD:
            raise ValueError('Changed CPU compensated verification')
        return r
    _,d,_,initial,_=frozen.load_context(p,t)
    step=chemistry.NumpyStepper(d,level['dt'],p['beta'],p['da'],p['db'])
    if step.count!=frozen.substeps(d,level['dt'],p['beta'],p['da'],p['db']):raise ValueError('Changed positivity substeps')
    u=initial.copy();correction=np.zeros_like(u);history=[u.copy()];clock=time.perf_counter()
    every=frozen._steps(p['interval'],level['dt'])
    for k in range(1,frozen._steps(p['duration'],level['dt'])+1):
        u,correction=step(u,correction)
        if k%every==0:history.append(u.copy())
    _validate_correction(u,correction)
    save_npz(file,times=np.arange(len(history))*p['interval'],trajectory=np.array(history),final_correction=correction)
    r=dict(trial=t,level=level,protocol_sha256=ph,paths_sha256=digest(file),method=chemistry.METHOD,wall_seconds=time.perf_counter()-clock)
    write_json(report,r);return r


def production(root,p,level,stop_after=None):
    """Mandatory concentration+correction checkpoints; optional test interrupt."""
    root=Path(root);ph=digest(root/'protocol.json');folder=root/'production'/level['level']
    file,report,checkpoint=folder/'paths.npz',folder/'result.json',folder/'checkpoint.npz'
    if report.exists():
        r=read(report)
        if r['protocol_sha256']!=ph or r['level']!=level or r['paths_sha256']!=digest(file) or r['checkpoint_sha256']!=digest(checkpoint) or r['method']!=chemistry.METHOD:
            raise ValueError('Changed completed compensated production')
        return r
    initial,delta,mapping,count=grouped.batch_layout(p);device=frozen.runtime_check(p)
    dt=level['dt'];total=frozen._steps(p['duration'],dt);every=frozen._steps(p['interval'],dt);save_every=frozen._steps(p['checkpoint_interval'],dt)
    if stop_after is not None and (not isinstance(stop_after,int) or stop_after<=0 or stop_after>=total or stop_after%save_every):
        raise ValueError('Test interruption must coincide with an interior checkpoint')
    completed,wall,history=0,0.,[frozen.unpack(initial,mapping,count)];state=initial.copy();correction=np.zeros_like(state)
    if checkpoint.exists():
        commit=read(folder/'checkpoint-commit.json')
        if commit['protocol_sha256']!=ph or commit['checkpoint_sha256']!=digest(checkpoint):raise ValueError('Changed/incomplete compensated checkpoint')
        with np.load(checkpoint) as z:
            meta=json.loads(str(z['metadata']));completed=meta['step'];wall=meta['wall_seconds']
            if 'correction' not in z:raise ValueError('Missing mandatory checkpoint correction; no reset')
            state,correction,history=z['state'].copy(),z['correction'].copy(),list(z['history'].copy())
        if (meta['protocol_sha256']!=ph or meta['level']!=level or meta['method']!=chemistry.METHOD
                or completed%save_every or not 0<=completed<=total or state.shape!=initial.shape
                or np.asarray(history).shape!=(completed//every+1,len(mapping),2,count)
                or not np.array_equal(history[0],frozen.unpack(initial,mapping,count))
                or not np.array_equal(history[-1],frozen.unpack(state,mapping,count))):
            raise ValueError('Compensated checkpoint clock/source/history mismatch')
        frozen.valid_state(state);frozen.valid_state(np.array(history));_validate_correction(state,correction)
    u=torch.as_tensor(state,dtype=torch.float64,device=device);e=torch.as_tensor(correction,dtype=torch.float64,device=device)
    matrix=torch.as_tensor(delta,dtype=torch.float64,device=device);stepper=chemistry.TorchStepper(matrix,dt,p['beta'],p['da'],p['db'])
    for c in p['contexts']:
        with np.load(c['source']) as z:
            if frozen.substeps(z['delta'],dt,p['beta'],p['da'],p['db'])!=stepper.count:raise ValueError('Batch changes positivity substeps')
    clock=time.perf_counter();resumed_step=completed
    for k in range(completed,total+1):
        if k%save_every==0:
            current=u.detach().cpu().numpy().copy();carry=e.detach().cpu().numpy().copy();_validate_correction(current,carry)
            meta=dict(protocol_sha256=ph,level=level,step=k,wall_seconds=wall+time.perf_counter()-clock,method=chemistry.METHOD)
            save_npz(checkpoint,state=current,correction=carry,history=np.array(history),metadata=np.array(json.dumps(meta)))
            write_json(folder/'checkpoint-commit.json',dict(protocol_sha256=ph,checkpoint_sha256=digest(checkpoint),level=level,step=k,method=chemistry.METHOD))
            write_json(folder/'status.json',dict(state='running',elapsed=k*dt,horizon=p['duration']))
            print('Compensated production:',level['level'],k*dt,flush=True)
        if stop_after is not None and k==stop_after:return dict(interrupted=True,step=k,checkpoint_sha256=digest(checkpoint))
        if k==total:break
        u,e=stepper(u,e)
        if (k+1)%every==0:history.append(frozen.unpack(u.detach().cpu().numpy(),mapping,count))
    paths=np.array(history);times=np.arange(len(paths))*p['interval']
    if paths.shape!=(len(times),len(mapping),2,count) or len(times)!=frozen._steps(p['duration'],p['interval'])+1:raise ValueError('Incomplete compensated paths')
    save_npz(file,times=times,trajectory=paths,trial_keys=np.array([m['trial'] for m in mapping]))
    r=dict(protocol_sha256=ph,level=level,paths_sha256=digest(file),checkpoint_sha256=digest(checkpoint),actual_paths=len(mapping),
        padding_paths=initial.shape[1]*initial.shape[-1]//count-len(mapping),wall_seconds=wall+time.perf_counter()-clock,
        method=chemistry.METHOD,resumed_from_step=resumed_step,correction_checkpointed=True)
    write_json(report,r);write_json(folder/'status.json',dict(state='completed',elapsed=p['duration']));return r


def run(root,workers=8):
    root=Path(root).resolve()
    if not 1<=workers<=8:raise ValueError('Use one to eight single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        if p['solver']['method']!=chemistry.METHOD:raise ValueError('Wrong compensated solver protocol')
        for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
        clock=time.perf_counter()
        try:
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                tasks=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,f in enumerate(as_completed(tasks),1):
                    if not f.result()['passed']:raise RuntimeError('Independent raw reference gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=len(tasks)))
                if not delayed.reference_qualification(root,p)['passed']:raise RuntimeError('Independent normalized reference gate failed')
                tasks=[pool.submit(cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for l in p['levels']:
                    for f in tasks:
                        if f.done():f.result()
                    write_json(root/'status.json',dict(state='running',stage='production_and_cpu_verification',level=l['level'],cpu_total=len(tasks)))
                    production(root,p,l)
                for i,f in enumerate(as_completed(tasks),1):
                    f.result();write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=len(tasks)))
            result=amount.assess(root);write_json(root/'timing.json',dict(coordinator_wall_seconds=time.perf_counter()-clock,workers=workers,cpu_verification_overlaps_gpu=True))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed'],total_production_paths=p['total_production_paths'],independent_histories=3,new_histories=0,method=chemistry.METHOD))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('action',choices=('prepare','run','assess'))
    a.add_argument('--output',type=Path,default=Path('output/studies/compensated-recipient-amount-delayed-response'));a.add_argument('--workers',type=int,default=8)
    args=a.parse_args();r=prepare(args.output) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else amount.assess(args.output)
    print({k:v for k,v in r.items() if k in ('passed','total_production_paths','unique_initial_value_problems')})
