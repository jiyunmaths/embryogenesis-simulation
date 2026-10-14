"""Qualify isolated compensated chemistry on exact existing delayed inputs."""
import argparse,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from multiprocessing import get_context
from pathlib import Path
import numpy as np
import torch
from . import compensated_chemistry as chemistry
from . import compensated_recipient_response as branch
from . import phase_carry_frozen_response as frozen
from . import phase_carry_recipient_state_response_grouped as grouped
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify

OLD=Path('output/studies/phase-carry-recipient-amount-delayed-response')
OLD_REVIEW=Path('archive/study-reviews/2026-10-08-recipient-amount-delayed-response/independent_review.json')
write_json=frozen.write_json
save_npz=frozen.save_npz


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    old=read(OLD/'protocol.json');verify(old);audit=read(OLD_REVIEW)
    if audit['passed'] or not audit['inputs_unchanged'] or not audit['input_audit']['passed'] or len([x for x in audit['orders'] if not x['passed']])!=1:
        raise ValueError('Preserved diagnosed original and qualified inputs required')
    for file,h in audit['inputs_sha256'].items():
        if digest(file)!=h:raise ValueError('Changed original numerical evidence')
    keys={c['key'] for c in old['contexts'] if c['matching_reference']=='S' and c['physics_label']=='S'}
    contexts=[c for c in old['contexts'] if c['key'] in keys];trials=[t for t in old['trials'] if t['context'] in keys]
    if len(contexts)!=12 or len(trials)!=36 or sum(t['target'] is None for t in trials)!=12:raise ValueError('Wrong numerical profile selection')
    files=[OLD/'protocol.json',OLD_REVIEW,Path('docs/compensated_chemistry.md'),Path('docs/compensated_recipient_response.md'),
        Path('archive/study-launches/2026-10-08-compensated-chemistry/kernel-tests.log'),Path('archive/study-launches/2026-10-08-compensated-chemistry/checkpoint-tests.log')]
    sources=[Path(__file__),Path(chemistry.__file__),Path(branch.__file__),Path('tests/test_compensated_chemistry.py'),Path('tests/test_compensated_recipient_response.py')]
    # Exact original sources and qualified independent references are read in place.
    for t in old['trials']:files += [Path(t['source']),OLD/'references'/t['key']/'paths.npz',OLD/'references'/t['key']/'result.json']
    files += [Path(c['source']) for c in old['contexts']]
    p={**old,'contexts':contexts,'trials':trials,'factorials':[],
        'source_sha256':{**old['source_sha256'],**frozen.pin(sources)},'input_sha256':{**old['input_sha256'],**frozen.pin(files)},
        'solver':dict(method=chemistry.METHOD),
        'scope':'Numerical qualification only,36 selected existing cases;72 qualified reference trajectories reused in place; all144 starts receive short GPU/CPU/restart checks. No scientific outcomes or new histories.'}
    initial,d,m,n=grouped.batch_layout(p);p['batching']=dict(physical_blocks=d.shape[0]//n,width=initial.shape[-1],padding_paths=initial.shape[1]*initial.shape[-1]//n-len(m))
    root.mkdir();write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',long_cases=36,all_start_prefix_cases=144))
    prefix=root/'prefix';prefix.mkdir();pref={**old,'duration':.6,'checkpoint_interval':.3,'solver':dict(method=chemistry.METHOD),
        'source_sha256':p['source_sha256'],'input_sha256':p['input_sha256'],'scope':'All144 actual delayed starts, numerical prefix/restart qualification only'}
    write_json(prefix/'protocol.json',pref);return p


def prefix(root):
    root=Path(root);folder=root/'prefix';p=read(folder/'protocol.json');verify(p);device=frozen.runtime_check(p);torch.set_num_threads(1)
    initial,d,m,n=grouped.batch_layout(p);rows=[]
    for level in p['levels']:
        dt=level['dt'];total=frozen._steps(.6,dt);every=frozen._steps(.15,dt)
        branch.production(folder,p,level,stop_after=frozen._steps(.3,dt));branch.production(folder,p,level)
        with np.load(folder/'production'/level['level']/'paths.npz') as z:gpu=z['trajectory'].copy()
        with np.load(folder/'production'/level['level']/'checkpoint.npz') as z:final=z['state'].copy();carry=z['correction'].copy()
        matrix=torch.as_tensor(d,dtype=torch.float64,device=device);step=chemistry.TorchStepper(matrix,dt,p['beta'],p['da'],p['db'])
        u=torch.as_tensor(initial.copy(),dtype=torch.float64,device=device);e=torch.zeros_like(u);shadow=[frozen.unpack(initial,m,n)]
        for k in range(1,total+1):
            u,e=step(u,e)
            if k%every==0:shadow.append(frozen.unpack(u.detach().cpu().numpy(),m,n))
        shadow=np.array(shadow);state=u.detach().cpu().numpy();ec=e.detach().cpu().numpy()
        exact=bool(np.array_equal(shadow,gpu) and np.array_equal(state,final) and np.array_equal(ec,carry))
        if not exact:raise ValueError('GPU checkpoint restart changes state/correction/history')
        frozen.save_npz(folder/'uninterrupted'/level['level']/'paths.npz',trajectory=shadow,state=state,correction=ec)
        cpu=[]
        for t in p['trials']:
            _,delta,_,values,_=frozen.load_context(p,t);step=chemistry.NumpyStepper(delta,dt,p['beta'],p['da'],p['db']);cc=np.zeros_like(values);history=[values.copy()]
            for k in range(1,total+1):
                values,cc=step(values,cc)
                if k%every==0:history.append(values.copy())
            cpu.append(history)
        cpu=np.array(cpu).transpose(1,0,2,3);frozen.save_npz(folder/'cpu-prefix'/level['level']/'paths.npz',trajectory=cpu)
        raw=float(np.max(abs(np.log(gpu/cpu))));normalized=[]
        for i,t in enumerate(p['trials']):
            if t['target'] is None:continue
            j=next(j for j,c in enumerate(p['trials']) if c['context']==t['context'] and c['target'] is None)
            normalized.append(float(np.max(abs(np.log(gpu[:,i]/gpu[:,j])-np.log(cpu[:,i]/cpu[:,j])))/t['initial_log_amplitude']))
        rows.append(dict(level=level,logical_cases=144,raw_log_error=raw,normalized_response_error=max(normalized),gpu_restart_bitwise_identical=exact,
            passed=bool(raw<=p['criteria']['cpu_torch_log_max'] and max(normalized)<=p['criteria']['normalized_response_max'])))
    report=dict(passed=bool(all(x['passed'] for x in rows)),method=chemistry.METHOD,levels=rows,actual_gpu=torch.cuda.get_device_name(device),protocol_sha256=digest(folder/'protocol.json'))
    write_json(folder/'summary.json',report);return report


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);old=read(OLD/'protocol.json');times=np.arange(401)*.15
    paths={};noise={};evidence={};accuracy=[];orders=[];rawpairs=[]
    for t in p['trials']:
        f=OLD/'references'/t['key']/'paths.npz';r=read(f.with_name('result.json'))
        if r['paths_sha256']!=digest(f) or r['protocol_sha256']!=digest(OLD/'protocol.json') or r['trial']!=t or not r['passed']:
            raise ValueError('Unqualified independent validation reference')
        with np.load(f) as z:
            if not np.array_equal(z['times'],times):raise ValueError('Changed reference clock')
            for label in ('dop853','radau'):
                paths[label,t['key']]=z[label].copy();frozen.valid_state(z[label])
                if not np.array_equal(z[label][0],frozen.load_context(p,t)[3]):raise ValueError('Reference start differs')
        noise[t['key']]=float(abs(np.log(paths['dop853',t['key']]/paths['radau',t['key']])).max());evidence.update(frozen.pin([f,f.with_name('result.json')]))
    for level in p['levels']:
        f=root/'production'/level['level']/'paths.npz';r=read(f.with_name('result.json'))
        if r['paths_sha256']!=digest(f) or r['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Changed numerical GPU paths')
        with np.load(f) as z:
            if not np.array_equal(z['times'],times) or list(z['trial_keys'])!=[t['key'] for t in p['trials']]:raise ValueError('Numerical GPU clock/keys differ')
            for key,values in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3)):paths[level['level'],str(key)]=values.copy();frozen.valid_state(values)
        evidence.update(frozen.pin([f,f.with_name('result.json'),f.with_name('checkpoint.npz'),f.with_name('checkpoint-commit.json')]))
        for t in p['trials']:
            f=root/'cpu-verification'/level['level']/t['key']/'paths.npz';r=read(f.with_name('result.json'))
            if r['paths_sha256']!=digest(f) or r['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Changed numerical CPU paths')
            with np.load(f) as z:
                if not np.array_equal(z['times'],times):raise ValueError('Numerical CPU clock differs')
                paths['cpu-'+level['level'],t['key']]=z['trajectory'].copy();frozen.valid_state(z['trajectory'])
            evidence.update(frozen.pin([f,f.with_name('result.json')]))
    for t in p['trials']:
        errors=[]
        for level in p['levels']:
            label=level['level'];error=max(float(np.max(abs(np.log(paths[label,t['key']]/paths[ref,t['key']])))) for ref in ('dop853','radau'))
            cpu=float(np.max(abs(np.log(paths[label,t['key']]/paths['cpu-'+label,t['key']]))));errors.append(error)
            accuracy.append(dict(trial=t['key'],level=label,raw_reference_log_error=error,cpu_torch_log_error=cpu,passed=bool(error<=.01 and cpu<=1e-8)))
        for a,b in [('coarse','fine'),('fine','diagnostic')]:
            error=float(np.max(abs(np.log(paths[a,t['key']]/paths[b,t['key']]))));rawpairs.append(dict(trial=t['key'],first=a,second=b,raw_log_error=error,passed=bool(error<=.01)))
        resolved=errors[-1]>100*noise[t['key']];values=np.log2(np.array(errors[:-1])/errors[1:]).tolist() if resolved and all(x>0 for x in errors) else None
        orders.append(dict(trial=t['key'],errors=errors,reference_noise=noise[t['key']],order_resolved=bool(resolved),orders=values,passed=bool(not resolved or values is not None and all(1.7<=x<=2.3 for x in values))))
    summaries={}
    for label in ('dop853','radau','coarse','fine','diagnostic'):
        metrics={}
        for t in p['trials']:
            if t['target'] is None:continue
            control=next(c for c in p['trials'] if c['context']==t['context'] and c['target'] is None)
            _,_,v,_,ids=frozen.load_context(p,t);i=ids.tolist().index(t['target'])
            metrics[t['key']]=frozen.response_metrics(paths[label,t['key']],paths[label,control['key']],times,v,i,t['factor'])
        summaries[label]=dict(metrics=metrics,effects=[])
    specs=[('dop853','radau'),('coarse','fine'),('fine','diagnostic')]+[(l['level'],r) for l in p['levels'] for r in ('dop853','radau')]
    comparisons=[frozen.compare_responses(paths,a,b,p,times,summaries,1e-8 if a=='dop853' else None) for a,b in specs]
    diagnosed=next(x for x in read(OLD/'summary.json')['orders'] if not x['passed'])['trial'];row=next(x for x in orders if x['trial']==diagnosed)
    control=paths['fine',diagnosed];drift=float(np.max(abs(np.log(control/control[0]))));diagnosis_pass=bool(not row['order_resolved'] and row['errors'][-1]<=100*row['reference_noise'] and drift>0)
    prefix_report=read(root/'prefix/summary.json')
    result=dict(passed=bool(prefix_report['passed'] and diagnosis_pass and all(x['passed'] for x in accuracy+orders+rawpairs+comparisons) and all(x<=1e-8 for x in noise.values())),
        method=chemistry.METHOD,protocol_sha256=digest(root/'protocol.json'),counts=dict(long_cases=36,controls=12,pulses=24,gpu_paths=108,cpu_paths=108,reused_reference_paths=72,actual_start_prefix_cases=144,independent_histories=3,new_histories=0),
        accuracy=accuracy,orders=orders,raw_pairs=rawpairs,response_comparisons=comparisons,diagnosed_control=dict(order_record=row,stored_max_log_relaxation=drift,passed=diagnosis_pass),prefix=prefix_report,evidence_sha256=evidence,criteria=p['criteria'],scope=p['scope'])
    write_json(root/'summary.json',result);return result


def run(root,workers=8):
    root=Path(root);p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1);clock=time.perf_counter()
    if not prefix(root)['passed']:raise RuntimeError('Actual GPU prefix/restart qualification failed')
    with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        tasks=[pool.submit(branch.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
        for l in p['levels']:
            write_json(root/'status.json',dict(state='running',stage='long_gpu_cpu_validation',level=l['level'],cpu_total=len(tasks)))
            branch.production(root,p,l)
        for i,f in enumerate(as_completed(tasks),1):
            f.result();write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=len(tasks)))
    result=assess(root);write_json(root/'timing.json',dict(seconds=time.perf_counter()-clock,workers=workers));write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',passed=result['passed']))
    return result


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('action',choices=('prepare','run','assess'));a.add_argument('--output',type=Path,default=Path('output/studies/compensated-chemistry-validation'));a.add_argument('--workers',type=int,default=8);args=a.parse_args()
    result=prepare(args.output) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k in ('passed','counts')})
