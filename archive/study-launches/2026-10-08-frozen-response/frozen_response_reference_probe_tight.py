import hashlib, json, time, os
for thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[thread_variable] = "1"
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from contextlib import nullcontext

ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
PARENT=ROOT/'outputs/phase-carry-redistribution-response'
OUT=Path('/tmp/frozen_response_reference_probe_tight')

def digest(file):
    h=hashlib.sha256()
    with Path(file).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''): h.update(block)
    return h.hexdigest()

def calculate(item):
    key,state,delta=item
    initial=np.array(state,dtype=float); matrix=np.array(delta,dtype=float);n=matrix.shape[0]
    def fun(t,y):
        a=y[:n];b=y[n:]
        return np.concatenate((np.square(a)/b-a+.02*(matrix@a),2*(np.square(a)-b)+.55*(matrix@b)))
    def jac(t,y):
        a=y[:n];b=y[n:]
        out=np.zeros((2*n,2*n));out[:n,:n]=.02*matrix;out[n:,n:]=.55*matrix
        ii=np.arange(n);out[ii,ii]+=2*a/b-1;out[ii,n+ii]=-np.square(a/b)
        out[n+ii,ii]=4*a;out[n+ii,n+ii]-=2
        return out
    times=np.arange(401,dtype=float)*.15;paths={}; rows={}
    with nullcontext():
        for method in ('DOP853','Radau'):
            started=time.perf_counter();kwargs={'jac':jac} if method=='Radau' else {}
            sol=solve_ivp(fun,(0.,60.),initial.ravel(),method=method,t_eval=times,rtol=1e-13,atol=1e-15,max_step=.15,**kwargs)
            path=sol.y.T.reshape(len(times),2,n)
            if not sol.success or not np.isfinite(path).all() or np.any(path<=0):raise RuntimeError(key+': invalid '+method)
            paths[method]=path;rows[method]=dict(wall_seconds=time.perf_counter()-started,nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,minimum=float(path.min()))
    raw=float(np.max(abs(np.log(paths['DOP853']/paths['Radau']))))
    if raw>1e-8: raise RuntimeError(key+': raw reference gate failed')
    return key,paths,dict(methods=rows,raw_solver_log_error=raw)

def recovery(times,deviation,threshold):
    ids=np.flatnonzero(np.logical_and.accumulate((deviation<=threshold)[::-1])[::-1] & (times<=36.))
    return float(times[ids[0]]) if len(ids) else None

if __name__=='__main__':
    p=json.loads((PARENT/'protocol.json').read_text()); jobs=[]; inputs={}; cases=[]; info={}
    for seed,target in ((8,20),(9,30)):
        for dose in (None,'amount'):
            job=next(j for j in p['jobs'] if j['seed']==seed and j['level']=='fine' and j['background']=='sham' and j['dose']==dose and (dose is None or j['target']==target))
            key=f'seed-{seed}_{"control" if dose is None else "amount-cell-"+str(target)}'
            folder=Path(p['reuse_parent'])/p['reused'][job['key']]['job']['key'] if job['key'] in p['reused'] else PARENT/job['key']
            history=folder/'history.json';row=json.loads(history.read_text())[0]
            cases.append((key,row['chemistry'],row['delta']));info[key]=dict(job=job,masses=row['volumes'],ids=row['ids'],history=str(history.resolve()))
            for file in (history,Path(job['source'])):inputs[str(file.resolve())]=digest(file)
    inputs[str((PARENT/'protocol.json').resolve())]=digest(PARENT/'protocol.json')
    started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=4,mp_context=get_context('spawn')) as pool: returned=list(pool.map(calculate,cases))
    records={key:record for key,paths,record in returned};trajectories={key:paths for key,paths,record in returned}
    times=np.arange(401)*.15; arrays={'times':times};pairs=[]
    for key,paths,_ in returned:
        for method,path in paths.items():arrays[key+'_'+method]=path
    for seed,target in ((8,20),(9,30)):
        ck=f'seed-{seed}_control';pk=f'seed-{seed}_amount-cell-{target}';meta=info[pk];f=meta['job']['factor'];amp=abs(np.log(f));cell=meta['ids'].index(target);m=np.array(meta['masses']);waves={};metrics={}
        for method in ('DOP853','Radau'):
            path=trajectories[pk][method];control=trajectories[ck][method];log=np.log(path/control);waves[method]=log/amp
            network=np.sqrt(np.sum((np.log(path)-np.log(control))**2*m,axis=(1,2))/(2*m.sum()))
            metrics[method]=dict(target_recovery_time=recovery(times,abs(log[:,:,cell]).max(axis=1),.1*amp),network_recovery_time=recovery(times,network,.1*network[0]),target_activator_log_auc_per_log_pulse=float(trapezoid(abs(log[:,0,cell]),times)/amp),target_activator_peak_gain=float(abs(log[:,0,cell]).max()/amp),target_inhibitor_peak_gain=float(abs(log[:,1,cell]).max()/amp))
        error=float(abs(waves['DOP853']-waves['Radau']).max());pairs.append(dict(seed=seed,target=target,factor=f,log_amplitude=float(amp),normalized_network_response_reference_error=error,metrics=metrics))
    arrays_file=OUT.with_suffix('.npz');np.savez_compressed(arrays_file,**arrays)
    report=dict(scope='Independent four-path frozen t510 ODE probe authorized by root; no project outputs or scientific sources changed.',rtol=1e-13,atol=1e-15,max_step=.15,times='0:.15:60, 401 samples',beta=2.,da=.02,db=.55,inputs_sha256=inputs,source_sha256={str(Path(__file__).resolve()):digest(__file__)},arrays_sha256=digest(arrays_file),total_wall_seconds=time.perf_counter()-started,records=records,pairs=pairs,metadata=info)
    OUT.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('inputs_sha256','source_sha256','metadata')},indent=2))
