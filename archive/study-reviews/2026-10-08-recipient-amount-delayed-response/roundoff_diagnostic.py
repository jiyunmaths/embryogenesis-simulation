"""Diagnose an unchanged failed order gate; no new acceptance rule or science."""
import hashlib,json,time
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.linalg import eigvals
HERE=Path(__file__).resolve().parent;ROOT=Path('output/studies/phase-carry-recipient-amount-delayed-response')
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
read=lambda f:json.loads(Path(f).read_text())
p=read(ROOT/'protocol.json');summary=read(ROOT/'summary.json');bad=[r for r in summary['orders'] if not r['passed']]
assert len(bad)==1;failed=bad[0];trial=next(t for t in p['trials'] if t['key']==failed['trial']);assert trial['target'] is None
context=next(c for c in p['contexts'] if c['key']==trial['context']);source=Path(context['source'])
with np.load(source) as z:u=z['initial'].copy();d=z['delta'].copy()
mp.mp.dps=80;n=u.shape[1];D=mp.matrix([[mp.mpf(float(x)) for x in row] for row in d]);Da,Db,beta=map(mp.mpf,(p['da'],p['db'],p['beta']))
x=mp.matrix([mp.mpf(float(v)) for v in u.ravel()])
def rhs(x):
 a=mp.matrix(x[:n]);b=mp.matrix(x[n:]);ad=D*a;bd=D*b
 return mp.matrix([a[i]**2/b[i]-a[i]+Da*ad[i] for i in range(n)]+[beta*(a[i]**2-b[i])+Db*bd[i] for i in range(n)])
def jac(x):
 J=mp.zeros(2*n)
 for i in range(n):
  for j in range(n):J[i,j]=Da*D[i,j];J[n+i,n+j]=Db*D[i,j]
  a,b=x[i],x[n+i];J[i,i]+=2*a/b-1;J[i,n+i]=-a*a/(b*b);J[n+i,i]=2*beta*a;J[n+i,n+i]-=beta
 return J
initial=x.copy();clock=time.perf_counter();iterations=[]
for iteration in range(8):
 f=rhs(x);residual=max(abs(q) for q in f);iterations.append(str(residual))
 if residual<mp.mpf('1e-65'):break
 x-=mp.lu_solve(jac(x),f)
assert max(abs(q) for q in rhs(x))<mp.mpf('1e-65') and all(q>0 for q in x)
equilibrium=np.array([np.longdouble(str(q)) for q in x]).reshape(u.shape)
J=np.array([[float(q) for q in row] for row in jac(x).tolist()]);rightmost=float(np.max(eigvals(J).real))
def ld_rhs(values):
 a,b=values;return np.stack([a*a/b-a+np.longdouble(p['da'])*(delta@a),np.longdouble(p['beta'])*(a*a-b)+np.longdouble(p['db'])*(delta@b)])
delta=d.astype(np.longdouble);records=[];files=[source,ROOT/'protocol.json',ROOT/'summary.json',ROOT/'status.json',Path(__file__)]
r0=ld_rhs(u.astype(np.longdouble));initial_rate=float(np.max(abs(r0/u)));paths={'times':np.arange(401)*.15,'equilibrium':equilibrium.astype(np.float64)}
for level in p['levels']:
 label=level['level'];file=ROOT/'production'/label/'paths.npz';files.append(file)
 with np.load(file) as z:saved=z['trajectory'][:,list(z['trial_keys']).index(trial['key'])].copy()
 dt=np.longdouble(level['dt']);values=u.astype(np.longdouble);history=[values.copy()];steps=round(60/level['dt']);every=round(.15/level['dt'])
 for k in range(1,steps+1):
  first=values+dt*ld_rhs(values);second=first+dt*ld_rhs(first);values=np.longdouble(.5)*values+np.longdouble(.5)*second
  assert np.isfinite(values).all() and np.all(values>0)
  if k%every==0:history.append(values.copy())
 history=np.array(history);assert history.shape==saved.shape
 drift=float(np.max(abs(np.log(saved/u))));ld_drift=float(np.max(abs(np.log(history/u))))
 records.append(dict(level=label,dt=float(dt),float64_saved_max_log_drift=drift,float64_saved_exactly_constant=bool(np.array_equal(saved,np.broadcast_to(u,saved.shape))),longdouble_max_log_drift=ld_drift,longdouble_endpoint_log_distance_to_equilibrium=float(np.max(abs(np.log(history[-1]/equilibrium)))),float64_endpoint_log_distance_to_equilibrium=float(np.max(abs(np.log(saved[-1].astype(np.longdouble)/equilibrium)))),first_euler_largest_increment_in_float64_ulps=float(np.max(abs(float(dt)*r0)/np.spacing(u)))))
 paths[label+'_longdouble']=history.astype(np.float64)
ref=ROOT/'references'/trial['key']/'paths.npz';files.append(ref)
with np.load(ref) as z:refs={label:float(np.max(abs(np.log(z[label][-1].astype(np.longdouble)/equilibrium)))) for label in ('dop853','radau')}
np.savez_compressed(HERE/'roundoff-diagnostic-paths.npz',**paths)
result=dict(diagnostic_complete=True,original_order_gate_unchanged=True,delayed_full_study_accepted=False,failed_order_record=failed,trial=trial['key'],mantissa_bits=dict(float64=np.finfo(np.float64).nmant,longdouble=np.finfo(np.longdouble).nmant),equilibrium_precision_digits=80,equilibrium_newton_residual_history=iterations,equilibrium_maximum_absolute_residual=str(max(abs(q) for q in rhs(x))),initial_maximum_relative_rate=initial_rate,initial_maximum_log_distance_to_equilibrium=float(np.max(abs(np.log(u.astype(np.longdouble)/equilibrium)))),reference_endpoint_log_distance_to_equilibrium=refs,approximate_rightmost_local_jacobian_eigenvalue=rightmost,levels=records,seconds=time.perf_counter()-clock,paths_sha256=sha(HERE/'roundoff-diagnostic-paths.npz'),files_sha256={str(f):sha(f) for f in files},scope='Numerical diagnosis only: an80-digit stationary solution and identical SSP-RK2 calculations with63-bit longdouble mantissa. No delayed control replacement, source projection, new scientific perturbation, acceptance-floor amendment or threshold change. The original failed order gate is retained.')
(HERE/'roundoff-diagnostic.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
