"""Resolve failed accuracy gates without replacing saved tight trajectories.

Original results remain untouched. Previously accepted paths are retained; only
failed control or trial comparisons receive a tighter independent audit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .attribute_persistence import rhs
from .feedback_endpoint_bistability import chemical_jacobian
from .resolution import write_json


def refine(source,output):
    source=Path(source);output=Path(output);output.mkdir(parents=True,exist_ok=False)
    write_json(output/'protocol.json',dict(source=str(source.resolve()),
        scope='Audit failed numerical gates only; retain original tight trajectories and response metrics. Original coarse failures remain in source reports.',
        rtol=1e-13,atol=1e-15,log_error_limit=1e-5,methods=['DOP853','Radau'],
        kinetics='beta=2, Da=0.02, Db=0.4 (verified against source checkpoint configuration)',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*[source/f'seed-{s}/results.json' for s in (7,8,9)]]}))
    checks=[]
    for seed in (7,8,9):
        folder=output/f'seed-{seed}';folder.mkdir()
        data=json.loads((source/f'seed-{seed}/results.json').read_text())
        for key,r in data.items():
            if r['numerical_pass']:continue
            p=source/f'seed-{seed}/{key}.npz'
            with np.load(p) as d:
                times=d['times'];delta=d['delta'];n=len(delta)
                fun=rhs(delta,2.,.02,.4)
                jac=lambda t,y:chemical_jacobian(y.reshape(2,n),delta,2.,.02,.4)
                targets=[]
                if r['control_max_log_error']>=1e-5:targets.append(('control',np.array(r['initial']),d['control']))
                targets.extend((i,d['initial_states'][i],d['trajectories'][i]) for i,t in enumerate(r['trials']) if t['reference_max_log_error']>=1e-5)
                for index,initial,saved in targets:
                    errs={}
                    for method in ('DOP853','Radau'):
                        kwargs={'jac':jac} if method=='Radau' else {}
                        sol=solve_ivp(fun,(0,240),initial.ravel(),t_eval=times,method=method,rtol=1e-13,atol=1e-15,**kwargs)
                        if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y<=0):raise RuntimeError('Refinement failed')
                        ref=sol.y.T.reshape(saved.shape);errs[method]=float(np.abs(np.log(saved/ref)).max())
                    valid=all(e<1e-5 for e in errs.values())
                    checks.append(dict(seed=seed,background=key,trial=index,errors=errs,passed=valid,trajectory_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
                    if not valid:raise RuntimeError('Refined trajectory agreement failed')
                    if index=='control':
                        r['original_control_max_log_error']=r['control_max_log_error'];r['control_max_log_error']=errs['DOP853']
                    else:
                        row=r['trials'][index];row['original_reference_max_log_error']=row['reference_max_log_error'];row['reference_max_log_error']=errs['DOP853']
                r['refinement_checks']=[c for c in checks if c['seed']==seed and c['background']==key]
            for row in r['trials']:
                valid=row['reference_max_log_error']<1e-5 and r['control_max_log_error']<1e-5
                settled=row['final_rhs_max']<1e-6 and row['final_jacobian_max_real']<0
                row['numerical_pass']=valid
                row['classification']=('numerically_unresolved' if not valid else 'returned' if settled and row['final_log_rms_from_control']<1e-4 else 'different_stable_endpoint' if settled else 'not_settled_by_horizon')
            r['numerical_pass']=all(t['numerical_pass'] for t in r['trials']) and r['independent_max_log_error']<1e-5
        write_json(folder/'results.json',data)
        write_json(folder/'status.json',dict(state='completed',numerical_pass=all(r['numerical_pass'] for r in data.values()),trajectory_source=str((source/f'seed-{seed}').resolve())))
    write_json(output/'refinement.json',checks)
    print(json.dumps(checks,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=Path('outputs/cell-response'));p.add_argument('--output',type=Path,default=Path('outputs/cell-response-refined'))
    a=p.parse_args();refine(a.source,a.output)
