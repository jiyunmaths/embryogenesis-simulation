"""Independent checks of formation and perturbed release in finite reservoirs."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .attribute_finite_reservoir import system


def verify(output):
    protocol=json.loads((output/'protocol.json').read_text());v=np.array(protocol['volumes']);n=len(v)
    rows=[]
    for arm in protocol['arms']:
        if arm['scale']!=4 or arm['kind'] not in ['formation','release']:
            continue
        if arm['kind']=='formation' and arm.get('seed') is None:
            continue
        data=np.load(output/(arm['name']+'.npz'));times=data['times'];start=data['initial'].copy()
        fun,jac=system(v,arm['size']*v.sum(),protocol['beta'],arm['ka'],arm['kb'])
        if arm['kind']=='release':
            rng=np.random.default_rng(7)
            changed=start[:,:n]*np.exp(1e-4*rng.normal(size=(2,n)))
            changed*=((start[:,:n]@v)/(changed@v))[:,None]
            start[:,:n]=changed
        sol=solve_ivp(fun,(0,times[-1]),start.ravel(),t_eval=times,method='DOP853',
                      rtol=1e-10,atol=1e-12,max_step=2.)
        if not sol.success or np.any(sol.y<=0):raise RuntimeError(sol.message)
        path=sol.y.T.reshape(len(times),2,n+1)
        if arm['kind']=='release':
            check=solve_ivp(fun,(0,times[-1]),start.ravel(),t_eval=times,method='Radau',jac=jac,
                            rtol=1e-10,atol=1e-12,max_step=2.)
            if not check.success or np.any(check.y<=0):raise RuntimeError(check.message)
            reference=check.y.T.reshape(path.shape)
        else:
            reference=data['trajectory']
        err=float(np.abs(np.log(path/reference)).max())
        spread=np.std(np.log(path[:,0,:n]),axis=1);late=times>=.8*times[-1]
        residual=float(np.abs(fun(times[-1],path[-1].ravel())).max())
        growth=float(np.linalg.eigvals(jac(times[-1],path[-1].ravel())).real.max())
        row=dict(name=arm['name'],intervention='amount-preserving log-noise 1e-4, seed 7' if arm['kind']=='release' else 'original seeded start',
                 independent_max_log_error=err,numerical_pass=err<1e-5,
                 late_min_log_activator_sd=float(spread[late].min()),final_rhs_max=residual,
                 final_jacobian_max_real=growth,persistent=bool(spread[late].min()>.1),
                 stable_stationary=bool(residual<1e-6 and growth<0),final_reservoir=path[-1,:,-1].tolist())
        rows.append(row)
        np.savez_compressed(output/('verified_'+arm['name']+'.npz'),times=times,trajectory=path,initial=start)
        (output/'verification.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('outputs/attribute-finite-reservoir'))
    verify(p.parse_args().output)
