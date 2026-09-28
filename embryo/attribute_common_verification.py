"""Additional tolerance refinement and independent single-cell Radau checks."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .attribute_common_environment import bath_rhs


def verify(output):
    protocol=json.loads((output/'protocol.json').read_text())
    results=json.loads((output/'results.json').read_text())
    data=np.load(output/'trajectories.npz');times=data['times']
    starts=data['initial_patterns'].transpose(1,0,2).reshape(2,-1)
    report={}
    for arm in protocol['arms']:
        name=arm['name'];row=results[name]
        fun=bath_rhs(protocol['config']['signal_beta'],row['ka'],row['kb'],row['bath'])
        sol=solve_ivp(fun,(0,240),starts.ravel(),method='DOP853',t_eval=times,rtol=1e-13,atol=1e-15)
        if not sol.success or np.any(sol.y<=0):
            raise RuntimeError('Refinement failed')
        path=sol.y.T.reshape(len(times),*starts.shape)
        errors=np.abs(np.log(path/data[name]))
        worst=int(np.unravel_index(np.argmax(errors),errors.shape)[2])
        # The worst affected initial state receives an independent implicit check.
        reference=solve_ivp(fun,(0,240),starts[:,worst],method='Radau',t_eval=times,rtol=1e-11,atol=1e-13)
        if not reference.success or np.any(reference.y<=0):
            raise RuntimeError('Radau reference failed')
        independent_error=float(np.abs(np.log(path[:,:,worst]/reference.y.T)).max())
        stable=np.array(row['stable_equilibria'])
        distances=np.sqrt(np.mean(np.log(path[-1].T[:,None,:]/stable[None,:,:])**2,axis=2))
        counts=np.bincount(distances.argmin(axis=1),minlength=len(stable)).tolist()
        max_error=float(errors.max())
        report[name]=dict(further_refinement_max_log_error=max_error,
            independent_worst_state_index=worst,independent_radau_max_log_error=independent_error,
            equilibrium_counts_unchanged=counts==row['nearest_stable_equilibrium_counts'],
            passed=bool(max_error<1e-5 and independent_error<1e-5 and counts==row['nearest_stable_equilibrium_counts']))
    (output/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('outputs/attribute-common-environment'))
    verify(p.parse_args().output)
