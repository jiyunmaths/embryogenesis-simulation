"""Joint stage-consistent regulator/fate integration in equilibrium deviations."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .causal_signaling import initial_signals, coefficients
from .model import Simulation
from .resolution import write_json, _steps


def rhs(state,delta,config,arm):
    u,v,z=state
    if arm=='no_self_activation':fu=(-u-v-u*v)/(1+v)
    elif arm=='no_inhibitor_action':fu=u*(1+u)
    else:fu=(1+u)*(u-v)/(1+v)
    fv=config.signal_beta*(-v if arm=='no_induced_inhibitor' else 2*u+u*u-v)
    da,dh=coefficients(arm,config.signal_da,config.signal_dh)
    gain=0. if arm=='no_signal_to_fate' else config.signal_fate_gain
    return np.stack([fu+da*(u@delta.T),fv+dh*(v@delta.T),config.fate_rate*(z-z**3+gain*u)])


def advance(state,delta,config,arm,dt):
    da,dh=coefficients(arm,config.signal_da,config.signal_dh)
    exit_rate=float(np.max(-np.diag(delta)))
    # Same signal positivity restriction as the existing integrator. Fate is
    # advanced jointly and monitored; no clipping or saturation modification.
    subdivisions=max(1,int(np.ceil(dt*max(1+da*exit_rate,config.signal_beta+dh*exit_rate)/.2)))
    value=np.array(state,float,copy=True);step=dt/subdivisions
    with np.errstate(over='raise',divide='raise',invalid='raise'):
        for _ in range(subdivisions):
            stage=value+step*rhs(value,delta,config,arm)
            if np.any(stage[0]<-1) or np.any(stage[1]<=-1):raise FloatingPointError('joint stage lost positivity')
            value=.5*value+.5*(stage+step*rhs(stage,delta,config,arm))
            if not np.isfinite(value).all() or np.any(value[0]<-1) or np.any(value[1]<=-1):raise FloatingPointError('joint solution lost positivity')
    return value


def trajectory(graph,config,seeds,arm,dt=.0075,until=60.,interval=.6):
    a,h=initial_signals(graph.masses,seeds,.001)
    state=np.stack([a-1,h-1,np.zeros_like(a)])
    every=_steps(interval,dt);steps=_steps(until,dt);history=[state.copy()]
    if steps%every:raise ValueError('observations must align')
    for i in range(1,steps+1):
        state=advance(state,graph.delta,config,arm,dt)
        if i%every==0:history.append(state.copy())
    return np.array(history)


def adaptive_reference(graph,config,seeds,arm,times,rtol=1e-10,atol=1e-13):
    a,h=initial_signals(graph.masses,seeds,.001);state=np.stack([a-1,h-1,np.zeros_like(a)])
    shape=state.shape
    solution=solve_ivp(lambda t,y:rhs(y.reshape(shape),graph.delta,config,arm).ravel(),
                       (0,float(times[-1])),state.ravel(),method='DOP853',t_eval=times,rtol=rtol,atol=atol)
    if not solution.success:raise RuntimeError(solution.message)
    return solution.y.T.reshape((len(times),)+shape)


def discrepancies(a,b,threshold=.55):
    return {'max_signal_rms':float(np.max(np.sqrt(np.mean((a[:,:2]-b[:,:2])**2,axis=-1)))),
            'max_fate_difference':float(np.max(abs(a[:,2]-b[:,2]))),
            'final_labels_match':bool(np.array_equal(np.where(a[-1,2]>threshold,1,np.where(a[-1,2]<-threshold,-1,0)),np.where(b[-1,2]>threshold,1,np.where(b[-1,2]<-threshold,-1,0))))}


def run(checkpoint,output):
    checkpoint=Path(checkpoint).resolve();output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh directory')
    sim=Simulation.restore(checkpoint);graph=sim.signaling_graph()
    p={'checkpoint':str(checkpoint),'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
       'arms':['full','no_self_activation','no_transport','equal_diffusion','no_signal_to_fate'],
       'seeds':list(range(20)),'dts':[.0075,.00375],'until':60.,'interval':.6,
       'reference_tolerances':[[1e-10,1e-13],[1e-12,1e-15]],
       'criteria':{'signal_rms_max':.01,'fate_difference_max':.05,'reference_fate_difference_max':.001},
       'scope':__doc__,'model_sha256':hashlib.sha256(Path(__file__).with_name('model.py').read_bytes()).hexdigest()}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'running'})
    try:
        rows=[];times=np.linspace(0,p['until'],round(p['until']/p['interval'])+1)
        for arm in p['arms']:
            coarse,fine=[trajectory(graph,sim.config,p['seeds'],arm,dt,p['until'],p['interval']) for dt in p['dts']]
            ref,tight=[adaptive_reference(graph,sim.config,p['seeds'],arm,times,*t) for t in p['reference_tolerances']]
            row={'arm':arm,'step_halving':discrepancies(coarse,fine),
                 'coarse_vs_reference':discrepancies(coarse,tight),'fine_vs_reference':discrepancies(fine,tight),
                 'reference_tightening':discrepancies(ref,tight)}
            np.savez_compressed(output/f'{arm}.npz',times=times,coarse=coarse,fine=fine,reference=tight)
            rows.append(row);write_json(output/'progress.json',rows);print(row,flush=True)
        c=p['criteria'];checks={
            'reference_independently_converged':all(r['reference_tightening']['max_fate_difference']<c['reference_fate_difference_max'] and r['reference_tightening']['final_labels_match'] for r in rows),
            'both_steps_match_reference_signals':all(r[key]['max_signal_rms']<c['signal_rms_max'] for r in rows for key in ('coarse_vs_reference','fine_vs_reference','step_halving')),
            'both_steps_match_reference_fates':all(r[key]['max_fate_difference']<c['fate_difference_max'] and r[key]['final_labels_match'] for r in rows for key in ('coarse_vs_reference','fine_vs_reference','step_halving')),
            'checkpoint_unchanged':hashlib.sha256(checkpoint.read_bytes()).hexdigest()==p['checkpoint_sha256']}
        result={'protocol':p,'rows':rows,'checks':checks,'passed':all(checks.values())}
        write_json(output/'comparison.json',result)
        lines=['# Joint fate integration validation','',f'All checks pass: **{result["passed"]}**','',
               '| Arm | Fate error dt=0.0075 vs reference | Fate error dt=0.00375 vs reference | Reference tightening |',
               '|---|---:|---:|---:|']
        lines += [f'| {r["arm"]} | {r["coarse_vs_reference"]["max_fate_difference"]:.6g} | {r["fine_vs_reference"]["max_fate_difference"]:.6g} | {r["reference_tightening"]["max_fate_difference"]:.6g} |' for r in rows]
        lines += ['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()]]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n');write_json(output/'status.json',{'state':'completed','passed':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--checkpoint',type=Path,default=Path('outputs/development-refinement/space-72/state-18.npz'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.checkpoint,args.output)
