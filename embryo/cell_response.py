"""Single-cell fractional-pulse response on validated frozen endpoint graphs.

Positive and negative instantaneous activator boluses are external interventions,
not amount-preserving redistribution. No prescribed cell identities are used.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, trapezoid

from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs, trajectory, distance
from .feedback_endpoint_bistability import chemical_jacobian
from .resolution import write_json

FACTORS = (.5, .75, .9, 1.1, 1.25, 1.5)
TIMES = np.unique(np.concatenate((np.arange(0, 10.01, .1),
                                  np.arange(10, 60.01, .5), np.arange(60, 240.01, 2.))))


def pulse(initial, cell, factor):
    values = np.array(initial, dtype=float, copy=True)
    if values.ndim != 2 or values.shape[0] != 2 or not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('Expected two positive finite species vectors')
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError('Pulse factor must be positive and finite')
    if not 0 <= cell < values.shape[1]:
        raise ValueError('Invalid target cell')
    values[0, cell] *= factor
    return values


def sustained_recovery(times, deviation, threshold, late_window=24.):
    """First sampled time after which all samples pass, with >=24 units observed."""
    passing = np.asarray(deviation) <= threshold
    suffix = np.logical_and.accumulate(passing[::-1])[::-1]
    indices = np.flatnonzero(suffix & (np.asarray(times) <= times[-1]-late_window))
    return float(times[indices[0]]) if len(indices) else None


def response_metrics(path, control, times, masses, cell, factor):
    log_change = np.log(path/control)
    amplitude = abs(np.log(factor))
    if amplitude == 0:
        raise ValueError('Response normalization requires a nonzero pulse')
    target = np.abs(log_change[:, :, cell]).max(axis=1)
    other = np.arange(len(masses)) != cell
    other_weights = masses[other]
    neighbor = np.sqrt(np.sum(log_change[:, :, other]**2*other_weights,axis=(1,2)) / (2*other_weights.sum()))
    local = np.abs(log_change[:, 0, cell])
    network = distance(path, control, masses)
    return dict(target_activator_peak_gain=float(local.max()/amplitude),
        target_inhibitor_peak_gain=float(np.abs(log_change[:,1,cell]).max()/amplitude),
        target_activator_log_auc_per_log_pulse=float(trapezoid(local,times)/amplitude),
        other_cells_peak_log_rms_per_log_pulse=float(neighbor.max()/amplitude),
        target_recovery_time=sustained_recovery(times,target,.1*amplitude),
        network_recovery_time=sustained_recovery(times,network,.1*network[0]),
        final_log_rms_from_control=float(network[-1]),
        late_max_log_rms_from_control=float(network[times>=times[-1]-24.].max()))


def run(seed, output):
    if seed not in (7,8,9):
        raise ValueError('This protocol uses the validated seeds 7, 8 and 9')
    output=Path(output)
    source=Path('outputs/feedback-endpoint-bistability'+('' if seed==7 else f'-seed-{seed}'))
    checkpoints=Path('outputs/feedback-survival' if seed==7 else f'outputs/feedback-survival-validation/seed-{seed}/survival')
    evidence=json.loads((source/'results.json').read_text())
    if not all(r['chemical_bistability_supported'] for r in evidence.values()):
        raise ValueError('Endpoint bistability validation must pass')
    prior=json.loads((source/'protocol.json').read_text())
    for filename,expected in prior['sha256'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Validated source/checkpoint changed: {filename}')
    inputs=[source/'results.json',source/'protocol.json']+[source/f'{arm}.npz' for arm in ('switch_on','keep_off')]
    inputs += [checkpoints/arm/'latest_state.npz' for arm in ('switch_on','keep_off')]
    inputs += [Path(__file__),Path(__file__).with_name('attribute_persistence.py'),Path(__file__).with_name('feedback_endpoint_bistability.py')]
    output.mkdir(parents=True,exist_ok=False)
    protocol=dict(seed=seed,factors=FACTORS,times=TIMES.tolist(),horizon=240,
        initial_states='Uniform and patterned equilibria from the validated 240-unit endpoint continuations; not raw moving t=150 chemistry.',
        intervention='At time zero multiply activator of one cell by factor, leaving every other concentration unchanged. External amount added or removed; no compensating redistribution.',
        controls='One unperturbed trajectory per equilibrium and graph, shared by all target-cell interventions.',
        recovery='First sampled time after which all remaining samples are <=10% of initial log perturbation, requiring >=24 units of subsequent observation. Unrecovered is right-censored, not proof of permanent memory.',
        criteria=dict(max_tolerance_log_error=1e-5,max_radau_log_error=1e-5,endpoint_rhs=1e-6,equilibrium_return_log_rms=1e-4),
        analysis='Continuous responses; no fate labels or clusters. Cell trials are nested interventions, not independent developmental replicates. Geometry/polarity frozen; no inference about active mechanical response.',
        sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    write_json(output/'protocol.json',protocol)
    results={}
    write_json(output/'status.json',dict(state='running',seed=seed,completed=0,total=384))
    try:
        for arm in ('switch_on','keep_off'):
            sim=AttributeSimulation.restore(checkpoints/arm/'latest_state.npz')
            g=sim.signaling_graph();c=sim.config;m=g.masses
            if len(m)!=16:raise ValueError("Protocol requires 16-cell endpoints")
            with np.load(source/f'{arm}.npz') as data:
                assert np.array_equal(sim.ids,data['ids'])
                assert np.allclose(g.delta,data['delta'],rtol=1e-12,atol=1e-12)
                assert np.allclose(m,data['volumes'],rtol=1e-12,atol=1e-12)
                equilibria={name:data['trajectories'][i,-1].copy() for i,name in enumerate(('pattern','uniform'))}
            fun=rhs(g.delta,c.signal_beta,c.signal_da,c.signal_dh)
            jac=lambda t,y:chemical_jacobian(y.reshape(2,-1),g.delta,c.signal_beta,c.signal_da,c.signal_dh)
            _,exposure=sim.contacts();centers=sim.centers();center=(centers*m[:,None]).sum(axis=0)/m.sum()
            context=dict(ids=sim.ids.tolist(),volume=m.tolist(),exposure=np.asarray(exposure).tolist(),
                weighted_degree=np.asarray(-np.diag(g.delta)*m).tolist(),radius=np.linalg.norm(centers-center,axis=1).tolist())
            for family,initial in equilibria.items():
                key=f'{arm}_{family}'
                initial_residual=float(np.abs(fun(0,initial.ravel())).max())
                initial_growth=float(np.linalg.eigvals(jac(0,initial.ravel())).real.max())
                if initial_residual>=1e-6 or initial_growth>=0:
                    raise ValueError('Initial state must be a validated stable chemical equilibrium')
                control=trajectory(fun,initial,TIMES,True)
                control_coarse=trajectory(fun,initial,TIMES)
                control_error=float(np.abs(np.log(control/control_coarse)).max())
                paths=[];starts=[];trials=[]
                for cell,cell_id in enumerate(sim.ids):
                    for factor in FACTORS:
                        start=pulse(initial,cell,factor)
                        coarse=trajectory(fun,start,TIMES);path=trajectory(fun,start,TIMES,True)
                        error=float(np.abs(np.log(path/coarse)).max())
                        residual=float(np.abs(fun(240,path[-1].ravel())).max())
                        growth=float(np.linalg.eigvals(jac(240,path[-1].ravel())).real.max())
                        metrics=response_metrics(path,control,TIMES,m,cell,factor)
                        returned=metrics['final_log_rms_from_control']<1e-4
                        numerical=error<1e-5 and control_error<1e-5
                        settled=residual<1e-6 and growth<0
                        classification=('numerically_unresolved' if not numerical else 'returned' if returned and settled
                                        else 'different_stable_endpoint' if settled else 'not_settled_by_horizon')
                        row=dict(cell_index=cell,cell_id=int(cell_id),factor=factor,
                            initial_activator=float(initial[0,cell]),initial_inhibitor=float(initial[1,cell]),
                            injected_activator_amount=float((factor-1)*initial[0,cell]*m[cell]),
                            reference_max_log_error=error,final_rhs_max=residual,final_jacobian_max_real=growth,
                            numerical_pass=numerical,classification=classification,**metrics)
                        trials.append(row);paths.append(path);starts.append(start)
                worst=int(np.argmax([r['reference_max_log_error'] for r in trials]))
                ref=solve_ivp(fun,(0,240),starts[worst].ravel(),method='Radau',jac=jac,t_eval=TIMES,rtol=1e-11,atol=1e-13)
                if not ref.success or not np.isfinite(ref.y).all() or np.any(ref.y<=0):
                    raise RuntimeError('Independent reference failed')
                err=float(np.abs(np.log(paths[worst]/ref.y.T.reshape(paths[worst].shape))).max())
                valid=control_error<1e-5 and err<1e-5 and all(r['numerical_pass'] for r in trials)
                results[key]=dict(context=context,initial=initial.tolist(),initial_residual=initial_residual,
                    initial_growth=initial_growth,control_max_log_error=control_error,
                    control_drift=float(distance(control[-1],initial,m)),independent_trial=worst,
                    independent_max_log_error=err,numerical_pass=valid,trials=trials)
                np.savez_compressed(output/f'{key}.npz',times=TIMES,control=control,trajectories=np.asarray(paths),
                    initial_states=np.asarray(starts),delta=g.delta,masses=m,cell_ids=sim.ids)
                write_json(output/'results.json',results)
                write_json(output/'status.json',dict(state='running',seed=seed,completed=sum(len(r['trials']) for r in results.values()),total=384))
                print(f'Seed {seed} {key}: {len(trials)} trials, numerical pass={valid}',flush=True)
        write_json(output/'status.json',dict(state='completed',seed=seed,completed=384,total=384,
                    numerical_pass=all(r['numerical_pass'] for r in results.values())))
    except Exception as exc:
        write_json(output/'status.json',dict(state='failed',seed=seed,error=str(exc)))
        raise
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.seed,a.output)
