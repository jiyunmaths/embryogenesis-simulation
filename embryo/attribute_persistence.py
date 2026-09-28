"""Frozen-context chemical recovery assay; no fate labels or clustering."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .attribute_development import AttributeSimulation
from .signaling import stability


def perturb(values, masses, amplitude, rng):
    """Positive multiplicative perturbation preserving each species' total amount."""
    changed = values * np.exp(amplitude * rng.normal(size=values.shape))
    return changed * ((values @ masses) / (changed @ masses))[:, None]


def rhs(delta, beta, da, db):
    n = len(delta)
    def evaluate(t, y):
        a, b = y[:n], y[n:]
        return np.concatenate((a*a/b-a+da*(delta@a),
                               beta*(a*a-b)+db*(delta@b)))
    return evaluate


def trajectory(fun, initial, times, tight=False):
    sol = solve_ivp(fun, (times[0], times[-1]), initial.ravel(),
                    method='DOP853', t_eval=times,
                    rtol=1e-11 if tight else 1e-8,
                    atol=1e-13 if tight else 1e-10)
    if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y <= 0):
        raise RuntimeError('Chemical integration failed positivity/finiteness check')
    return sol.y.T.reshape(len(times), *initial.shape)


def distance(x, y, masses):
    """Volume-weighted log RMS, retaining cell identities and species equally."""
    return np.sqrt(np.sum((np.log(x)-np.log(y))**2*masses, axis=(-1,-2))
                   / (2*masses.sum()))


def run(source, output):
    output.mkdir(parents=True, exist_ok=False)
    times = np.linspace(0, 120, 241)
    protocol = dict(source=str(source.resolve()), duration=120, sample_interval=.5,
                    amplitudes=[.01,.05], seeds=list(range(20)),
                    recovery_ratio_limit=.1, reference_log_rms_limit=1e-5,
                    scope='Frozen final geometry; chemical recovery relative to matched unperturbed continuation. No fate dynamics, clustering, or mechanical recovery.',
                    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [Path(__file__),Path(__file__).with_name('attribute_development.py'),
                                Path(__file__).with_name('model.py'),Path(__file__).with_name('transport.py'),
                                Path(__file__).with_name('signaling.py')]},
                    checkpoint_sha256={mode:hashlib.sha256((source/mode/'final_state.npz').read_bytes()).hexdigest()
                                       for mode in ['direct','no_feedback']})
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    result = {}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(10,4))
    for ax, mode in zip(axes,['direct','no_feedback']):
        sim = AttributeSimulation.restore(source/mode/'final_state.npz')
        graph = sim.signaling_graph(); masses = graph.masses
        c = sim.config
        fun = rhs(graph.delta,c.signal_beta,c.signal_da,c.signal_dh)
        initial = np.array([sim.activator,sim.inhibitor])
        control = trajectory(fun,initial,times)
        reference = trajectory(fun,initial,times,True)
        max_error = float(distance(control,reference,masses).max())
        trials=[]; paths=[]
        for amplitude in protocol['amplitudes']:
            curves=[]
            for seed in protocol['seeds']:
                start=perturb(initial,masses,amplitude,np.random.default_rng(seed))
                path=trajectory(fun,start,times)
                ref=trajectory(fun,start,times,True)
                max_error=max(max_error,float(distance(path,ref,masses).max()))
                d=distance(path,control,masses)
                ratio=float(d[times>=100].max()/d[0])
                trials.append(dict(amplitude=amplitude,seed=seed,initial_distance=float(d[0]),
                                   late_max_distance=float(d[times>=100].max()),
                                   recovery_ratio=ratio,recovered=ratio<.1))
                paths.append(path);curves.append(d)
            curves=np.array(curves)
            ax.plot(times,np.median(curves,axis=0),label=f'noise scale {amplitude}')
            ax.fill_between(times,curves.min(axis=0),curves.max(axis=0),alpha=.2)
        y=control[-1].ravel(); n=len(masses); a,b=control[-1]
        jac=np.block([[np.diag(2*a/b-1)+c.signal_da*graph.delta,np.diag(-a*a/b**2)],
                      [np.diag(2*c.signal_beta*a),-c.signal_beta*np.eye(n)+c.signal_dh*graph.delta]])
        result[mode]=dict(trials=trials,reference_max_log_rms=max_error,
            reference_pass=max_error<1e-5,
            control_drift_from_start=float(distance(control[-1],initial,masses)),
            control_late_drift=float(distance(control[-1],control[times==100][0],masses)),
            final_rhs_max=float(np.abs(fun(120,y)).max()),
            final_jacobian_max_real=float(np.linalg.eigvals(jac).real.max()),
            homogeneous_stability=stability(graph,c.signal_beta,c.signal_da,c.signal_dh))
        np.savez_compressed(output/f'{mode}.npz',times=times,control=control,
                            perturbed=np.array(paths),masses=masses,delta=graph.delta,
                            initial=initial,cell_ids=sim.ids)
        ax.set(title=mode,xlabel='Time since geometry frozen',ylabel='Log chemical distance from control',yscale='log');ax.legend()
    fig.tight_layout();fig.savefig(output/'recovery.png',dpi=160);plt.close(fig)
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('outputs/attribute-development'))
    parser.add_argument('--output',type=Path,default=Path('outputs/attribute-persistence'))
    args=parser.parse_args();run(args.source,args.output)
