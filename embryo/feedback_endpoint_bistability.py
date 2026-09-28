"""Test coexistence of chemical attractors on post-switch endpoint geometries."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs, trajectory, perturb, distance
from .resolution import write_json
from .signaling import stability


def chemical_jacobian(state,delta,beta,da,db):
    a,b=state;n=len(a)
    return np.block([[np.diag(2*a/b-1)+da*delta,np.diag(-a*a/b**2)],
                     [np.diag(2*beta*a),-beta*np.eye(n)+db*delta]])


def run(source,output):
    source=Path(source);output=Path(output)
    status=json.loads((source/'status.json').read_text())
    comparison=json.loads((source/'comparison.json').read_text())
    if status['state']!='completed' or not all(r['quality_pass'] for r in comparison['arms'].values()):
        raise ValueError('Requires completed, quality-passing moving experiment')
    if output.exists():raise FileExistsError('Choose a fresh output directory')
    output.mkdir(parents=True)
    files=[Path(__file__),Path(__file__).with_name('attribute_development.py'),Path(__file__).with_name('attribute_persistence.py'),
           Path(__file__).with_name('model.py'),Path(__file__).with_name('transport.py'),Path(__file__).with_name('signaling.py')]
    checkpoints={arm:source/arm/'latest_state.npz' for arm in ['switch_on','keep_off']}
    protocol=dict(duration=240.,interval=1.,seeds=list(range(20)),near_uniform_log_noise=.001,pattern_log_noise=.01,
        criteria=dict(max_reference_log_error=1e-5,final_rhs_max=1e-6,pattern_log_sd_min=.1,
                      uniform_max_log_deviation=1e-4,return_log_rms_max=1e-4),
        scope='Frozen t=150 endpoint graphs. Local chemical attractor coexistence, not stability of the complete moving system. One developed embryo; perturbations are not independent developmental replicates.',
        sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in files+list(checkpoints.values())})
    write_json(output/'protocol.json',protocol)
    times=np.arange(241.);results={}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,(arm,checkpoint) in zip(axes,checkpoints.items()):
        sim=AttributeSimulation.restore(checkpoint)
        if abs(sim.time-150)>1e-9:raise ValueError('Requires t=150 endpoint')
        graph=sim.signaling_graph();c=sim.config;v=graph.masses
        fun=rhs(graph.delta,c.signal_beta,c.signal_da,c.signal_dh)
        jac=lambda t,y:chemical_jacobian(y.reshape(2,-1),graph.delta,c.signal_beta,c.signal_da,c.signal_dh)
        pattern=np.array([sim.activator,sim.inhibitor]);unit=np.ones_like(pattern)
        starts=[pattern,unit];labels=['pattern_control','uniform_control']
        for group,base,scale in [('pattern_perturbed',pattern,.01),('uniform_perturbed',unit,.001)]:
            for seed in range(20):
                starts.append(perturb(base,v,scale,np.random.default_rng(seed)));labels.append(f'{group}_{seed}')
        paths=[];errors=[];rows=[]
        for label,start in zip(labels,starts):
            coarse=trajectory(fun,start,times);path=trajectory(fun,start,times,True)
            error=float(np.abs(np.log(coarse/path)).max());errors.append(error);paths.append(path)
            end=path[-1];residual=float(np.abs(fun(240,end.ravel())).max())
            growth=float(np.linalg.eigvals(jac(240,end.ravel())).real.max())
            rows.append(dict(label=label,reference_max_log_error=error,final_rhs_max=residual,
                final_log_activator_sd=float(np.std(np.log(end[0]))),max_log_deviation_from_uniform=float(np.abs(np.log(end)).max()),
                final_jacobian_max_real=growth,
                distance_to_pattern_control=float(distance(end,paths[0][-1],v))))
        # Independent Radau reference for the most tolerance-sensitive trajectory.
        worst=int(np.argmax(errors));reference=solve_ivp(fun,(0,240),starts[worst].ravel(),method='Radau',jac=jac,
            t_eval=times,rtol=1e-11,atol=1e-13)
        if not reference.success or np.any(reference.y<=0):raise RuntimeError('Independent reference failed')
        independent_error=float(np.abs(np.log(paths[worst]/reference.y.T.reshape(paths[worst].shape))).max())
        pattern_rows=[r for r in rows if r['label'].startswith('pattern')]
        uniform_rows=[r for r in rows if r['label'].startswith('uniform')]
        numerical=all(r['reference_max_log_error']<1e-5 and r['final_rhs_max']<1e-6 for r in rows) and independent_error<1e-5
        pattern_pass=all(r['final_log_activator_sd']>.1 and r['final_jacobian_max_real']<0 and r['distance_to_pattern_control']<1e-4 for r in pattern_rows)
        uniform_pass=all(r['max_log_deviation_from_uniform']<1e-4 and r['final_jacobian_max_real']<0 for r in uniform_rows)
        results[arm]=dict(numerical_pass=numerical,pattern_basin_pass=pattern_pass,uniform_basin_pass=uniform_pass,
            chemical_bistability_supported=bool(numerical and pattern_pass and uniform_pass),
            independent_reference_label=labels[worst],independent_max_log_error=independent_error,
            homogeneous_spectrum=stability(graph,c.signal_beta,c.signal_da,c.signal_dh),trials=rows)
        np.savez_compressed(output/(arm+'.npz'),times=times,trajectories=np.array(paths),initial_states=np.array(starts),
                            labels=np.array(labels),volumes=v,delta=graph.delta,ids=sim.ids)
        for label,path in zip(labels,paths):
            patterned=label.startswith('pattern')
            ax.plot(times,np.std(np.log(path[:,0]),axis=1),color='tab:blue' if patterned else 'tab:orange',
                alpha=.35 if 'perturbed' in label else 1.,label='Developed pattern + perturbations' if label=='pattern_control' else 'Uniform + perturbations' if label=='uniform_control' else None)
        ax.set(title=arm,xlabel='Time on frozen endpoint geometry',ylabel='Across-cell SD of log activator');ax.legend(fontsize=8)
        write_json(output/'results.json',results)
    fig.tight_layout();fig.savefig(output/'bistability.png',dpi=160);plt.close(fig)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/feedback-survival'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-endpoint-bistability'))
    args=p.parse_args();run(args.source,args.output)
