"""Snapshot-based signaling replay; fidelity is tested before causal attribution."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .joint_fate import advance, discrepancies
from .moving_causal import JointMovingSimulation
from .resolution import write_json, _steps

ARMS = ('frozen', 'replay', 'no_dilution', 'frozen_transport')


def geometry(times, conductance, volumes, t):
    index = int(np.clip(np.searchsorted(times, t, side='right') - 1, 0, len(times)-2))
    fraction = np.clip((t-times[index])/(times[index+1]-times[index]), 0, 1)
    g = (1-fraction)*conductance[index] + fraction*conductance[index+1]
    v = (1-fraction)*volumes[index] + fraction*volumes[index+1]
    k = np.diag(g.sum(axis=1))-g
    return -k/v[:, None], v


def trajectory(times, conductance, volumes, initial, config, arm, dt, observation_times):
    if arm not in ARMS:
        raise ValueError('unknown replay arm')
    state = initial.copy()
    frozen, _ = geometry(times, conductance, volumes, times[0])
    every = _steps(observation_times[1]-observation_times[0], dt)
    steps = _steps(times[-1]-times[0], dt)
    result = [state.copy()]
    for step in range(steps):
        t = times[0]+step*dt
        delta, before = geometry(times, conductance, volumes, t)
        _, after = geometry(times, conductance, volumes, t+dt)
        if arm in ('frozen', 'frozen_transport'):
            delta = frozen
        state = advance(state, delta, config, 'full', dt)
        if arm in ('replay', 'frozen_transport'):
            for species in (0, 1):
                state[species] = state[species]*(before/after)+(before-after)/after
        if not np.isfinite(state).all() or np.any(state[:2] <= -1):
            raise FloatingPointError('replay lost finite positive regulators')
        if (step+1) % every == 0:
            result.append(state.copy())
    return np.asarray(result)


def spectral_history(times, conductance, volumes, config):
    rows = []
    jacobian = np.array([[1., -1.], [2*config.signal_beta, -config.signal_beta]])
    for t, g, v in zip(times, conductance, volumes):
        k = np.diag(g.sum(axis=1))-g
        values = np.linalg.eigvalsh(k/np.sqrt(v[:,None]*v[None,:]))
        rates = [float(np.linalg.eigvals(jacobian-l*np.diag([config.signal_da,config.signal_dh])).real.max()) for l in values]
        spatial = [r for l,r in zip(values,rates) if l>1e-10]
        rows.append({'time':float(t),'eigenvalues':values.tolist(),'growth_rates':rates,
                     'unstable_modes':sum(r>0 for r in spatial), 'max_spatial_rate':max(spatial) if spatial else None})
    rates = [r['max_spatial_rate'] for r in rows]
    integral = cumulative_trapezoid(rates, times, initial=0).tolist() if all(r is not None for r in rates) else None
    return {'rows':rows, 'integrated_instantaneous_max_rate':integral,
            'interpretation':'Trapezoidal integral of frozen-equilibrium maximum rates only. Not actual amplification: excludes mode rotation, dilution, finite departures from equilibrium, and changing eigenvectors.'}


def run(source, output):
    source=Path(source).resolve(); output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    source_status=json.loads((source/'status.json').read_text())
    if source_status['state']!='completed' or not source_status['all_quality_pass']:
        raise ValueError('requires completed moving controls with passing quality')
    paths=sorted((source/'full').glob('state-*.npz'),key=lambda p:float(p.stem.split('-')[1]))
    if len(paths)!=11:raise ValueError('requires eleven t=18..78 source checkpoints')
    p={'source':str(source),'arms':list(ARMS),'dts':[.0075,.00375], 'snapshot_strides':[1,2],
       'criteria':{'signal_rms_max':.01,'fate_difference_max':.05,'persistent_contrast_min':.1},
       'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[source/'full/history.json',source/'comparison.json']},
       'source_sha256':{str(Path(__file__).with_name(f).resolve()):hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in ['geometry_replay.py','joint_fate.py','model.py','moving_causal.py','transport.py']},
       'scope':'One seed and geometry history. Piecewise linear symmetric conductances and positive measured volumes from six-unit checkpoints; compare twelve-unit decimation. This is approximate prescribed replay, not exact per-step geometry. No-dilution and frozen-transport/variable-volume arms are intentionally nonconservative diagnostic interventions.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'running'})
    try:
        simulations=[JointMovingSimulation.restore(f) for f in paths]
        first=simulations[0];times=np.array([s.time for s in simulations])
        np.testing.assert_allclose(times,np.arange(18,79,6),atol=1e-9)
        for sim in simulations:
            np.testing.assert_array_equal(sim.ids,first.ids)
            if sim.divisions:raise ValueError('replay requires fixed cell identities')
        conductance=np.array([sim.signaling_graph().weights for sim in simulations])
        volumes=np.array([sim.volumes() for sim in simulations])
        initial=first.joint_state.copy();config=first.config
        del simulations
        history=json.loads((source/'full/history.json').read_text())
        observations=np.array([r['time'] for r in history])
        np.testing.assert_allclose(observations,np.linspace(18,78,101),atol=1e-9)
        recorded=np.array([[np.array(r['activator'])-1,np.array(r['inhibitor'])-1,r['fate']] for r in history])
        np.savez_compressed(output/'geometry.npz',times=times,conductance=conductance,volumes=volumes,initial=initial,ids=first.ids)
        spectra=spectral_history(times,conductance,volumes,config);write_json(output/'spectra.json',spectra)
        rows=[];summaries={};runs={}
        for arm in ARMS:
            for stride in p['snapshot_strides']:
                for dt in p['dts']:
                    state=trajectory(times[::stride],conductance[::stride],volumes[::stride],initial,config,arm,dt,observations)
                    runs[arm,stride,dt]=state
                    np.savez_compressed(output/f'{arm}-spacing-{6*stride}-dt-{dt}.npz',times=observations,state=state)
            coarse=runs[arm,1,p['dts'][0]];fine=runs[arm,1,p['dts'][1]]
            row={'arm':arm,'step_halving':discrepancies(coarse,fine), 'snapshot_decimation':discrepancies(fine,runs[arm,2,p['dts'][1]])}
            rows.append(row)
            contrast=np.std(fine[:,0],axis=1);late=observations>=63-1e-10
            labels=np.where(fine[-1,2]>.55,1,np.where(fine[-1,2]<-.55,-1,0))
            capacities=np.repeat(volumes[:1],len(observations),axis=0) if arm=='frozen' else np.array([geometry(times,conductance,volumes,t)[1] for t in observations])
            amounts=np.sum((fine[:,:2]+1)*capacities[:,None,:],axis=2)
            summaries[arm]={'late_contrast_min':float(contrast[late].min()),'late_contrast_max':float(contrast[late].max()),
                            'final_a':int(sum(labels==1)),'final_b':int(sum(labels==-1)),'final_uncommitted':int(sum(labels==0)),
                            'final_amounts':amounts[-1].tolist()}
            write_json(output/'progress.json',{'rows':rows,'summary':summaries});print(arm,row,flush=True)
        c=p['criteria']
        acceptable=lambda d:d['max_signal_rms']<c['signal_rms_max'] and d['max_fate_difference']<c['fate_difference_max'] and d['final_labels_match']
        fidelity=discrepancies(runs['replay',1,p['dts'][1]],recorded)
        checks={'all_chemical_step_checks':all(acceptable(r['step_halving']) for r in rows),
                'all_snapshot_decimation_checks':all(acceptable(r['snapshot_decimation']) for r in rows),
                'replay_matches_recorded_full':acceptable(fidelity)}
        result={'protocol':p,'checks':checks,'attribution_ready':all(checks.values()),'rows':rows,'replay_fidelity':fidelity,'summary':summaries,
                'limitation':'Decimation agreement alone is not proof of sufficient snapshot density; fidelity to recorded full chemistry/fate is also required. Frozen operator with changing volumes need not conserve volume-weighted amount. Reactions change total amounts in every arm.'}
        write_json(output/'comparison.json',result)
        lines=['# Geometry replay diagnostic','',f'All attribution prerequisites pass: **{result["attribution_ready"]}**','',
               '| Arm | Late contrast min–max | Final A / B / uncommitted |','|---|---:|---:|']
        for arm,r in summaries.items():lines.append(f'| {arm} | {r["late_contrast_min"]:.6g}–{r["late_contrast_max"]:.6g} | {r["final_a"]} / {r["final_b"]} / {r["final_uncommitted"]} |')
        lines+=['',*[f'- {k}: {v}' for k,v in checks.items()], '',f'Replay versus recorded trajectory: {fidelity}', '',result['limitation']]
        lines+=['','If prerequisites fail, do not attribute the observed live suppression to a replay intervention. Capture denser geometry and repeat before drawing a mechanism conclusion.']
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
        for arm in ARMS:ax.plot(observations,np.std(runs[arm,1,p['dts'][1]][:,0],axis=1),label=arm)
        ax.plot(observations,np.std(recorded[:,0],axis=1),'k--',label='recorded full')
        ax.axhline(.1,color='gray',ls=':');ax.axvspan(63,78,alpha=.08,color='gray')
        ax.set(xlabel='Developmental time',ylabel='Cell-to-cell activator standard deviation',title='Approximate geometry replay; interpret only after fidelity checks');ax.legend();fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)
        write_json(output/'status.json',{'state':'completed','attribution_ready':result['attribution_ready']})
        return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,default=Path('outputs/moving-causal'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.source,args.output)
