"""Matched single-cell pulse pilot with coevolving chemistry and mechanics."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs, trajectory
from .cell_response import pulse, response_metrics
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .resolution import write_json, _steps


def initialize(checkpoint, initial, target, factor):
    sim=AttributeSimulation.restore(checkpoint)
    if sim.divisions or len(sim.phi)!=sim.config.max_cells:
        raise ValueError('Requires a mature checkpoint without further division')
    values=np.asarray(initial,dtype=float)
    if values.shape!=(2,len(sim.ids)) or not np.isfinite(values).all() or np.any(values<=0):
        raise ValueError('Invalid initial chemical state')
    if target is None:
        if factor!=1.:raise ValueError('A control cannot have a pulse')
        values=values.copy()
    else:
        index=list(sim.ids).index(target);values=pulse(values,index,factor)
    sim.activator,sim.inhibitor=values.copy()
    return sim


def name(job):
    return job['family']+'_'+('control' if job['target'] is None else f"cell-{job['target']}_factor-{job['factor']:g}")


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('Choose a fresh output directory')
    checkpoint=Path('outputs/feedback-survival/switch_on/latest_state.npz').resolve()
    endpoint=Path('outputs/feedback-endpoint-bistability')
    old=json.loads((endpoint/'protocol.json').read_text())
    for path,value in old['sha256'].items():
        if digest(path)!=value:raise ValueError('Frozen source changed: '+path)
    evidence=json.loads((endpoint/'results.json').read_text())['switch_on']
    if not evidence['chemical_bistability_supported']:raise ValueError('Requires accepted endpoint equilibria')
    sim=AttributeSimulation.restore(checkpoint);c=sim.config;g=sim.signaling_graph()
    if abs(sim.time-150)>1e-9 or len(sim.ids)!=16 or sim.attribute_mode!='direct' or not c.feedback or c.signal_transport!='conservative':
        raise ValueError('Requires the full-coupling conservative 16-cell t=150 checkpoint')
    with np.load(endpoint/'switch_on.npz') as d:
        if not np.array_equal(d['ids'],sim.ids) or not np.allclose(g.delta,d['delta'],rtol=1e-12,atol=1e-12):
            raise ValueError('Endpoint graph mismatch')
        initial={k:d['trajectories'][i,-1].copy() for i,k in enumerate(('pattern','uniform'))}
    # Target selection uses initial activator only, never observed response amplitude.
    order=sorted(range(len(sim.ids)),key=lambda i:(initial['pattern'][0,i],int(sim.ids[i])))
    targets=[int(sim.ids[order[0]]),int(sim.ids[order[-1]])]
    jobs=[dict(family=family,target=target,factor=factor) for family in ('pattern','uniform')
          for target,factor in [(None,1.),*[(target,factor) for target in targets for factor in (.9,1.1)]]]
    jobs.sort(key=lambda job: job['target'] is not None)
    duration=60.;interval=.15;every=_steps(interval,c.dt);steps=_steps(duration,c.dt)
    if steps%every:raise ValueError('Align sampling with integration steps')
    output.mkdir(parents=True)
    np.savez_compressed(output/'initial_states.npz',**initial,ids=sim.ids,masses=g.masses)
    times=np.arange(steps//every+1)*interval;fun=rhs(g.delta,c.signal_beta,c.signal_da,c.signal_dh)
    frozen={};frozen_errors={}
    for job in jobs:
        i=None if job['target'] is None else list(sim.ids).index(job['target'])
        x=initial[job['family']] if i is None else pulse(initial[job['family']],i,job['factor'])
        tight=trajectory(fun,x,times,True)
        ref=solve_ivp(fun,(0,duration),x.ravel(),method='Radau',t_eval=times,rtol=1e-11,atol=1e-13)
        if not ref.success or not np.isfinite(ref.y).all() or np.any(ref.y<=0):raise RuntimeError('Frozen reference failed')
        error=float(np.abs(np.log(tight/ref.y.T.reshape(tight.shape))).max())
        if error>=1e-5:raise ValueError('Frozen pilot reference accuracy failed')
        frozen[name(job)]=tight;frozen_errors[name(job)]=error
    np.savez_compressed(output/'frozen_reference.npz',times=times,**frozen)
    files=['cell_response_moving.py','cell_response.py','attribute_development.py','attribute_persistence.py','feedback_long.py',
           'model.py','transport.py','signaling.py','polarity.py','resolution.py','mechanics.py']
    paths=[Path(__file__).with_name(f) for f in files if Path(__file__).with_name(f).exists()]
    inputs=[checkpoint,endpoint/'protocol.json',endpoint/'results.json',endpoint/'switch_on.npz',output/'initial_states.npz',output/'frozen_reference.npz']
    p=dict(seed=7,checkpoint=str(checkpoint),start=sim.time,dt=c.dt,duration=duration,interval=interval,jobs=jobs,
        targets=targets,target_selection='Minimum and maximum initial patterned activator, ties by cell ID; same cells on uniform background.',
        source_config=asdict(c),frozen_reference_errors=frozen_errors,
        criteria=dict(volume_max=.05,radius_min=4.,clipping_max=0.,boundary_max=.01,recovery_ratio=.1,recovery_hold=24.),
        scope='Pilot: one seed-7 full-coupling endpoint geometry, two chemical equilibria, two prespecified cells, +/-10% activator boluses, no further cleavage. Equilibrated chemistry is transplanted onto its source t=150 geometry/polarity. Each family has a matched evolving unperturbed control. Not a fresh zygote or isolated polarity intervention.',
        endpoints='Pulse-control log responses, sampled sustained recovery, final chemical and geometric differences. No static-equilibrium or moving-system stability claim.',
        time_refinement='Main pilot uses dt=0.0075; moving pulse-response timestep refinement is not yet performed.',
        source_sha256={str(f.resolve()):digest(f) for f in paths},input_sha256={str(f.resolve()):digest(f) for f in inputs})
    write_json(output/'protocol.json',p)
    write_json(output/'status.json',dict(state='prepared',protocol_sha256=digest(output/'protocol.json'),completed=0,total=len(jobs)))
    return p


def observe(sim,elapsed):
    g=sim.signaling_graph();m=sim.metrics()
    return dict(elapsed=elapsed,time=sim.time,ids=sim.ids.tolist(),chemistry=np.array([sim.activator,sim.inhibitor]).tolist(),
        volumes=g.masses.tolist(),centers=sim.centers().tolist(),polarity=sim.polarity.tolist(),delta=g.delta.tolist(),
        axis_ratio=m['axis_ratio'],boundary_occupancy=m['boundary_occupancy'],log_activator_sd=float(np.std(np.log(sim.activator))))


def worker(args):
    root,job=args;root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    folder=root/name(job);folder.mkdir(exist_ok=True)
    if (folder/'result.json').exists():
        result=json.loads((folder/'result.json').read_text())
        if result['protocol_sha256']!=ph:raise ValueError('Result protocol mismatch')
        return result
    checkpoint=folder/'latest_state.npz'
    if checkpoint.exists():sim,audit,history=restore_checkpoint(checkpoint,ph,job)
    else:
        with np.load(root/'initial_states.npz') as d:initial=d[job['family']]
        sim=initialize(p['checkpoint'],initial,job['target'],job['factor']);history=[]
        audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    step0=round(p['start']/p['dt']);stop=step0+_steps(p['duration'],p['dt']);every=_steps(p['interval'],p['dt'])
    checkpoint_every=_steps(3.,p['dt']);old=audit['elapsed_seconds'];started=time.monotonic()
    last=history[-1]['time'] if history else -1.
    try:
        while sim.step_number<=stop:
            v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(np.max(abs(v/sim.target-1))))
            audit['min_radius']=min(audit['min_radius'],float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx))
            audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
            if not np.isfinite(v).all() or not np.isfinite(sim.activator).all() or not np.isfinite(sim.inhibitor).all() or np.any(sim.activator<=0) or np.any(sim.inhibitor<=0):
                raise RuntimeError('Nonfinite or nonpositive state')
            if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:
                raise RuntimeError('Numerical quality limit exceeded')
            if (sim.step_number-step0)%every==0 and sim.time>last+1e-9:
                row=observe(sim,sim.time-p['start']);history.append(row);last=sim.time
                audit['max_sampled_boundary']=max(audit['max_sampled_boundary'],row['boundary_occupancy'])
                if audit['max_sampled_boundary']>=.01:raise RuntimeError('Boundary quality limit exceeded')
                audit['elapsed_seconds']=old+time.monotonic()-started
                write_json(folder/'history.json',history)
                write_json(folder/'status.json',dict(state='running',elapsed=row['elapsed'],until=p['duration'],audit=audit))
                if (sim.step_number-step0)%checkpoint_every==0 or sim.step_number==stop:
                    save_checkpoint(sim,checkpoint,audit,history,ph,job)
                    print(f'{name(job)}: elapsed model time {row["elapsed"]:.2f}; wall {audit["elapsed_seconds"]:.1f}s',flush=True)
            if sim.step_number==stop:break
            sim.step()
        result=dict(job=job,quality_pass=True,protocol_sha256=ph,audit=audit)
        write_json(folder/'result.json',result);write_json(folder/'status.json',dict(state='completed',elapsed=p['duration'],audit=audit))
        return result
    except Exception as exc:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(exc),audit=audit));raise


def compare(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    histories={}
    for job in p['jobs']:
        folder=root/name(job);result=json.loads((folder/'result.json').read_text())
        if not result['quality_pass'] or result['protocol_sha256']!=ph:raise ValueError('Unvalidated input')
        histories[name(job)]=json.loads((folder/'history.json').read_text())
    with np.load(root/'initial_states.npz') as d:masses=d['masses'];ids=d['ids']
    rows=[]
    with np.load(root/'frozen_reference.npz') as frozen:
        for job in p['jobs']:
            if job['target'] is None:continue
            key=name(job);control_key=job['family']+'_control';h=histories[key];control=histories[control_key]
            times=np.array([r['elapsed'] for r in h]);t0=np.array([r['elapsed'] for r in control])
            if not np.allclose(times,t0,rtol=0,atol=1e-9):raise ValueError('Control time mismatch')
            if any(not np.array_equal(r['ids'],ids) for r in h+control):raise ValueError('Cell order mismatch')
            cell=list(ids).index(job['target']);path=np.array([r['chemistry'] for r in h]);base=np.array([r['chemistry'] for r in control])
            moving=response_metrics(path,base,times,masses,cell,job['factor'])
            fixed=response_metrics(frozen[key],frozen[control_key],frozen['times'],masses,cell,job['factor'])
            end=h[-1];ref=control[-1]
            geometry=dict(final_centroid_rms=float(np.sqrt(np.sum(masses*np.sum((np.array(end['centers'])-ref['centers'])**2,axis=1))/masses.sum())),
                final_relative_axis_ratio_difference=float(abs(end['axis_ratio']/ref['axis_ratio']-1)),
                final_polarity_rms=float(np.sqrt(np.mean((np.array(end['polarity'])-ref['polarity'])**2))),
                final_relative_operator_difference=float(np.linalg.norm(np.array(end['delta'])-ref['delta'])/np.linalg.norm(ref['delta'])))
            rows.append(dict(job=job,moving=moving,frozen=fixed,geometry=geometry))
    report=dict(scope=p['scope'],trials=rows,time_refinement=p['time_refinement'],
                interpretation='Response relative to a matched moving control; residual differences do not establish another equilibrium or cell identity.')
    write_json(root/'comparison.json',report)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for row in rows:
        label=f"{row['job']['family']} / cell {row['job']['target']} / {row['job']['factor']}"
        for ax,metric in zip(axes,('target_activator_log_auc_per_log_pulse','other_cells_peak_log_rms_per_log_pulse')):
            ax.scatter(row['frozen'][metric],row['moving'][metric],label=label)
    for ax in axes:
        lim=max(ax.get_xlim()[1],ax.get_ylim()[1]);ax.plot([0,lim],[0,lim],':',color='gray');ax.set(xlabel='Frozen geometry',ylabel='Moving geometry')
    axes[0].set_title('Normalized target response integral');axes[1].set_title('Normalized other-cell response');axes[1].legend(fontsize=6)
    fig.tight_layout();fig.savefig(root/'comparison.png',dpi=180);plt.close(fig)
    return report


def run(root,workers=2):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    if json.loads((root/'status.json').read_text())['protocol_sha256']!=ph:raise ValueError('Protocol changed')
    for path,value in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(path)!=value:raise ValueError('Frozen input changed: '+path)
    failures=[];done=0
    write_json(root/'status.json',dict(state='running',protocol_sha256=ph,completed=0,total=len(p['jobs']),workers=workers))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,(str(root),job)):job for job in p['jobs']}
        for future in as_completed(futures):
            try:future.result();done+=1
            except Exception as exc:failures.append(dict(job=futures[future],error=str(exc)))
            write_json(root/'status.json',dict(state='running',protocol_sha256=ph,completed=done,total=len(p['jobs']),failures=failures,workers=workers))
    if not failures:
        try:compare(root)
        except Exception as exc:failures.append(dict(stage='comparison',error=str(exc)))
    write_json(root/'status.json',dict(state='failed' if failures else 'completed',protocol_sha256=ph,completed=done,total=len(p['jobs']),failures=failures))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['prepare','run','compare']);p.add_argument('--output',type=Path,default=Path('outputs/cell-response-moving'));p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    if a.workers<1 or a.workers>4:raise ValueError('Use 1–4 workers')
    if a.command=='prepare':prepare(a.output)
    elif a.command=='run':run(a.output,a.workers)
    else:compare(a.output)
