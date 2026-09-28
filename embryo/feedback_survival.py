"""Switch direct mechanics on at developmental t=90 and test pattern survival."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from .attribute_development import AttributeSimulation, attributes
from .attribute_persistence import rhs, trajectory, distance
from .feedback_spectrum import summarize
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .resolution import write_json, _steps

ARMS=('switch_on','keep_off')


def initialize(checkpoint,arm):
    if arm not in ARMS:raise ValueError('Unknown feedback-switch arm')
    sim=AttributeSimulation.restore(checkpoint)
    if sim.attribute_mode!='no_feedback' or sim.config.feedback or sim.divisions or len(sim.phi)!=sim.config.max_cells:
        raise ValueError('Requires mature no-feedback attribute checkpoint')
    if arm=='switch_on':
        sim.attribute_mode='direct';sim.config.feedback=True
    # No other state, coefficient, random stream, or clock is changed.
    return sim


def similarity(current,initial,volumes):
    """Descriptive cell-ID matched chemical memory, with fixed reference weights."""
    weights=volumes/volumes.sum();x=np.log(current[0]);y=np.log(initial[0])
    x=x-np.dot(weights,x);y=y-np.dot(weights,y)
    denom=np.sqrt(np.dot(weights,x*x)*np.dot(weights,y*y))
    corr=float(np.dot(weights,x*y)/denom) if denom>1e-12 else None
    return dict(initial_log_activator_correlation=corr,
                initial_chemical_log_rms=float(distance(current,initial,volumes)))


def prepare(checkpoint,output,duration=60.,interval=.6):
    checkpoint=Path(checkpoint).resolve();output=Path(output)
    if output.exists():raise FileExistsError('Choose a fresh output directory')
    sim=initialize(checkpoint,'keep_off');c=sim.config
    if abs(sim.time-90)>1e-9 or len(sim.phi)!=16 or c.signal_transport!='conservative':
        raise ValueError('This protocol requires the sixteen-cell conservative t=90 state')
    if _steps(duration,c.dt)%_steps(interval,c.dt):raise ValueError('Sample times must align')
    initial=np.array([sim.activator,sim.inhibitor]);volumes=sim.volumes()
    spread=float(np.std(np.log(initial[0])))
    if spread<=.1:raise ValueError('Source has insufficient chemical contrast')
    graph=sim.signaling_graph();fun=rhs(graph.delta,c.signal_beta,c.signal_da,c.signal_dh)
    times=np.linspace(0,duration,round(duration/interval)+1)
    coarse=trajectory(fun,initial,times);reference=trajectory(fun,initial,times,True)
    error=float(np.abs(np.log(coarse/reference)).max())
    if error>=1e-5:raise ValueError('Frozen reference tolerance check failed')
    output.mkdir(parents=True)
    np.savez_compressed(output/'frozen_reference.npz',times=times+sim.time,trajectory=reference,
        initial=initial,volumes=volumes,delta=graph.delta,ids=sim.ids)
    files=['feedback_survival.py','feedback_long.py','attribute_development.py','attribute_persistence.py',
           'feedback_spectrum.py','model.py','transport.py','signaling.py','polarity.py','resolution.py']
    p=dict(checkpoint=str(checkpoint),checkpoint_sha256=digest(checkpoint),start=sim.time,
        duration=duration,interval=interval,dt=c.dt,late_duration=min(15.,duration/2),arms=list(ARMS),
        source_config=asdict(c),initial_log_activator_sd=spread,
        initial_spectrum=summarize(sim)[0],frozen_reference_max_log_error=error,
        frozen_reference_sha256=digest(output/'frozen_reference.npz'),
        criteria=dict(late_min_log_activator_sd=.1,late_min_correlation=.8,
                      volume_max=.05,radius_min=4.,boundary_max=.01,clipping_max=0.),
        scope='Exact no-feedback developmental t=90 state. Only feedback flag/mode changes in switch arm. No chemical reset, frozen preconditioning, new noise, or coefficient sweep. Single developed embryo; sixty-unit survival test, not full developmental causality or proof of cell identity.',
        source_sha256={str(Path(__file__).with_name(f).resolve()):digest(Path(__file__).with_name(f)) for f in files})
    write_json(output/'protocol.json',p);write_json(output/'status.json',dict(state='prepared',protocol_sha256=digest(output/'protocol.json')))
    return p


def worker(args):
    root,arm=args;root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    folder=root/arm;folder.mkdir(exist_ok=True);job=dict(experiment='feedback_survival',arm=arm)
    if (folder/'result.json').exists():return json.loads((folder/'result.json').read_text())
    if (folder/'latest_state.npz').exists():
        sim,audit,history=restore_checkpoint(folder/'latest_state.npz',ph,job)
    else:
        sim=initialize(p['checkpoint'],arm);history=[]
        audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    with np.load(root/'frozen_reference.npz') as f:initial=f['initial'];weights=f['volumes']
    last=history[-1]['metrics']['time'] if history else -1.
    started=time.monotonic();old_elapsed=audit['elapsed_seconds'];step0=round(p['start']/p['dt'])
    stop=round((p['start']+p['duration'])/p['dt']);every=round(p['interval']/p['dt'])
    try:
        while sim.step_number<=stop:
            v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(np.max(abs(v/sim.target-1))))
            audit['min_radius']=min(audit['min_radius'],float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx))
            audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
            if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:
                raise RuntimeError('Numerical quality limit exceeded')
            if (sim.step_number-step0)%every==0 and sim.time>last+1e-9:
                row=attributes(sim);row['spectrum']=summarize(sim)[0]
                row['similarity']=similarity(np.array([sim.activator,sim.inhibitor]),initial,weights)
                history.append(row);last=sim.time
                audit['max_sampled_boundary']=max(audit['max_sampled_boundary'],row['metrics']['boundary_occupancy'])
                if audit['max_sampled_boundary']>=.01:raise RuntimeError('Boundary quality limit exceeded')
                audit['elapsed_seconds']=old_elapsed+time.monotonic()-started
                write_json(folder/'history.json',history)
                write_json(folder/'status.json',dict(state='running',time=sim.time,until=p['start']+p['duration'],audit=audit))
                if (sim.step_number-step0)%(10*every)==0 or sim.step_number==stop:
                    save_checkpoint(sim,folder/'latest_state.npz',audit,history,ph,job)
                    print(f'{arm}: t={sim.time:g}, elapsed={audit["elapsed_seconds"]:.1f}s',flush=True)
            if sim.step_number==stop:break
            sim.step()
        late=[r for r in history if r['metrics']['time']>=p['start']+p['duration']-p['late_duration']-1e-9]
        spreads=[r['attribute_std'][0] for r in late];corr=[r['similarity']['initial_log_activator_correlation'] for r in late]
        survives=min(spreads)>.1;ordering=all(x is not None and x>=.8 for x in corr)
        result=dict(arm=arm,quality_pass=True,audit=audit,late_min_log_activator_sd=min(spreads),
            late_min_initial_correlation=min(corr) if all(x is not None for x in corr) else None,
            contrast_survives=survives,initial_ordering_retained=ordering,
            outcome='contrast_and_ordering_retained' if survives and ordering else 'contrast_survives_but_reorganized' if survives else 'contrast_not_sustained',
            final_spectrum=history[-1]['spectrum'])
        write_json(folder/'result.json',result);write_json(folder/'status.json',dict(state='completed',time=sim.time,audit=audit));return result
    except Exception as error:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(error),audit=audit));raise


def compare(root,results):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());results={r['arm']:r for r in results}
    with np.load(root/'frozen_reference.npz') as f:
        t=f['times'];spread=np.std(np.log(f['trajectory'][:,0]),axis=1)
    late=t>=p['start']+p['duration']-p['late_duration']-1e-9
    a=results['switch_on']['contrast_survives'];b=results['keep_off']['contrast_survives']
    interpretation=('Switch and moving no-feedback control both retain chemical contrast.' if a and b else
                    'Switch loses contrast while the moving no-feedback control retains it.' if not a and b else
                    'Both moving arms lose sustained contrast; loss is not specific to switching feedback.' if not a and not b else
                    'Switch retains contrast while the moving no-feedback control loses it.')
    report=dict(arms=results,frozen_late_min_log_activator_sd=float(spread[late].min()),interpretation=interpretation,
                scope='Finite-horizon chemical survival and cell-ID association; not proof of all cell attributes or identity.')
    write_json(root/'comparison.json',report)
    lines=['# Moving-geometry feedback-switch survival','',''+interpretation,'',
           '| Arm | Late minimum log-activator SD | Minimum initial-state correlation | Outcome |','|---|---:|---:|---|']
    lines += [f'| {k} | {v["late_min_log_activator_sd"]:.6g} | {v["late_min_initial_correlation"]} | {v["outcome"]} |' for k,v in results.items()]
    (root/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for arm in ARMS:
        h=json.loads((root/arm/'history.json').read_text());times=[r['metrics']['time'] for r in h]
        axes[0].plot(times,[r['attribute_std'][0] for r in h],label=arm)
        axes[1].plot(times,[r['similarity']['initial_log_activator_correlation'] for r in h],label=arm)
        axes[2].plot(times,[r['spectrum']['maximum_spatial_growth'] for r in h],label=arm)
    axes[0].plot(t,spread,'--',label='frozen reference');axes[0].axhline(.1,color='grey',ls=':')
    axes[1].axhline(.8,color='grey',ls=':');axes[2].axhline(0,color='grey',ls=':')
    for ax,label in zip(axes,['SD of log activator','Correlation with initial cell states','Frozen homogeneous growth rate']):
        ax.set(xlabel='Model time',ylabel=label);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(root/'comparison.png',dpi=160);plt.close(fig)


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());status=json.loads((root/'status.json').read_text());ph=digest(root/'protocol.json')
    if ph!=status['protocol_sha256']:raise ValueError('Protocol changed')
    for path,value in {**p['source_sha256'],p['checkpoint']:p['checkpoint_sha256'],str(root/'frozen_reference.npz'):p['frozen_reference_sha256']}.items():
        if digest(path)!=value:raise ValueError('Frozen input changed: '+path)
    write_json(root/'status.json',dict(state='running',protocol_sha256=ph,completed=0,total=2))
    results=[];failures=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs={pool.submit(worker,(str(root),arm)):arm for arm in ARMS}
        for future in as_completed(jobs):
            try:results.append(future.result())
            except Exception as error:failures.append(dict(arm=jobs[future],error=str(error)))
            write_json(root/'status.json',dict(state='running',protocol_sha256=ph,completed=len(results),total=2,failures=failures))
    if not failures:compare(root,results)
    write_json(root/'status.json',dict(state='failed' if failures else 'completed',protocol_sha256=ph,completed=len(results),total=2,failures=failures))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['prepare','run'])
    p.add_argument('--checkpoint',type=Path,default=Path('outputs/attribute-development/no_feedback/state-90.npz'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-survival'))
    p.add_argument('--duration',type=float,default=60.);p.add_argument('--interval',type=float,default=.6)
    args=p.parse_args()
    if args.command=='prepare':prepare(args.checkpoint,args.output,args.duration,args.interval)
    else:run(args.output)
