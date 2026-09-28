"""Long matched mechanical-component controls for initiation and persistence."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from .attribute_development import AttributeSimulation, attributes
from .attribute_persistence import rhs, trajectory
from .causal_signaling import initial_signals
from .feedback_spectrum import summarize
from .resolution import write_json, _steps

ARMS={'baseline':(0.,0.,0.),'tension':(.25,0.,0.),'adhesion':(0.,.35,0.),
      'polarity':(0.,0.,.35),'full':(.25,.35,.35)}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def initialize(checkpoint,signals,arm):
    sim=AttributeSimulation.restore(checkpoint)
    if arm not in ARMS or sim.divisions or len(sim.phi)!=sim.config.max_cells:
        raise ValueError('Requires a known arm and completed cell divisions')
    sim.attribute_mode='direct';sim.config.feedback=True
    sim.config.fate_tension,sim.config.fate_adhesion,sim.config.polarity_tension=ARMS[arm]
    sim.activator,sim.inhibitor=np.asarray(signals,dtype=float).copy()
    # All arms retain identical source polarity; only its mechanical contribution varies.
    sim.config.validate()
    return sim


def prepare(checkpoint,output,duration=60.,interval=.6):
    checkpoint=Path(checkpoint).resolve();output=Path(output)
    if output.exists():raise FileExistsError('Choose a fresh output directory')
    sim=AttributeSimulation.restore(checkpoint);c=sim.config
    if sim.divisions or len(sim.phi)!=c.max_cells or c.signal_transport!='conservative':
        raise ValueError('Requires mature conservative attribute checkpoint')
    steps=_steps(duration,c.dt);every=_steps(interval,c.dt)
    if steps%every or duration<interval:raise ValueError('Align duration and sample interval')
    a,b=initial_signals(sim.volumes(),[7],.001);near=np.stack([a[0],b[0]])
    graph=sim.signaling_graph();fun=rhs(graph.delta,c.signal_beta,c.signal_da,c.signal_dh)
    times=np.linspace(0,240,481)
    coarse=trajectory(fun,near,times);tight=trajectory(fun,near,times,True)
    error=float(np.abs(np.log(coarse/tight)).max());residual=float(np.abs(fun(240,tight[-1].ravel())).max())
    spread=float(np.std(np.log(tight[-1,0])))
    if error>=1e-5 or residual>=1e-6 or spread<=.1:
        raise ValueError(f'Frozen patterned-start gate failed: error={error}, residual={residual}, spread={spread}')
    spectrum,_=summarize(sim)
    if spectrum['maximum_spatial_growth']<=0:raise ValueError('Initiation geometry must support a growing mode')
    output.mkdir(parents=True)
    np.savez_compressed(output/'initial_states.npz',formation=near,persistence=tight[-1],
                        precondition_times=times,precondition_trajectory=tight,cell_ids=sim.ids)
    files=['feedback_long.py','attribute_development.py','attribute_persistence.py','feedback_spectrum.py',
           'causal_signaling.py','model.py','transport.py','signaling.py','polarity.py','resolution.py']
    p=dict(checkpoint=str(checkpoint),checkpoint_sha256=digest(checkpoint),start=sim.time,
        duration=duration,interval=interval,dt=c.dt,late_duration=min(15.,duration/2),seed=7,
        arms=ARMS,jobs=[dict(state=state,arm=arm) for arm in ['baseline','full','tension','adhesion','polarity'] for state in ['formation','persistence']],
        initial_states_sha256=digest(output/'initial_states.npz'),initial_spectrum=spectrum,
        preconditioning=dict(duration=240.,max_log_error=error,final_rhs_max=residual,final_log_activator_sd=spread),
        criteria=dict(volume_max=.05,radius_min=4.,boundary_max=.01,clipping_max=0.,late_log_activator_sd_min=.1),
        scope='One matched 16-cell geometry and perturbation seed; no cleavage or new zygote ensemble. Formation and persistence differ only in chemical initial state, with equal starting geometry/polarity/lineage. Persistence is preconditioned on the frozen graph, not developed in the moving system. Original linear constitutive law, no imposed fate switch.',
        source_sha256={str(Path(__file__).with_name(f).resolve()):digest(Path(__file__).with_name(f)) for f in files})
    write_json(output/'protocol.json',p);write_json(output/'status.json',dict(state='prepared',protocol_sha256=digest(output/'protocol.json')))
    return p


def save_checkpoint(sim,path,audit,history,protocol_hash,job):
    raw=path.with_name(path.stem+'.raw.npz');tmp=path.with_name(path.stem+'.tmp.npz')
    sim.checkpoint(raw)
    with np.load(raw,allow_pickle=False) as data:payload={k:data[k].copy() for k in data.files}
    payload['long_experiment']=np.array(json.dumps(dict(audit=audit,history=history,protocol_hash=protocol_hash,job=job)))
    np.savez_compressed(tmp,**payload);tmp.replace(path);raw.unlink()


def restore_checkpoint(path,protocol_hash,job):
    with np.load(path,allow_pickle=False) as data:meta=json.loads(str(data['long_experiment']))
    if meta['protocol_hash']!=protocol_hash or meta['job']!=job:raise ValueError('Checkpoint protocol/job mismatch')
    return AttributeSimulation.restore(path),meta['audit'],meta['history']


def worker(args):
    output,job=args;output=Path(output);p=json.loads((output/'protocol.json').read_text());ph=digest(output/'protocol.json')
    name=job['state']+'_'+job['arm'];folder=output/name;folder.mkdir(exist_ok=True)
    result_file=folder/'result.json'
    if result_file.exists():return json.loads(result_file.read_text())
    checkpoint=folder/'latest_state.npz'
    if checkpoint.exists():sim,audit,history=restore_checkpoint(checkpoint,ph,job)
    else:
        with np.load(output/'initial_states.npz') as data:signals=data[job['state']]
        sim=initialize(p['checkpoint'],signals,job['arm']);history=[]
        audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    start_clock=time.monotonic();base_elapsed=audit['elapsed_seconds'];last_saved=history[-1]['metrics']['time'] if history else -1.
    stop=round((p['start']+p['duration'])/p['dt']);every=round(p['interval']/p['dt']);initial_step=round(p['start']/p['dt'])
    try:
        while sim.step_number<=stop:
            v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(np.max(abs(v/sim.target-1))))
            audit['min_radius']=min(audit['min_radius'],float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx))
            audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
            if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:
                raise RuntimeError('Numerical quality limit exceeded; inspect audit')
            if (sim.step_number-initial_step)%every==0 and sim.time>last_saved+1e-9:
                row=attributes(sim);row['spectrum']=summarize(sim)[0];history.append(row);last_saved=sim.time
                audit['max_sampled_boundary']=max(audit['max_sampled_boundary'],row['metrics']['boundary_occupancy'])
                if audit['max_sampled_boundary']>=.01:raise RuntimeError('Boundary quality limit exceeded')
                audit['elapsed_seconds']=base_elapsed+time.monotonic()-start_clock
                write_json(folder/'history.json',history)
                write_json(folder/'status.json',dict(state='running',time=sim.time,until=p['start']+p['duration'],audit=audit))
                # An atomic state bundle holds its own matching history and audit for restart.
                if (sim.step_number-initial_step)%(10*every)==0 or sim.step_number==stop:
                    save_checkpoint(sim,checkpoint,audit,history,ph,job)
                    print(f'{name}: t={sim.time:g}, elapsed={audit["elapsed_seconds"]:.1f}s',flush=True)
            if sim.step_number==stop:break
            sim.step()
        late=[r for r in history if r['metrics']['time']>=p['start']+p['duration']-p['late_duration']-1e-9]
        spreads=[r['attribute_std'][0] for r in late]
        result=dict(job=job,audit=audit,completed=True,quality_pass=True,
            late_min_log_activator_sd=min(spreads),late_mean_log_activator_sd=float(np.mean(spreads)),
            persistent=min(spreads)>.1,late_mean_axis_ratio=float(np.mean([r['metrics']['axis_ratio'] for r in late])),
            final_spectrum=history[-1]['spectrum'])
        write_json(result_file,result);write_json(folder/'status.json',dict(state='completed',time=sim.time,audit=audit))
        return result
    except Exception as error:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(error),audit=audit));raise


def compare(output,rows):
    write_json(output/'comparison.json',dict(results=rows,all_quality_pass=all(r['quality_pass'] for r in rows)))
    lines=['# Long matched feedback controls','','Single-geometry pilot; formation and persistence use matched geometry and different initial chemistry.','',
           '| Initial chemistry | Mechanical arm | Late minimum log-activator SD | Persistent | Late mean axis ratio |',
           '|---|---|---:|---|---:|']
    lines += [f'| {r["job"]["state"]} | {r["job"]["arm"]} | {r["late_min_log_activator_sd"]:.6g} | {r["persistent"]} | {r["late_mean_axis_ratio"]:.6f} |' for r in rows]
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for col,state in enumerate(['formation','persistence']):
        for arm in ARMS:
            h=json.loads((output/(state+'_'+arm)/'history.json').read_text());t=[r['metrics']['time'] for r in h]
            axes[0,col].plot(t,[r['attribute_std'][0] for r in h],label=arm)
            axes[1,col].plot(t,[r['spectrum']['maximum_spatial_growth'] for r in h],label=arm)
        axes[0,col].set(title=state,ylabel='Across-cell SD of log activator');axes[0,col].axhline(.1,color='grey',ls='--')
        axes[1,col].set(xlabel='Model time',ylabel='Frozen homogeneous growth rate');axes[1,col].axhline(0,color='grey',ls='--')
        axes[0,col].legend(fontsize=8)
    fig.tight_layout();fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)


def run(output,workers=2):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());status=json.loads((output/'status.json').read_text());ph=digest(output/'protocol.json')
    if ph!=status['protocol_sha256']:raise ValueError('Protocol changed')
    for path,expected in {**p['source_sha256'],p['checkpoint']:p['checkpoint_sha256'],str(output/'initial_states.npz'):p['initial_states_sha256']}.items():
        if digest(path)!=expected:raise ValueError(f'Frozen input changed: {path}')
    completed=[];failures=[]
    write_json(output/'status.json',dict(state='running',protocol_sha256=ph,total=len(p['jobs']),completed=0))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,(str(output),job)):job for job in p['jobs']}
        for future in as_completed(futures):
            try:completed.append(future.result())
            except Exception as error:failures.append(dict(job=futures[future],error=str(error)))
            write_json(output/'status.json',dict(state='running',protocol_sha256=ph,total=len(p['jobs']),completed=len(completed),failures=failures))
    if not failures:compare(output,completed)
    write_json(output/'status.json',dict(state='failed' if failures else 'completed',protocol_sha256=ph,total=len(p['jobs']),completed=len(completed),failures=failures))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['prepare','run'])
    parser.add_argument('--output',type=Path,default=Path('outputs/feedback-long'))
    parser.add_argument('--checkpoint',type=Path,default=Path('outputs/attribute-development/no_feedback/state-18.npz'))
    parser.add_argument('--duration',type=float,default=60.);parser.add_argument('--interval',type=float,default=.6)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    if args.command=='prepare':prepare(args.checkpoint,args.output,args.duration,args.interval)
    else:run(args.output,args.workers)
