"""Paired frozen-geometry signal withdrawal and intrinsic-memory controls."""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from .causal_signaling import initial_signals
from .joint_fate import rhs, discrepancies
from .model import Simulation
from .resolution import write_json, _steps


def neutral_fate(initial, elapsed, rate, bistable=True):
    """Exact unforced flow; no integration error is added after withdrawal."""
    initial=np.asarray(initial,float)
    decay=np.exp(-rate*np.asarray(elapsed,float))
    if not bistable:return initial*decay
    return initial/np.sqrt(decay**2+(1-decay**2)*initial**2)


def joint_rhs(state, delta, config):
    base=rhs(state[:3],delta,config,'full')
    relaxing=config.fate_rate*(-state[3]+config.signal_fate_gain*state[0])
    return np.concatenate([base,relaxing[None]],axis=0)


def driven_history(graph, config, seeds, dt, until=120., interval=.6):
    a,b=initial_signals(graph.masses,seeds,.001)
    state=np.stack([a-1,b-1,np.zeros_like(a),np.zeros_like(a)])
    history=[state.copy()];every=_steps(interval,dt)
    exit_rate=float(np.max(-np.diag(graph.delta)))
    substeps=max(1,int(np.ceil(dt*max(1+config.signal_da*exit_rate,config.signal_beta+config.signal_dh*exit_rate)/.2)))
    step=dt/substeps
    for i in range(1,_steps(until,dt)+1):
        for _ in range(substeps):
            stage=state+step*joint_rhs(state,graph.delta,config)
            if np.any(stage[:2]<=-1):raise FloatingPointError('nonpositive regulator stage')
            state=.5*state+.5*(stage+step*joint_rhs(stage,graph.delta,config))
            if not np.isfinite(state).all() or np.any(state[:2]<=-1):raise FloatingPointError('invalid joint state')
        if i%every==0:history.append(state.copy())
    return np.array(history)


def withdrawal_history(driven,times,withdrawal,rate,bistable):
    index=int(np.argmin(abs(times-withdrawal)))
    if abs(times[index]-withdrawal)>1e-9:raise ValueError('withdrawal must align with observations')
    result=driven[:,2 if bistable else 3].copy()
    result[index:]=neutral_fate(result[index],(times[index:]-withdrawal)[:,None,None],rate,bistable)
    return result


def labels(fate,threshold=.55):
    return np.where(fate>threshold,1,np.where(fate<-threshold,-1,0))


def summarize(fate,times,withdrawal,threshold):
    final=labels(fate[-1],threshold);at=labels(fate[np.argmin(abs(times-withdrawal))],threshold)
    committed=at!=0
    first=np.argmax(abs(fate)>threshold,axis=0)
    crossed=np.any(abs(fate)>threshold,axis=0)
    crossing=np.where(crossed,times[first],np.nan)
    return {'both_fates_trials':int(np.sum(np.any(final==1,axis=1)&np.any(final==-1,axis=1))),
            'final_committed_cells':int(np.sum(final!=0)), 'committed_at_withdrawal':int(np.sum(committed)),
            'retained_committed_labels':int(np.sum(committed&(final==at))),
            'newly_committed_cells_after_withdrawal':int(np.sum(~committed&(final!=0))),
            'mean_final_absolute_fate':float(np.mean(abs(fate[-1]))),
            'sampled_first_crossing_time_range':([float(np.nanmin(crossing)),float(np.nanmax(crossing))] if crossed.any() else None),
            'final_counts_per_seed':[{'a':int(np.sum(v==1)),'b':int(np.sum(v==-1)),'uncommitted':int(np.sum(v==0))} for v in final]}


def run(checkpoint,output):
    checkpoint=Path(checkpoint).resolve();output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    sim=Simulation.restore(checkpoint);graph=sim.signaling_graph()
    if sim.divisions or sim.config.signal_transport!='conservative':raise ValueError('requires mature conservative source')
    p={'checkpoint':str(checkpoint),'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
       'seeds':list(range(20)),'withdrawal_times':[0.,.6,3.,6.,12.,24.,48.,60.],
       'until':120.,'interval':.6,'dts':[.0075,.00375],'threshold':sim.config.fate_threshold,
       'criteria':{'signal_rms_max':.01,'fate_difference_max':.05,'both_fates_trials_min':16},
       'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['fate_memory.py','joint_fate.py','causal_signaling.py','model.py','transport.py']},
       'scope':'One frozen mature graph, 20 paired chemical perturbations. Both fate laws receive identical endogenous activator until withdrawal; thereafter input is neutralized to a=1 only in the fate equation. Chemistry continues unchanged. Exact analytic post-withdrawal flow. Same gain/rate and zero initial fate. Continuous exposure is an additional control. Exposure durations are sampled, not estimates of a universal minimum.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p);write_json(output/'status.json',{'state':'running'})
    try:
        times=np.linspace(0,p['until'],_steps(p['until'],p['interval'])+1)
        driven=[]
        for dt in p['dts']:
            driven.append(driven_history(graph,sim.config,p['seeds'],dt,p['until'],p['interval']))
            print(f'driven dt={dt} complete',flush=True)
        np.savez_compressed(output/'driven.npz',times=times,coarse=driven[0],fine=driven[1])
        rows=[];checks=[];fine_fates={}
        for law,bistable in [('bistable',True),('relaxing',False)]:
            for withdrawal in p['withdrawal_times']+[None]:
                trajectories=[d[:,2 if bistable else 3] if withdrawal is None else withdrawal_history(d,times,withdrawal,sim.config.fate_rate,bistable) for d in driven]
                coarse,fine=trajectories
                difference=float(np.max(abs(coarse-fine)));same=bool(np.array_equal(labels(coarse[-1],p['threshold']),labels(fine[-1],p['threshold'])))
                row={'law':law,'withdrawal':withdrawal,'max_step_fate_difference':difference,'final_labels_match':same,
                     **summarize(fine,times,p['until'] if withdrawal is None else withdrawal,p['threshold'])}
                rows.append(row);checks.append(difference<p['criteria']['fate_difference_max'] and same)
                name=f'{law}-'+('continuous' if withdrawal is None else f'withdraw-{withdrawal:g}')
                fine_fates[name]=fine
                np.savez_compressed(output/f'{name}.npz',times=times,coarse=coarse,fine=fine)
        signal=float(np.max(np.sqrt(np.mean((driven[0][:,:2]-driven[1][:,:2])**2,axis=-1))))
        checks_dict={'signal_step_refinement':signal<p['criteria']['signal_rms_max'],'all_fate_step_and_label_checks':all(checks),
                     'neutral_initial_fate_remains_zero':all(r['final_committed_cells']==0 for r in rows if r['withdrawal']==0),
                     'source_checkpoint_unchanged':hashlib.sha256(checkpoint.read_bytes()).hexdigest()==p['checkpoint_sha256']}
        result={'protocol':p,'checks':checks_dict,'passed':all(checks_dict.values()),'max_signal_rms_difference':signal,'rows':rows}
        write_json(output/'comparison.json',result)
        lines=['# Signal withdrawal and fate memory','',f'Numerical and neutral-control checks pass: **{result["passed"]}**','',
               '| Fate law | Withdraw at | Both labels / 20 | Committed at withdrawal | Retained labels | New commitments | Mean final absolute fate |',
               '|---|---:|---:|---:|---:|---:|---:|']
        for r in rows:lines.append(f'| {r["law"]} | {r["withdrawal"] if r["withdrawal"] is not None else "continuous"} | {r["both_fates_trials"]} | {r["committed_at_withdrawal"]} | {r["retained_committed_labels"]} | {r["newly_committed_cells_after_withdrawal"]} | {r["mean_final_absolute_fate"]:.6g} |')
        lines+=['',f'Maximum signal timestep discrepancy: {signal:.6g}.',f'Maximum fate timestep discrepancy: {max(r["max_step_fate_difference"] for r in rows):.6g}.','',
                'Neutralization removes activator forcing of fate, not the underlying chemical pattern. Persistent identity is expected from the assumed bistable law; the relaxing control has exact exponential decay after withdrawal. Zero initial fate stays zero without input. Any nonzero bistable fate eventually approaches its signed stable state in exact arithmetic, so apparent exposure thresholds depend on the observation horizon. This is not evidence that sustained Turing patterning is necessary, irreversible biological commitment, or robustness to noise or reversal.',
                'Continuous-exposure rows evaluate commitment at the final observation, so retention/new-commitment counts there do not constitute a withdrawal assay. First threshold crossings are sampled and need not be permanent.']
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        for ax,law in zip(axes,['bistable','relaxing']):
            for duration in [0.,.6,6.,24.,60.]:
                fate=fine_fates[f'{law}-withdraw-{duration:g}']
                ax.plot(times,np.mean(abs(fate)>p['threshold'],axis=(1,2)),label=f'withdraw {duration:g}')
            ax.plot(times,np.mean(abs(fine_fates[f'{law}-continuous'])>p['threshold'],axis=(1,2)),'k--',label='continuous')
            ax.set(title=law,xlabel='Time since signal reset',ylabel='Committed fraction (320 cells)',ylim=(-.03,1.03));ax.grid(alpha=.2)
        axes[1].legend(fontsize=8);fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)
        write_json(output/'status.json',{'state':'completed','passed':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{'state':'failed','error':str(error)});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--checkpoint',type=Path,default=Path('outputs/development-refinement/space-72/state-18.npz'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.checkpoint,args.output)
