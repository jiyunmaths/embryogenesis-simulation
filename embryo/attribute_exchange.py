"""Exchange chemical states on a fixed developed graph, without identity labels."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from .attribute_persistence import distance, perturb, rhs, trajectory


def exchange(values, masses, i, j):
    """Swap pair concentrations, then rescale only the pair to conserve amounts."""
    out = values.copy()
    pair = [i, j]
    swapped = values[:, [j, i]]
    scale = (values[:, pair] @ masses[pair]) / (swapped @ masses[pair])
    out[:, pair] = swapped * scale[:, None]
    return out


def classify(destination_ratio, donor_ratio, stationary, informative):
    if not informative:
        return 'uninformative'
    if not stationary:
        return 'unsettled'
    if destination_ratio < .1:
        return 'destination'
    if donor_ratio < .1:
        return 'transferred'
    return 'reorganized'


def run(source, development, output):
    data = np.load(source/'no_feedback.npz')
    initial = data['control'][-1]
    masses, delta, ids = data['masses'], data['delta'], data['cell_ids']
    config = json.loads((development/'protocol.json').read_text())['config']
    history = json.loads((development/'no_feedback/history.json').read_text())
    # The exported history is a list of attribute snapshots.
    cells = {cell['id']:cell for cell in history[-1]['cells']}
    exposure = np.array([cells[int(i)]['context']['exposure'] for i in ids])
    context = np.column_stack([exposure, -np.diag(delta), masses])
    normalized = (context-context.mean(axis=0))/np.maximum(context.std(axis=0),1e-12)
    pairs = list(combinations(range(len(masses)),2))
    context_distances = np.array([np.linalg.norm(normalized[i]-normalized[j]) for i,j in pairs])
    cutoff = float(np.median(context_distances))
    files = [Path(__file__),Path(__file__).with_name('attribute_persistence.py'),
             source/'no_feedback.npz',development/'protocol.json',development/'no_feedback/history.json']
    protocol = dict(duration=1200,late_start=1000,sample_interval=2,
        pair_selection='All unordered pairs; contrasting context means at least median standardized context distance.',
        context_features=['exposure','transport_exit_rate','volume'],context_distance_cutoff=cutoff,
        informative_pair_distance_min=.1,return_ratio_max=.1,stationarity_rhs_max=1e-6,
        reference_log_rms_max=1e-5,reset_noise=.001,reset_seeds=list(range(20)),
        initial='End of 120-unit frozen-context no-feedback continuation, following developmental t=90.',
        exchange='Swap both concentrations; species-specific pair rescaling preserves pair amounts. Other cells unchanged.',
        config=config,sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    output.mkdir(parents=True,exist_ok=False)
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    times = np.arange(0,1201,2.)
    fun = rhs(delta,config['signal_beta'],config['signal_da'],config['signal_dh'])
    max_error = 0.
    def solve(start):
        nonlocal max_error
        path = trajectory(fun,start,times)
        ref = trajectory(fun,start,times,True)
        max_error = max(max_error,float(distance(path,ref,masses).max()))
        return ref  # Analyze the tightened trajectory.
    control = solve(initial)
    late = times >= 1000
    rows=[]; paths=[]
    for (i,j),cd in zip(pairs,context_distances):
        start=exchange(initial,masses,i,j)
        path=solve(start)
        pair=[i,j]; pair_masses=masses[pair]
        separation=float(distance(start[:,pair],initial[:,pair],pair_masses))
        donor=np.array([exchange(frame,masses,i,j) for frame in control])
        dest_distance=distance(path[:,:,pair],control[:,:,pair],pair_masses)
        donor_distance=distance(path[:,:,pair],donor[:,:,pair],pair_masses)
        dest_ratio=float(dest_distance[late].max()/max(separation,1e-30))
        donor_ratio=float(donor_distance[late].max()/max(separation,1e-30))
        residual=float(np.max(np.abs(fun(times[-1],path[-1].ravel()))))
        informative=separation>.1
        rows.append(dict(cell_ids=[int(ids[i]),int(ids[j])],indices=pair,
            context_distance=float(cd),contrasting_context=bool(cd>=cutoff),
            initial_pair_distance=separation,informative=informative,
            destination_ratio=dest_ratio,transferred_ratio=donor_ratio,
            global_destination_distance=float(distance(path[-1],control[-1],masses)),
            final_rhs_max=residual,
            outcome=classify(dest_ratio,donor_ratio,residual<1e-6,informative)))
        paths.append(path)
    uniform=np.repeat(((initial@masses)/masses.sum())[:,None],len(masses),axis=1)
    reset=solve(uniform)
    resets=[];reset_rows=[]
    for seed in range(20):
        path=solve(perturb(uniform,masses,.001,np.random.default_rng(seed)))
        resets.append(path)
        reset_rows.append(dict(seed=seed,final_log_activator_sd=float(np.std(np.log(path[-1,0]))),
            destination_distance=float(distance(path[-1],control[-1],masses)),
            final_rhs_max=float(np.max(np.abs(fun(times[-1],path[-1].ravel()))))))
    def counts(selected):
        return {key:sum(r['outcome']==key for r in selected)
                for key in ['destination','transferred','reorganized','unsettled','uninformative']}
    result=dict(all_pairs=counts(rows),contrasting_pairs=counts([r for r in rows if r['contrasting_context']]),
        reference_max_log_rms=max_error,reference_pass=max_error<1e-5,
        control_drift=float(distance(control[-1],initial,masses)),
        exact_uniform_final_log_sd=float(np.std(np.log(reset[-1,0]))),
        pairs=rows,seeded_uniform_resets=reset_rows)
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(output/'trajectories.npz',times=times,control=control,exchanges=np.array(paths),
        uniform_reset=reset,seeded_resets=np.array(resets),pairs=np.array(pairs),
        initial=initial,masses=masses,delta=delta,cell_ids=ids,context=context)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for key,color in [('destination','tab:blue'),('transferred','tab:orange'),('reorganized','tab:green'),('unsettled','tab:red')]:
        selected=[r for r in rows if r['outcome']==key]
        if selected:
            axes[0].scatter([r['destination_ratio'] for r in selected],
                            [r['transferred_ratio'] for r in selected],label=f'{key} ({len(selected)})',color=color)
    axes[0].axvline(.1,color='grey',ls='--');axes[0].axhline(.1,color='grey',ls='--')
    axes[0].set(xlabel='Distance to destination / initial separation',ylabel='Distance to transferred state / initial separation',xscale='log',yscale='log',title='Informative pair exchanges');axes[0].legend()
    axes[1].plot(times,np.std(np.log(control[:,0]),axis=1),color='black',label='Unmodified pattern')
    axes[1].plot(times,np.std(np.log(reset[:,0]),axis=1),ls='--',label='Uniform reset (numerical symmetry)')
    for k,path in enumerate(resets):
        axes[1].plot(times,np.std(np.log(path[:,0]),axis=1),alpha=.4,color='tab:green',label='Seeded uniform resets' if k==0 else None)
    axes[1].set(xlabel='Time after intervention',ylabel='Across-cell SD of log activator',title='Pattern regeneration controls');axes[1].legend()
    fig.tight_layout();fig.savefig(output/'exchange.png',dpi=160);plt.close(fig)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('outputs/attribute-persistence'))
    parser.add_argument('--development',type=Path,default=Path('outputs/attribute-development'))
    parser.add_argument('--output',type=Path,default=Path('outputs/attribute-exchange'))
    args=parser.parse_args();run(args.source,args.development,args.output)
