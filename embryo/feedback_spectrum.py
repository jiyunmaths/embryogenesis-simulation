"""Developmental transport spectra and geometric counterfactuals."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .attribute_development import AttributeSimulation
from .signaling import stability
from .transport import conservative_transport, transport_graph


def cv(x):
    x=np.asarray(x)
    return float(np.std(x)/np.mean(x)) if len(x) and np.mean(x)>0 else None


def summarize(sim):
    graph=sim.signaling_graph();g=graph.weights;v=graph.masses
    edges=g[np.triu_indices(len(v),1)];edges=edges[edges>0]
    report=stability(graph,sim.config.signal_beta,sim.config.signal_da,sim.config.signal_dh)
    report.update(time=sim.time,edge_conductance_cv=cv(edges),volume_cv=cv(v),
                  weighted_degree_cv=cv(g.sum(axis=1)),mean_positive_conductance=float(edges.mean()) if len(edges) else 0.,
                  total_conductance=float(edges.sum()),log_activator_sd=float(np.std(np.log(sim.activator))))
    return report,graph


def run(source,output):
    output.mkdir(parents=True,exist_ok=False)
    snapshots={mode:sorted((source/mode).glob('state-*.npz'),key=lambda x:float(x.stem.split('-')[1])) for mode in ['direct','no_feedback']}
    protocol=dict(source=str(source),scope='Frozen instantaneous spectra do not alone establish stability of a time-dependent, dividing, diluting system.',
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for paths in snapshots.values() for p in paths})
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    results={};graphs={}
    for mode,paths in snapshots.items():
        results[mode]=[]
        for p in paths:
            sim=AttributeSimulation.restore(p);row,graph=summarize(sim)
            results[mode].append(row);graphs[(mode,int(round(sim.time)))]=graph
        (output/'spectra.json').write_text(json.dumps(results,indent=2)+'\n')
    # At t=90 both graph arrays follow the same ID order (checked explicitly).
    sims=[AttributeSimulation.restore(source/mode/'state-90.npz') for mode in ['direct','no_feedback']]
    if not np.array_equal(sims[0].ids,sims[1].ids):raise ValueError('Cannot mix matrices with different cell orders')
    gd=graphs[('direct',90)];gn=graphs[('no_feedback',90)];cfg=sims[0].config
    counter={}
    for label,g,v in [('direct_G_direct_V',gd.weights,gd.masses),('direct_G_no_feedback_V',gd.weights,gn.masses),
                      ('no_feedback_G_direct_V',gn.weights,gd.masses),('no_feedback_G_no_feedback_V',gn.weights,gn.masses)]:
        counter[label]=stability(transport_graph(conservative_transport(g,v)),cfg.signal_beta,cfg.signal_da,cfg.signal_dh)
    # Keep graph support and masses, replace each positive conductance by its mean.
    for label,graph in [('direct',gd),('no_feedback',gn)]:
        mask=graph.weights>0;g=np.where(mask,graph.weights[mask].mean(),0.)
        counter[label+'_equal_positive_edges']=stability(transport_graph(conservative_transport(g,graph.masses)),cfg.signal_beta,cfg.signal_da,cfg.signal_dh)
    scale=gd.weights.sum()/gn.weights.sum()
    counter['no_feedback_G_scaled_to_direct_total']=stability(transport_graph(conservative_transport(gn.weights*scale,gn.masses)),cfg.signal_beta,cfg.signal_da,cfg.signal_dh)
    (output/'counterfactuals.json').write_text(json.dumps(dict(conductance_scale=float(scale),cases=counter),indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for mode,color in [('direct','tab:blue'),('no_feedback','tab:orange')]:
        rows=results[mode];times=[r['time'] for r in rows]
        for row in rows:
            axes[0,0].scatter([row['time']]*len(row['eigenvalues'][1:]),row['eigenvalues'][1:],s=10,color=color,alpha=.65)
        axes[0,1].plot(times,[r['maximum_spatial_growth'] for r in rows],label=mode,color=color)
        axes[1,0].plot(times,[r['edge_conductance_cv'] for r in rows],label=mode,color=color)
        axes[1,1].plot(times,[r['total_conductance'] for r in rows],label=mode,color=color)
    band=results['direct'][-1]['continuous_lambda_band'];axes[0,0].axhspan(*band,color='green',alpha=.12,label='Unstable band')
    axes[0,0].set(ylabel='Eigenvalues of −ΔV',title='Sampled developmental spectra');axes[0,0].legend()
    axes[0,1].axhline(0,color='grey',ls='--');axes[0,1].set(ylabel='Maximum spatial growth rate',title='Frozen homogeneous stability');axes[0,1].legend()
    axes[1,0].set(ylabel='Positive-edge coefficient of variation',title='Conductance heterogeneity');axes[1,0].legend()
    axes[1,1].set(ylabel='Sum of unique-edge conductances',title='Total conductance');axes[1,1].legend()
    for ax in axes.ravel():ax.set_xlabel('Developmental time')
    fig.tight_layout();fig.savefig(output/'spectra.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/attribute-development'))
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-spectrum'))
    args=p.parse_args();run(args.source,args.output)
