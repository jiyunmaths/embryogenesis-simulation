"""Assess the matched response screen against its mechanical baseline."""
import argparse
import json
from pathlib import Path
import numpy as np
from .attribute_development import AttributeSimulation


def summarize(output):
    rows=json.loads((output/'results.json').read_text());baseline=next(r for r in rows if r['arm']['name']=='baseline')
    for row in rows:
        row['growth_change_relative_to_baseline']=row['growth_change']-baseline['growth_change']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for law,color in [('linear','tab:blue'),('exponential','tab:orange')]:
        selected=sorted([r for r in rows if r['arm']['law']==law and not r['arm']['polarity']],key=lambda r:r['arm']['gamma'])
        axes[0].plot([r['arm']['gamma'] for r in selected],[r['growth_change_relative_to_baseline'] for r in selected],'-o',label=law,color=color)
    axes[0].axhline(0,color='grey',ls='--');axes[0].set(xlabel='Tension coupling strength',ylabel='Growth-rate change relative to baseline',title='Spectral response after mechanics evolves');axes[0].legend()
    for name in ['baseline','polarity_only','linear_0.25','historical_full','exponential_2']:
        h=json.loads((output/name/'history.json').read_text())
        axes[1].plot([r['elapsed'] for r in h],[r['maximum_spatial_growth'] for r in h],label=name)
    axes[1].axhline(0,color='grey',ls='--');axes[1].set(xlabel='Time after intervention',ylabel='Instantaneous homogeneous growth rate',title='Matched moving-geometry interventions');axes[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(output/'response.png',dpi=160);plt.close(fig)
    sim=AttributeSimulation.restore(json.loads((output/'protocol.json').read_text())['checkpoint'])
    g0=sim.signaling_graph();v0=sim.volumes();a0=sim.activator.copy();x0=sim.centers()
    distance0=np.linalg.norm(x0[:,None]-x0[None,:],axis=-1)
    edges=np.triu(g0.weights>0,1)
    contact0=sim.contacts()[0]
    geometry=[]
    for row in rows:
        name=row['arm']['name'];d=np.load(output/name/'response_state.npz')
        sim.phi=d['phi'];sim.activator=d['activator'];sim.inhibitor=d['inhibitor']
        graph=sim.signaling_graph();x=sim.centers();v=sim.volumes()
        dist=np.linalg.norm(x[:,None]-x[None,:],axis=-1)
        a,b=sim.activator,sim.inhibitor;n=len(a);c=sim.config
        jac=np.block([[np.diag(2*a/b-1)+c.signal_da*graph.delta,np.diag(-a*a/b**2)],
                      [np.diag(2*c.signal_beta*a),-c.signal_beta*np.eye(n)+c.signal_dh*graph.delta]])
        geometry.append(dict(name=name,mean_original_edge_distance_change=float(np.mean(dist[edges]-distance0[edges])),
            original_edge_overlap_ratio=float(sim.contacts()[0][edges].sum()/contact0[edges].sum()),
            correlation_initial_activity_volume_change=float(np.corrcoef(a0,v/v0-1)[0,1]),
            final_chemical_jacobian_max_real=float(np.linalg.eigvals(jac).real.max())))
    (output/'assessment.json').write_text(json.dumps(dict(rows=rows,geometry=geometry),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('outputs/feedback-response'))
    summarize(p.parse_args().output)
