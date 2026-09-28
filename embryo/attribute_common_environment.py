"""Separate autonomous chemical persistence from shared-reservoir persistence."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .attribute_persistence import trajectory, rhs


def bath_rhs(beta, ka, kb, bath):
    def evaluate(t,y):
        a,b=y.reshape(2,-1)
        return np.concatenate((a*a/b-a+ka*(bath[0]-a),
                               beta*(a*a-b)+kb*(bath[1]-b)))
    return evaluate


def equilibria(beta,ka,kb,bath):
    if ka==0 and kb==0:
        roots=[1.]
    else:
        roots=np.roots([-beta*(1+ka),beta*ka*bath[0]+beta+kb,
                       -kb*bath[1]*(1+ka),ka*bath[0]*kb*bath[1]])
        roots=sorted(float(a.real) for a in roots if abs(a.imag)<1e-8 and a.real>0)
    rows=[]
    for a in roots:
        b=(beta*a*a+kb*bath[1])/(beta+kb)
        jac=np.array([[2*a/b-1-ka,-a*a/b**2],[2*beta*a,-beta-kb]])
        growth=float(np.linalg.eigvals(jac).real.max())
        rows.append(dict(activator=a,inhibitor=b,max_real_eigenvalue=growth,locally_stable=growth<0))
    return rows


def run(source,development,output):
    data=np.load(source/'trajectories.npz')
    # Original stationary pattern plus all exchange and seeded-reset endpoints.
    patterns=np.concatenate([data['control'][-1][None],data['exchanges'][:,-1],data['seeded_resets'][:,-1]])
    starts=patterns.transpose(1,0,2).reshape(2,-1)
    config=json.loads((development/'protocol.json').read_text())['config']
    beta=config['signal_beta'];da=config['signal_da'];db=config['signal_dh']
    rate=float(np.median(-np.diag(data['delta'])))
    mean=patterns[0]@data['masses']/data['masses'].sum()
    arms=[dict(name='isolated',bath=[1.,1.],scale=0.)]
    for name,bath in [('unit',[1.,1.]),('mean',mean.tolist())]:
        for scale in [.25,1.,4.]:
            arms.append(dict(name=f'{name}_{scale:g}',bath=bath,scale=scale))
    files=[Path(__file__),Path(__file__).with_name('attribute_persistence.py'),source/'trajectories.npz',development/'protocol.json']
    protocol=dict(duration=240,sample_interval=1,late_start=200,stationarity_limit=1e-6,
        reference_max_log_error_limit=1e-5,equilibrium_log_distance_limit=1e-5,
        median_graph_exit_rate=rate,arms=arms,config=config,
        initial_state_order='Original control endpoint, 120 exchange endpoints, 20 seeded-reset endpoints; each contains 16 cells.',
        scope='2256 chemical initial states drawn from one developed geometry; duplicates retained, not independent replicates. Frozen infinite reservoirs; no mechanics, polarity, fate dynamics or cell-cell transport in assay arms.',
        sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    output.mkdir(parents=True,exist_ok=False)
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    times=np.arange(241.)
    result={};saved=dict(times=times,initial_patterns=patterns,cell_ids=data['cell_ids'])
    # Same-state network continuation is a positive control for spatial contrast.
    network_fun=rhs(data['delta'],beta,da,db)
    network=trajectory(network_fun,patterns[0],times,True)
    saved['network_control']=network
    result['network_control']=dict(final_log_activator_sd=float(np.std(np.log(network[-1,0]))),
        final_rhs_max=float(np.abs(network_fun(240,network[-1].ravel())).max()))
    for arm in arms:
        ka=arm['scale']*rate*da;kb=arm['scale']*rate*db
        fun=bath_rhs(beta,ka,kb,arm['bath'])
        coarse=trajectory(fun,starts,times)
        path=trajectory(fun,starts,times,True)
        error=float(np.abs(np.log(coarse/path)).max())
        eq=equilibria(beta,ka,kb,arm['bath'])
        stable=np.array([[r['activator'],r['inhibitor']] for r in eq if r['locally_stable']])
        residual=float(np.abs(fun(240,path[-1].ravel())).max())
        d=np.sqrt(np.mean((np.log(path[-1].T[:,None,:]/stable[None,:,:]))**2,axis=2))
        closest=np.argmin(d,axis=1)
        counts=np.bincount(closest,minlength=len(stable))
        late=path[times>=200]
        result[arm['name']]=dict(ka=ka,kb=kb,bath=arm['bath'],equilibria=eq,
            nearest_stable_equilibrium_counts=counts.tolist(),stable_equilibria=stable.tolist(),
            max_final_equilibrium_log_distance=float(d.min(axis=1).max()),
            max_final_rhs=residual,max_reference_log_error=error,
            final_log_activator_range=[float(np.log(path[-1,0]).min()),float(np.log(path[-1,0]).max())],
            late_max_log_change=float(np.abs(np.log(late/late[-1])).max()),
            numerical_pass=bool(error<1e-5 and residual<1e-6 and d.min(axis=1).max()<1e-5))
        saved[arm['name']]=path
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(output/'trajectories.npz',**saved)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].plot(times,np.std(np.log(network[:,0]),axis=1),label='Original network')
    for name in ['isolated','unit_1','mean_1']:
        original=saved[name][:,0,:16]
        axes[0].plot(times,np.std(np.log(original),axis=1),label=name)
    axes[0].set(xlabel='Time in new environment',ylabel='Across-cell SD of log activator',title='Original pattern under identical conditions');axes[0].legend()
    for k,arm in enumerate(arms):
        name=arm['name'];vals=saved[name][-1,0]
        axes[1].scatter(np.full(len(vals),k),vals,s=10,alpha=.1,color='tab:blue')
        for eq in result[name]['equilibria']:
            axes[1].scatter(k,eq['activator'],marker='_' if eq['locally_stable'] else 'x',s=120,color='black' if eq['locally_stable'] else 'tab:red')
    axes[1].set_xticks(range(len(arms)),[a['name'] for a in arms],rotation=45,ha='right')
    axes[1].set(ylabel='Final activator',title='Initial-state outcomes and analytic equilibria')
    fig.tight_layout();fig.savefig(output/'common-environment.png',dpi=160);plt.close(fig)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/attribute-exchange'))
    p.add_argument('--development',type=Path,default=Path('outputs/attribute-development'))
    p.add_argument('--output',type=Path,default=Path('outputs/attribute-common-environment'))
    args=p.parse_args();run(args.source,args.development,args.output)
