"""Amount-conservative finite shared reservoir with intracellular GM reactions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp


def system(volumes, reservoir_volume, beta, ka, kb, reactions=True, clamped=False):
    volumes=np.asarray(volumes,dtype=float)
    if np.any(volumes<=0) or reservoir_volume<=0:
        raise ValueError('Compartment volumes must be positive')
    n=len(volumes);m=n+1;ratio=volumes/reservoir_volume
    def evaluate(t,y):
        a,b=y.reshape(2,m);x,z=a[:-1],b[:-1]
        fa=ka*(a[-1]-x);fb=kb*(b[-1]-z)
        ra=-np.dot(ratio,fa) if not clamped else 0.
        rb=-np.dot(ratio,fb) if not clamped else 0.
        if reactions:
            fa=fa+x*x/z-x;fb=fb+beta*(x*x-z)
        return np.concatenate([fa,[ra],fb,[rb]])
    def jacobian(t,y):
        a,b=y.reshape(2,m);x,z=a[:-1],b[:-1]
        jac=np.zeros((2*m,2*m));i=np.arange(n)
        jac[i,i]=-ka;jac[i,m-1]=ka
        jac[m+i,m+i]=-kb;jac[m+i,2*m-1]=kb
        if not clamped:
            jac[m-1,i]=ratio*ka;jac[m-1,m-1]=-ka*ratio.sum()
            jac[2*m-1,m+i]=ratio*kb;jac[2*m-1,2*m-1]=-kb*ratio.sum()
        if reactions:
            jac[i,i]+=2*x/z-1;jac[i,m+i]=-x*x/z**2
            jac[m+i,i]=2*beta*x;jac[m+i,m+i]-=beta
        return jac
    return evaluate,jacobian


def integrate(fun,jac,initial,times,tight=False,method='Radau'):
    sol=solve_ivp(fun,(0,times[-1]),initial.ravel(),method=method,jac=jac,
        t_eval=times,rtol=1e-10 if tight else 1e-8,atol=1e-12 if tight else 1e-10)
    if not sol.success or np.any(sol.y<=0) or not np.isfinite(sol.y).all():
        raise RuntimeError('Reservoir integration failed: '+sol.message)
    return sol.y.T.reshape(len(times),*initial.shape)


def run(source,exchange,output):
    prior=np.load(source/'trajectories.npz');geo=np.load(exchange/'trajectories.npz')
    p=json.loads((source/'protocol.json').read_text());v=geo['masses'];n=len(v)
    beta=p['config']['signal_beta'];rate=p['median_graph_exit_rate']
    arms=[]
    for scale in [.25,1.,4.]:
        ka=scale*rate*p['config']['signal_da'];kb=scale*rate*p['config']['signal_dh']
        for size in [.1,1.,10.,100.]:
            # Twenty slow-species bath exchange times, at least 1200 model units.
            duration=max(1200.,20*size/ka)
            for bath_name in ['unit','mean']:
                key=f'{bath_name}_{scale:g}'
                bath=next(a['bath'] for a in p['arms'] if a['name']==key)
                start=np.column_stack([prior[key][-1,:,:n],bath])
                arms.append(dict(name=f'release_{key}_R{size:g}',kind='release',scale=scale,size=size,
                                 ka=ka,kb=kb,duration=duration,initial=start.tolist()))
            for seed in [None,0,1,2,3,4]:
                start=np.ones((2,n+1))
                if seed is not None:
                    rng=np.random.default_rng(seed)
                    start[:,:n]=np.exp(.001*rng.normal(size=(2,n)))
                    # Neutral reservoir and exact total-cell amount for each species.
                    start[:,:n]*=(v.sum()/(start[:,:n]@v))[:,None]
                arms.append(dict(name=f'formation_s{scale:g}_R{size:g}_seed{seed}',kind='formation',seed=seed,
                                 scale=scale,size=size,ka=ka,kb=kb,duration=duration,initial=start.tolist()))
    # Matched clamped controls for release, plus isolated original pattern.
    for arm in p['arms']:
        if arm['scale']==0:continue
        key=arm['name'];scale=arm['scale']
        arms.append(dict(name='clamped_'+key,kind='clamped',scale=scale,size=1.,
            ka=scale*rate*p['config']['signal_da'],kb=scale*rate*p['config']['signal_dh'],duration=1200.,
            initial=np.column_stack([prior[key][-1,:,:n],arm['bath']]).tolist()))
    arms.append(dict(name='isolated',kind='isolated',scale=0.,size=1.,ka=0.,kb=0.,duration=1200.,
                     initial=np.column_stack([geo['control'][-1],np.ones(2)]).tolist()))
    files=[Path(__file__),source/'trajectories.npz',source/'protocol.json',exchange/'trajectories.npz']
    protocol=dict(arms=arms,volumes=v.tolist(),beta=beta,late_fraction=.8,
        criteria=dict(persistent_min_log_activator_sd=.1,stationary_max_rhs=1e-6,
                      max_trajectory_log_error=1e-5,max_amount_balance_residual=1e-10),
        scope='One embryo geometry; 24 reservoir-release arms, 60 seeded formation arms, 12 exact-uniform controls, 6 clamped controls, 1 isolation control. Uniform exchange coefficients; no direct cell contacts, reservoir reactions, external replenishment, mechanics, or fate switch.',
        sha256={str(f.resolve()):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    output.mkdir(parents=True,exist_ok=False)
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    results=[]
    for index,arm in enumerate(arms):
        start=np.array(arm['initial']);rv=arm['size']*v.sum()
        fun,jac=system(v,rv,beta,arm['ka'],arm['kb'],clamped=arm['kind']=='clamped')
        # Dense early sampling plus full-horizon sampling, including the late window.
        times=np.unique(np.concatenate([np.linspace(0,min(240,arm['duration']),241),np.linspace(0,arm['duration'],301)]))
        coarse=integrate(fun,jac,start,times)
        path=integrate(fun,jac,start,times,True)
        err=float(np.max(np.abs(np.log(coarse/path))))
        late=times>=.8*arm['duration'];spread=np.std(np.log(path[:,0,:n]),axis=1)
        residual=float(np.abs(fun(times[-1],path[-1].ravel())).max())
        balance=0.
        if arm['kind']!='clamped':
            for frame in path:
                derivative=fun(0,frame.ravel()).reshape(2,n+1)
                a,b=frame[:,:n]
                reaction=np.array([a*a/b-a,beta*(a*a-b)])
                balance=max(balance,float(np.abs(derivative@np.append(v,rv)-reaction@v).max()))
        eig=np.linalg.eigvals(jac(0,path[-1].ravel()))
        result=dict(name=arm['name'],kind=arm['kind'],size=arm['size'],scale=arm['scale'],seed=arm.get('seed'),
            duration=arm['duration'],late_min_log_activator_sd=float(spread[late].min()),
            final_log_activator_sd=float(spread[-1]),final_rhs_max=residual,
            final_reservoir=path[-1,:,-1].tolist(),reference_max_log_error=err,
            amount_balance_residual=balance if arm['kind']!='clamped' else None,
            final_jacobian_max_real=float(eig.real.max()),
            persistent=bool(spread[late].min()>.1),stationary=bool(residual<1e-6),
            numerical_pass=bool(err<1e-5 and balance<1e-10))
        results.append(result)
        np.savez_compressed(output/(arm['name']+'.npz'),times=times,trajectory=path,initial=start)
        (output/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        (output/'status.json').write_text(json.dumps(dict(completed=index+1,total=len(arms),state='running'))+'\n')
    (output/'status.json').write_text(json.dumps(dict(completed=len(arms),total=len(arms),state='completed'))+'\n')
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('outputs/attribute-common-environment'))
    parser.add_argument('--exchange',type=Path,default=Path('outputs/attribute-exchange'))
    parser.add_argument('--output',type=Path,default=Path('outputs/attribute-finite-reservoir'))
    args=parser.parse_args();run(args.source,args.exchange,args.output)
