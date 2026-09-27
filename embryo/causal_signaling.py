"""Matched activator-inhibitor ablations on one frozen resolved embryo graph.

Twenty perturbation seeds isolate biochemical and signal-to-fate causality.
Geometry is inherited and frozen: this does not test shape causality.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .model import Simulation
from .resolution import write_json, _steps


ARMS=('full','no_self_activation','no_induced_inhibitor','no_inhibitor_action',
      'no_transport','equal_diffusion','no_signal_to_fate')


def reaction(a,h,beta,arm):
    production=1/h if arm=='no_self_activation' else a*a if arm=='no_inhibitor_action' else a*a/h
    inhibitor_source=np.ones_like(a) if arm=='no_induced_inhibitor' else a*a
    return production-a,beta*(inhibitor_source-h)


def jacobian(beta,arm):
    return np.array([[-1 if arm=='no_self_activation' else 1,0 if arm=='no_inhibitor_action' else -1],
                     [0 if arm=='no_induced_inhibitor' else 2*beta,-beta]],float)


def coefficients(arm,da,dh):
    return (0.,0.) if arm=='no_transport' else (da,da) if arm=='equal_diffusion' else (da,dh)


def preflight(eigenvalues,beta,da,dh,arm):
    da,dh=coefficients(arm,da,dh);j=jacobian(beta,arm)
    rates=np.array([np.max(np.linalg.eigvals(j-l*np.diag([da,dh])).real) for l in eigenvalues])
    local=float(np.max(np.linalg.eigvals(j).real))
    return {'jacobian':j.tolist(),'diffusivities':[da,dh],'local_max_growth':local,
            'local_stable':local<0,'growth_rates':rates.tolist(),
            'unstable_spatial_modes':np.flatnonzero((eigenvalues>1e-10)&(rates>1e-10)).tolist(),
            'interpretation':'Diffusion-driven amplification requires local stability; locally unstable knockouts are dysregulation controls, not Turing-pattern evidence.'}


def initial_signals(volumes,seeds,amplitude):
    weights=volumes/volumes.sum();a=[];h=[]
    for seed in seeds:
        draw=np.random.default_rng(seed).normal(size=(2,len(volumes)))
        draw-=np.sum(draw*weights,axis=1)[:,None]
        draw*=amplitude/np.sqrt(np.sum(draw*draw*weights,axis=1))[:,None]
        a.append(1+draw[0]);h.append(1+draw[1])
    return np.array(a),np.array(h)


def simulate(graph,config,seeds,arm,dt=.01,until=60.,interval=.5,amplitude=.001,ceiling=10.):
    if arm not in ARMS:raise ValueError('unknown causal arm')
    count=_steps(until,dt);every=_steps(interval,dt)
    if count%every:raise ValueError('observation times must align')
    weights=graph.masses/graph.masses.sum()
    a,h=initial_signals(graph.masses,seeds,amplitude);fate=np.zeros_like(a)
    da,dh=coefficients(arm,config.signal_da,config.signal_dh)
    substeps=max(1,int(np.ceil(dt*max(1+da*np.max(-np.diag(graph.delta)),config.signal_beta+dh*np.max(-np.diag(graph.delta)))/.2)))
    step=dt/substeps;active=np.ones(len(seeds),bool);stopped=[None]*len(seeds);history=[]
    def record(time):
        means=a@weights;contrast=np.sqrt(((a-means[:,None])**2)@weights)
        return {'time':time,'active':active.tolist(),'activator':a.tolist(),'inhibitor':h.tolist(),'fate':fate.tolist(),
                'activator_contrast':contrast.tolist(),'activator_mean':means.tolist(),
                'fate_a':np.sum(fate>config.fate_threshold,axis=1).tolist(),
                'fate_b':np.sum(fate < -config.fate_threshold,axis=1).tolist()}
    history.append(record(0.))
    with np.errstate(over='raise',divide='raise',invalid='raise'):
        for number in range(1,count+1):
            for sub in range(substeps):
                indices=np.flatnonzero(active)
                if not len(indices):break
                x,y=a[indices],h[indices]
                def rhs(x,y):
                    fa,fh=reaction(x,y,config.signal_beta,arm)
                    return fa+da*(x@graph.delta.T),fh+dh*(y@graph.delta.T)
                fa,fh=rhs(x,y);xa,ha=x+step*fa,y+step*fh
                fa,fh=rhs(xa,ha)
                a[indices]=.5*x+.5*(xa+step*fa);h[indices]=.5*y+.5*(ha+step*fh)
                if not np.isfinite(a).all() or not np.isfinite(h).all() or np.any(a<0) or np.any(h<=0):
                    raise FloatingPointError('kinetics lost finite positivity')
                # Record finite exits, never clip or treat runaway as organization.
                for i in indices[np.maximum(a[indices].max(axis=1),h[indices].max(axis=1))>=ceiling]:
                    active[i]=False;stopped[i]={'time':(number-1)*dt+(sub+1)*step,'reason':'regulator reached predeclared ceiling','maximum':float(max(a[i].max(),h[i].max()))}
            indices=np.flatnonzero(active)
            gain=0. if arm=='no_signal_to_fate' else config.signal_fate_gain
            z=fate[indices];signal=a[indices]
            def drift(value):return config.fate_rate*(value-value**3+gain*(signal-1))
            first=z+dt*drift(z);fate[indices]=.5*z+.5*(first+dt*drift(first))
            if not np.isfinite(fate).all():raise FloatingPointError('nonfinite fate')
            if number%every==0:history.append(record(number*dt))
    outcomes=[]
    late=[r for r in history if r['time']>=until-20.]
    for i,seed in enumerate(seeds):
        outcomes.append({'seed':seed,'completed':bool(active[i]),'stopped':stopped[i],
                         'final_contrast':history[-1]['activator_contrast'][i],
                         'late_minimum_contrast':min(r['activator_contrast'][i] for r in late),
                         'persistent_contrast_above_0_1':bool(active[i] and min(r['activator_contrast'][i] for r in late)>=.1),
                         'late_contrast_below_0_01':bool(active[i] and max(r['activator_contrast'][i] for r in late)<.01),
                         'both_fates_at_end':bool(active[i] and history[-1]['fate_a'][i]>0 and history[-1]['fate_b'][i]>0),
                         'final_fate_a':history[-1]['fate_a'][i],'final_fate_b':history[-1]['fate_b'][i]})
    return {'arm':arm,'dt':dt,'seeds':list(seeds),'preflight':preflight(graph.eigenvalues,config.signal_beta,config.signal_da,config.signal_dh,arm),
            'outcomes':outcomes,'history':history}


def prepare(checkpoint,output):
    checkpoint=Path(checkpoint).resolve();output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    sim=Simulation.restore(checkpoint);metrics=sim.metrics()
    if sim.config.signal_transport!='conservative' or sim.divisions or len(sim.phi)!=16 or metrics['min_radius_grid_cells']<4 or metrics['boundary_occupancy']>=.01:
        raise ValueError('requires a resolved, boundary-clear, post-cleavage 16-cell conservative state')
    p={'checkpoint':str(checkpoint),'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
       'arms':list(ARMS),'seeds':list(range(20)),'dt':.01,'until':60.,'interval':.5,'perturbation_rms':.001,'regulator_ceiling':10.,
       'criteria':{'full_persistence_minimum_seeds':16,'suppression_minimum_seeds':16,'temporal_signal_rms_max':.01,'temporal_fate_absolute_max':.05},
       'temporal_confirmation_arms':['full','no_self_activation','equal_diffusion'],'confirmation_dt':.005,
       'source_geometry_metrics':metrics,'scope':__doc__,
       'interpretation':'Signals reset to homogeneous equilibrium plus matched unbiased perturbations; fate reset to zero; geometry and conductances fixed. No inherited labels or imposed spatial axis. All reaction ablations retain equilibrium (1,1); inhibition knockouts can destroy its local stability. Regulator ceilings are finite scientific exits, not clipped successful patterns. No shape outcome is tested. Twenty chemical perturbation seeds on one geometry are not twenty developmental replicates.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',p)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return p


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:raise ValueError('requires unchanged protocol')
    if hashlib.sha256(Path(p['checkpoint']).read_bytes()).hexdigest()!=p['checkpoint_sha256']:raise ValueError('checkpoint changed')
    sim=Simulation.restore(p['checkpoint']);graph=sim.signaling_graph();completed=[]
    try:
        reports={}
        for arm in p['arms']:
            write_json(output/'status.json',{**state,'state':'running','arm':arm,'completed_arms':completed})
            report=simulate(graph,sim.config,p['seeds'],arm,p['dt'],p['until'],p['interval'],p['perturbation_rms'],p['regulator_ceiling'])
            write_json(output/(arm+'.json'),report);reports[arm]=report;completed.append(arm)
            print(arm,[(key,sum(r[key] for r in report['outcomes'])) for key in ('completed','persistent_contrast_above_0_1','both_fates_at_end')],flush=True)
        temporal=[]
        for arm in p['temporal_confirmation_arms']:
            fine=simulate(graph,sim.config,p['seeds'],arm,p['confirmation_dt'],p['until'],p['interval'],p['perturbation_rms'],p['regulator_ceiling'])
            write_json(output/(arm+'-half-dt.json'),fine)
            coarse=reports[arm]
            signal=max(float(np.max(np.sqrt(np.mean((np.array(a[key])-b[key])**2,axis=1)))) for a,b in zip(coarse['history'],fine['history']) for key in ('activator','inhibitor'))
            fate=max(float(np.max(abs(np.array(a['fate'])-b['fate']))) for a,b in zip(coarse['history'],fine['history']))
            temporal.append({'arm':arm,'max_signal_rms_difference':signal,'max_fate_difference':fate,
                             'all_trajectories_complete':all(x['completed'] for x in coarse['outcomes']+fine['outcomes']),
                             'final_fate_counts_match':all(a['final_fate_a']==b['final_fate_a'] and a['final_fate_b']==b['final_fate_b'] for a,b in zip(coarse['outcomes'],fine['outcomes']))})
        c=p['criteria'];evidence={
            'full_loop_persistent_in_at_least_16_seeds':sum(x['persistent_contrast_above_0_1'] for x in reports['full']['outcomes'])>=c['full_persistence_minimum_seeds'],
            'no_self_activation_suppresses_contrast':sum(x['late_contrast_below_0_01'] for x in reports['no_self_activation']['outcomes'])>=c['suppression_minimum_seeds'],
            'no_transport_suppresses_contrast':sum(x['late_contrast_below_0_01'] for x in reports['no_transport']['outcomes'])>=c['suppression_minimum_seeds'],
            'equal_diffusion_suppresses_contrast':sum(x['late_contrast_below_0_01'] for x in reports['equal_diffusion']['outcomes'])>=c['suppression_minimum_seeds'],
            'full_loop_generates_both_fates_in_at_least_16_seeds':sum(x['both_fates_at_end'] for x in reports['full']['outcomes'])>=c['full_persistence_minimum_seeds'],
            'removing_signal_to_fate_keeps_all_fates_uncommitted':all(x['completed'] and x['final_fate_a']==x['final_fate_b']==0 for x in reports['no_signal_to_fate']['outcomes'])}
        quality={'time_refinement_signals_agree':all(t['max_signal_rms_difference']<c['temporal_signal_rms_max'] and t['all_trajectories_complete'] for t in temporal),
                 'time_refinement_fates_agree':all(t['max_fate_difference']<c['temporal_fate_absolute_max'] and t['final_fate_counts_match'] for t in temporal),
                 'checkpoint_unchanged':hashlib.sha256(Path(p['checkpoint']).read_bytes()).hexdigest()==p['checkpoint_sha256']}
        summary=[{'arm':arm,'locally_stable':r['preflight']['local_stable'],'unstable_modes':len(r['preflight']['unstable_spatial_modes']),
                  **{key:sum(x[key] for x in r['outcomes']) for key in ('completed','persistent_contrast_above_0_1','both_fates_at_end')}} for arm,r in reports.items()]
        result={'protocol':p,'summary':summary,'causal_evidence':evidence,'numerical_checks':quality,'temporal_confirmation':temporal,
                'hypothesis_supported_in_this_screen':all(evidence.values()) and all(quality.values()),'shape_causality_tested':False}
        write_json(output/'comparison.json',result)
        lines=['# Fixed-geometry causal activator-inhibitor screen','',f'All declared causal and numerical checks supported: **{result["hypothesis_supported_in_this_screen"]}**','',p['interpretation'],'',
               '| Arm | Local stability | Completed / 20 | Persistent contrast / 20 | Both fates / 20 |','|---|---|---:|---:|---:|']
        lines += [f'| {r["arm"]} | {r["locally_stable"]} | {r["completed"]} | {r["persistent_contrast_above_0_1"]} | {r["both_fates_at_end"]} |' for r in summary]
        lines+=['',*[f'- {k}: {v}' for k,v in {**evidence,**quality}.items()]]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json',{**state,'state':'completed','completed_arms':completed,'hypothesis_supported':result['hypothesis_supported_in_this_screen']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','completed_arms':completed,'error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,default=Path('outputs/development-refinement/space-72/state-18.npz'))
    args=parser.parse_args();prepare(args.checkpoint,args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
