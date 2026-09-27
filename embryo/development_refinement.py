"""Full zygote-to-aggregate refinement with independent voxel/time changes."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.ndimage import label
from scipy.special import expit
from scipy.stats import wasserstein_distance

from .model import Config, Simulation, occupancy
from .resolution import write_json, _steps
from .shape import observe


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def initial(config, simulation_type=Simulation):
    sim=simulation_type(config)
    # Simulation samples the radius-0.8 analytic zygote directly. Give every
    # resolution the same continuous target, rather than a voxel-sum target.
    width=config.interface_width
    volume=quad(lambda r:4*np.pi*r*r*occupancy(expit(np.sqrt(2)*(.8-r)/width)),
                0,.8+30*width,epsabs=1e-12,epsrel=1e-12)[0]
    sim.target[:]=volume
    sim.initial_volume=volume
    return sim


def prepare(output, grids=(56,72,88), until=90., dts=(.015,.0075,.00375), interval=.6):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    if len(grids)!=3 or list(grids)!=sorted(set(grids)) or len(dts)!=3 or list(dts)!=sorted(set(dts),reverse=True):
        raise ValueError('three increasing grids and decreasing time steps required')
    cases={}
    for grid in grids:
        cases[f'space-{grid}']=asdict(Config(grid=grid,extent=2.24,dt=dts[1],steps=_steps(until,dts[1]),save_every=_steps(interval,dts[1])))
    for dt in (dts[0],dts[2]):
        cases[f'time-{dt:g}']=asdict(Config(grid=grids[1],extent=2.24,dt=dt,steps=_steps(until,dt),save_every=_steps(interval,dt)))
    for values in cases.values():
        Config(**values).validate()
        if values['steps']%values['save_every']:raise ValueError('end time must align with observation interval')
    protocol={'cases':cases,'space_cases':[f'space-{g}' for g in grids],
              'time_cases':[f'time-{dts[0]:g}',f'space-{grids[1]}',f'time-{dts[2]:g}'],
              'until':until,'sample_interval':interval,'checkpoint_every_samples':10,
              'initialization':'Analytic radius-0.8 zygote; common continuous radial volume target; seed 7; independently sampled voxels; standard shape-selected divisions and unaltered noise.',
              'criteria':{'axis_ratio_relative_max':.01,'signal_wasserstein_max':.01,'signal_contrast_max':.005,
                          'fate_fraction_max':1/16,'maturity_time_difference_max':.6,'cell_count_difference_max':1,
                          'volume_error_max':.05,'boundary_max':.01,'radius_min':4.,'abscission_amount_jump_max':1e-6},
              'code_sha256':{str(Path(__file__).with_name(n).resolve()):sha(Path(__file__).with_name(n)) for n in ('model.py','signaling.py','transport.py','polarity.py','development_refinement.py')},
              'scope':'Single-seed full development at fixed interface width and domain. Same random seed does not force identical event order, axes, or cellwise perturbations. Distribution comparisons do not prove identical spatial signaling patterns. No imposed axes, resampled mature states, or retuned physics. Historical coarse runs are not reused because the common analytic target changes initialization slightly. Earlier cleavage calibration failure remains unresolved.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':sha(output/'protocol.json')})
    return protocol


class DevelopmentSimulation(Simulation):
    def _complete_division(self,*args):
        volumes=self.volumes()
        before=np.array([volumes.sum(),volumes@self.activator,volumes@self.inhibitor])
        super()._complete_division(*args)
        volumes=self.volumes()
        after=np.array([volumes.sum(),volumes@self.activator,volumes@self.inhibitor])
        self.event_audits.append({'time':self.time,'relative_jumps':abs(after/before-1).tolist()})


def record(sim):
    metrics=observe(sim)
    # A pinching mother can have disconnected lobes before abscission. Only
    # non-dividing cells are screened as intact compartments.
    metrics['nondividing_cell_components']=[int(label(field>=.5)[1]) for i,field in enumerate(sim.phi) if int(sim.ids[i]) not in sim.divisions]
    return {'metrics':metrics,'ids':sim.ids.tolist(),'parents':sim.parents.tolist(),
            'volumes':sim.volumes().tolist(),'centers':sim.centers().tolist(),
            **{name:getattr(sim,name).tolist() for name in ('activator','inhibitor','fate','polarity')}}


def run_case(output,name,p,state,completed):
    path=output/name;path.mkdir()
    config=Config(**p['cases'][name])
    # Same initializer, instrumented only at the abscission method boundary.
    sim=initial(config,DevelopmentSimulation);sim.event_audits=[]
    history=[];frames=[];maximum_volume=0.;minimum_radius=float('inf');clipping=0.;started=time.monotonic()
    for step in range(config.steps+1):
        volumes=sim.volumes()
        maximum_volume=max(maximum_volume,float(np.max(abs(volumes/sim.target-1))))
        minimum_radius=min(minimum_radius,float(np.min((3*volumes/(4*np.pi))**(1/3))/sim.dx))
        clipping=max(clipping,sim.clipped_fraction)
        if step%config.save_every==0:
            row=record(sim);history.append(row)
            write_json(path/'history.json',history)
            write_json(path/'progress.json',{'time':sim.time,'max_volume_error':maximum_volume,'min_radius_grid_cells':minimum_radius,'max_clipped_fraction':clipping,'event_audits':sim.event_audits,'lineage':sim.lineage})
            write_json(output/'status.json',{**state,'state':'running','case':name,'time':sim.time,'completed_cases':completed})
            if step%(config.save_every*p['checkpoint_every_samples'])==0 or step==config.steps:
                checkpoint=path/f'state-{sim.time:g}.npz'
                sim.checkpoint(checkpoint)
                frames.append({'metrics':row['metrics'],'cells':sim.surfaces(max_points=200),'graph':sim.graph_snapshot()})
                print(f'{name}: t={sim.time:g}, cells={len(sim.phi)}, radius/dx={minimum_radius:.3f}, elapsed={time.monotonic()-started:.1f}s',flush=True)
        if step==config.steps:break
        sim.step()
    sim.checkpoint(path/'final_state.npz')
    result={'config':asdict(config),'history':history,'lineage':sim.lineage,'event_audits':sim.event_audits,
            'max_volume_error':maximum_volume,'min_radius_grid_cells':minimum_radius,
            'max_clipped_fraction':clipping,'elapsed_seconds':time.monotonic()-started}
    write_json(path/'analysis.json',result)
    from .surface import viewer_template
    payload=json.dumps({'config':asdict(config),'frames':frames},allow_nan=False,separators=(',',':'))
    (path/'trajectory.json').write_text(payload)
    (path/'viewer.html').write_text(viewer_template().replace('__SIMULATION_DATA__',payload))
    return result


def pair_metrics(left,right):
    a,b=left['history'],right['history']
    if len(a)!=len(b) or any(not np.isclose(x['metrics']['time'],y['metrics']['time'],rtol=0,atol=1e-10) for x,y in zip(a,b)):
        raise ValueError('comparisons require matching observation times')
    result={}
    for species in ('activator','inhibitor'):
        result[species+'_max_wasserstein']=max(float(wasserstein_distance(x[species],y[species],u_weights=x['volumes'],v_weights=y['volumes'])) for x,y in zip(a,b))
        result[species+'_max_contrast_difference']=max(abs(x['metrics'][species+'_std']-y['metrics'][species+'_std']) for x,y in zip(a,b))
    result['max_axis_ratio_relative_difference']=max(abs(x['metrics']['axis_ratio']/y['metrics']['axis_ratio']-1) for x,y in zip(a,b))
    result['max_cell_count_difference']=max(abs(x['metrics']['cells']-y['metrics']['cells']) for x,y in zip(a,b))
    result['max_fate_fraction_difference']=max(abs(x['metrics'][key]/x['metrics']['cells']-y['metrics'][key]/y['metrics']['cells']) for x,y in zip(a,b) for key in ('fate_a','fate_b','uncommitted'))
    # Exact last-abscission time, not a sample-rounded maturity time.
    def mature(report):
        if report['history'][-1]['metrics']['cells']!=report['config']['max_cells']:return None
        times=[r['division'] for r in report['lineage'] if r['division'] is not None]
        return max(times) if times else None
    times=[mature(left),mature(right)]
    result['maturity_times']=times
    result['maturity_time_difference']=abs(times[0]-times[1]) if all(t is not None for t in times) else None
    result['final_cells']=[a[-1]['metrics']['cells'],b[-1]['metrics']['cells']]
    return result


def compare(output,p):
    reports={n:json.loads((output/n/'analysis.json').read_text()) for n in p['cases']}
    pairs={family:[{'cases':[a,b],**pair_metrics(reports[a],reports[b])} for a,b in zip(p[family+'_cases'][:-1],p[family+'_cases'][1:])] for family in ('space','time')}
    checks={};limits=p['criteria']
    for family,rows in pairs.items():
        fine=rows[-1]
        checks[family+'_shape_within_1_percent']=fine['max_axis_ratio_relative_difference']<limits['axis_ratio_relative_max']
        checks[family+'_signal_distributions_within_0_01']=max(fine[k+'_max_wasserstein'] for k in ('activator','inhibitor'))<limits['signal_wasserstein_max']
        checks[family+'_signal_contrast_within_0_005']=max(fine[k+'_max_contrast_difference'] for k in ('activator','inhibitor'))<limits['signal_contrast_max']
        checks[family+'_fate_fractions_within_one_sixteenth']=fine['max_fate_fraction_difference']<=limits['fate_fraction_max']
        checks[family+'_cell_counts_within_one']=fine['max_cell_count_difference']<=limits['cell_count_difference_max']
        checks[family+'_maturity_timing_agrees']=fine['maturity_time_difference'] is not None and fine['maturity_time_difference']<=limits['maturity_time_difference_max']
        for key in ('max_axis_ratio_relative_difference','activator_max_wasserstein','inhibitor_max_wasserstein'):
            checks[family+'_'+key+'_nonincreasing']=fine[key]<=rows[0][key]+1e-10
    quality={}
    for name,r in reports.items():
        quality[name]={'all_divisions_complete':r['history'][-1]['metrics']['cells']==r['config']['max_cells'] and r['history'][-1]['metrics']['dividing_cells']==0,
                       'volumes_within_5_percent':r['max_volume_error']<limits['volume_error_max'],
                       'cells_resolved':r['min_radius_grid_cells']>=limits['radius_min'],
                       'no_clipping':r['max_clipped_fraction']==0,
                       'sampled_boundaries_clear':all(h['metrics']['boundary_occupancy']<limits['boundary_max'] for h in r['history']),
                       'sampled_nondividing_cells_connected':all(all(c==1 for c in h['metrics']['nondividing_cell_components']) for h in r['history']),
                       'abscission_conserves_amount_and_volume':bool(r['event_audits']) and all(max(e['relative_jumps'])<limits['abscission_amount_jump_max'] for e in r['event_audits'])}
    relevant=set(p['space_cases'][-2:]+p['time_cases'][-2:])
    # Coarsest run is a reference, not a candidate resolved production grid.
    result={'protocol':p,'pairs':pairs,'comparison_checks':checks,'quality_by_case':quality,
            'comparison_checks_pass':all(checks.values()),
            'candidate_quality_pass':all(all(quality[n].values()) for n in relevant),
            'candidate_cases':sorted(relevant)}
    result['passed']=result['comparison_checks_pass'] and result['candidate_quality_pass']
    write_json(output/'comparison.json',result)
    lines=['# Full developmental refinement','',f'All declared candidate checks pass: **{result["passed"]}**','',p['scope'],'']
    lines += [f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()]
    lines += ['','Quality failures (including coarse reference):']+[f'- {n}: {k}' for n,q in quality.items() for k,v in q.items() if not v]
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for name,r in reports.items():
        times=[h['metrics']['time'] for h in r['history']]
        for ax,key in zip(axes.flat,('cells','activator_std','fate_a','axis_ratio')):
            ax.plot(times,[h['metrics'][key] for h in r['history']],label=name)
            ax.set(xlabel='Model time',ylabel=key.replace('_',' '))
    axes[0,0].legend(fontsize=8);fig.savefig(output/'comparison.png',dpi=140);plt.close(fig)
    return result


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or sha(output/'protocol.json')!=state['protocol_sha256']:raise ValueError('requires unchanged prepared protocol')
    if any(sha(path)!=value for path,value in p['code_sha256'].items()):raise ValueError('simulation code changed since preparation')
    completed=[]
    try:
        for name in p['cases']:
            run_case(output,name,p,state,completed);completed.append(name)
        result=compare(output,p)
        write_json(output/'status.json',{**state,'state':'completed','completed_cases':completed,'all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','completed_cases':completed,'error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
