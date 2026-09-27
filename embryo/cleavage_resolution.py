"""Controlled small-cell cleavage: voxel/time refinement and axis sensitivity."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import hashlib

import numpy as np
from scipy.integrate import quad
from scipy.ndimage import label
from scipy.special import expit

from .model import Config, Simulation, occupancy
from .projection_audit import project, relative_error
from .resolution import write_json, _steps


class AuditedSimulation(Simulation):
    def _complete_division(self, *args):
        volumes = self.volumes()
        before = np.array([volumes.sum(), volumes@self.activator, volumes@self.inhibitor])
        field = occupancy(self.phi).sum(axis=0)
        super()._complete_division(*args)
        volumes = self.volumes()
        after = np.array([volumes.sum(), volumes@self.activator, volumes@self.inhibitor])
        self.abscission_audit = {'time':self.time, 'relative_jumps':(abs(after/before-1)).tolist(),
                                'aggregate_occupancy_relative_jump':relative_error(occupancy(self.phi).sum(axis=0),field)}


def initial(config, direction, radius=.4):
    sim = AuditedSimulation(config)
    distance = np.sqrt(np.sum(sim.xyz**2,axis=0))
    sim.phi = expit(np.sqrt(2)*(radius-distance)/config.interface_width)[None].astype(np.float32)
    volume = quad(lambda r:4*np.pi*r*r*occupancy(expit(np.sqrt(2)*(radius-r)/config.interface_width)),
                  0,radius+30*config.interface_width,epsabs=1e-12,epsrel=1e-12)[0]
    sim.target[:]=volume
    sim.initial_volume=volume
    sim.abscission_audit=None
    sim.divide(0,direction)
    return sim


def prepare(output, grids=(56,72,88), dts=(.015,.0075,.00375), until=3.):
    output=Path(output)
    if output.exists(): raise FileExistsError('choose a fresh output directory')
    if len(grids)!=3 or list(grids)!=sorted(set(grids)) or len(dts)!=3 or list(dts)!=sorted(set(dts),reverse=True):
        raise ValueError('three increasing grids and decreasing time steps required')
    cases={}
    for axis,direction in [('axial',[1,0,0]),('oblique',[1,2,3])]:
        for grid in grids:
            config=Config(grid=grid,extent=2.24,max_cells=2,dt=dts[-1],steps=_steps(until,dts[-1]),signal_partition_noise=0.)
            config.validate()
            cases[f'{axis}-{grid}']={'config':asdict(config),'direction':direction}
    for dt in dts[:-1]:
        config=Config(**{**cases[f'oblique-{grids[1]}']['config'],'dt':dt,'steps':_steps(until,dt)})
        config.validate()
        cases[f'time-{dt:g}']={'config':asdict(config),'direction':[1,2,3]}
    protocol={'cases':cases,'grids':list(grids),'dts':list(dts),'until':until,'radius':.4,
              'sample_interval':.15,'probe_grids':[grids[-1],grids[-1]+24],'interpolation_order':3,
              'time_cases':[f'time-{dt:g}' for dt in dts[:-1]]+[f'oblique-{grids[1]}'],
              'model_sha256':hashlib.sha256(Path(__file__).with_name('model.py').read_bytes()).hexdigest(),
              'criteria':{'field_relative_max':.01,'initial_calibration_max':.0025,
                          'timing_difference_max':.05,'orientation_timing_difference_max':.1,
                          'axis_ratio_relative_max':.01,'conservation_jump_max':1e-6,
                          'volume_error_max':.05,'boundary_max':.01,'finest_radius_min':4.},
              'scope':'One prescribed division of an analytic radius-0.4 cell; common continuous volume target; fixed interface width and domain. All standard coupling enabled, no partition noise. With two cells below fate competence this does not test differentiation. Prescribed axes test lattice sensitivity, not spontaneous orientation. No full developmental or sharp-interface convergence claim.'}
    for case in cases.values():_steps(protocol['sample_interval'],case['config']['dt'])
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run_case(output,name,p):
    case=p['cases'][name]
    sim=initial(Config(**case['config']),case['direction'],p['radius'])
    path=output/name;path.mkdir()
    history=[]; frames=[]; max_volume=0.; min_radius=float('inf'); clipping=0.
    interval=_steps(p['sample_interval'],sim.config.dt)
    started=time.monotonic()
    for step in range(sim.config.steps+1):
        volumes=sim.volumes()
        max_volume=max(max_volume,float(np.max(abs(volumes/sim.target-1))))
        min_radius=min(min_radius,float(np.min((3*volumes/(4*np.pi))**(1/3))/sim.dx))
        clipping=max(clipping,sim.clipped_fraction)
        if step%interval==0 or step==sim.config.steps:
            row=sim.metrics()
            row['daughter_components']=[label(field>=.5)[1] for field in sim.phi] if len(sim.phi)==2 else None
            history.append(row)
            if step==0 or step==sim.config.steps:
                frames.append({'metrics':row,'cells':sim.surfaces(max_points=200),'graph':sim.graph_snapshot()})
            write_json(path/'history.json',history)
            print(f'{name}: t={sim.time:.3f}, cells={len(sim.phi)}',flush=True)
        if step==sim.config.steps:break
        sim.step()
        if sim.abscission_audit and 'daughter_components' not in sim.abscission_audit:
            sim.abscission_audit['daughter_components']=[label(field>=.5)[1] for field in sim.phi]
    sim.checkpoint(path/'final_state.npz')
    report={'history':history,'abscission':sim.abscission_audit,'max_volume_error':max_volume,
            'min_radius_grid_cells':min_radius,'max_clipped_fraction':clipping,'elapsed_seconds':time.monotonic()-started}
    write_json(path/'analysis.json',report)
    from .surface import viewer_template
    payload=json.dumps({'config':asdict(sim.config),'frames':frames},allow_nan=False)
    (path/'viewer.html').write_text(viewer_template().replace('A freely evolving 3D aggregate of deformable cells.','Prescribed small-cell cleavage: initial and final surfaces.').replace('__SIMULATION_DATA__',payload))
    return report


def compare(output,p):
    reports={n:json.loads((output/n/'analysis.json').read_text()) for n in p['cases']}
    limits=p['criteria']; checks={}; pairs=[]; calibration=[]
    def pair(left,right,probe):
        a,b=(Simulation.restore(output/n/'final_state.npz') for n in (left,right))
        matched=np.array_equal(a.ids,b.ids) and len(a.ids)==2
        ra,rb=reports[left],reports[right]
        return {'cases':[left,right],'probe_grid':probe,'matched_daughters':matched,
                'field_relative_l2':relative_error(project(a,probe,3),project(b,probe,3)) if matched else None,
                'abscission_time_difference':abs(ra['abscission']['time']-rb['abscission']['time']) if ra['abscission'] and rb['abscission'] else None,
                'axis_ratio_relative_difference':abs(ra['history'][-1]['axis_ratio']/rb['history'][-1]['axis_ratio']-1)}
    for axis in ('axial','oblique'):
        names=[f'{axis}-{g}' for g in p['grids']]
        for probe in p['probe_grids']:
            rows=[pair(a,b,probe) for a,b in zip(names[:-1],names[1:])]; pairs+=rows
            fine=rows[-1]
            checks[f'{axis}_{probe}_field_below_1_percent']=fine['field_relative_l2'] is not None and fine['field_relative_l2']<limits['field_relative_max']
            checks[f'{axis}_{probe}_field_decreases']=all(r['field_relative_l2'] is not None for r in rows) and rows[-1]['field_relative_l2']<rows[0]['field_relative_l2']
        checks[f'{axis}_timing_agrees']=fine['abscission_time_difference'] is not None and fine['abscission_time_difference']<limits['timing_difference_max']
        checks[f'{axis}_shape_agrees']=fine['axis_ratio_relative_difference']<limits['axis_ratio_relative_max']
    time_rows=[pair(a,b,p['grids'][1]) for a,b in zip(p['time_cases'][:-1],p['time_cases'][1:])];pairs+=time_rows
    checks['time_field_below_1_percent']=time_rows[-1]['field_relative_l2'] is not None and time_rows[-1]['field_relative_l2']<limits['field_relative_max']
    checks['time_field_decreases']=all(r['field_relative_l2'] is not None for r in time_rows) and time_rows[-1]['field_relative_l2']<time_rows[0]['field_relative_l2']
    checks['time_abscission_agrees']=time_rows[-1]['abscission_time_difference'] is not None and time_rows[-1]['abscission_time_difference']<limits['timing_difference_max']
    finest=[reports[f'{axis}-{p["grids"][-1]}'] for axis in ('axial','oblique')]
    checks['orientation_timing_agrees']=all(r['abscission'] for r in finest) and abs(finest[0]['abscission']['time']-finest[1]['abscission']['time'])<limits['orientation_timing_difference_max']
    checks['orientation_shape_agrees']=abs(finest[0]['history'][-1]['axis_ratio']/finest[1]['history'][-1]['axis_ratio']-1)<limits['axis_ratio_relative_max']
    for probe in p['probe_grids']:
        config=Config(**{**p['cases'][f'axial-{p["grids"][0]}']['config'],'grid':probe})
        exact=occupancy(initial(config,[1,0,0],p['radius']).phi)
        for grid in p['grids']:
            sim=initial(Config(**p['cases'][f'axial-{grid}']['config']),[1,0,0],p['radius'])
            calibration.append({'grid':grid,'probe':probe,'relative_l2':relative_error(project(sim,probe,3),exact)})
    checks['initial_reconstruction_calibrated']=max(r['relative_l2'] for r in calibration)<limits['initial_calibration_max']
    checks['all_divisions_complete']=all(r['abscission'] for r in reports.values())
    checks['abscission_conserves_volume_amount_and_occupancy']=all(r['abscission'] and max(r['abscission']['relative_jumps']+[r['abscission']['aggregate_occupancy_relative_jump']])<limits['conservation_jump_max'] for r in reports.values())
    checks['daughters_connected']=all(r['abscission'] and r['abscission']['daughter_components']==[1,1] and all(h['daughter_components'] in (None,[1,1]) for h in r['history']) for r in reports.values())
    checks['volumes_within_5_percent']=all(r['max_volume_error']<limits['volume_error_max'] for r in reports.values())
    checks['boundaries_clear']=all(h['boundary_occupancy']<limits['boundary_max'] for r in reports.values() for h in r['history'])
    checks['no_clipping']=all(r['max_clipped_fraction']==0 for r in reports.values())
    checks['finest_daughters_resolved']=all(r['min_radius_grid_cells']>=limits['finest_radius_min'] for r in finest)
    result={'protocol':p,'checks':checks,'passed':bool(all(checks.values())),'pairs':pairs,'initial_calibration':calibration}
    write_json(output/'comparison.json',result)
    (output/'RESULTS.md').write_text('# Controlled cleavage resolution\n\nAll declared checks pass: **'+str(result['passed'])+'**\n\n'+'\n'.join(f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items())+'\n\n'+p['scope']+'\n')
    return result


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());status=json.loads((output/'status.json').read_text())
    if status['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=status['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    if hashlib.sha256(Path(__file__).with_name('model.py').read_bytes()).hexdigest()!=p['model_sha256']:
        raise ValueError('model changed after preparation')
    completed=[]
    try:
        for name in p['cases']:
            write_json(output/'status.json',{**status,'state':'running','case':name,'completed_cases':completed})
            run_case(output,name,p);completed.append(name)
        result=compare(output,p)
        write_json(output/'status.json',{**status,'state':'completed','completed_cases':completed,'all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**status,'state':'failed','completed_cases':completed,'error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
