"""Calibrated native-phase reconstruction for archived cleavage measurements.

Reconstruct phi, then evaluate h(phi). No clipping, solver changes, or overwrite
of the original occupancy-interpolation result. Analytic holdouts calibrate the
measurement; evolved-field cross-checks do not provide an exact error bound.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import expit

from .model import Simulation, occupancy
from .projection_audit import relative_error
from .resolution import write_json


def reconstruct(phi, extent, probe, order=3):
    if order not in (3,5) or type(probe) is not int or probe<2:
        raise ValueError('requires cubic/quintic order and integer probe >=2')
    phi=np.asarray(phi)
    if phi.ndim!=4 or len(set(phi.shape[1:]))!=1:
        raise ValueError('requires cell-by-cubic-grid phase fields')
    grid=phi.shape[1]
    if grid==probe:return phi.copy()
    axis=(np.arange(probe)+.5)*2*extent/probe-extent
    xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))
    indices=(xyz+extent)/(2*extent/grid)-.5
    return np.stack([map_coordinates(f,indices,order=order,mode='nearest',prefilter=True) for f in phi])


def analytic(grid,extent,width,radius,center,stretch,angle):
    axis=(np.arange(grid)+.5)*2*extent/grid-extent
    xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))-np.array(center)[:,None,None,None]
    c,s=np.cos(angle),np.sin(angle)
    rotation=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
    local=np.einsum('ij,jxyz->ixyz',rotation,xyz)/np.array(stretch)[:,None,None,None]
    distance=np.sqrt(np.sum(local**2,axis=0))
    return expit(np.sqrt(2)*(radius-distance)/width)[None].astype(np.float32)


def prepare(source,output):
    source,output=Path(source).resolve(),Path(output)
    if output.exists():raise FileExistsError('choose a fresh validation directory')
    old=json.loads((source/'protocol.json').read_text())
    names=[f'{axis}-{g}' for axis in ('axial','oblique') for g in old['grids']]
    paths=[source/'protocol.json',source/'comparison.json']+[source/n/'final_state.npz' for n in names]
    configs=[old['cases'][n]['config'] for n in names]
    if len({(c['extent'],c['interface_width']) for c in configs})!=1:
        raise ValueError('requires common physical domain and interface width')
    protocol={'source':str(source),'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'grids':old['grids'],'extent':configs[0]['extent'],'interface_width':configs[0]['interface_width'],
              'orders':[3,5],'probes':[88,112,128],
              'criteria':{'analytic_error_max':.0025,'finest_field_error_max':.01,'diagnostic_spread_max':.001},
              'analytic_geometries':[
                  {'name':'original_sphere','radius':.4,'center':[0,0,0],'stretch':[1,1,1],'angle':0},
                  {'name':'translated_smaller_holdout','radius':.36,'center':[.017,-.023,.031],'stretch':[1,1,1],'angle':0},
                  {'name':'rotated_ellipsoid_holdout','radius':.44,'center':[-.029,.013,.019],'stretch':[.9,1.05,1.1],'angle':.43}],
              'scope':__doc__,'selection':'Exploratory original-sphere calibration suggested native-phi reconstruction. Translated/deformed holdouts and archived evolved-field checks are specified here before their evaluation. Cubic is primary; quintic and three probes are required cross-checks. Original historical failure remains unchanged.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());status=json.loads((output/'status.json').read_text())
    if status['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=status['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    def unchanged():return all(hashlib.sha256(Path(path).read_bytes()).hexdigest()==value for path,value in p['source_sha256'].items())
    if not unchanged():raise ValueError('archived source changed')
    try:
        write_json(output/'status.json',{**status,'state':'running','stage':'analytic_holdouts'})
        calibration=[]
        for geometry in p['analytic_geometries']:
            kwargs={k:v for k,v in geometry.items() if k!='name'}
            for probe in p['probes']:
                exact=occupancy(analytic(probe,p['extent'],p['interface_width'],**kwargs))
                for grid in p['grids']:
                    field=analytic(grid,p['extent'],p['interface_width'],**kwargs)
                    for order in p['orders']:
                        projected=reconstruct(field,p['extent'],probe,order)
                        calibration.append({'geometry':geometry['name'],'grid':grid,'probe':probe,'order':order,
                                            'relative_l2':relative_error(occupancy(projected),exact),
                                            'phase_range':[float(projected.min()),float(projected.max())]})
            print(f'calibrated {geometry["name"]}',flush=True)
        write_json(output/'calibration.json',calibration)
        write_json(output/'status.json',{**status,'state':'running','stage':'archived_cleavage_comparison'})
        rows=[]
        for axis in ('axial','oblique'):
            names=[f'{axis}-{g}' for g in p['grids']]
            states=[Simulation.restore(Path(p['source'])/n/'final_state.npz') for n in names]
            if any(len(s.ids)!=2 or not np.array_equal(s.ids,states[-1].ids) or not np.isclose(s.time,states[-1].time) for s in states):
                raise ValueError('requires matched completed daughter states')
            for probe in p['probes']:
                for order in p['orders']:
                    fields=[occupancy(reconstruct(s.phi,p['extent'],probe,order)) for s in states]
                    for i in (0,1):
                        rows.append({'axis':axis,'cases':names[i:i+2],'probe':probe,'order':order,
                                     'relative_l2':relative_error(fields[i],fields[i+1]),
                                     'occupancy_range':[float(fields[i].min()),float(fields[i].max())]})
            print(f'compared {axis} daughter fields',flush=True)
        limits=p['criteria'];checks={}
        for order in p['orders']:
            checks[f'order_{order}_all_analytic_calibrations_below_0_25_percent']=max(r['relative_l2'] for r in calibration if r['order']==order)<limits['analytic_error_max']
        for axis in ('axial','oblique'):
            fine=[r for r in rows if r['axis']==axis and r['cases'][-1]==f'{axis}-{p["grids"][-1]}']
            coarse={(r['probe'],r['order']):r for r in rows if r['axis']==axis and r not in fine}
            values=[r['relative_l2'] for r in fine]
            checks[axis+'_all_finest_fields_below_1_percent']=max(values)<limits['finest_field_error_max']
            checks[axis+'_all_fields_decrease']=all(r['relative_l2']<coarse[r['probe'],r['order']]['relative_l2'] for r in fine)
            checks[axis+'_diagnostic_spread_below_0_1_percentage_point']=max(values)-min(values)<limits['diagnostic_spread_max']
        checks['archived_inputs_unchanged']=unchanged()
        original=json.loads((Path(p['source'])/'comparison.json').read_text())
        checks['other_original_checks_pass']=all(value for key,value in original['checks'].items() if key!='initial_reconstruction_calibrated')
        result={'protocol':p,'calibration':calibration,'pairs':rows,'checks':checks,'passed':all(checks.values()),
                'original_screen_passed':original['passed'],
                'interpretation':'Separate calibrated measurement validation on archived dynamics. No original pass/fail is overwritten; analytic calibration is not a rigorous bound on evolved geometry, and no new physical-model validation is inferred.'}
        write_json(output/'comparison.json',result)
        lines=['# Cleavage measurement validation','',f'All revised measurement checks pass: **{result["passed"]}**','',result['interpretation'],'',
               '| Geometry | Maximum cubic analytic error | Maximum quintic analytic error |','|---|---:|---:|']
        for g in p['analytic_geometries']:
            values=[max(r['relative_l2'] for r in calibration if r['geometry']==g['name'] and r['order']==order) for order in p['orders']]
            lines.append(f'| {g["name"]} | {values[0]:.4%} | {values[1]:.4%} |')
        lines+=['','| Orientation | Minimum finest-pair error | Maximum finest-pair error |','|---|---:|---:|']
        for axis in ('axial','oblique'):
            values=[r['relative_l2'] for r in rows if r['axis']==axis and r['cases'][-1]==f'{axis}-{p["grids"][-1]}']
            lines.append(f'| {axis} | {min(values):.4%} | {max(values):.4%} |')
        lines+=['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()]]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json',{**status,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**status,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run'])
    parser.add_argument('--source',type=Path,default=Path('outputs/cleavage-resolution'));parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.source,args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
