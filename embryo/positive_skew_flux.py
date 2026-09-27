"""Positive conservative wider-stencil diffusion for affine shears |s| <= 1.

The tensor is decomposed into nonnegative directional second differences.
Diagonal exchanges are routed equally along two grid paths to reconstruct
face fluxes. This is not an arbitrary-mesh or live-cell transport replacement.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .resolution import write_json
from .skew_flux import metric, manufactured


def directions(shear):
    if not np.isfinite(shear) or abs(shear)>1:
        raise ValueError('this stencil requires finite |shear| <= 1')
    t=abs(float(shear));r=-1 if shear>=0 else 1
    return [((1,0,0),1-t),((0,1,0),1+shear*shear-t),((0,0,1),1.),((1,r,0),t)]


def neighbor(field,direction):
    return np.roll(field,tuple(-v for v in direction),axis=(0,1,2))


def validate(field):
    field=np.asarray(field,float)
    if field.ndim!=3 or len(set(field.shape))!=1 or field.shape[0]<4 or not np.isfinite(field).all():
        raise ValueError('requires finite cubic field, grid >=4')
    return field


def divergence(field,shear):
    field=validate(field);n=field.shape[0]
    return n*n*sum(w*(neighbor(field,d)+neighbor(field,tuple(-x for x in d))-2*field) for d,w in directions(shear))


def face_fluxes(field,shear):
    field=validate(field);h=1/field.shape[0];terms=directions(shear)
    q=[h*terms[a][1]*(np.roll(field,-1,a)-field) for a in range(3)]
    d,w=terms[-1];r=d[1];diagonal=h*w*(neighbor(field,d)-field)
    # Each diagonal exchange is half-routed x-then-y and half y-then-x.
    q[0]+=.5*(diagonal+np.roll(diagonal,r,axis=1))
    if r==1:
        q[1]+=.5*(diagonal+np.roll(diagonal,1,axis=0))
    else:
        shifted=np.roll(diagonal,-1,axis=1)
        q[1]-=.5*(shifted+np.roll(shifted,1,axis=0))
    return q


def symbol(n,shear):
    if type(n) is not int or n<4:raise ValueError('requires integer grid >=4')
    angles=np.meshgrid(*([2*np.pi*np.fft.fftfreq(n)]*3),indexing='ij')
    return -4*n*n*sum(w*np.sin(sum(d[a]*angles[a] for a in range(3))/2)**2 for d,w in directions(shear))


def exact_evolve(field,shear,time):
    field=validate(field)
    if not np.isfinite(time) or time<0:raise ValueError('time must be nonnegative')
    return np.fft.ifftn(np.fft.fftn(field)*np.exp(time*symbol(field.shape[0],shear))).real


def positivity_step(n,shear):
    return 1/(2*n*n*sum(w for _,w in directions(shear)))


def ssprk2(field,shear,time,steps):
    """Convex-combination Euler stages; no clipping or renormalization."""
    field=validate(field).copy()
    if type(steps) is not int or steps<1 or not np.isfinite(time) or time<=0:
        raise ValueError('positive time and integer step count required')
    n=field.shape[0];dt=time/steps
    if dt>positivity_step(n,shear)*(1+1e-14):raise ValueError('time step exceeds positivity bound')
    terms=directions(shear);center=1-2*dt*n*n*sum(w for _,w in terms)
    if center<0:raise ValueError('time step has negative central weight')
    def euler(value):
        return center*value+dt*n*n*sum(w*(neighbor(value,d)+neighbor(value,tuple(-x for x in d))) for d,w in terms)
    for _ in range(steps):
        field=.5*field+.5*euler(euler(field))
    return field


def linear_patch(shear,n=8):
    deformation,_=metric(shear);h=1/n
    axis=(np.arange(n)+.5)*h
    xyz=np.einsum('ij,jxyz->ixyz',deformation,np.stack(np.meshgrid(axis,axis,axis,indexing='ij')))
    normals=np.linalg.inv(deformation).T*h*h;rows=[]
    for gradient in ([1,0,0],[0,1,0],[1,1,0],[0,0,1]):
        values=np.einsum('i,ixyz->xyz',gradient,xyz)
        for a,q in enumerate(face_fluxes(values,shear)):
            exact=float(np.array(gradient)@normals[:,a]);measured=float(q[3,3,3])
            rows.append({'gradient':gradient,'face_axis':a,'exact_flux':exact,'flux':measured,
                         'error_per_area':abs(measured-exact)/np.linalg.norm(normals[:,a])})
    return rows


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    protocol={'grids':[8,16,32,64],'shears':[0.,.25,.5,1.], 'time':.02,
              'pulse_grid':16,'pulse_time_over_h_squared':.01,'time_refinement_duration':.01,
              'criteria':{'linear_flux_per_area_max':1e-12,'mass_drift_max':1e-12,
                          'finest_relative_l2_max':.001,'last_order_min':1.8,
                          'positive_eigenvalue_max':1e-10,'negative_concentration_tolerance':1e-12,
                          'operator_agreement_max':1e-10},
              'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ('positive_skew_flux.py','skew_flux.py')},
              'scope':__doc__,'limits':'Fixed periodic affine meshes only; diagonal neighbors are available here. Independent shear range limited to |s|<=1. No general unstructured geometry, reflecting walls, evolving topology, reactions, or phase-field face extraction validation.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    if any(hashlib.sha256(Path(k).read_bytes()).hexdigest()!=v for k,v in p['source_sha256'].items()):raise ValueError('prototype sources changed')
    write_json(output/'status.json',{**state,'state':'running'})
    try:
        cases=[]
        for shear in p['shears']:
            refinement=[]
            for n in p['grids']:
                initial,exact=manufactured(n,shear,p['time']);final=exact_evolve(initial,shear,p['time'])
                refinement.append({'grid':n,'relative_l2':float(np.linalg.norm((final-exact).ravel())/np.linalg.norm(exact.ravel())),
                                   'mass_drift':float(abs(final.mean()-initial.mean())),'maximum_eigenvalue':float(symbol(n,shear).max())})
            errors=np.array([r['relative_l2'] for r in refinement]);orders=np.log2(errors[:-1]/errors[1:]).tolist()
            n=p['pulse_grid'];pulse=np.zeros((n,n,n));pulse[n//2,n//2,n//2]=1.
            duration=p['pulse_time_over_h_squared']/n**2
            exact_pulse=exact_evolve(pulse,shear,duration)
            bound=positivity_step(n,shear)
            numerical_pulse=ssprk2(pulse,shear,10*.9*bound,10)
            test=np.random.default_rng(7).random((n,n,n));rhs=divergence(test,shear)
            flux=face_fluxes(test,shear);face_rhs=n**3*sum(q-np.roll(q,1,a) for a,q in enumerate(flux))
            spectral_rhs=np.fft.ifftn(np.fft.fftn(test)*symbol(n,shear)).real
            initial,_=manufactured(n,shear,p['time_refinement_duration'])
            reference=exact_evolve(initial,shear,p['time_refinement_duration'])
            base=int(np.ceil(p['time_refinement_duration']/(.9*bound)));time_rows=[]
            for steps in (base,2*base,4*base):
                value=ssprk2(initial,shear,p['time_refinement_duration'],steps)
                time_rows.append({'steps':steps,'relative_l2':float(np.linalg.norm((value-reference).ravel())/np.linalg.norm(reference.ravel())),
                                  'mass_drift':float(abs(value.mean()-initial.mean())),'minimum':float(value.min())})
            time_errors=np.array([r['relative_l2'] for r in time_rows])
            cases.append({'shear':shear,'direction_weights':[{'direction':d,'weight':w} for d,w in directions(shear)],
                          'linear_patch':linear_patch(shear),'refinement':refinement,'orders':orders,
                          'exact_pulse_minimum':float(exact_pulse.min()),'rk2_pulse_minimum':float(numerical_pulse.min()),
                          'pulse_mass_drift':float(max(abs(exact_pulse.mean()-pulse.mean()),abs(numerical_pulse.mean()-pulse.mean()))),
                          'face_operator_max_difference':float(np.max(abs(rhs-face_rhs))),
                          'symbol_operator_max_difference':float(np.max(abs(rhs-spectral_rhs))),
                          'time_refinement':time_rows,'time_orders':np.log2(time_errors[:-1]/time_errors[1:]).tolist()})
        c=p['criteria']
        checks={'all_rates_nonnegative':all(r['weight']>=0 for case in cases for r in case['direction_weights']),
                'all_linear_face_fluxes_exact':all(max(r['error_per_area'] for r in case['linear_patch'])<c['linear_flux_per_area_max'] for case in cases),
                'face_and_symbol_operators_agree':all(max(case['face_operator_max_difference'],case['symbol_operator_max_difference'])<c['operator_agreement_max'] for case in cases),
                'all_mass_drifts_below_1e_12':all(max([r['mass_drift'] for r in case['refinement']+case['time_refinement']]+[case['pulse_mass_drift']])<c['mass_drift_max'] for case in cases),
                'spectra_nonpositive':all(r['maximum_eigenvalue']<=c['positive_eigenvalue_max'] for case in cases for r in case['refinement']),
                'all_spatial_errors_decrease':all(all(a['relative_l2']>b['relative_l2'] for a,b in zip(case['refinement'][:-1],case['refinement'][1:])) for case in cases),
                'spatial_last_orders_above_1_8':all(case['orders'][-1]>c['last_order_min'] for case in cases),
                'finest_spatial_errors_below_0_1_percent':all(case['refinement'][-1]['relative_l2']<c['finest_relative_l2_max'] for case in cases),
                'exact_pulses_nonnegative':all(case['exact_pulse_minimum']>=-c['negative_concentration_tolerance'] for case in cases),
                'rk2_pulses_nonnegative':all(case['rk2_pulse_minimum']>=-c['negative_concentration_tolerance'] for case in cases),
                'temporal_last_orders_above_1_8':all(case['time_orders'][-1]>c['last_order_min'] for case in cases)}
        result={'protocol':p,'cases':cases,'checks':checks,'passed':all(checks.values()),'live_integration_ready':False}
        write_json(output/'comparison.json',result)
        lines=['# Positive wider-stencil skew diffusion','',f'All declared checks pass: **{result["passed"]}**','',
               '| Shear | Finest spatial error | Spatial order | Temporal order | Exact pulse minimum | RK2 pulse minimum |',
               '|---:|---:|---:|---:|---:|---:|']
        lines += [f'| {r["shear"]} | {r["refinement"][-1]["relative_l2"]:.6%} | {r["orders"][-1]:.4f} | {r["time_orders"][-1]:.4f} | {r["exact_pulse_minimum"]:.3g} | {r["rk2_pulse_minimum"]:.3g} |' for r in cases]
        lines += ['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()],'',p['limits']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n');write_json(output/'status.json',{**state,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
