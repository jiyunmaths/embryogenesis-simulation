"""Face-normal-aware conservative diffusion on an affine sheared periodic mesh.

A standalone centered multipoint prototype, not live embryo transport. Exact
semidiscrete Fourier evolution isolates spatial error from time integration.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .resolution import write_json


def metric(shear):
    if not np.isfinite(shear):raise ValueError('shear must be finite')
    deformation=np.array([[1.,0,0],[shear,1.,0],[0,0,1.]])
    inverse=np.linalg.inv(deformation)
    return deformation,inverse@inverse.T


def face_fluxes(field,shear):
    """D=1; positive flux transfers from the forward neighbor into this cell.

    Coordinate map x=xi, y=s*xi+eta, z=zeta; det(F)=1. Flux through
    logical-a face = h^2 (F^-1 F^-T grad_logical(c))_a. Normal derivative
    is two-point; tangential derivatives are centered and face-averaged.
    """
    field=np.asarray(field,float)
    if field.ndim!=3 or len(set(field.shape))!=1 or field.shape[0]<4:
        raise ValueError('requires cubic cell array with grid >=4')
    _,tensor=metric(shear);h=1/field.shape[0]
    gradients=[(np.roll(field,-1,b)-np.roll(field,1,b))/(2*h) for b in range(3)]
    flux=[]
    for a in range(3):
        value=tensor[a,a]*(np.roll(field,-1,a)-field)/h
        for b in range(3):
            if b!=a:
                value+=tensor[a,b]*(gradients[b]+np.roll(gradients[b],-1,a))/2
        flux.append(h*h*value)
    return flux


def divergence(field,shear):
    h=1/field.shape[0]
    flux=face_fluxes(field,shear)
    return sum(q-np.roll(q,1,a) for a,q in enumerate(flux))/h**3


def symbol(n,shear,legacy=False):
    if type(n) is not int or n<4:raise ValueError('grid must be an integer >=4')
    _,tensor=metric(shear)
    theta=2*np.pi*np.fft.fftfreq(n)
    k=np.meshgrid(theta,theta,theta,indexing='ij');h=1/n
    if legacy:
        # Exact face area / Euclidean center distance on the same mesh.
        coefficients=[1/np.sqrt(1+shear*shear),np.sqrt(1+shear*shear),1.]
        return -sum(4*coefficients[a]*np.sin(k[a]/2)**2/h**2 for a in range(3))
    result=-sum(4*tensor[a,a]*np.sin(k[a]/2)**2/h**2 for a in range(3))
    for a in range(3):
        for b in range(a+1,3):
            result-=2*tensor[a,b]*np.sin(k[a])*np.sin(k[b])/h**2
    return result


def evolve(field,shear,time,legacy=False):
    if not np.isfinite(time) or time<0:raise ValueError('time must be finite and nonnegative')
    return np.fft.ifftn(np.fft.fftn(field)*np.exp(time*symbol(field.shape[0],shear,legacy))).real


def linear_patch(shear,n=8):
    deformation,_=metric(shear);h=1/n
    axis=(np.arange(n)+.5)*h
    logical=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))
    xyz=np.einsum('ij,jxyz->ixyz',deformation,logical)
    normals=np.linalg.inv(deformation).T*h*h # columns = outward face area vectors
    rows=[]
    for gradient in ([1,0,0],[0,1,0],[1,1,0],[0,0,1]):
        field=np.einsum('i,ixyz->xyz',gradient,xyz)
        for a,flux in enumerate(face_fluxes(field,shear)):
            measured=float(flux[2,2,2]) # stencil wholly inside patch, no periodic seam
            exact=float(np.asarray(gradient)@normals[:,a])
            rows.append({'gradient':gradient,'face_axis':a,'exact_flux':exact,'flux':measured,
                         'absolute_error_per_area':abs(measured-exact)/np.linalg.norm(normals[:,a])})
    return rows


def manufactured(n,shear,time):
    _,tensor=metric(shear)
    axis=(np.arange(n)+.5)/n
    xyz=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'))
    field=np.ones((n,n,n));exact=field.copy()
    # Exact cell averages, not point samples. These modes are periodic in the
    # logical parallelepiped and include mixed, normal, and z variation.
    for mode,amplitude in [([1,1,0],.12),([2,1,1],.08),([0,1,1],.06)]:
        m=np.array(mode)
        wave=amplitude*np.prod(np.sinc(m/n))*np.cos(2*np.pi*np.einsum('i,ixyz->xyz',m,xyz))
        field+=wave;exact+=wave*np.exp(-4*np.pi**2*(m@tensor@m)*time)
    return field,exact


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    protocol={'grids':[8,16,32,64],'shears':[0.,.25,.5,1.], 'time':.02,
              'pulse_grid':16,'pulse_time_over_h_squared':.01,
              'criteria':{'linear_flux_per_area_max':1e-12,'mass_drift_max':1e-12,
                          'finest_relative_l2_max':.001,'last_order_min':1.8,
                          'positive_eigenvalue_max':1e-10,'negative_concentration_tolerance':1e-12,
                          'face_symbol_agreement_max':1e-10},
              'prototype_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':__doc__,'limitations':'Periodic affine mesh with exact faces, normals, volumes and centroids. Does not validate free boundaries, irregular unstructured meshes, phase-field geometry extraction, or evolving topology. Linear affine tests exclude the periodic seam. Fourier integration is exact for this semidiscrete operator, not a production time solver.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());status=json.loads((output/'status.json').read_text())
    if status['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=status['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=p['prototype_sha256']:raise ValueError('prototype changed since preparation')
    write_json(output/'status.json',{**status,'state':'running'})
    try:
        cases=[]
        for shear in p['shears']:
            rows=[]
            for n in p['grids']:
                start,exact=manufactured(n,shear,p['time'])
                final=evolve(start,shear,p['time']);old=evolve(start,shear,p['time'],legacy=True)
                eigen=symbol(n,shear)
                rows.append({'grid':n,'relative_l2':float(np.linalg.norm((final-exact).ravel())/np.linalg.norm(exact.ravel())),
                             'legacy_relative_l2':float(np.linalg.norm((old-exact).ravel())/np.linalg.norm(exact.ravel())),
                             'mass_drift':float(abs(final.mean()-start.mean())),
                             'max_eigenvalue':float(eigen.max()),'min_eigenvalue':float(eigen.min()),
                             'smooth_minimum':float(final.min())})
            errors=[r['relative_l2'] for r in rows]
            orders=np.log2(np.array(errors[:-1])/errors[1:]).tolist()
            n=p['pulse_grid'];pulse=np.zeros((n,n,n));pulse[n//2,n//2,n//2]=1.
            dt=p['pulse_time_over_h_squared']/n**2
            advanced=evolve(pulse,shear,dt)
            derivative=divergence(pulse,shear)
            off_diagonal=derivative.copy();off_diagonal[n//2,n//2,n//2]=0
            rng=np.random.default_rng(7);random=rng.normal(size=(n,n,n))
            spectral_derivative=np.fft.ifftn(np.fft.fftn(random)*symbol(n,shear)).real
            case={'shear':shear,'refinement':rows,'orders':orders,'linear_patch':linear_patch(shear),
                  'pulse':{'time':dt,'minimum':float(advanced.min()),'mass_drift':float(abs(advanced.mean()-pulse.mean())),
                           'minimum_off_diagonal_rate':float(off_diagonal.min()),
                           'negative_entries':int(np.count_nonzero(advanced < -1e-12))},
                  'face_symbol_max_difference':float(np.max(abs(divergence(random,shear)-spectral_derivative)))}
            cases.append(case)
        c=p['criteria']
        checks={'all_linear_flux_patches_exact':all(max(r['absolute_error_per_area'] for r in case['linear_patch'])<c['linear_flux_per_area_max'] for case in cases),
                'face_divergence_matches_fourier_operator':all(case['face_symbol_max_difference']<c['face_symbol_agreement_max'] for case in cases),
                'all_masses_preserved':all(max([r['mass_drift'] for r in case['refinement']]+[case['pulse']['mass_drift']])<c['mass_drift_max'] for case in cases),
                'all_spectra_nonpositive':all(r['max_eigenvalue']<=c['positive_eigenvalue_max'] for case in cases for r in case['refinement']),
                'all_spatial_errors_decrease':all(all(a['relative_l2']>b['relative_l2'] for a,b in zip(case['refinement'][:-1],case['refinement'][1:])) for case in cases),
                'all_last_orders_above_1_8':all(case['orders'][-1]>c['last_order_min'] for case in cases),
                'all_finest_errors_below_0_1_percent':all(case['refinement'][-1]['relative_l2']<c['finest_relative_l2_max'] for case in cases),
                'all_positive_pulses_remain_nonnegative':all(case['pulse']['minimum']>=-c['negative_concentration_tolerance'] for case in cases)}
        result={'protocol':p,'cases':cases,'checks':checks,'passed':all(checks.values()),
                'live_integration_ready':False,
                'decision':'Centered face-normal correction is a diagnostic prototype. A negative exact-semigroup pulse is a spatial positivity failure; reducing the time step cannot make it a positivity-preserving generator. No live replacement is authorized by smooth accuracy alone.'}
        write_json(output/'comparison.json',result)
        lines=['# Face-normal-aware skew-mesh prototype','',f'All declared checks pass: **{result["passed"]}**','',
               '| Shear | Finest corrected L2 | Finest A/ell L2 | Last order | Positive-pulse minimum |',
               '|---:|---:|---:|---:|---:|']
        lines += [f'| {case["shear"]} | {case["refinement"][-1]["relative_l2"]:.6%} | {case["refinement"][-1]["legacy_relative_l2"]:.6%} | {case["orders"][-1]:.4f} | {case["pulse"]["minimum"]:.6g} |' for case in cases]
        lines+=['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()],'',result['decision'],'',p['limitations']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json',{**status,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**status,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
