"""Positive wall-tangential conductance correction for a uniform skew mesh.

Derived from the leading boundary-row error and the conormal condition.
Known affine geometry only; neither an empirical parameter fit nor a live fix.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .skew_boundary import closed_transport, average, P, C, manufactured
from .transport import conservative_transport
from .resolution import write_json


def corrected_transport(n,shear):
    base=closed_transport(n,shear);t=abs(shear);h=1/n
    gamma_x=.5*t*(1-t)
    gamma_y=.5*t*(1-t/(1+shear*shear))
    ids=np.arange(n**3).reshape((n,n,n));rows=[];columns=[];values=[]
    for wall in (0,n-1):
        for a,b,coefficient in ((ids[wall,:-1,:],ids[wall,1:,:],gamma_x),
                                 (ids[:-1,wall,:],ids[1:,wall,:],gamma_y)):
            if coefficient==0:continue
            a,b=a.ravel(),b.ravel()
            rows.extend([a,b]);columns.extend([b,a]);values.extend([np.full(len(a),h*coefficient)]*2)
    if rows:
        extra=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(columns))),shape=base.conductance.shape).tocsr()
        return conservative_transport(base.conductance+extra,base.volumes)
    return base


def reference(n,shear,family):
    if family=='xi':return manufactured(n,shear)
    if family!='eta':raise ValueError('unknown boundary family')
    k=shear/(1+shear*shear)
    x=[average(P.deriv(j),n)[:,None,None] for j in range(4)]
    y=[average(C.deriv(j),n)[None,:,None] for j in range(3)]
    z=(np.cos(2*np.pi*(np.arange(n)+.5)/n)*np.sinc(1/n))[None,None,:]
    f=64*(x[0]+k*y[0]*x[1])*z
    laplace=64*(x[2]+k*y[0]*x[3]-2*shear*k*y[1]*x[2]+(1+shear*shear)*k*y[2]*x[1])*z-4*np.pi**2*f
    return f.ravel(),laplace.ravel()


def boundary_residual(shear,family):
    v=np.linspace(0,1,37);errors=[]
    k=shear if family=='xi' else shear/(1+shear*shear)
    for wall in (0.,1.):
        if family=='xi':
            gx=k*C.deriv()(wall)*P.deriv()(v);gy=P.deriv()(v)+k*C(wall)*P.deriv(2)(v)
            errors.extend(abs(gx-shear*gy))
            gx=k*C.deriv()(v)*P.deriv()(wall);gy=P.deriv()(wall)+k*C(v)*P.deriv(2)(wall)
            errors.extend(abs(-shear*gx+(1+shear*shear)*gy))
        else:
            gx=P.deriv()(v)+k*C(wall)*P.deriv(2)(v);gy=k*C.deriv()(wall)*P.deriv()(v)
            errors.extend(abs(-shear*gx+(1+shear*shear)*gy))
            gx=P.deriv()(wall)+k*C(v)*P.deriv(2)(wall);gy=k*C.deriv()(v)*P.deriv()(wall)
            errors.extend(abs(gx-shear*gy))
    return float(max(errors))


def evolve(transport,n,shear,family,time):
    f,laplace=reference(n,shear,family);initial=1+.2*f;source=.2*(-f-laplace)
    operator=sparse.bmat([[transport.delta,sparse.csr_matrix(source[:,None])],
                         [sparse.csr_matrix((1,n**3)),sparse.csr_matrix([[-1.]])]],format='csr')
    final=expm_multiply(time*operator,np.r_[initial,1.])[:-1]
    return initial,final,1+.2*np.exp(-time)*f,source


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    protocol={'grids':[8,16,32],'confirmation_grid':64,'confirmation_shears':[.25,.5],
              'shears':[0.,.25,.5,1.],'families':['xi','eta'],'time':.02,
              'criteria':{'error_max':.001,'order_min':1.8,'mass_drift_max':1e-12,
                          'positive_tolerance':1e-12,'invariant_residual_max':1e-9},
              'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ('skew_boundary_correction.py','skew_boundary.py','positive_skew_flux.py','transport.py')},
              'rule':'Add h*|s|(1-|s|)/2 to eta-direction edges in xi-wall layers, and h*|s|(1-|s|/(1+s^2))/2 to xi-direction edges in eta-wall layers. No coefficient fitting.',
              'scope':__doc__,'limits':'Uniform-volume static affine mesh, |s|<=1, flat reflecting walls. Does not establish arbitrary-corner, unstructured geometry, or moving-boundary consistency. Original uncorrected report retained.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:raise ValueError('requires unchanged protocol')
    if any(hashlib.sha256(Path(k).read_bytes()).hexdigest()!=v for k,v in p['source_sha256'].items()):raise ValueError('sources changed')
    write_json(output/'status.json',{**state,'state':'running'})
    try:
        cases=[]
        for shear in p['shears']:
            for family in p['families']:
                grids=p['grids']+([p['confirmation_grid']] if shear in p['confirmation_shears'] else [])
                rows=[]
                for n in grids:
                    transport=corrected_transport(n,shear)
                    initial,final,exact,source=evolve(transport,n,shear,family,p['time'])
                    ijk=np.indices((n,n,n));wall=np.any((ijk==0)|(ijk==n-1),axis=0).ravel();error=final-exact
                    rows.append({'grid':n,'relative_l2':float(np.linalg.norm(error)/np.linalg.norm(exact)),
                                 'boundary_relative_l2':float(np.linalg.norm(error[wall])/np.linalg.norm(exact[wall])),
                                 'mass_drift':float(abs(transport.volumes@(final-initial))),
                                 'source_integral':float(transport.volumes@source),
                                 'constant_residual':float(np.max(abs(transport.delta@np.ones(n**3)))),
                                 'mass_residual':float(np.max(abs(transport.volumes@transport.delta)))})
                    print(f'{family}, shear={shear}, n={n}: {rows[-1]["relative_l2"]:.5%}',flush=True)
                errors=np.array([r['relative_l2'] for r in rows]);orders=np.log2(errors[:-1]/errors[1:]).tolist()
                cases.append({'shear':shear,'family':family,'boundary_reference_residual':boundary_residual(shear,family),'refinement':rows,'orders':orders})
        pulses=[]
        for shear in p['shears']:
            n=16;t=corrected_transport(n,shear)
            off=t.delta.copy();off.setdiag(0);off.eliminate_zeros()
            for location in ((0,8,8),(8,0,8),(0,0,0)):
                start=np.zeros(n**3);start.reshape((n,n,n))[location]=1
                final=expm_multiply(t.delta*(.01/n**2),start)
                pulses.append({'shear':shear,'location':location,'minimum':float(final.min()),'maximum':float(final.max()),
                               'mass_drift':float(abs(t.volumes@(final-start))),
                               'minimum_off_diagonal_rate':float(off.data.min())})
        c=p['criteria'];checks={
            'both_manufactured_wall_families_exact':all(x['boundary_reference_residual']<1e-12 for x in cases),
            'mass_and_constants_preserved':all(max(r['constant_residual'],r['mass_residual'])<c['invariant_residual_max'] for x in cases for r in x['refinement']),
            'source_integral_and_mass_drift_below_1e_12':all(max(abs(r['source_integral']),r['mass_drift'])<c['mass_drift_max'] for x in cases for r in x['refinement']),
            'all_refinement_errors_decrease':all(all(a['relative_l2']>b['relative_l2'] for a,b in zip(x['refinement'][:-1],x['refinement'][1:])) for x in cases),
            'all_original_16_to_32_orders_above_1_8':all(x['orders'][1]>c['order_min'] for x in cases),
            'all_64_confirmation_orders_above_1_8':all(x['orders'][-1]>c['order_min'] for x in cases if x['shear'] in p['confirmation_shears']),
            'all_32_errors_below_0_1_percent':all(x['refinement'][2]['relative_l2']<c['error_max'] for x in cases),
            'rates_nonnegative':all(x['minimum_off_diagonal_rate']>=0 for x in pulses),
            'face_and_corner_pulses_positive_conservative':all(x['minimum']>=-c['positive_tolerance'] and x['maximum']<=1+c['positive_tolerance'] and x['mass_drift']<c['mass_drift_max'] for x in pulses)}
        result={'protocol':p,'cases':cases,'pulses':pulses,'checks':checks,'passed':all(checks.values()),'live_integration_ready':False}
        write_json(output/'comparison.json',result)
        lines=['# Positive reflecting-wall correction','',f'All declared checks pass: **{result["passed"]}**','',
               '| Shear | Wall family | 32³ error | 16–32 order | 64³ error | 32–64 order |',
               '|---:|---|---:|---:|---:|---:|']
        for x in cases:
            fine=f'{x["refinement"][3]["relative_l2"]:.6%}' if len(x['refinement'])==4 else '—'
            order=f'{x["orders"][2]:.4f}' if len(x['orders'])==3 else '—'
            lines.append(f'| {x["shear"]} | {x["family"]} | {x["refinement"][2]["relative_l2"]:.6%} | {x["orders"][1]:.4f} | {fine} | {order} |')
        lines+=['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()],'',p['limits']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n');write_json(output/'status.json',{**state,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
