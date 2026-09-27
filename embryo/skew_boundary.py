"""Reflecting-wall validation of the positive directional skew-mesh stencil.

Deleting links outside the domain preserves conservation/positivity; accuracy
of the resulting boundary closure is tested independently, not assumed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .positive_skew_flux import directions
from .resolution import write_json
from .transport import conservative_transport


P=Polynomial([0,0,0,1,-3,3,-1])
C=Polynomial([0,1,-3,2])


def closed_transport(n,shear):
    if type(n) is not int or n<4:raise ValueError('grid must be an integer >=4')
    h=1/n;index=np.arange(n**3).reshape((n,n,n))
    rows=[];cols=[];values=[]
    for direction,weight in directions(shear):
        if weight==0:continue
        left=[];right=[]
        for d in direction:
            left.append(slice(max(0,-d),min(n,n-d)))
            right.append(slice(max(0,d),min(n,n+d)))
        a=index[tuple(left)].ravel();b=index[tuple(right)].ravel()
        rows.extend([a,b]);cols.extend([b,a]);values.extend([np.full(len(a),h*weight)]*2)
    conductance=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(n**3,n**3)).tocsr()
    return conservative_transport(conductance,np.full(n**3,h**3))


def average(poly,n):
    edges=np.linspace(0,1,n+1);integral=poly.integ()
    return np.diff(integral(edges))*n


def manufactured(n,shear):
    """Exact cell averages for f=64 [P(eta)+s C(xi) P'(eta)] cos(2pi z).

    At xi walls: f_xi=s f_eta, hence (B grad f)_xi=0. At eta
    walls P=P'=P''=0, and at z walls sin(2pi z)=0. Tangential
    gradients on xi walls are generally nonzero. f has zero volume mean.
    Exact solution c=1+0.2 exp(-t) f, source=0.2 exp(-t)(-f-Lf).
    """
    x=[average(C.deriv(k),n)[:,None,None] for k in range(3)]
    y=[average(P.deriv(k),n)[None,:,None] for k in range(4)]
    z=(np.cos(2*np.pi*(np.arange(n)+.5)/n)*np.sinc(1/n))[None,None,:]
    f=64*(y[0]+shear*x[0]*y[1])*z
    laplace=64*(shear*x[2]*y[1]-2*shear**2*x[1]*y[2]+(1+shear**2)*(y[2]+shear*x[0]*y[3]))*z-4*np.pi**2*f
    return f.ravel(),laplace.ravel()


def boundary_reference(shear):
    values=np.linspace(0,1,37);errors=[];tangent=[]
    for x in (0.,1.):
        gx=shear*C.deriv()(x)*P.deriv()(values)
        gy=P.deriv()(values)+shear*C(x)*P.deriv(2)(values)
        errors.extend(abs(gx-shear*gy));tangent.extend(abs(gy))
    for y in (0.,1.):
        gx=shear*C.deriv()(values)*P.deriv()(y)
        gy=P.deriv()(y)+shear*C(values)*P.deriv(2)(y)
        errors.extend(abs(-shear*gx+(1+shear**2)*gy))
    return {'max_conormal_residual':float(max(errors)),
            'max_logical_tangential_derivative_on_xi_wall':float(max(tangent))}


def evolve_manufactured(transport,n,shear,time):
    f,laplace=manufactured(n,shear);initial=1+.2*f
    source=.2*(-f-laplace)
    augmented=sparse.bmat([[transport.delta,sparse.csr_matrix(source[:,None])],
                          [sparse.csr_matrix((1,n**3)),sparse.csr_matrix([[-1.]])]],format='csr')
    final=expm_multiply(time*augmented,np.r_[initial,1.])[:-1]
    exact=1+.2*np.exp(-time)*f
    return initial,final,exact,source


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    protocol={'grids':[8,16,32],'shears':[0.,.25,.5,1.],'time':.02,
              'pulse_grid':16,'pulse_time_over_h_squared':.01,
              'criteria':{'finest_relative_l2_max':.001,'last_order_min':1.8,
                          'mass_drift_max':1e-12,'conormal_reference_max':1e-12,
                          'negative_concentration_tolerance':1e-12,'invariant_residual_max':1e-9},
              'source_sha256':{str(Path(__file__).with_name(n).resolve()):hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ('skew_boundary.py','positive_skew_flux.py','transport.py')},
              'scope':__doc__,
              'boundary_rule':'Omit every directional link with an endpoint outside the logical cube. Retained exchanges use the same nonnegative weights as the periodic prototype; no periodic wrap or flux clipping.',
              'reference':'Exact polynomial/trigonometric cell averages with homogeneous conormal Neumann walls and zero-mean manufactured source. Augmented matrix exponential isolates spatial error from time stepping. Positivity tested separately without source.',
              'limits':'Static affine parallelepiped, constant coefficients, uniform volumes. Boundary-cell norms include edges/corners. Does not validate irregular faces, variable capacities, live phase-field geometry or evolving topology.'}
    output.mkdir(parents=True);write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());state=json.loads((output/'status.json').read_text())
    if state['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=state['protocol_sha256']:raise ValueError('requires unchanged prepared protocol')
    if any(hashlib.sha256(Path(k).read_bytes()).hexdigest()!=v for k,v in p['source_sha256'].items()):raise ValueError('prototype sources changed')
    write_json(output/'status.json',{**state,'state':'running'})
    try:
        cases=[]
        for shear in p['shears']:
            rows=[]
            for n in p['grids']:
                transport=closed_transport(n,shear)
                initial,final,exact,source=evolve_manufactured(transport,n,shear,p['time'])
                ijk=np.indices((n,n,n));wall=np.any((ijk==0)|(ijk==n-1),axis=0).ravel()
                error=final-exact
                rows.append({'grid':n,'relative_l2':float(np.linalg.norm(error)/np.linalg.norm(exact)),
                             'boundary_relative_l2':float(np.linalg.norm(error[wall])/np.linalg.norm(exact[wall])),
                             'interior_relative_l2':float(np.linalg.norm(error[~wall])/np.linalg.norm(exact[~wall])),
                             'mass_drift':float(abs(transport.volumes@(final-initial))),
                             'source_integral':float(transport.volumes@source),
                             'minimum_concentration':float(final.min()),
                             'constant_residual':float(np.max(abs(transport.delta@np.ones(n**3)))),
                             'mass_residual':float(np.max(abs(transport.volumes@transport.delta)))})
            errors=np.array([r['relative_l2'] for r in rows]);orders=np.log2(errors[:-1]/errors[1:]).tolist()
            n=p['pulse_grid'];t=closed_transport(n,shear);pulse=np.zeros(n**3);pulse.reshape((n,n,n))[0,n//2,n//2]=1
            final=expm_multiply(t.delta*(p['pulse_time_over_h_squared']/n**2),pulse)
            off=t.delta.copy();off.setdiag(0);off.eliminate_zeros()
            case={'shear':shear,'reference_boundary':boundary_reference(shear),'refinement':rows,'orders':orders,
                  'pulse_minimum':float(final.min()),'pulse_mass_drift':float(abs(t.volumes@(final-pulse))),
                  'minimum_off_diagonal_rate':float(off.data.min()) if off.nnz else 0.,
                  'pulse_maximum':float(final.max())}
            cases.append(case)
            print(f'shear={shear}: error={errors[-1]:.4%}, order={orders[-1]:.4f}',flush=True)
        c=p['criteria'];checks={
            'manufactured_conormal_boundary_exact':all(x['reference_boundary']['max_conormal_residual']<c['conormal_reference_max'] for x in cases),
            'nonzero_tangential_gradient_exercised':all(x['reference_boundary']['max_logical_tangential_derivative_on_xi_wall']>.001 for x in cases),
            'mass_and_constants_preserved':all(max(r['constant_residual'],r['mass_residual'])<c['invariant_residual_max'] for x in cases for r in x['refinement']),
            'zero_source_integral_and_mass_drift':all(max(abs(r['source_integral']),r['mass_drift'])<c['mass_drift_max'] for x in cases for r in x['refinement']),
            'rates_nonnegative':all(x['minimum_off_diagonal_rate']>=0 for x in cases),
            'wall_pulses_positive_and_conservative':all(x['pulse_minimum']>=-c['negative_concentration_tolerance'] and x['pulse_maximum']<=1+c['negative_concentration_tolerance'] and x['pulse_mass_drift']<c['mass_drift_max'] for x in cases),
            'all_spatial_errors_decrease':all(all(a['relative_l2']>b['relative_l2'] for a,b in zip(x['refinement'][:-1],x['refinement'][1:])) for x in cases),
            'all_last_orders_above_1_8':all(x['orders'][-1]>c['last_order_min'] for x in cases),
            'all_finest_errors_below_0_1_percent':all(x['refinement'][-1]['relative_l2']<c['finest_relative_l2_max'] for x in cases)}
        result={'protocol':p,'cases':cases,'checks':checks,'passed':all(checks.values()),'live_integration_ready':False}
        write_json(output/'comparison.json',result)
        lines=['# Reflecting-wall skew-mesh validation','',f'All declared checks pass: **{result["passed"]}**','',
               '| Shear | Finest overall error | Finest boundary error | Last overall order | Wall pulse minimum |',
               '|---:|---:|---:|---:|---:|']
        lines += [f'| {x["shear"]} | {x["refinement"][-1]["relative_l2"]:.6%} | {x["refinement"][-1]["boundary_relative_l2"]:.6%} | {x["orders"][-1]:.4f} | {x["pulse_minimum"]:.3g} |' for x in cases]
        lines+=['',*[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in checks.items()],'',p['reference'],'',p['limits']]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json',{**state,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**state,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
