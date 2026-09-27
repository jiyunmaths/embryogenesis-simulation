"""Independent geometric tests of the live diffuse-contact conductance closure.

Known reference areas and linear-field fluxes distinguish quadrature error from
closure error. No embryo parameters or running simulations are modified.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from .resolution import write_json
from .transport import contact_transport


def shells(s, epsilon, gap=0.):
    a=expit(-np.sqrt(2)*(s+gap/2)/epsilon)
    b=expit(np.sqrt(2)*(s-gap/2)/epsilon)
    return a*(1-a)*b*(1-b)


def flat_overlap(epsilon,gap=0.):
    # Symmetry reduction of a planar 3D interface with unit transverse area.
    bound=abs(gap)/2+30*epsilon
    return quad(lambda s:shells(s,epsilon,gap),-bound,bound,epsabs=1e-13,epsrel=1e-12)[0]


def adapter(overlap,epsilon,centers=None):
    centers=np.array([[-.5,0,0],[.5,0,0]]) if centers is None else np.asarray(centers)
    return contact_transport(np.array([[0.,overlap],[overlap,0.]]),np.ones(2),centers,epsilon,.02)


def spherical_area(radius,epsilon):
    # The outer cell is complementary to the inner one. This tests area only:
    # their actual centroids coincide, so the two-point conductance is undefined.
    upper=radius+30*epsilon
    overlap=quad(lambda r:4*np.pi*r*r*shells(r-radius,epsilon),0,upper,epsabs=1e-13,epsrel=1e-11)[0]
    return 6*np.sqrt(2)*overlap/epsilon


def spherical_voxel_area(grid,radius,epsilon,extent=1.):
    # Integrate the same 3D shell product on voxel centers, in z slices.
    dx=2*extent/grid
    axis=(np.arange(grid)+.5)*dx-extent
    xy=axis[:,None]**2+axis[None,:]**2
    total=0.
    for z in axis:
        signed=np.sqrt(xy+z*z)-radius
        total+=float(np.sum(shells(signed,epsilon),dtype=np.float64))
    return 6*np.sqrt(2)*total*dx**3/epsilon


def flux_patch(shear,gradient,epsilon=.05):
    # Exact centroids of sharp sheared prisms: x in [-1,0] / [0,1],
    # y=shear*x+u, u,z in [-1/2,1/2]. Shared face x=0 has area 1.
    d=np.array([1.,shear,0.]);centers=np.stack([-d/2,d/2])
    transport=adapter(flat_overlap(epsilon),epsilon,centers)
    gradient=np.asarray(gradient,float)
    estimated=float(transport.conductance[0,1]*(gradient@d))
    exact=float(gradient[0]) # D=1, signed transfer into left cell
    projected_distance_estimate=float(gradient@d) # A/d_normal shortcut
    return {'shear':shear,'gradient':gradient.tolist(),'exact_flux':exact,'estimated_flux':estimated,
            'absolute_error':abs(estimated-exact),'projected_distance_flux':projected_distance_estimate,
            'mass_residual':float(np.max(abs(transport.volumes@transport.delta.toarray()))),
            'constant_residual':float(np.max(abs(transport.delta@np.ones(2))))}


def prepare(output):
    output=Path(output)
    if output.exists():raise FileExistsError('choose a fresh output directory')
    protocol={'flat_epsilons':[.02,.05,.1], 'gap_over_epsilon':[0.,.25,.5,1.,2.,4.,8.],
              'fixed_gap':.1,'fixed_gap_epsilons':[.1,.05,.025,.0125],
              'radius':.5,'curvature_epsilon_over_radius':[.05,.1,.2,.3],
              'voxel_grids':[64,96,128],'voxel_extent':1.,'shears':[0.,.25,.5,1.],
              'gradients':[[1,0,0],[0,1,0],[1,1,0]],
              'criteria':{'flat_area_relative_max':1e-8,'finest_voxel_quadrature_relative_max':.001,
                          'curved_area_relative_max':.01,'positive_gap_area_relative_max':.01,
                          'linear_flux_absolute_max':.01,'conservation_residual_max':1e-12},
              'transport_sha256':hashlib.sha256(Path(__file__).with_name('transport.py').read_bytes()).hexdigest(),
              'scope':__doc__,
              'interpretation_rules':'Sharp direct-contact reference: positive physical gap has zero common face; diffuse exchange across a gap is leakage under this interpretation, not proof of invalid extracellular signaling. Curvature tests area, not nested-cell conductance; concentric centers are inadmissible. Nonorthogonal patch uses exact sharp centroids to isolate A/ell rather than add diffuse-centroid bias. Independent reference geometry is manufactured, not biological calibration.'}
    output.mkdir(parents=True)
    write_json(output/'protocol.json',protocol)
    write_json(output/'status.json',{'state':'prepared','protocol_sha256':hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()})
    return protocol


def run(output):
    output=Path(output);p=json.loads((output/'protocol.json').read_text());status=json.loads((output/'status.json').read_text())
    if status['state']!='prepared' or hashlib.sha256((output/'protocol.json').read_bytes()).hexdigest()!=status['protocol_sha256']:
        raise ValueError('requires unchanged prepared protocol')
    if hashlib.sha256(Path(__file__).with_name('transport.py').read_bytes()).hexdigest()!=p['transport_sha256']:
        raise ValueError('transport changed after preparation')
    write_json(output/'status.json',{**status,'state':'running'})
    try:
        flat=[]
        for epsilon in p['flat_epsilons']:
            estimate=float(adapter(flat_overlap(epsilon),epsilon).conductance[0,1])
            flat.append({'epsilon':epsilon,'estimated_area':estimate,'relative_error':abs(estimate-1)})
        gaps=[]
        for ratio in p['gap_over_epsilon']:
            epsilon=.05
            t=adapter(flat_overlap(epsilon,ratio*epsilon),epsilon)
            gaps.append({'gap_over_epsilon':ratio,'estimated_area_over_unit_face':float(t.conductance[0,1]),
                         'edge_retained_at_default_cutoff':bool(t.conductance.nnz)})
        fixed=[]
        for epsilon in p['fixed_gap_epsilons']:
            t=adapter(flat_overlap(epsilon,p['fixed_gap']),epsilon)
            fixed.append({'epsilon':epsilon,'gap':p['fixed_gap'],'estimated_area_over_unit_face':float(t.conductance[0,1])})
        curved=[]
        exact=4*np.pi*p['radius']**2
        for ratio in p['curvature_epsilon_over_radius']:
            epsilon=ratio*p['radius'];reference=spherical_area(p['radius'],epsilon)
            voxel=[{'grid':grid,'area':spherical_voxel_area(grid,p['radius'],epsilon,p['voxel_extent'])} for grid in p['voxel_grids']]
            for row in voxel:row['relative_error_to_diffuse_integral']=abs(row['area']/reference-1)
            curved.append({'epsilon_over_radius':ratio,'exact_sharp_area':exact,'diffuse_integral_area':reference,
                           'relative_closure_bias':abs(reference/exact-1),'voxel_quadrature':voxel})
        patches=[flux_patch(shear,gradient) for shear in p['shears'] for gradient in p['gradients']]
        limits=p['criteria']
        numerical={'flat_calibration_exact':all(r['relative_error']<limits['flat_area_relative_max'] for r in flat),
                   'finest_3d_quadrature_matches_diffuse_integral':all(r['voxel_quadrature'][-1]['relative_error_to_diffuse_integral']<limits['finest_voxel_quadrature_relative_max'] for r in curved),
                   'mass_and_constants_preserved':all(max(r['mass_residual'],r['constant_residual'])<limits['conservation_residual_max'] for r in patches),
                   'orthogonal_linear_flux_exact':all(r['absolute_error']<limits['linear_flux_absolute_max'] for r in patches if r['shear']==0)}
        physical={'all_curved_areas_within_1_percent':all(r['relative_closure_bias']<limits['curved_area_relative_max'] for r in curved),
                  'all_positive_gaps_below_1_percent_of_face_area':all(r['estimated_area_over_unit_face']<limits['positive_gap_area_relative_max'] for r in gaps if r['gap_over_epsilon']>0),
                  'all_nonorthogonal_linear_fluxes_within_0_01':all(r['absolute_error']<limits['linear_flux_absolute_max'] for r in patches if r['shear']>0)}
        trends={'curvature_bias_decreases_with_interface_width':all(a['relative_closure_bias']<b['relative_closure_bias'] for a,b in zip(curved[:-1],curved[1:])),
                'fixed_gap_leakage_decreases_with_interface_width':all(a['estimated_area_over_unit_face']>b['estimated_area_over_unit_face'] for a,b in zip(fixed[:-1],fixed[1:]))}
        result={'protocol':p,'flat':flat,'gaps':gaps,'fixed_gap':fixed,'curvature':curved,'linear_flux_patches':patches,
                'numerical_checks':numerical,'physical_closure_checks':physical,'trend_checks':trends,
                'passed':all(numerical.values()) and all(physical.values()) and all(trends.values())}
        write_json(output/'comparison.json',result)
        lines=['# Geometric transport approximation validation','',f'All checks pass: **{result["passed"]}**','',p['interpretation_rules'],'',
               '| Width / curvature radius | Area bias vs sharp sphere | Finest voxel error vs diffuse integral |','|---:|---:|---:|']
        lines += [f'| {r["epsilon_over_radius"]} | {r["relative_closure_bias"]:.4%} | {r["voxel_quadrature"][-1]["relative_error_to_diffuse_integral"]:.4%} |' for r in curved]
        lines += ['','| Gap / width | Inferred area / unit planar area |','|---:|---:|']+[f'| {r["gap_over_epsilon"]} | {r["estimated_area_over_unit_face"]:.6f} |' for r in gaps]
        lines += ['','| Shear | Gradient | Exact flux | Estimated flux |','|---:|---|---:|---:|']+[f'| {r["shear"]} | {r["gradient"]} | {r["exact_flux"]:.6g} | {r["estimated_flux"]:.6g} |' for r in patches]
        lines += ['']+[f'- {k}: {"PASS" if v else "FAIL"}' for k,v in {**numerical,**physical,**trends}.items()]
        (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
        write_json(output/'status.json',{**status,'state':'completed','all_checks_pass':result['passed']})
        return result
    except Exception as error:
        write_json(output/'status.json',{**status,'state':'failed','error':str(error)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.output) if args.action=='prepare' else run(args.output)


if __name__=='__main__':main()
