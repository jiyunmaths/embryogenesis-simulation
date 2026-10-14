import json, math
from pathlib import Path
import numpy as np
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
parent=ROOT/'outputs/phase-carry-redistribution-response'
p=json.loads((parent/'protocol.json').read_text())
def folder(j):
    return Path(p['reuse_parent'])/p['reused'][j['key']]['job']['key'] if j['key'] in p['reused'] else parent/j['key']
def first(j): return json.loads((folder(j)/'history.json').read_text())[0]
controls=[j for j in p['jobs'] if j['level']=='fine' and j['dose'] is None]
records=[]; all_initial_equal=True; source_equal=True; dose_error=0; source_ids_equal=True
for j in controls:
    row=first(j); m=np.array(row['volumes']); d=np.array(row['delta']); x=np.array(row['chemistry']); ids=np.array(row['ids'])
    g=m[:,None]*d; off=g.copy(); np.fill_diagonal(off,0)
    scale=np.max(np.sum(abs(g),axis=1)); symmetric=-np.sqrt(m)[:,None]*d/np.sqrt(m)[None,:]
    a,b=x; n=len(a)
    jac=np.block([[np.diag(2*a/b-1)+.02*d,np.diag(-a*a/b**2)],[np.diag(4*a),-2*np.eye(n)+.55*d]])
    eigen=np.linalg.eigvals(jac); q=float(max(-np.diag(d)))
    for other in p['jobs']:
        if other['seed']!=j['seed'] or other['background']!=j['background']: continue
        r=first(other); xx=np.array(r['chemistry']); expected=x.copy()
        if other['dose']:
            index=list(ids).index(other['target']); expected[0,index]*=other['factor']
            actual=(x[0,index]-xx[0,index])*m[index]
            dose_error=max(dose_error,abs(actual/other['requested_removed_amount']-1))
        all_initial_equal &= np.array_equal(xx,expected) and np.array_equal(r['delta'],d) and np.array_equal(r['volumes'],m) and np.array_equal(r['ids'],ids)
        with np.load(other['source']) as z:
            source_equal &= np.array_equal(z['activator'],xx[0]) and np.array_equal(z['inhibitor'],xx[1])
            source_ids_equal &= np.array_equal(z['ids'],ids)
    records.append(dict(seed=j['seed'],background=j['background'],time=row['time'],cells=n,ids=ids.tolist(),min_concentration=float(x.min()),max_concentration=float(x.max()),volume_min=float(m.min()),volume_max=float(m.max()),volume_sum=float(m.sum()),amounts=(x@m).tolist(),max_exit=q,laplacian_lambda_max=float(np.linalg.eigvalsh(symmetric).max()),initial_jacobian_real_min=float(eigen.real.min()),initial_jacobian_real_max=float(eigen.real.max()),initial_jacobian_imag_abs_max=float(abs(eigen.imag).max()),initial_rhs_abs_max=float(np.max(abs(np.array([a*a/b-a+.02*d@a,2*(a*a-b)+.55*d@b])))),row_sum_relative=float(abs(d.sum(axis=1)).max()/q),mass_conservation_relative=float(abs(m@d).max()/scale),reciprocity_relative=float(abs(g-g.T).max()/scale),offdiagonal_min=float(off.min()),substeps={str(dt):max(1,math.ceil(dt*max(1+.02*q,2+.55*q)/.2)) for dt in [.00375,.001875,.0009375]}))
pulses=[j for j in p['jobs'] if j['dose']]
minimum=min(pulses,key=lambda j:abs(math.log(j['factor'])))
out=dict(scope='Read-only audit of archived starts and matrices; no trajectories integrated.',contexts=records,logical_paths=len(p['jobs']),unique_fine_paths=len([j for j in p['jobs'] if j['level']=='fine']),all_initial_conditions_exact=bool(all_initial_equal),all_source_chemistry_exact=bool(source_equal),all_source_ids_exact=bool(source_ids_equal),max_requested_pulse_amount_relative_error=dose_error,min_pulse_log_amplitude=abs(math.log(minimum['factor'])),min_pulse_job=minimum['key'],reference_raw_log_limit=1e-8,conservative_response_reference_error_bound=2e-8/abs(math.log(minimum['factor'])),conservative_effect_reference_error_bound=4e-8/abs(math.log(minimum['factor'])))
Path('/tmp/frozen_response_design_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='contexts'},indent=2))
for r in records: print(r['seed'],r['background'],'q',r['max_exit'],'lambda_max',r['laplacian_lambda_max'],'initial_jac_real',r['initial_jacobian_real_min'],r['initial_jacobian_real_max'],'minchem',r['min_concentration'],'rhsmax',r['initial_rhs_abs_max'],'substeps',r['substeps'])
