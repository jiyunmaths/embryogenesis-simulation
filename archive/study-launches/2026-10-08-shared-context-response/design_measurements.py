import json,hashlib
from pathlib import Path
import numpy as np
root=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/output/studies/phase-carry-frozen-response')
p=json.loads((root/'protocol.json').read_text());records=[];pins={str(root/'protocol.json'):hashlib.sha256((root/'protocol.json').read_bytes()).hexdigest()}
contexts={c['key']:c for c in p['contexts']}
def get(seed,bg):
    c=next(c for c in p['contexts'] if c['seed']==seed and c['background']==bg);f=Path(c['source']);pins[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
    with np.load(f) as z:return {k:z[k].copy() for k in z.files}
for seed in p['histories']:
    s=get(seed,'sham')
    for cell in p['selections'][str(seed)]['selected_ids']:
        r=get(seed,f'redistributed-cell-{cell}');assert np.array_equal(s['ids'],r['ids']);i=s['ids'].tolist().index(cell)
        amounts={c+v:(chem['initial']@phys['volumes']) for c,chem in [('S',s),('R',r)] for v,phys in [('S',s),('R',r)]}
        amin=min(s['initial'][0,i],r['initial'][0,i]);factors={c:float(1.-.1*amin/z['initial'][0,i]) for c,z in [('S',s),('R',r)]}
        doses={v:float(.1*min(s['initial'][0,i]*phys['volumes'][i],r['initial'][0,i]*phys['volumes'][i])) for v,phys in [('S',s),('R',r)]}
        errors=[]
        for c,chem in [('S',s),('R',r)]:
            for v,phys in [('S',s),('R',r)]:errors.append(abs((chem['initial'][0,i]-chem['initial'][0,i]*factors[c])*phys['volumes'][i]/doses[v]-1))
        records.append(dict(seed=seed,cell=cell,recipient_concentrations={c:z['initial'][:,i].tolist() for c,z in [('S',s),('R',r)]},recipient_volumes={v:z['volumes'][i] for v,z in [('S',s),('R',r)]},amounts={key:value.tolist() for key,value in amounts.items()},chemical_amount_ratio_red_sham_at_sham_physics=(amounts['RS']/amounts['SS']).tolist(),chemical_amount_ratio_red_sham_at_red_physics=(amounts['RR']/amounts['SR']).tolist(),physical_amount_ratio_red_sham_at_sham_chemistry=(amounts['SR']/amounts['SS']).tolist(),physical_amount_ratio_red_sham_at_red_chemistry=(amounts['RR']/amounts['RS']).tolist(),common_factors=factors,common_amount_by_physics=doses,largest_relative_common_dose_roundoff=max(errors)))
report=dict(scope='Static source calculations only; no ODE integration.',records=records,inputs_sha256=pins)
Path('/tmp/shared_context_design_measurements.json').write_text(json.dumps(report,indent=2)+'\n')
for r in records:print(r['seed'],r['cell'],'chemical amount R/S @S,@R',r['chemical_amount_ratio_red_sham_at_sham_physics'],r['chemical_amount_ratio_red_sham_at_red_physics'],'physical amounts R/S at C=S,R',r['physical_amount_ratio_red_sham_at_sham_chemistry'],r['physical_amount_ratio_red_sham_at_red_chemistry'],'factors',r['common_factors'],'doseerr',r['largest_relative_common_dose_roundoff'])
