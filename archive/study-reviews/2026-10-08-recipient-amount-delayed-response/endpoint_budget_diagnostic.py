"""Post-run descriptive comparison, without new acceptance thresholds or solves."""
from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
root=Path('output/studies/phase-carry-recipient-amount-response')
pfile=root/'protocol.json';p=json.loads(pfile.read_text());oldroot=Path(p['parent']);ofile=oldroot/'protocol.json';old=json.loads(ofile.read_text())
files=[pfile,ofile,root/'production/fine/paths.npz',oldroot/'production/fine/paths.npz']
def load(f):
 with np.load(f) as z:return z['times'].copy(),{str(k):u for k,u in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
time,new=load(files[2]);other_time,prior=load(files[3]);np.testing.assert_array_equal(time,other_time)
rows=[]
for c in p['contexts']:
 first=next(t for t in p['trials'] if t['context']==c['key'] and t['target'] is None)
 second=next(t for t in old['trials'] if t['context']==c['parent_context'] and t['target'] is None)
 a,b=new[first['key']],prior[second['key']]
 rows.append(dict(context=c['key'],seed=c['seed'],cell=c['targets'][0],matching_reference=c['matching_reference'],physics=c['physics_label'],chemistry=c['chemical_origin'],maximum_initial_log_discrepancy=float(np.max(abs(np.log(a[0]/b[0])))),maximum_late_log_discrepancy=float(np.max(abs(np.log(a[time>=36]/b[time>=36])))),maximum_endpoint_log_discrepancy=float(np.max(abs(np.log(a[-1]/b[-1]))))))
out=dict(post_run_descriptive=True,no_new_acceptance_threshold=True,no_new_integrations=True,contexts=len(rows),maximum_endpoint_log_discrepancy=max(r['maximum_endpoint_log_discrepancy'] for r in rows),maximum_late_log_discrepancy=max(r['maximum_late_log_discrepancy'] for r in rows),rows=rows,files_sha256={str(f):sha(f) for f in [*files,Path(__file__)]},scope='Initial budget-matched vs accepted unmatched-budget unpulsed trajectories under the same local-reference/chemical-origin/physical-context labels. Endpoints near each other do not prove identical basins or global stability.')
(HERE/'endpoint-budget-diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print('Post-run endpoint discrepancy',out['maximum_endpoint_log_discrepancy'])
