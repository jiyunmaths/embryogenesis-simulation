"""Post-run descriptive control endpoint comparison; no integration or gates."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
FOLDER=Path(__file__).resolve().parent
SOURCE=ROOT/'output/studies/phase-carry-recipient-state-response-grouped'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    files=[SOURCE/'protocol.json',SOURCE/'production/fine/paths.npz'];pins={str(p):sha(p) for p in files}
    p=json.loads(files[0].read_text())
    with np.load(files[1],allow_pickle=False) as z:
        paths={str(k):v.copy() for k,v in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))};times=z['times']
    rows=[]
    for seed in p['histories']:
        for cell in p['selections'][str(seed)]['selected_ids']:
            fs={L:next(f for f in p['factorials'] if (f['seed'],f['cell'],f['matching_reference'])==(seed,cell,L)) for L in ('S','R')}
            for C in ('S','R'):
                for P in ('S','R'):
                    controls={L:next(t['key'] for t in p['trials'] if t['context']==fs[L]['contexts'][C+'_'+P] and t['target'] is None) for L in ('S','R')}
                    difference=abs(np.log(paths[controls['S']]/paths[controls['R']]))
                    rows.append(dict(seed=seed,cell=cell,surrounding=C,physics=P,
                                     endpoint_max_raw_log_reference_difference=float(difference[-1].max()),
                                     late_max_raw_log_reference_difference=float(difference[times>=36].max())))
    assert all(sha(Path(f))==h for f,h in pins.items())
    r=dict(post_run_descriptive=True,no_new_acceptance_threshold=True,rows=rows,
           maximum_endpoint_raw_log_reference_difference=max(r['endpoint_max_raw_log_reference_difference'] for r in rows),
           maximum_late_raw_log_reference_difference=max(r['late_max_raw_log_reference_difference'] for r in rows),
           inputs_sha256=pins,inputs_unchanged=True,script_sha256=sha(Path(__file__)),
           scope='Alternative initial recipient reference pairs at fixed surrounding start and physical context. Descriptive endpoint/late comparison only; no equivalence, attractor uniqueness or new scientific replicate.')
    (FOLDER/'reference-endpoint-diagnostic.json').write_text(json.dumps(r,indent=2)+'\n')
    print({k:v for k,v in r.items() if k.startswith('maximum')})
if __name__=='__main__':main()
