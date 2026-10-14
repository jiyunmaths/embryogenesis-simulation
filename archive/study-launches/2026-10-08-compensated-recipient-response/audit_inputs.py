"""Verify all new starts and their exact all-start GPU prefix coverage."""
import importlib.util,json,hashlib
from pathlib import Path
import numpy as np

SOURCE=Path('archive/study-launches/2026-10-08-recipient-amount-delayed-response/audit_inputs.py')
spec=importlib.util.spec_from_file_location('accepted_input_audit',SOURCE)
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ROOT=Path('output/studies/compensated-recipient-amount-delayed-response')
METHOD='float64 SSP-RK2 increment with Kahan accumulation v1'

def main():
    result=r.audit(ROOT);p=r.read(ROOT/'protocol.json');old=r.read(Path(p['historical_failed_study'])/'protocol.json');q=r.read(p['backend_qualification'])
    assert p['solver']['method']==q['method']==METHOD and q['passed']
    for file,h in q['files_sha256'].items():r.pin(file,h)
    prefix=Path('output/studies/compensated-chemistry-validation/prefix');pref=r.read(prefix/'protocol.json');pre=r.read(prefix/'summary.json')
    assert pref['trials']==old['trials'] and pref['contexts']==old['contexts'] and pre['passed']
    assert pref['levels']==p['levels'] and all(x['logical_cases']==144 and x['gpu_restart_bitwise_identical'] and x['passed'] for x in pre['levels'])
    assert pref['solver']['method']==p['solver']['method']
    for a,b in zip(p['contexts'],old['contexts']):
        assert a['key']==b['key'];x,y=r.arrays(a['source']),r.arrays(b['source'])
        for k in ('initial','delta','volumes','ids'):np.testing.assert_array_equal(x[k],y[k])
    for a,b in zip(p['trials'],old['trials']):
        assert {k:v for k,v in a.items() if k!='source'}=={k:v for k,v in b.items() if k!='source'}
        np.testing.assert_array_equal(r.arrays(a['source'])['initial'],r.arrays(b['source'])['initial'])
    for file,h in r.WATCHED.items():assert r.digest(file)==h
    result.update(method=METHOD,exact_original_failed_preparations=True,actual_start_prefix_coverage=144,prefix_reuse_justification='All144 inputs, contexts, clocks, levels and method are exact copies of the independently qualified prefix; numerical evidence only, not new scientific trajectories.',prefix_summary=str(prefix/'summary.json'),backend_qualification_verified=True,inputs_sha256=r.WATCHED,script_sha256=r.digest(__file__))
    Path(__file__).with_name('preparation-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','logical_cases','unique_initial_value_problems','unique_physical_initial_value_problems','actual_start_prefix_coverage')},indent=2))

if __name__=='__main__':main()
