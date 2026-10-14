"""Launch only this independently qualified frozen GPU study; detach safely."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
FOLDER=Path(__file__).resolve().parent
STUDY=ROOT/'output/studies/pulse-time-amount-response'
UUID='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e'
sys.path.insert(0,str(ROOT))

from embryo import phase_carry_frozen_response as frozen
from embryo.parameter_robustness import verify


def sha(file):return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def read(file):return json.loads(Path(file).read_text())


def main():
    os.chdir(ROOT);p=read(STUDY/'protocol.json');verify(p)
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==UUID
    frozen.runtime_check(p)
    pre=read(FOLDER/'preparation-audit.json');pref=read(STUDY/'prefix/independent-audit.json')
    assert pre['passed'] and pre['protocol_sha256']==sha(STUDY/'protocol.json')
    assert pref['passed'] and pref['actual_context_timestep_cases']==1008
    assert pref['protocol_sha256']==sha(STUDY/'prefix/protocol.json')
    assert pref['summary_sha256']==sha(STUDY/'prefix/summary.json')
    for file,h in pref['evidence_sha256'].items():assert sha(file)==h,file
    assert read(STUDY/'status.json')['state']=='qualified'
    assert not (FOLDER/'process-launch.json').exists()
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    command=[sys.executable,'-m','embryo.pulse_time_amount_response','run','--workers','8']
    with (FOLDER/'run.log').open('xb') as log:
        child=subprocess.Popen(command,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    record=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,command=command,
        gpu_uuid=UUID,physical_gpu=1,qualified_gpu=p['runtime']['expected_gpu_name'],workers=8,blas_threads=1,
        protocol_sha256=sha(STUDY/'protocol.json'),script_sha256=sha(__file__),
        preparation_audit_sha256=sha(FOLDER/'preparation-audit.json'),prefix_audit_sha256=sha(STUDY/'prefix/independent-audit.json'),
        logical_starts=336,distinct_physical_problems=p['unique_physical_initial_value_problems'],
        production_paths=1008,cpu_paths=1008,reference_trajectories=672,independent_histories=3,new_histories=0,
        scientific_results_pending=True,premature_audit_note='Initial read-only prefix audit was attempted before qualification summary existed; failed without accepting results. Its log is preserved. The completed qualification and separate rerun both pass.')
    (FOLDER/'process-launch.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
