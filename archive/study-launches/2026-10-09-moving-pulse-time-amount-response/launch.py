"""Detach only the pinned moving matched-start test; every gate runs first."""
from datetime import datetime,timezone
import hashlib,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
FOLDER=Path(__file__).resolve().parent
STUDY=ROOT/'output/studies/moving-pulse-time-amount-response'
UUID='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e'
sys.path.insert(0,str(ROOT))
from embryo import phase_carry_frozen_response as frozen
from embryo import phase_carry_exchange_response as moving
from embryo.parameter_robustness import verify

def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def read(f):return json.loads(Path(f).read_text())
def main():
    os.chdir(ROOT);p=read(STUDY/'protocol.json');verify(p)
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==UUID
    frozen.runtime_check(p);moving.require_device(p['cuda_validation'])
    pre=read(FOLDER/'preparation-audit.json');assert pre['passed'] and pre['protocol_sha256']==sha(STUDY/'protocol.json')
    for file,h in pre['inputs_sha256'].items():assert sha(file)==h,file
    assert read(STUDY/'status.json')['state']=='prepared' and len(p['jobs'])==48
    assert not (FOLDER/'process-launch.json').exists()
    env=os.environ.copy()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='1'
    command=[sys.executable,'-m','embryo.moving_pulse_time_amount_response','run','--workers','4']
    with (FOLDER/'run.log').open('xb') as log:
        child=subprocess.Popen(command,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    record=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,command=command,gpu_uuid=UUID,physical_gpu=1,qualified_gpu=p['runtime']['expected_gpu_name'],native_prefix_workers=4,native_threads_per_worker=2,gpu_workers=1,blas_threads=1,protocol_sha256=sha(STUDY/'protocol.json'),script_sha256=sha(__file__),preparation_audit_sha256=sha(FOLDER/'preparation-audit.json'),qualification_audit_required_before_production=p['qualification_audit'],all_context_pilot_required_before_full=True,total_paths=48,logical_starts=24,distinct_initial_value_problems=p['distinct_initial_value_problems'],independent_histories=3,new_histories=0,start=570.,end=630.,estimated_hours=p['estimated_moving_seconds']/3600,scientific_results_pending=True)
    (FOLDER/'process-launch.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
