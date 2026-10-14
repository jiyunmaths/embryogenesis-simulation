"""Queue the full study behind the existing qualifier, never bypass its gates."""
import fcntl,json,subprocess,sys,hashlib
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path('output/studies/moving-recipient-amount-response')
HERE=Path(__file__).parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    # This helper runs asynchronously; the UI is not waiting on this lock.
    with (ROOT/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        status=json.loads((ROOT/'status.json').read_text())
        if status['state']!='qualified_wait_starts' or status['qualified_jobs']!=24:
            raise RuntimeError('Full study not started: actual-start qualification incomplete/failed')
        with (HERE/'qualification-audit.log').open('w') as log:
            subprocess.run([sys.executable,str(HERE/'audit_qualification.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
        audit=json.loads((HERE/'qualification-audit.json').read_text());assert audit['passed'] and audit['actual_start_gates']==24
        report=dict(passed=True,state='launching_full_study',utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=sha(ROOT/'protocol.json'),method=audit['method'],actual_start_gates=24,qualification_audit_sha256=sha(HERE/'qualification-audit.json'),physical_gpu=1,native_prefix_workers=4,native_threads_per_worker=2,wait_paths=24,response_paths=72,total_new_paths=96,independent_histories=3,new_histories=0,scientific_assessment_pending=True,helper_sha256=sha(__file__))
        (HERE/'launch-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    subprocess.run([sys.executable,'-m','embryo.moving_recipient_amount_response','run','--workers','4'],check=True)

if __name__=='__main__':main()
