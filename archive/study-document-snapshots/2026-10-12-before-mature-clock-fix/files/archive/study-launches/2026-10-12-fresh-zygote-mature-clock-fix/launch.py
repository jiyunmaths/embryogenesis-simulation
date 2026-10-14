"""Launch the authorized short mature-clock qualification after independent audit."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=Path.cwd()
STUDY=ROOT/'output/studies/fresh-zygote-mature-clock-fix'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


if __name__=='__main__':
    assert not (HERE/'process-launch.json').exists(),'Already launched'
    p=json.loads((STUDY/'protocol.json').read_text());audit=json.loads((HERE/'prelaunch-audit.json').read_text())
    assert audit['passed'] and audit['protocol_sha256']==sha(STUDY/'protocol.json')
    assert audit['test_log_sha256']==sha(HERE/'tests.log')
    assert audit['test_verification_sha256']==sha(HERE/'test-verification.json')
    assert json.loads((STUDY/'status.json').read_text())==dict(state='prepared',protocol_sha256=audit['protocol_sha256'])
    assert shutil.disk_usage(STUDY).free>20*2**30,'Insufficient data-drive headroom'
    availability=subprocess.check_output(['nvidia-smi','--id='+p['gpu_uuid'],
        '--query-gpu=uuid,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
    uuid,used,util=[x.strip() for x in availability.split(',')]
    assert uuid==p['gpu_uuid'] and int(used)<100 and int(util)==0,'Dedicated GPU is busy'
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=p['gpu_uuid'],OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONUNBUFFERED='1')
    migration=ROOT/'archive/output-relocation/2026-10-10-output/studies-manifest.json'
    argv=[sys.executable,str(ROOT/'scripts/with_study_output_path_alias.py'),
        '--record',str(migration),'--record-sha256',audit['migration_record_sha256'],
        '-m','embryo.fresh_zygote_mature_clock_fix','run',str(STUDY)]
    with (HERE/'run.log').open('x') as log:
        child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    record=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,argv=argv,cwd=str(ROOT),
        gpu_uuid=p['gpu_uuid'],histories=p['histories'],branches=p['branches'],new_scientific_histories=0,
        new_numerical_paths=24,reused_prefixes=1,original_failure_retained=True,physical_window=[18,24],long_scientific_continuation_started=False,
        protocol_sha256=audit['protocol_sha256'],prelaunch_audit_sha256=sha(HERE/'prelaunch-audit.json'),
        launcher_sha256=sha(__file__),availability_before_launch=availability,
        external_free_bytes=shutil.disk_usage(STUDY).free)
    (HERE/'process-launch.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
