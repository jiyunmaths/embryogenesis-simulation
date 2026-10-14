"""Launch one serial GPU timestep ladder after the pinned instrumentation audit."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess
import sys

if __name__=='__main__':
    folder=Path(__file__).resolve().parent;root=Path.cwd();study=root/'output/studies/development-time-diagnosis'
    if (folder/'process-launch.json').exists():raise RuntimeError('Already launched')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    audit=json.loads((folder/'prelaunch-audit.json').read_text())
    assert audit['passed'] and audit['protocol_sha256']==sha(study/'protocol.json')
    assert audit['instrumentation_gate_sha256']==sha(study/'instrumentation-gate.json')
    assert audit['test_log_sha256']==sha(folder/'tests.log')
    assert json.loads((study/'status.json').read_text())['state']=='prepared'
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',
        OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    argv=[sys.executable,'-m','embryo.development_time_diagnosis','run','--output',str(study)]
    with (folder/'run.log').open('w') as log:
        child=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    result=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,argv=argv,cwd=str(root),gpu_uuid=env['CUDA_VISIBLE_DEVICES'],
        cpu_full_ladder_workers=0,total_gpu_paths=4,protocol_sha256=audit['protocol_sha256'],
        prelaunch_audit_sha256=sha(folder/'prelaunch-audit.json'),instrumentation_gate_sha256=audit['instrumentation_gate_sha256'],
        launcher_sha256=sha(Path(__file__)))
    (folder/'process-launch.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
