"""Launch four serial GPU paths after the pinned analytic/calibration audit."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess
import sys


if __name__ == '__main__':
    folder = Path(__file__).resolve().parent; root = Path.cwd()
    study = root/'output/studies/development-space-refinement'
    if (folder/'process-launch.json').exists(): raise RuntimeError('Already launched')
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    audit = json.loads((folder/'prelaunch-audit.json').read_text())
    protocol = json.loads((study/'protocol.json').read_text())
    assert audit['passed'] and audit['protocol_sha256'] == sha(study/'protocol.json')
    assert audit['calibration_sha256'] == sha(study/'calibration.json')
    assert audit['test_log_sha256'] == sha(folder/'tests.log')
    assert audit['inherited_instrumentation_gate_sha256'] == sha(protocol['instrumentation_gate'])
    assert json.loads((study/'status.json').read_text())['state'] == 'prepared'
    env = os.environ.copy(); env.update(CUDA_VISIBLE_DEVICES=protocol['gpu_uuid'],
        OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    argv = [sys.executable,'-m','embryo.development_space_refinement','run','--output',str(study)]
    with (folder/'run.log').open('w') as log:
        child = subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    result = dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,argv=argv,cwd=str(root),
        gpu_uuid=env['CUDA_VISIBLE_DEVICES'],new_gpu_paths=4,exact_reused_paths=1,cpu_full_ladder_workers=0,
        protocol_sha256=audit['protocol_sha256'],prelaunch_audit_sha256=sha(folder/'prelaunch-audit.json'),
        calibration_sha256=audit['calibration_sha256'],launcher_sha256=sha(__file__),new_scientific_histories=0)
    (folder/'process-launch.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
