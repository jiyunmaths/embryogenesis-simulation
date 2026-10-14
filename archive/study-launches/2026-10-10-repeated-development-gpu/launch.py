"""Launch the fixed repeated-development study with bounded reference workers."""
from pathlib import Path
import os,json,hashlib,subprocess,sys
from datetime import datetime,timezone

if __name__=='__main__':
    folder=Path(__file__).resolve().parent;root=Path.cwd();study=root/'output/studies/repeated-development-gpu-qualification-v2'
    if (folder/'process-launch.json').exists():raise RuntimeError('Already launched')
    audit=json.loads((folder/'prelaunch-audit.json').read_text())
    assert audit['passed'] and audit['protocol_sha256']==hashlib.sha256((study/'protocol.json').read_bytes()).hexdigest()
    assert json.loads((study/'status.json').read_text())['state']=='prepared'
    assert '35 passed' in (folder/'tests.log').read_text()
    assert json.loads((study/'calibration.json').read_text())['passed']
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',
        OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    argv=[sys.executable,'-m','embryo.qualify_repeated_development','run','--output',str(study),'--workers','4']
    with (folder/'run.log').open('w') as log:
        child=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    record=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,argv=argv,cwd=str(root),cpu_workers=4,
        gpu_uuid=env['CUDA_VISIBLE_DEVICES'],protocol_sha256=audit['protocol_sha256'],
        calibration_sha256=hashlib.sha256((study/'calibration.json').read_bytes()).hexdigest(),
        prelaunch_audit_sha256=hashlib.sha256((folder/'prelaunch-audit.json').read_bytes()).hexdigest(),
        launcher_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (folder/'process-launch.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
