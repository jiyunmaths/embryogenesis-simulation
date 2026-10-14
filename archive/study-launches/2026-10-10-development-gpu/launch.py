"""Launch only the fixed, prepared cleavage qualification."""
import hashlib,json,os,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

if __name__=='__main__':
    folder=Path(__file__).resolve().parent;root=Path.cwd();study=root/'output/studies/development-gpu-cleavage-qualification-v2'
    if (folder/'process-launch.json').exists():raise RuntimeError('Already launched')
    assert json.loads((study/'status.json').read_text())['state']=='prepared'
    assert '28 passed' in (folder/'tests.log').read_text()
    assert '5 passed' in (folder/'gpu-regression-tests.log').read_text()
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',
        OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    argv=[sys.executable,'-m','embryo.qualify_development_gpu','run','--output',str(study),'--workers','4']
    with (folder/'run.log').open('w') as log:process=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    record=dict(utc=datetime.now(timezone.utc).isoformat(),pid=process.pid,argv=argv,cwd=str(root),workers=4,gpu_uuid=env['CUDA_VISIBLE_DEVICES'],
        protocol_sha256=hashlib.sha256((study/'protocol.json').read_bytes()).hexdigest(),launcher_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (folder/'process-launch.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
