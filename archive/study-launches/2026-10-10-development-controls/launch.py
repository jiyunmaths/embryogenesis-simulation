"""Launch the authorized serial control-refinement and live-restart study."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,subprocess,sys

if __name__=='__main__':
    here=Path(__file__).resolve().parent;root=Path.cwd()
    study=root/'output/studies/development-control-qualification'
    if (here/'process-launch.json').exists():raise RuntimeError('Already launched')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    p=json.loads((study/'protocol.json').read_text());audit=json.loads((here/'prelaunch-audit.json').read_text())
    assert audit['passed'] and audit['protocol_sha256']==sha(study/'protocol.json')
    assert audit['test_log_sha256']==sha(here/'tests.log')
    assert audit['calibration_sha256']==sha(p['calibration_reference'])
    assert audit['inherited_instrumentation_sha256']==sha(p['instrumentation_gate'])
    assert json.loads((study/'status.json').read_text())['state']=='prepared'
    availability=subprocess.check_output(['nvidia-smi','--id='+p['gpu_uuid'],
        '--query-gpu=uuid,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
    uuid,used,util=[x.strip() for x in availability.split(',')]
    assert uuid==p['gpu_uuid'] and int(used)<100 and int(util)==0,'Designated GPU became busy'
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=p['gpu_uuid'],OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    argv=[sys.executable,'-m','embryo.development_control_qualification','run','--output',str(study)]
    with (here/'run.log').open('x') as log:
        child=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    result=dict(utc=datetime.now(timezone.utc).isoformat(),pid=child.pid,argv=argv,cwd=str(root),
        gpu_uuid=p['gpu_uuid'],new_gpu_paths=12,exact_reused_paths=1,full_restart_replays=12,
        cpu_full_ladder_workers=0,protocol_sha256=audit['protocol_sha256'],
        prelaunch_audit_sha256=sha(here/'prelaunch-audit.json'),launcher_sha256=sha(__file__),
        new_scientific_histories=0,availability_before_launch=availability)
    (here/'process-launch.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
