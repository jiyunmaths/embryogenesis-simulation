"""Resume the unchanged coordinator after verified storage migration.

Completed paths are verified and reused. Only the interrupted incomplete path
is restarted from its unchanged initial checkpoint; its entire attempt is kept.
The historical source folder/protocol keys remain logical through a narrow alias.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys

ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation');HERE=Path(__file__).resolve().parent
STUDY=ROOT/'output/studies/development-control-qualification'
spec=importlib.util.spec_from_file_location('relocated_study_alias',ROOT/'scripts/with_study_output_path_alias.py')
alias=importlib.util.module_from_spec(spec);spec.loader.exec_module(alias)

def sha(p):
    digest=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for b in iter(lambda:stream.read(1024*1024),b''):digest.update(b)
    return digest.hexdigest()

def write(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')

def main(action):
    manifest=HERE/'studies-manifest.json';manifest_sha=sha(manifest)
    with alias.base.validated_alias(manifest,manifest_sha,source_root=alias.SOURCE_ROOT,target_root=alias.TARGET_ROOT):
        sys.path.insert(0,str(ROOT));from embryo import development_control_qualification as c
        p=c.read(STUDY/'protocol.json');interrupted=c.read(HERE/'interruption.json')
        assert interrupted['worker_exited'] and not Path(f'/proc/{interrupted["old_pid"]}').exists()
        c.verify_preserved(p)
        assert c.read(HERE/'retirement.json')['passed']
        if action=='launch':
            assert not (HERE/'resume.json').exists() and not (STUDY/'summary.json').exists()
            assert interrupted['status_before']['stage']=='production'
            assert not list(STUDY.glob('*-refinement.json')) and not (STUDY/'replays').exists()
            case=interrupted['status_before']['case'];job=next(j for j in p['jobs'] if j['key']==case)
            folder=STUDY/case;assert not job['reused']
            completed=[]
            for item in p['jobs']:
                if (STUDY/item['key']/'gpu-result.json').exists():
                    assert c.case_report(STUDY,item,p)['passed'];completed.append(item['key'])
            assert case not in completed,'Interruption case completed; review before changing paths'
            archived=STUDY/'storage-interrupted-attempts/2026-10-10'/case
            assert not archived.exists();archived.parent.mkdir(parents=True,exist_ok=True)
            files=c.pin(f for f in folder.rglob('*') if f.is_file())
            folder.rename(archived);folder.mkdir();shutil.copy2(archived/'source.npz',folder/'source.npz')
            assert sha(folder/'source.npz')==p['input_sha256'][str(folder/'source.npz')]
            relocated={str((archived/Path(f).relative_to(folder)).resolve()):h for f,h in files.items()}
            for f,h in relocated.items():assert sha(f)==h
            write(HERE/'interrupted-attempt.json',dict(case=case,original_folder=str(folder),
                archived_folder=str(archived),files_sha256=relocated,source_sha256=sha(folder/'source.npz')))
            (HERE/'status-before-resume.json').write_bytes((STUDY/'status.json').read_bytes())
            c.write(STUDY/'status.json',dict(state='prepared',protocol_sha256=c.digest(STUDY/'protocol.json')))
            env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=p['gpu_uuid'],OPENBLAS_NUM_THREADS='1',
                OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
            availability=subprocess.check_output(['nvidia-smi','--id='+p['gpu_uuid'],
                '--query-gpu=uuid,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
            uuid,used,util=[x.strip() for x in availability.split(',')]
            assert uuid==p['gpu_uuid'] and int(used)<100 and int(util)==0
            argv=[sys.executable,str(HERE/'resume_worker.py'),'worker']
            with (HERE/'resume.log').open('x') as stream:
                child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            previous_seconds=sum(c.case_report(STUDY,j,p)['seconds'] for j in p['jobs']
                if j['key'] in completed and not j['reused'])
            write(HERE/'resume.json',dict(pid=child.pid,utc=datetime.now(timezone.utc).isoformat(),
                argv=argv,gpu_uuid=p['gpu_uuid'],completed_paths_reused=completed,
                interrupted_case_restarted=case,prior_completed_production_seconds=previous_seconds,
                manifest_sha256=manifest_sha,protocol_sha256=c.digest(STUDY/'protocol.json'),
                script_sha256=sha(__file__),source_and_model_unchanged=True,new_scientific_histories=0))
            print(c.read(HERE/'resume.json'),flush=True)
        else:
            plan=c.read(HERE/'resume.json');assert plan['manifest_sha256']==manifest_sha
            assert plan['script_sha256']==sha(__file__)
            original=c.execute_path;case=plan['interrupted_case_restarted']
            evidence=c.read(HERE/'interrupted-attempt.json');archived=Path(evidence['archived_folder'])
            def execute(root,job):
                report=Path(root)/job['key']/'gpu-result.json'
                if report.exists():
                    assert job['key'] in plan['completed_paths_reused']
                    result=c.case_report(root,job,p);assert result['passed']
                    print(job['key'],'verified completed path reused after migration',flush=True);return result
                result=original(root,job)
                if job['key']==case and result['passed']:
                    import numpy as np
                    old_history=c.read(archived/'gpu-history.json')
                    new_history=c.read(Path(root)/case/'gpu-history.json')
                    assert new_history[:len(old_history)]==old_history,'Interrupted history prefix changed'
                    checked=[]
                    for old in archived.glob('gpu-t[0-9][0-9].npz'):
                        with np.load(old,allow_pickle=False) as a,np.load(Path(root)/case/old.name,allow_pickle=False) as b:
                            assert set(a.files)==set(b.files)
                            for key in a.files:assert np.array_equal(a[key],b[key]),(old.name,key)
                        checked.append(old.name)
                    write(HERE/'prefix-reproduction.json',dict(passed=True,case=case,
                        exact_history_rows=len(old_history),exact_checkpoint_arrays=checked,
                        original_attempt_hashes_unchanged=all(sha(f)==h for f,h in evidence['files_sha256'].items())))
                return result
            c.execute_path=execute
            c.run(STUDY)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('launch','worker'))
    main(parser.parse_args().action)
