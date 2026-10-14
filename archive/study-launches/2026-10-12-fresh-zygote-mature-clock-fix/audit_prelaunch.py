"""Independent alias/provenance audit for the zero-checkpoint scheduler fix."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=Path('output/studies/fresh-zygote-mature-clock-fix')
OLD=Path('output/studies/fresh-zygote-mature-qualification')
OLD_SHA='f6c099718c3d368646912660fb91fc2ce1280ce11130fb105cf07790f9811131'
PRIOR_AUDIT=Path('archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/prelaunch-audit.json')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    output=HERE/'prelaunch-audit.json';assert not output.exists()
    p=read(ROOT/'protocol.json');old=read(OLD/'protocol.json');ph=sha(ROOT/'protocol.json')
    assert sha(OLD/'protocol.json')==OLD_SHA
    prior=read(PRIOR_AUDIT);assert prior['passed'] and prior['protocol_sha256']==OLD_SHA
    assert prior['all24_complete_starts_preserved'] and prior['losslessly_checked_start_arrays']==432
    assert read(ROOT/'status.json')==dict(state='prepared',protocol_sha256=ph)
    for name,expected in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(name)==expected,name
    excluded={'jobs','study_root','source_sha256','input_sha256','stop_policy'}
    for k,v in old.items():
        if k not in excluded:assert p[k]==v,k
    for mapping in ('source_sha256','input_sha256'):
        for k,v in old[mapping].items():assert p[mapping][k]==v,k
    assert p['jobs']==[dict(j,source=str(Path(p['study_root'])/j['key']/'source.npz')) for j in old['jobs']]
    for job,original in zip(p['jobs'],old['jobs']):
        source=Path(job['source']);assert source.is_symlink() and source.resolve()==Path(original['source']).resolve()
        assert not (ROOT/job['key']/'gpu-result.json').exists()
    assert p['new_numerical_paths']==24 and p['new_scientific_histories']==0 and len(p['time_pairs'])==12
    assert read(OLD/'status.json')['state']=='failed'
    first=OLD/'h1101_off_coarse';failed=read(first/'gpu-result.json');prefix=read(first/'actual-start-prefix.json')
    assert failed['final_time']==18. and failed['failure']['elapsed']==0. and failed['failure']['message']=='duration and dt must be positive and finite'
    assert failed['lineage_cycle_rng_unchanged'] and failed['audit']['max_clipping']==0
    assert prefix['passed'] and prefix['steps']==8 and prefix['exact_gpu_restart'] and prefix['protocol_sha256']==OLD_SHA
    for k,v in prefix['errors'].items():assert v<=p['prefix_limits'][k],k
    for name,expected in prefix['evidence_sha256'].items():assert sha(name)==expected,name
    assert p['reused_prefix_case']=='h1101_off_coarse'
    assert p['scheduling_fix']=='Elapsed zero maps to integer step zero; positive durations retain the original resolver.'
    tests=read(HERE/'test-verification.json');assert tests['exit_code']==0 and sha(HERE/'tests.log')==tests['log_sha256']
    migration=Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration)=='8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        original_protocol_sha256=OLD_SHA,original_failure_retained=True,
        aliased_unchanged_start_files=24,original_array_audit_sha256=sha(PRIOR_AUDIT),
        reused_qualified_prefix_count=1,reused_prefix_sha256=sha(first/'actual-start-prefix.json'),
        model_integrator_and_thresholds_unchanged=True,new_numerical_paths=24,new_scientific_histories=0,
        test_log_sha256=sha(HERE/'tests.log'),test_verification_sha256=sha(HERE/'test-verification.json'),
        migration_record_sha256=sha(migration),runtime_binary_sha256=p['runtime_binary_sha256'],audit_script_sha256=sha(__file__))
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
