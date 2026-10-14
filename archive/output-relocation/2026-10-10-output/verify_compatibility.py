"""Check all scientific byte pins and completed control provenance after cutover."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,importlib.util,json,sys

ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation');HERE=Path(__file__).parent
spec=importlib.util.spec_from_file_location('output_migration_alias',ROOT/'scripts/with_study_output_path_alias.py')
alias=importlib.util.module_from_spec(spec);spec.loader.exec_module(alias)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

if __name__=='__main__':
    manifest=HERE/'manifest.json';digest=sha(manifest)
    assert json.loads((HERE/'cutover.json').read_text())['manifest_sha256']==digest
    with alias.base.validated_alias(manifest,digest,source_root=ROOT/'output',target_root=Path('/mnt/data/Documents/Python/embryogenesis-simulation/output')):
        sys.path.insert(0,str(ROOT))
        from embryo import development_control_qualification as c
        study=ROOT/'output/studies/development-control-qualification';p=c.read(study/'protocol.json')
        c.verify_preserved(p)
        interrupted=json.loads((HERE/'interruption.json').read_text());completed=[]
        for job in p['jobs']:
            report=study/job['key']/'gpu-result.json'
            if report.exists():
                saved=c.case_report(study,job,p);assert saved['passed'];completed.append(job['key'])
        assert set(interrupted['status_before']['completed'])<=set(completed)
        assert 'off120' in completed
        launch=json.loads((ROOT/'archive/study-launches/2026-10-10-development-controls/process-launch.json').read_text())
        assert c.digest(study/'protocol.json')==launch['protocol_sha256']
        old_proof=json.loads((ROOT/'archive/output-relocation/2026-10-07/cutover.json').read_text())
        assert sha(ROOT/'archive/output-relocation/2026-10-07/manifest.json')==old_proof['manifest_sha256']
        assert (ROOT/'outputs').is_symlink() and (ROOT/'outputs').resolve()==Path('/mnt/data/Documents/Python/embryogenesis-simulation/outputs')
        testlog=(HERE/'path-alias-tests.log').read_text();assert '11 passed' in testlog
        result=dict(passed=True,manifest_sha256=digest,checked_utc=datetime.now(timezone.utc).isoformat(),
            source_pins=len(p['source_sha256']),input_pins=len(p['input_sha256']),
            preserved_pins=len(c.read(p['preserved_manifest'])['preserved_sha256']),
            completed_cases=completed,protocol_sha256=launch['protocol_sha256'],
            baseline_reuse_original_provenance_preserved=True,older_plural_output_link_unchanged=True,
            alias_tests_passed=11,scientific_sources_unchanged=True,new_scientific_histories=0,
            scope='Storage byte/path compatibility; no new numerical accuracy or scientific acceptance.')
    target=HERE/'compatibility-verification.json';assert not target.exists()
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
