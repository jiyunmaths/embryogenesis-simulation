"""Keep small tracked deliverables physical; link large output subdirectories."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,shutil

ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
SOURCE=ROOT/'output';TARGET=Path('/mnt/data/Documents/Python/embryogenesis-simulation/output')
HERE=Path(__file__).parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
if __name__=='__main__':
    manifest=HERE/'manifest.json';full=json.loads(manifest.read_text())
    assert json.loads((HERE/'retirement.json').read_text())['passed']
    assert SOURCE.is_symlink() and SOURCE.resolve()==TARGET
    assert not (HERE/'layout.json').exists()
    before=shutil.disk_usage(ROOT).free
    SOURCE.unlink();SOURCE.mkdir()
    for folder in ('figures','presentations'):
        shutil.copytree(TARGET/folder,SOURCE/folder,symlinks=True,copy_function=shutil.copy2)
        for name,entry in full['files'].items():
            if name.startswith(folder+'/'):assert sha(SOURCE/name)==entry['sha256'],name
    for folder in ('studies','kernel-cache','pdf'):
        (SOURCE/folder).symlink_to(TARGET/folder,target_is_directory=True)
    shutil.copystat(TARGET,SOURCE)
    prefix='studies/'
    files={name[len(prefix):]:entry for name,entry in full['files'].items() if name.startswith(prefix)}
    links={name[len(prefix):]:entry for name,entry in full['symlinks'].items() if name.startswith(prefix)}
    record=dict(passed=True,source_root=str(SOURCE/'studies'),target_root=str(TARGET/'studies'),
        files=files,symlinks=links,regular_files=len(files),
        parent_manifest=str(manifest.resolve()),parent_manifest_sha256=sha(manifest),
        scope='Exact verified study subtree of the full migration; filenames retain the original logical paths.')
    (HERE/'studies-manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),
        local_tracked_folders=['output/figures','output/presentations'],
        links={str(SOURCE/name):str(TARGET/name) for name in ('studies','kernel-cache','pdf')},
        restored_local_bytes=before-shutil.disk_usage(ROOT).free,
        studies_manifest_sha256=sha(HERE/'studies-manifest.json'),
        full_manifest_sha256=sha(manifest))
    (HERE/'layout.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
