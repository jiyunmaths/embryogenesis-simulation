"""Private smoke fixtures for the launch archive; never uses project outputs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import zipfile

SCRIPT = Path(__file__).resolve().parents[1] / 'archive/output-compression/2026-10-07-redistribution-launch/manage_archive.py'
spec = importlib.util.spec_from_file_location('archive_tool', SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture():
    root = Path(tempfile.mkdtemp(prefix='redistribution-archive-smoke-'))
    m.ROOT = root
    m.RECORD = root / 'archive/output-compression/test'
    m.RECORD.mkdir(parents=True)
    study = root / 'outputs/complete'
    study.mkdir(parents=True)
    status = study / 'status.json'
    status.write_text('{"state":"completed"}\n')
    keeper = study / 'immutable.npz'
    with zipfile.ZipFile(keeper, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('phi.npy', b'float-data-' * 10000)
    duplicate = study / 'duplicate.npz'
    duplicate.write_bytes(keeper.read_bytes())
    raw = study / 'trajectory.json'
    raw.write_text('[{"time":1,"cells":[]},' + '{"time":1,"cells":[]},' * 10000 + '{}]\n')
    raw.chmod(0o640)
    os.utime(raw, ns=(123400000000, 123500000000))
    raw_copy = study / 'trajectory-copy.json'
    raw_copy.write_bytes(raw.read_bytes())
    none = study / 'pinned-none.txt'
    none.write_text('must retain inode and stat')
    source = root / 'source.py'
    source.write_text('source = 1\n')
    audit = {
        'protected_paths_sha256': {str(keeper): m.sha256(keeper), str(none): None, str(source): m.sha256(source)},
        'guard_files_sha256': {str(status): m.sha256(status)},
        'completed_roots': {str(study): {'status': str(status), 'state': 'completed', 'sha256': m.sha256(status)}},
        'json_zstd10_candidates': [{'path': str(p), 'sha256': m.sha256(p), 'original_bytes': p.stat().st_size,
                                   'stat_guard': {k: v for k, v in m.metadata(p).items() if k in ('device','inode','mtime_ns')},
                                   'completed_root': str(study)} for p in (raw, raw_copy)],
        'duplicate_npz_groups': [{'keeper': str(keeper), 'replace_candidates': [str(duplicate)],
                                 'sha256': m.sha256(keeper), 'bytes': keeper.stat().st_size,
                                 'guards': {str(p): {k:v for k,v in m.metadata(p).items() if k in ('device','inode','nlink')} for p in (keeper, duplicate)}}]
    }
    audit_path = root / 'audit.json'
    audit_path.write_text(json.dumps(audit))
    return root, raw, keeper, duplicate, none, source, audit_path


def cli(*args):
    sys.argv = [str(SCRIPT), *args]
    m.main()


checks = []
root, raw, keeper, duplicate, none, source, audit = fixture()
old = {str(p): (p.read_bytes(), m.metadata(p)) for p in (raw, keeper, duplicate, none, source)}
m.pack(audit)
index = json.loads((m.RECORD/'manifest.json').read_text())
assert index['state'] == 'completed' and len(index['files']) == 3
assert len({r['blob'] for r in index['files'] if r['codec']=='zstd'}) == 1
assert not raw.exists() and not duplicate.exists()
assert sum(r['codec']=='reference' for r in index['files']) == 1
for p in (keeper, none, source): assert m.metadata(p) == old[str(p)][1]
cli('verify')
cli('restore', 'outputs/complete/trajectory.json')
cli('restore', '--all')
for p in (raw, duplicate):
    assert p.read_bytes() == old[str(p)][0]
    assert (p.stat().st_mode, p.stat().st_mtime_ns) == (old[str(p)][1]['mode'], old[str(p)][1]['mtime_ns'])
assert duplicate.stat().st_ino != keeper.stat().st_ino
assert duplicate.stat().st_nlink == raw.stat().st_nlink == 1
cli('restore', '--all')
checks.append('pack, same-pack duplicate blob reuse, blob+reference verification, prefix/all restore, exact bytes/mode/mtime, independent inode, idempotence')
raw.write_bytes(b'changed')
try: cli('restore', 'outputs/complete/trajectory.json')
except FileExistsError: checks.append('changed original cannot be overwritten')
else: raise AssertionError('overwrite accepted')
try: m.safe('outputs/../secret', True)
except ValueError: checks.append('path traversal rejected')
else: raise AssertionError('traversal accepted')
link = root / 'outputs/link'
link.symlink_to(root / 'outputs/complete', target_is_directory=True)
try: m.safe('outputs/link/immutable.npz', True)
except ValueError: checks.append('symlink ancestor rejected')
else: raise AssertionError('symlink accepted')
keeper.write_bytes(b'changed immutable source')
try: cli('verify')
except ValueError: checks.append('modified protected reference rejected')
else: raise AssertionError('modified reference accepted')
root, raw, keeper, duplicate, none, source, audit = fixture()
raw.write_text('changed candidate')
try: m.pack(audit)
except ValueError: checks.append('stale candidate rejected before manifest/blob/unlink')
else: raise AssertionError('stale accepted')
assert not (m.RECORD/'manifest.json').exists() and duplicate.exists()
root, raw, keeper, duplicate, none, source, audit = fixture()
audit_data = json.loads(audit.read_text())
duplicate.unlink()
os.link(keeper, duplicate)
audit_data['duplicate_npz_groups'][0]['guards'][str(duplicate)] = {k:v for k,v in m.metadata(duplicate).items() if k in ('device','inode','nlink')}
audit.write_text(json.dumps(audit_data))
try: m.pack(audit)
except ValueError: checks.append('protected shared inode rejected')
else: raise AssertionError('shared pinned inode accepted')
assert raw.exists() and duplicate.exists() and not (m.RECORD/'manifest.json').exists()
print(json.dumps({'passed': True, 'checks':checks}, indent=2))
