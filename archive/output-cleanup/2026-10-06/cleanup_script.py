"""Remove only the individually verified redundant files from the dated audit."""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
AUDIT = Path('/tmp/embryo-outputs-cleanup-audit.json')
OUT = ROOT / 'archive/output-cleanup/2026-10-06'
audit = json.loads(AUDIT.read_text())
active_root = ROOT / audit['active_root']
protocol = active_root / 'protocol.json'

def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def rel(p):
    return str(p.relative_to(ROOT))

assert digest(protocol) == audit['active_protocol_sha256']
active_paths = set(audit['active_protocol_output_paths'])
plan = []

def add(p, kind, retained=None, **details):
    name = rel(p)
    assert name in audit['files']
    assert name not in active_paths
    assert not p.is_relative_to(active_root)
    assert not audit['references'].get(name), (name, audit['references'].get(name))
    assert p.is_file() and not p.is_symlink()
    stat = p.stat()
    previous = audit['files'][name]
    assert (stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev) == (
        previous['bytes'], previous['mtime_ns'], previous['inode'], previous['device'])
    entry = {'path': name, 'kind': kind, 'bytes': stat.st_size,
        'allocated_bytes': stat.st_blocks * 512, 'sha256': digest(p), **details}
    if retained:
        entry['retained_path'] = rel(retained)
        entry['retained_sha256'] = digest(retained)
    plan.append(entry)

for name in ('source.npz', 'initial_states.npz', 'frozen_reference.npz'):
    p = ROOT / 'outputs/cell-exchange-moving-incomplete-preparation' / name
    retained = ROOT / 'outputs/cell-exchange-moving' / name
    assert digest(p) == digest(retained)
    add(p, 'unused_abandoned_preparation_exact_duplicate', retained)

for level in ('coarse', 'fine'):
    name = f'partial-seed-7_{level}_fixed-conductances.json'
    p = ROOT / 'outputs/moving-initiation-controls-interim/2026-10-05T023407Z' / name
    retained = ROOT / f'outputs/moving-initiation-controls/seed-7_{level}_fixed-conductances/history.json'
    partial = json.loads(p.read_text())
    complete = json.loads(retained.read_text())
    assert isinstance(partial, list) and len(partial) < len(complete)
    assert partial == complete[:len(partial)]
    add(p, 'orphan_working_snapshot_complete_history_prefix', retained,
        retained_prefix_rows=len(partial), complete_rows=len(complete))

p = ROOT / 'outputs/cytokinesis-baseline/__pycache__/model.cpython-312.pyc'
assert p.with_name(p.name.split('.')[0]).name == 'model'
source = p.parent.parent / 'model.py'
assert source.is_file()
add(p, 'regenerable_python_bytecode', source)

active_existence_before = {name: (ROOT / name).exists() for name in active_paths}
free_before = shutil.disk_usage(ROOT).free
OUT.mkdir(parents=True, exist_ok=False)
with gzip.open(OUT / 'inventory-before.json.gz', 'wb') as f:
    f.write(AUDIT.read_bytes())
shutil.copy2('/tmp/embryo-outputs-cleanup-audit.py', OUT / 'audit_script.py')
shutil.copy2(__file__, OUT / 'cleanup_script.py')

manifest = {'captured_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
    'state': 'verified_plan_before_cleanup', 'active_root': audit['active_root'],
    'active_protocol_sha256': audit['active_protocol_sha256'],
    'scanned_files': len(audit['files']), 'scanned_text_files': audit['scanned_text_files'],
    'scan_errors': audit['scan_errors'], 'removed': plan,
    'removed_directories': [], 'deleted_bytes': 0, 'free_bytes_before': free_before,
    'scope': 'Only unused exact preparation copies, orphan working snapshots recoverable from retained full histories, and bytecode. Scientific evidence and active inputs remain at their original paths.',
    'historical_inventory_note': 'The dated results ledger and its archives are unchanged; they continue to describe their original directory inventories.'}
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

for entry in plan:
    p = ROOT / entry['path']
    assert digest(p) == entry['sha256']
    p.unlink()
    manifest['deleted_bytes'] += entry['bytes']
    entry['removed'] = True
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

# rmdir only: never recursively remove a folder that has gained a file.
for name in (
    'outputs/cell-exchange-moving-incomplete-preparation',
    'outputs/moving-initiation-controls-interim/2026-10-05T023407Z',
    'outputs/cytokinesis-baseline/__pycache__',
    'outputs/live-transport-validation',
):
    p = ROOT / name
    if p.is_dir() and not any(p.iterdir()):
        p.rmdir()
        manifest['removed_directories'].append(name)

assert digest(protocol) == audit['active_protocol_sha256']
for name, existed in active_existence_before.items():
    assert not existed or (ROOT / name).exists(), name
for entry in plan:
    assert not (ROOT / entry['path']).exists()
    retained = ROOT / entry['retained_path']
    assert digest(retained) == entry['retained_sha256']
manifest.update({'state': 'completed', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
    'active_protocol_unchanged': True, 'active_protocol_existing_paths_preserved':
        sum(active_existence_before.values()), 'retained_replacements_verified': True,
    'deleted_allocated_bytes': sum(p['allocated_bytes'] for p in plan),
    'free_bytes_after': shutil.disk_usage(ROOT).free})
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({k: v for k, v in manifest.items() if k != 'removed'}, indent=2))
