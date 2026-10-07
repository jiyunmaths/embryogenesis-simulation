from __future__ import annotations

import collections
import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
OUTPUTS = ROOT / 'outputs'
ACTIVE = OUTPUTS / 'phase-carry-network-context'
REPORT = Path('/tmp/embryo-outputs-cleanup-audit.json')
TEXT_SUFFIXES = {'.json', '.md', '.py', '.tex', '.bib', '.toml', '.yaml', '.yml', '.sh', '.txt', '.html'}
SKIP_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', 'node_modules'}
PATH_PATTERN = re.compile(r'outputs/[A-Za-z0-9_.+*/\-]+')

def relative(p: Path) -> str:
    return str(p.relative_to(ROOT))

def strings(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from strings(v)
    elif isinstance(obj, str):
        yield obj

files = {}
for base, dirs, names in os.walk(OUTPUTS):
    for name in names:
        p = Path(base) / name
        if not p.is_symlink() and p.is_file():
            s = p.stat()
            files[relative(p)] = {'bytes': s.st_size, 'allocated_bytes': s.st_blocks * 512,
                'mtime_ns': s.st_mtime_ns, 'inode': s.st_ino, 'device': s.st_dev}

references = collections.defaultdict(set)
scanned_files = []
scan_errors = []
for base, dirs, names in os.walk(ROOT):
    dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not (Path(base) / d).is_symlink()]
    for name in names:
        p = Path(base) / name
        if p.suffix not in TEXT_SUFFIXES or p.is_symlink():
            continue
        referrer = relative(p)
        try:
            # Scan streaming text to avoid parsing large numerical histories or embedded viewers.
            with p.open('r', encoding='utf-8', errors='replace') as f:
                tail = ''
                while chunk := f.read(1024 * 1024):
                    block = tail + chunk
                    for match in PATH_PATTERN.finditer(block):
                        raw = match.group(0)
                        idx = raw.find('outputs/')
                        target = raw[idx:].rstrip('.')
                        references[target].add(referrer)
                    tail = block[-4096:]
            scanned_files.append(referrer)
            if len(scanned_files) % 1000 == 0:
                print('Scanned', len(scanned_files), 'text files', flush=True)
        except OSError as exc:
            scan_errors.append({'file': referrer, 'error': str(exc)})

active_protocol = json.loads((ACTIVE / 'protocol.json').read_text())
active_paths = set()
for value in strings(active_protocol):
    for match in PATH_PATTERN.finditer(value):
        raw = match.group(0)
        active_paths.add(raw[raw.find('outputs/'):].rstrip('.'))

directory_inventory = []
for folder in sorted(OUTPUTS.iterdir()):
    if not folder.is_dir():
        continue
    prefix = relative(folder) + '/'
    rows = [(k, v) for k, v in files.items() if k.startswith(prefix)]
    direct = references.get(prefix.rstrip('/'), set())
    file_refs = set().union(*(references.get(k, set()) for k, _ in rows)) if rows else set()
    external = sorted(r for r in direct | file_refs if not r.startswith(prefix))
    pinned = sorted(k for k, _ in rows if k in active_paths)
    directory_inventory.append({
        'path': prefix.rstrip('/'), 'file_count': len(rows),
        'bytes': sum(v['bytes'] for _, v in rows),
        'allocated_bytes': sum(v['allocated_bytes'] for _, v in rows),
        'active': folder == ACTIVE, 'active_protocol_paths': pinned,
        'external_referrers': external,
        'unreferenced_files': [k for k, _ in rows if not references.get(k)],
    })

size_groups = collections.defaultdict(list)
for name, stat in files.items():
    if stat['bytes'] > 1024 * 1024 and not name.startswith(relative(ACTIVE) + '/'):
        size_groups[stat['bytes']].append(name)
duplicate_candidates = [{'bytes_each': size, 'paths': names}
    for size, names in size_groups.items() if len(names) > 1]

audit = {'captured_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
    'scope': 'Read-only outputs audit; absence of a literal path reference does not prove scientific data are unused.',
    'active_root': relative(ACTIVE),
    'active_protocol_sha256': hashlib.sha256((ACTIVE / 'protocol.json').read_bytes()).hexdigest(),
    'files': files, 'references': {k: sorted(v) for k, v in references.items()},
    'directories': directory_inventory,
    'active_protocol_output_paths': sorted(active_paths),
    'scanned_text_files': len(scanned_files), 'scan_errors': scan_errors,
    'same_size_candidates': duplicate_candidates}
REPORT.write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps({'report': str(REPORT), 'files': len(files),
    'logical_GiB': sum(v['bytes'] for v in files.values()) / 2**30,
    'scanned_text_files': len(scanned_files), 'references': len(references),
    'active_protocol_output_paths': len(active_paths), 'scan_errors': scan_errors}, indent=2))
print('Folders without external literal references other than the ledger inventory:')
for d in directory_inventory:
    refs = [r for r in d['external_referrers'] if r not in {'docs/results_ledger.json'}]
    if not refs:
        print(d['path'], round(d['bytes']/2**20, 3), 'MiB', 'files', d['file_count'])
print('Temporary bytecode references:', references.get('outputs/cytokinesis-baseline/__pycache__/model.cpython-312.pyc', set()))
