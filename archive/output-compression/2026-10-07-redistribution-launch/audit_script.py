"""Read-only, bounded-memory storage candidates; never change old evidence."""
import concurrent.futures as futures
import collections
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'
DEST = ROOT / 'tmp/redistribution-storage-audit.json'
REVIEWS = [
    OUT / 'phase-carry-context-dosage-review/2026-10-07T211547Z-final/assessment.json',
    OUT / 'phase-carry-delayed-response-review/2026-10-07T105523Z-final/assessment.json',
]
PROTOCOLS = [OUT / f'{n}/protocol.json' for n in (
    'phase-carry-context-dosage', 'phase-carry-delayed-response')]

def digest(p):
    before = p.stat()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(4 << 20), b''):
            h.update(chunk)
    after = p.stat()
    assert (before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_ino, after.st_size, after.st_mtime_ns), p
    return h.hexdigest()

protected = {}
guard_files = {}
for p in REVIEWS + PROTOCOLS:
    d = json.loads(p.read_text())
    guard_files[str(p)] = digest(p)
    for key in ('evidence_sha256', 'source_sha256', 'input_sha256'):
        for filename, expected in d.get(key, {}).items():
            filename = str(Path(filename).resolve())
            assert filename not in protected or protected[filename] == expected, filename
            protected[filename] = expected
for p in (OUT / 'phase-carry-context-dosage').rglob('*'):
    if p.is_file():
        protected.setdefault(str(p.resolve()), None)
# Protect every review artifact, regardless of whether its assessment hashes itself.
for root in ('phase-carry-context-dosage-review', 'phase-carry-delayed-response-review'):
    for p in (OUT / root).rglob('*'):
        if p.is_file(): protected.setdefault(str(p.resolve()), None)

complete_states = {'completed', 'completed_with_unresolved_checks'}
roots = {}
for p in OUT.iterdir():
    s = p / 'status.json'
    if p.is_dir() and not p.name.startswith('.') and s.exists():
        try:
            data = json.loads(s.read_text())
            if data.get('state') in complete_states:
                roots[p] = {'status': str(s), 'state': data['state'], 'sha256': digest(s)}
        except (ValueError, OSError): pass

def eligible(p):
    root = OUT / p.relative_to(OUT).parts[0]
    return root in roots and str(p.resolve()) not in protected and not p.is_symlink()

all_npz = []
json_candidates = []
denied = ('protocol', 'status', 'summary', 'results', 'assessment', 'report', 'provenance')
for root in roots:
    for p in root.rglob('*'):
        if not p.is_file() or p.is_symlink(): continue
        n = p.stat().st_size
        if p.suffix == '.npz' and n > 8 << 20:
            all_npz.append(p)
        elif p.suffix == '.json' and n > 1 << 20 and eligible(p) and not any(x in p.name.lower() for x in denied):
            # Only raw trajectory/mesh/viewer records, not scientific decision reports.
            d = json.loads(p.read_text())
            raw_list = isinstance(d, list) and d and isinstance(d[0], dict) and any(k in d[0] for k in ('time', 'elapsed', 'step', 'cells'))
            frames = isinstance(d, dict) and isinstance(d.get('frames'), list)
            if raw_list or frames:
                json_candidates.append(p)
                del d

size_groups = collections.defaultdict(list)
for p in all_npz: size_groups[p.stat().st_size].append(p)
potential = [p for group in size_groups.values() if len(group) > 1 for p in group]
print('Protected paths', len(protected), 'completed roots', len(roots), flush=True)
print('JSON candidates', len(json_candidates), 'GiB', sum(p.stat().st_size for p in json_candidates) / 2**30, flush=True)
print('NPZ size-prefilter files', len(potential), 'groups', sum(len(g)>1 for g in size_groups.values()), flush=True)

hashes = {}
with futures.ThreadPoolExecutor(max_workers=2) as pool:
    for i, (p, sha) in enumerate(pool.map(lambda p: (p, digest(p)), potential), 1):
        hashes[p] = sha
        if i % 100 == 0: print('Hashed NPZ', i, '/', len(potential), flush=True)
duplicates = collections.defaultdict(list)
for p, sha in hashes.items(): duplicates[(p.stat().st_size, sha)].append(p)
duplicate_groups = []
for (size, sha), group in duplicates.items():
    if len(group) < 2: continue
    pinned = [p for p in group if not eligible(p)]
    keeper = sorted(pinned or group, key=str)[0]
    replace = [p for p in group if p != keeper and eligible(p)]
    if not replace: continue
    stats = {str(p): {'device': p.stat().st_dev, 'inode': p.stat().st_ino, 'nlink': p.stat().st_nlink} for p in group}
    # Hardlinked files do not yield additional bytes. Only unique nonkeeper inodes count.
    keeper_inode = (keeper.stat().st_dev, keeper.stat().st_ino)
    inode_paths = collections.defaultdict(list)
    for p in group: inode_paths[(p.stat().st_dev,p.stat().st_ino)].append(p)
    saved = 0
    for inode, members in inode_paths.items():
        if inode != keeper_inode and all(p in replace for p in members) and members[0].stat().st_nlink <= len(members): saved += size
    duplicate_groups.append({'sha256': sha, 'bytes': size, 'keeper': str(keeper), 'replace_candidates': list(map(str, replace)), 'same_device': len({p.stat().st_dev for p in group}) == 1, 'saving_bytes': saved, 'guards': stats})
print('NPZ duplicate groups',len(duplicate_groups),'GiB savings',sum(g['saving_bytes'] for g in duplicate_groups)/2**30,flush=True)

def measured_zstd(p):
    before=p.stat()
    h=hashlib.sha256(); count=0
    # SHA before compression plus stat afterwards detects accidental mutation.
    sha=digest(p)
    process=subprocess.Popen(['zstd','-q','-10','-T1','--stdout',str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    while True:
        chunk=process.stdout.read(4<<20)
        if not chunk: break
        count+=len(chunk); h.update(chunk)
    stderr=process.stderr.read().decode()
    assert process.wait()==0,(p,stderr)
    after=p.stat()
    assert (before.st_ino,before.st_size,before.st_mtime_ns)==(after.st_ino,after.st_size,after.st_mtime_ns),p
    return {'path':str(p),'sha256':sha,'original_bytes':before.st_size,'zstd10_bytes':count,'compressed_sha256':h.hexdigest(),'saving_bytes':before.st_size-count,'stat_guard':{'device':before.st_dev,'inode':before.st_ino,'mtime_ns':before.st_mtime_ns},'completed_root':str(OUT/p.relative_to(OUT).parts[0])}

compressions=[]
with futures.ThreadPoolExecutor(max_workers=2) as pool:
    for i, value in enumerate(pool.map(measured_zstd,json_candidates),1):
        compressions.append(value)
        if i%50==0: print('Measured JSON compression',i,'/',len(json_candidates),'saved GiB',sum(x['saving_bytes'] for x in compressions)/2**30,flush=True)
assert all(digest(Path(p))==sha for p,sha in guard_files.items())
disk=shutil.disk_usage(ROOT)
report={
    'created_unix':time.time(),'read_only':True,'no_old_files_modified':True,
    'free_bytes':disk.free,'requested_additional_bytes':int(4.09*2**30),
    'guard_files_sha256':guard_files,'protected_paths_sha256':protected,
    'protected_paths_count':len(protected),'completed_roots':{str(k):v for k,v in roots.items()},
    'duplicate_npz_groups':duplicate_groups,'dedup_saving_bytes':sum(g['saving_bytes'] for g in duplicate_groups),
    'json_zstd10_candidates':compressions,'json_zstd10_saving_bytes':sum(x['saving_bytes'] for x in compressions),
    'npz_size_prefilter_files':len(potential),'json_candidates':len(compressions),
    'limitations':['No deletion, hardlinking, archival, or compression was performed. Compression stream outputs were counted and discarded.', 'Candidates are outside both final immutable review maps, both source/input protocol registries, the entire current context-dosage study, and both review roots.', 'Hardlink savings apply only when all links to each replaced inode are accounted for. Archives would need manifest and round-trip verification before removal.']}
report['total_potential_saving_bytes']=report['dedup_saving_bytes']+report['json_zstd10_saving_bytes']
report['duplicate_keeper_breakdown'] = {}
for protected_keeper in (True, False):
    groups = [g for g in duplicate_groups if (g['keeper'] in protected) == protected_keeper]
    report['duplicate_keeper_breakdown']['protected' if protected_keeper else 'eligible'] = {
        'groups': len(groups), 'replace_files': sum(len(g['replace_candidates']) for g in groups),
        'saving_bytes_with_keeper_preserved': sum(g['saving_bytes'] for g in groups),
        'restore_method': 'Reference the protected existing keeper by path + SHA; restore to an independent file.' if protected_keeper else 'Keep eligible keeper in place, or archive a single verified blob for the entire identical-content group.'}
DEST.write_text(json.dumps(report,indent=2)+'\n')
print('Audit saved',DEST,'measured total GiB',report['total_potential_saving_bytes']/2**30,flush=True)
