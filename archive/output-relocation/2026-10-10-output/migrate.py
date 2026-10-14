"""Verify and cut over a copied outputs tree without changing scientific bytes."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
from zoneinfo import ZoneInfo

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
SOURCE = ROOT/'output'
TARGET = Path('/mnt/data/Documents/Python/embryogenesis-simulation/output')
ARCHIVE = ROOT/'archive/output-relocation/2026-10-10-output'
MANIFEST = ARCHIVE/'manifest.json'
BACKUP = ROOT/'output.relocation-original-2026-10-10'


def read(path):
    return json.loads(Path(path).read_text())


def metadata(path):
    s = Path(path).lstat()
    return dict(dev=s.st_dev, ino=s.st_ino, mode=s.st_mode, uid=s.st_uid,
        gid=s.st_gid, nlink=s.st_nlink, size=s.st_size, mtime_ns=s.st_mtime_ns)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def inventory(root):
    files, links, directories = [], {}, []
    for folder, dirs, names in os.walk(root, followlinks=False):
        base = Path(folder)
        directories.append(str(base.relative_to(root)))
        for name in list(dirs):
            path = base/name
            if path.is_symlink():
                links[str(path.relative_to(root))] = os.readlink(path)
                dirs.remove(name)
        for name in names:
            path = base/name
            if path.is_symlink():
                links[str(path.relative_to(root))] = os.readlink(path)
            else:
                assert stat.S_ISREG(path.lstat().st_mode), f'Unsupported special file: {path}'
                files.append(str(path.relative_to(root)))
    return sorted(files), links, sorted(directories)


def verify():
    assert SOURCE.is_dir() and not SOURCE.is_symlink()
    assert TARGET.is_dir() and not TARGET.is_symlink()
    initial = read(ARCHIVE/'initial-copy.json')
    assert initial['source'] == str(SOURCE) and initial['destination'] == str(TARGET)
    assert initial['created_destination_empty'] and initial['initial_rsync_exit_code'] in (0, 24)
    stopped = read(ARCHIVE/'interruption.json')
    assert stopped['worker_exited'] is True
    assert not Path(f'/proc/{stopped["old_pid"]}').exists(), 'Old worker still exists'
    assert not MANIFEST.exists(), 'Already verified; do not overwrite migration proof'
    source_files, source_links, source_dirs = inventory(SOURCE)
    assert (source_files, source_links, source_dirs) == inventory(TARGET), 'Tree membership mismatch'
    compatible = ('mode', 'uid', 'gid', 'nlink', 'size', 'mtime_ns')

    def check(relative):
        a, b = SOURCE/relative, TARGET/relative
        before_a, before_b = metadata(a), metadata(b)
        assert all(before_a[k] == before_b[k] for k in compatible), f'Metadata mismatch: {relative}'
        with a.open('rb') as left, b.open('rb') as right:
            digest_a, digest_b = hashlib.sha256(), hashlib.sha256()
            while True:
                x, y = left.read(2**20), right.read(2**20)
                assert len(x) == len(y), relative
                if not x:
                    break
                digest_a.update(x)
                digest_b.update(y)
        assert digest_a.digest() == digest_b.digest(), f'Content mismatch: {relative}'
        assert metadata(a) == before_a and metadata(b) == before_b, f'File changed while verified: {relative}'
        return relative, dict(sha256=digest_a.hexdigest(), source_stat=before_a, destination_stat=before_b)

    checked = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for i, (relative, value) in enumerate(pool.map(check, source_files), 1):
            checked[relative] = value
            if i % 256 == 0 or i == len(source_files):
                print(f'Verified {i}/{len(source_files)} files', flush=True)
    symlinks = {}
    for relative, target in source_links.items():
        a, b = metadata(SOURCE/relative), metadata(TARGET/relative)
        assert all(a[k] == b[k] for k in compatible), relative
        symlinks[relative] = dict(target=target, source_stat=a, destination_stat=b)
    directories = {}
    for relative in source_dirs:
        a, b = metadata(SOURCE/relative), metadata(TARGET/relative)
        assert all(a[k] == b[k] for k in ('mode', 'uid', 'gid', 'mtime_ns')), relative
        directories[relative] = dict(source_stat=a, destination_stat=b)
    assert (source_files, source_links, source_dirs) == inventory(SOURCE) == inventory(TARGET)
    now = datetime.now(timezone.utc)
    result = dict(passed=True, source_root=str(SOURCE), target_root=str(TARGET),
        verified_utc=now.isoformat(), verified_local=now.astimezone(ZoneInfo('America/New_York')).isoformat(),
        regular_files=len(checked), symlink_count=len(symlinks), directory_count=len(directories),
        logical_bytes=sum(v['source_stat']['size'] for v in checked.values()),
        files=checked, symlinks=symlinks, directories=directories,
        initial_copy_sha256=sha(ARCHIVE/'initial-copy.json'),
        interruption_sha256=sha(ARCHIVE/'interruption.json'), verifier_sha256=sha(__file__),
        scope='Exact stable copied bytes, links, permissions and modification times before cutover; device/inode mapping preserves historical archive guards.')
    write(MANIFEST, result)
    print({k:v for k,v in result.items() if k not in ('files', 'symlinks', 'directories')}, flush=True)
    print('Manifest SHA256:', sha(MANIFEST), flush=True)


def cutover():
    record = read(MANIFEST)
    assert record['passed'] and record['source_root'] == str(SOURCE) and record['target_root'] == str(TARGET)
    assert not SOURCE.is_symlink() and SOURCE.is_dir() and not BACKUP.exists()
    assert TARGET.is_dir() and not TARGET.is_symlink()
    assert not Path(f'/proc/{read(ARCHIVE/"interruption.json")["old_pid"]}').exists()
    files, links, dirs = inventory(SOURCE)
    assert set(files) == set(record['files']) and set(links) == set(record['symlinks'])
    assert set(dirs) == set(record['directories'])
    assert (files, links, dirs) == inventory(TARGET)
    for relative, item in record['files'].items():
        assert metadata(SOURCE/relative) == item['source_stat'], relative
        assert metadata(TARGET/relative) == item['destination_stat'], relative
    for relative, item in record['symlinks'].items():
        assert os.readlink(SOURCE/relative) == os.readlink(TARGET/relative) == item['target']
        assert metadata(SOURCE/relative) == item['source_stat'], relative
        assert metadata(TARGET/relative) == item['destination_stat'], relative
    for relative, item in record['directories'].items():
        assert metadata(SOURCE/relative) == item['source_stat'], relative
        assert metadata(TARGET/relative) == item['destination_stat'], relative
    temporary = ROOT/'.output-relocation-link-2026-10-10'
    assert not temporary.exists() and not temporary.is_symlink()
    temporary.symlink_to(TARGET, target_is_directory=True)
    SOURCE.rename(BACKUP)
    try:
        temporary.rename(SOURCE)
    except BaseException:
        BACKUP.rename(SOURCE)
        temporary.unlink(missing_ok=True)
        raise
    assert SOURCE.is_symlink() and os.readlink(SOURCE) == str(TARGET)
    write(ARCHIVE/'cutover.json', dict(passed=True, source_root=str(SOURCE), target_root=str(TARGET),
        original_copy_retained=str(BACKUP), manifest_sha256=sha(MANIFEST),
        completed_utc=datetime.now(timezone.utc).isoformat()))
    print('Original outputs path now points to the verified new drive; original copy retained for validation.')


def retire():
    proof = read(ARCHIVE/'compatibility-verification.json')
    record = read(MANIFEST)
    assert record['passed'] and record['source_root'] == str(SOURCE) and record['target_root'] == str(TARGET)
    assert read(ARCHIVE/'cutover.json')['manifest_sha256'] == sha(MANIFEST)
    assert proof['passed'] and proof['manifest_sha256'] == sha(MANIFEST)
    assert SOURCE.is_symlink() and os.readlink(SOURCE) == str(TARGET)
    assert BACKUP.is_dir() and not BACKUP.is_symlink()
    assert TARGET.is_dir() and not TARGET.is_symlink() and TARGET.stat().st_dev != BACKUP.stat().st_dev
    files, links, dirs = inventory(BACKUP)
    assert set(files) == set(record['files']) and set(links) == set(record['symlinks'])
    assert set(dirs) == set(record['directories'])
    assert (files, links, dirs) == inventory(TARGET)
    for relative, item in record['files'].items():
        assert metadata(BACKUP/relative) == item['source_stat'], relative
        assert metadata(TARGET/relative) == item['destination_stat'], relative
    for relative, item in record['symlinks'].items():
        assert os.readlink(BACKUP/relative) == os.readlink(TARGET/relative) == item['target'], relative
        assert metadata(BACKUP/relative) == item['source_stat'], relative
        assert metadata(TARGET/relative) == item['destination_stat'], relative
    before = shutil.disk_usage(ROOT).free
    shutil.rmtree(BACKUP)
    after = shutil.disk_usage(ROOT).free
    write(ARCHIVE/'retirement.json', dict(passed=True, manifest_sha256=sha(MANIFEST),
        removed_original_copy=str(BACKUP), measured_freed_bytes=after-before,
        old_disk_free_bytes=after, new_disk_free_bytes=shutil.disk_usage(TARGET).free,
        completed_utc=datetime.now(timezone.utc).isoformat()))
    print('Removed only the verified original copy; measured freed GiB:', (after-before)/2**30)


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O; migration checks must remain enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'cutover', 'retire'))
    args = parser.parse_args()
    {'verify': verify, 'cutover': cutover, 'retire': retire}[args.phase]()
