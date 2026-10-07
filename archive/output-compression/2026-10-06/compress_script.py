"""One-time operation for audit.json. Every original is verified before unlinking.

Do not edit simulation data or their hash records. Restoration recreates their
exact bytes; identical originals share a blob, never a writable hard link.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from manage_archive import RECORD, ROOT, sha256, verify_blob, zstd

MANIFEST = RECORD / "manifest.json"
BLOBS = ROOT / "outputs/.compressed-archive/2026-10-06/blobs"


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x") as stream:
        stream.write(json.dumps(value, indent=2) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def signature(path):
    s = path.stat()
    return dict(bytes=s.st_size, allocated_bytes=s.st_blocks * 512,
                mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns,
                inode=s.st_ino, device=s.st_dev, mode=s.st_mode, nlink=s.st_nlink)


def main():
    if MANIFEST.exists() or BLOBS.exists():
        raise FileExistsError("Dated archive already exists; use manage_archive.py to inspect/restore")
    audit = json.loads((RECORD / "audit.json").read_text())
    active = ROOT / audit["active"]
    protocol = active / "protocol.json"
    if sha256(protocol) != audit["protocol_sha256"]:
        raise ValueError("Live protocol changed since the audit")
    current = json.loads(protocol.read_text())
    required = {**current["source_sha256"], **current["input_sha256"]}
    before = {name: signature(Path(name)) for name in required}
    atomic_json(RECORD / "active-input-stat-before.json", before)
    selected = set(audit["text_candidates"])
    expected = {}
    for group in audit["duplicate_npz_groups"]:
        selected.update(group["paths"])
        expected.update({name: group["sha256"] for name in group["paths"]})
    selected = sorted(selected, key=lambda name: (ROOT / name).stat().st_size, reverse=True)
    for name in selected:
        folder = "/".join(name.split("/")[:2])
        if folder in audit["protected_roots"] or str(ROOT / name) in required:
            raise ValueError(f"Protected path in selection: {name}")
        if json.loads((ROOT / folder / "status.json").read_text()).get("state") != "completed":
            raise ValueError(f"Study is not completed: {folder}")
        path = ROOT / name
        if path.is_symlink() or path.stat().st_nlink != 1:
            raise ValueError(f"Refusing a shared file: {path}")
    manifest = dict(state="packing", started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                    active=audit["active"], active_protocol_sha256=audit["protocol_sha256"],
                    protected_roots=audit["protected_roots"],
                    codec="zstd 1.5.7, level 10 for text; level 3 for NPZ; one thread",
                    disk_free_before=shutil.disk_usage(ROOT).free, files=[])
    BLOBS.mkdir(parents=True)
    atomic_json(MANIFEST, manifest)
    blobs = {}
    started = time.monotonic()
    for i, name in enumerate(selected):
        path = ROOT / name
        original = signature(path)
        access = path.stat().st_atime_ns
        digest = sha256(path)
        if name in expected and digest != expected[name]:
            raise ValueError(f"Duplicate content changed: {name}")
        if signature(path) != original:
            raise ValueError(f"Original changed while hashing: {name}")
        blob = BLOBS / f"{digest}.zst"
        record = dict(path=name, sha256=digest, atime_ns=access, **original,
                      blob=str(blob.relative_to(ROOT)), original_removed=False)
        if digest not in blobs:
            temporary = blob.with_suffix(".zst.tmp")
            with temporary.open("xb") as stream:
                level = "-3" if path.suffix == ".npz" else "-10"
                subprocess.run([zstd(), "-q", "-T1", level, "-c", str(path)],
                               stdout=stream, check=True)
                stream.flush()
                os.fsync(stream.fileno())
            record.update(compressed_sha256=sha256(temporary),
                          compressed_bytes=temporary.stat().st_size,
                          compressed_allocated_bytes=temporary.stat().st_blocks * 512)
            os.replace(temporary, blob)
            verify_blob(record)
            blobs[digest] = {key: record[key] for key in
                ("compressed_sha256", "compressed_bytes", "compressed_allocated_bytes")}
        else:
            record.update(blobs[digest])
        # A verified blob and a durable recovery record exist before unlinking.
        manifest["files"].append(record)
        atomic_json(MANIFEST, manifest)
        if signature(path) != original:
            raise ValueError(f"Original changed while compressing: {name}")
        if sha256(protocol) != audit["protocol_sha256"]:
            raise ValueError("Live protocol changed; originals are still recoverable")
        path.unlink()
        record["original_removed"] = True
        atomic_json(MANIFEST, manifest)
        if (i + 1) % 20 == 0 or i + 1 == len(selected):
            print(f"Archived {i+1}/{len(selected)} files, {len(blobs)} unique blobs; "
                  f"{time.monotonic()-started:.1f} seconds", flush=True)
    unchanged = all(signature(Path(name)) == info for name, info in before.items())
    if not unchanged or sha256(protocol) != audit["protocol_sha256"]:
        raise ValueError("Protected input stat or protocol changed")
    # Verify every pinned scientific source and input against the live protocol,
    # not just its existence or modification time.
    for i, (name, expected_digest) in enumerate(required.items()):
        if sha256(Path(name)) != expected_digest:
            raise ValueError(f"Pinned input hash mismatch: {name}")
        if (i + 1) % 200 == 0:
            print(f"Verified {i+1}/{len(required)} live pinned source/input hashes", flush=True)
    original_bytes = sum(r["bytes"] for r in manifest["files"])
    original_allocated = sum(r["allocated_bytes"] for r in manifest["files"])
    compressed_bytes = sum(r["compressed_bytes"] for r in blobs.values())
    compressed_allocated = sum(r["compressed_allocated_bytes"] for r in blobs.values())
    manifest.update(state="complete", finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                    unique_blobs=len(blobs), original_bytes=original_bytes,
                    compressed_bytes=compressed_bytes, saved_bytes=original_bytes-compressed_bytes,
                    original_allocated_bytes=original_allocated,
                    compressed_allocated_bytes=compressed_allocated,
                    saved_allocated_bytes=original_allocated-compressed_allocated,
                    verified_pinned_files=len(required), protected_input_stats_unchanged=unchanged,
                    disk_free_after=shutil.disk_usage(ROOT).free,
                    active_status_after=json.loads((active / "status.json").read_text()))
    atomic_json(MANIFEST, manifest)
    print(json.dumps({key: manifest[key] for key in
                     ("state", "unique_blobs", "original_bytes", "compressed_bytes", "saved_bytes",
                      "saved_allocated_bytes", "verified_pinned_files", "disk_free_before", "disk_free_after")}, indent=2))


if __name__ == "__main__":
    main()
