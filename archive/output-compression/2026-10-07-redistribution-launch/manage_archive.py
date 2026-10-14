"""Reversible, byte-exact completed-output archive; never overwrites evidence.

pack --audit tmp/redistribution-storage-audit.json
verify [outputs/original-file-or-prefix] / restore <prefix> / restore --all
All restored files have independent inodes, original modes, and original mtimes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import time
import zipfile

RECORD = Path(__file__).resolve().parent
ROOT = RECORD.parents[2]
CHUNK = 4 << 20


def safe(name, output=False):
    raw = Path(name)
    if ".." in raw.parts:
        raise ValueError(f"Path traversal: {name}")
    path = raw if raw.is_absolute() else ROOT / raw
    relative = path.relative_to(ROOT)
    if output and (not relative.parts or relative.parts[0] != "outputs"):
        raise ValueError(f"Not an output path: {name}")
    for part in (path, *path.parents):
        if part == ROOT.parent:
            break
        if part.is_symlink():
            raise ValueError(f"Symbolic link: {part}")
    path.resolve().relative_to(ROOT.resolve())
    return path


def sha256(path):
    before = metadata(path)
    with Path(path).open("rb") as stream:
        value = hashlib.file_digest(stream, "sha256").hexdigest()
    if metadata(path) != before:
        raise ValueError(f"File changed while hashing: {path}")
    return value


def metadata(path):
    value = Path(path).stat()
    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"Not a regular file: {path}")
    return {"device": value.st_dev, "inode": value.st_ino,
            "nlink": value.st_nlink, "bytes": value.st_size,
            "mtime_ns": value.st_mtime_ns, "mode": value.st_mode}


def fsync_dir(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(value, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
        fsync_dir(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def zstd():
    executable = shutil.which("zstd")
    if executable is None:
        raise RuntimeError("zstd must be available on PATH")
    return executable


def stream_source(record, destination=None):
    digest = hashlib.sha256()
    count = 0
    if record["codec"] == "reference":
        source = safe(record["reference"], output=True)
        if metadata(source) != record["reference_guard"]:
            raise ValueError(f"Reference changed: {source}")
        process = None
        stream = source.open("rb")
    else:
        source = safe(record["blob"], output=True)
        if sha256(source) != record["compressed_sha256"]:
            raise ValueError(f"Blob hash mismatch: {source}")
        process = subprocess.Popen([zstd(), "-q", "-d", "-c", str(source)],
                                   stdout=subprocess.PIPE)
        stream = process.stdout
    try:
        while block := stream.read(CHUNK):
            digest.update(block)
            count += len(block)
            if destination is not None:
                destination.write(block)
    finally:
        stream.close()
        if process is not None and process.wait() != 0:
            raise RuntimeError(f"Decompression failed: {source}")
    if record["codec"] == "reference" and metadata(source) != record["reference_guard"]:
        raise ValueError(f"Reference changed while reading: {source}")
    if (digest.hexdigest(), count) != (record["sha256"], record["bytes"]):
        raise ValueError(f"Recovery hash/size mismatch: {source}")


def guards(audit, original=None):
    """Hash every pinned source/input; None-valued pins retain exact inode/stat."""
    expected = dict(audit["protected_paths_sha256"])
    for name, digest in audit["guard_files_sha256"].items():
        if name in expected and expected[name] not in (None, digest):
            raise ValueError(f"Conflicting protected hashes: {name}")
        expected[name] = digest
    for item in audit["completed_roots"].values():
        expected[item["status"]] = item["sha256"]
    result = {}
    for name, digest in expected.items():
        path = safe(name)
        result[str(path)] = metadata(path)
        if digest is not None and sha256(path) != digest:
            raise ValueError(f"Protected hash mismatch: {path}")
    if original is not None and result != original:
        raise ValueError("Protected file inode/stat changed")
    return result


def completed(audit, path):
    root = ROOT / "outputs" / path.relative_to(ROOT / "outputs").parts[0]
    item = audit["completed_roots"].get(str(root))
    if item is None:
        raise ValueError(f"Not a completed study: {root}")
    status = safe(item["status"], output=True)
    if sha256(status) != item["sha256"]:
        raise ValueError(f"Completed status changed: {status}")
    state = json.loads(status.read_text()).get("state")
    if state != item["state"] or state not in ("completed", "completed_with_unresolved_checks"):
        raise ValueError(f"Study no longer completed: {status}")


def checked_original(audit, name, digest, size, guard, pins):
    path = safe(name, output=True)
    if str(path) in pins:
        raise ValueError(f"Candidate is protected: {path}")
    completed(audit, path)
    value = metadata(path)
    if value["bytes"] != size or any(value.get(key) != val for key, val in guard.items()):
        raise ValueError(f"Candidate stat mismatch: {path}")
    if sha256(path) != digest:
        raise ValueError(f"Candidate hash mismatch: {path}")
    return path, value


def create_blob(path, digest, size, level):
    directory = safe("outputs/.compressed-archive/2026-10-07-redistribution-launch/blobs", True)
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / f"{digest}.zst"
    if destination.exists():
        raise FileExistsError(f"Archive collision: {destination}")
    with tempfile.NamedTemporaryFile(dir=directory, prefix=".packing-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            result = subprocess.run([zstd(), "-q", f"-{level}", "-T1", "--stdout", str(path)],
                                    stdout=stream, check=True)
            stream.flush()
            os.fsync(stream.fileno())
            os.link(temporary, destination)  # Publish without overwriting any concurrent file.
            fsync_dir(directory)
        finally:
            temporary.unlink(missing_ok=True)
    record = {"codec": "zstd", "blob": str(destination.relative_to(ROOT)),
              "sha256": digest, "bytes": size, "level": level,
              "compressed_sha256": sha256(destination),
              "compressed_bytes": destination.stat().st_size}
    stream_source(record)
    return record


def pack(audit_path):
    audit_path = safe(audit_path)
    audit_sha = sha256(audit_path)
    audit = json.loads(audit_path.read_text())
    manifest_path = RECORD / "manifest.json"
    if manifest_path.exists():
        raise FileExistsError("Manifest already exists; refuse a second pack")
    pins = guards(audit)
    protected_inodes = {(m["device"], m["inode"]) for m in pins.values()}
    measured_json = {item["path"]: item for item in audit["json_zstd10_candidates"]}
    groups = []
    seen = set()
    for item in audit["json_zstd10_candidates"]:
        groups.append((item["sha256"], item["original_bytes"], 10,
                       [(item["path"], item["stat_guard"])], None))
    for item in audit["duplicate_npz_groups"]:
        keeper = safe(item["keeper"], True)
        # This launch deliberately excludes unprotected-keeper groups: JSON
        # compression plus immutable references supplies the required headroom.
        if str(keeper) not in audit["protected_paths_sha256"]:
            continue
        reference = str(keeper)
        names = list(item["replace_candidates"])
        groups.append((item["sha256"], item["bytes"], 3,
                       [(name, item["guards"][name]) for name in names], reference))
    # Validate the complete candidate set before creating blobs or removing files.
    ready = []
    for digest, size, level, candidates, reference in groups:
        entries = []
        for name, guard in candidates:
            path, value = checked_original(audit, name, digest, size, guard, pins)
            if str(path) in seen:
                raise ValueError(f"Duplicate candidate: {path}")
            if (value["device"], value["inode"]) in protected_inodes:
                raise ValueError(f"Candidate shares a protected inode: {path}")
            seen.add(str(path))
            entry = {"path": str(path.relative_to(ROOT)), "sha256": digest,
                     "bytes": size, "original_guard": value, "atime_ns": path.stat().st_atime_ns,
                     "mode": value["mode"], "mtime_ns": value["mtime_ns"], "removed": False}
            entries.append(entry)
        if level == 3:
            with zipfile.ZipFile(ROOT / entries[0]["path"]) as archive:
                if any(member.compress_type != zipfile.ZIP_DEFLATED for member in archive.infolist()):
                    raise ValueError("NPZ must already use DEFLATE compression")
        if reference and (sha256(safe(reference)) != digest or metadata(safe(reference))["bytes"] != size):
            raise ValueError("Immutable reference hash/size mismatch")
        ready.append((digest, size, level, entries, reference))
    index = {"version": 1, "state": "packing", "started_unix": time.time(),
             "audit": str(audit_path.relative_to(ROOT)), "audit_sha256": audit_sha,
             "protected_stat": pins, "files": []}
    atomic_json(manifest_path, index)
    created_blobs = {}
    for digest, size, level, entries, reference in ready:
        if reference:
            recovery = {"codec": "reference", "reference": str(safe(reference).relative_to(ROOT)),
                        "reference_guard": pins[reference], "sha256": digest, "bytes": size,
                        "compressed_bytes": 0}
            stream_source(recovery)
        else:
            if digest in created_blobs:
                recovery = created_blobs[digest]
                if (recovery["bytes"], recovery["level"]) != (size, level):
                    raise ValueError("Inconsistent duplicate blob metadata")
                stream_source(recovery)
            else:
                recovery = create_blob(ROOT / entries[0]["path"], digest, size, level)
                created_blobs[digest] = recovery
            if level == 10:
                measured = measured_json[str(ROOT / entries[0]["path"])]
                if "compressed_sha256" in measured and (
                    recovery["compressed_sha256"], recovery["compressed_bytes"]
                ) != (measured["compressed_sha256"], measured["zstd10_bytes"]):
                    raise ValueError("Compressed bytes differ from read-only measurement")
        for entry in entries:
            entry.update(recovery)
        index["files"].extend(entries)
        atomic_json(manifest_path, index)  # Durable recovery metadata precedes every unlink.
        for entry in entries:
            path = safe(entry["path"], True)
            current = metadata(path)
            # Unlinking sibling hardlinks can change nlink; inode/bytes/mtime/mode must remain fixed.
            expected = entry["original_guard"]
            if any(current[k] != expected[k] for k in expected if k != "nlink"):
                raise ValueError(f"Original changed before unlink: {path}")
            completed(audit, path)
            if sha256(path) != digest:
                raise ValueError(f"Original hash changed before unlink: {path}")
            path.unlink()
            fsync_dir(path.parent)
            entry["removed"] = True
            atomic_json(manifest_path, index)
        print(f"Archived {len(entries)} file(s), {size / 2**20:.2f} MiB each", flush=True)
    guards(audit, pins)
    if sha256(audit_path) != audit_sha:
        raise ValueError("Audit changed during pack")
    index["state"] = "completed"
    index["completed_unix"] = time.time()
    index["removed_bytes"] = sum(e["bytes"] for e in index["files"] if e["removed"])
    blobs = {e["blob"]: e["compressed_bytes"] for e in index["files"] if e["codec"] == "zstd"}
    index["archive_bytes"] = sum(blobs.values())
    index["logical_saving_bytes"] = index["removed_bytes"] - index["archive_bytes"]
    atomic_json(manifest_path, index)
    print(f"Completed: {len(index['files'])} originals, {len(blobs)} blobs; all protected guards unchanged")


def restore_file(record):
    destination = safe(record["path"], True)
    if destination.exists():
        if sha256(destination) != record["sha256"]:
            raise FileExistsError(f"Refusing to overwrite changed original: {destination}")
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".restoring-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream_source(record, stream)
            stream.flush()
            os.fchmod(stream.fileno(), stat.S_IMODE(record["mode"]))
            os.utime(temporary, ns=(record["atime_ns"], record["mtime_ns"]))
            os.fsync(stream.fileno())
            # This links only a newly created temporary inode, never an archived/reference file.
            os.link(temporary, destination)
            fsync_dir(destination.parent)
        finally:
            temporary.unlink(missing_ok=True)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("pack", "verify", "restore"))
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--audit", default="tmp/redistribution-storage-audit.json")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if args.action == "pack":
        if args.paths or args.all:
            parser.error("pack uses all audited candidates; paths/--all are unsupported")
        pack(args.audit)
        return
    index = json.loads((RECORD / "manifest.json").read_text())
    audit_path = safe(index["audit"])
    if sha256(audit_path) != index["audit_sha256"]:
        raise ValueError("Audit no longer matches manifest")
    audit = json.loads(audit_path.read_text())
    guards(audit, index["protected_stat"])
    rows = index["files"]
    if args.paths:
        prefixes = [str(safe(p, True).relative_to(ROOT)).rstrip("/") for p in args.paths]
        rows = [r for r in rows if any(r["path"] == p or r["path"].startswith(p + "/") for p in prefixes)]
        if not rows:
            parser.error("No matching archived files")
    elif args.action == "restore" and not args.all:
        parser.error("restore requires original paths or --all")
    checked = set()
    count = 0
    for row in rows:
        if args.action == "verify":
            key = row.get("blob", row.get("reference"))
            if key not in checked:
                stream_source(row)
                checked.add(key)
        else:
            count += restore_file(row)
    print(f"Verified {len(checked)} recoveries for {len(rows)} files" if args.action == "verify"
          else f"Restored {count} independent files; recoveries remain available")
    guards(audit, index["protected_stat"])


if __name__ == "__main__":
    main()
