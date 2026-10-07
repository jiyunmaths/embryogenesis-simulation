"""Lossless historical-output archive; run from any working directory.

list / verify / restore <original-file-or-directory> / restore --all
Archive payloads stay in ignored outputs; this script and the index are small.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import subprocess
import tempfile
from pathlib import Path

RECORD = Path(__file__).resolve().parent
ROOT = RECORD.parents[2]
MANIFEST = RECORD / "manifest.json"
CHUNK = 1024 * 1024


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def inside_root(name):
    path = Path(name)
    if path.is_absolute() or ".." in path.parts or path.parts[0] != "outputs":
        raise ValueError(f"Invalid archive path: {name}")
    destination = ROOT / path
    destination.resolve().relative_to(ROOT)
    if destination.is_symlink():
        raise ValueError(f"Refusing a symbolic link: {destination}")
    return destination


def zstd():
    executable = shutil.which("zstd")
    if executable is None:
        raise RuntimeError("zstd must be on PATH to read this archive")
    return executable


def read_blob(blob, destination=None):
    """Stream a complete decompression and validate the decoder's exit status."""
    digest = hashlib.sha256()
    size = 0
    with subprocess.Popen([zstd(), "-q", "-d", "-c", str(blob)],
                          stdout=subprocess.PIPE) as decoder:
        while block := decoder.stdout.read(CHUNK):
            digest.update(block)
            size += len(block)
            if destination is not None:
                destination.write(block)
        if decoder.wait() != 0:
            raise RuntimeError(f"Decompression failed: {blob}")
    return digest.hexdigest(), size


def verify_blob(record):
    blob = inside_root(record["blob"])
    if sha256(blob) != record["compressed_sha256"]:
        raise ValueError(f"Compressed file hash mismatch: {blob}")
    digest, size = read_blob(blob)
    if digest != record["sha256"] or size != record["bytes"]:
        raise ValueError(f"Original content hash/size mismatch: {blob}")


def restore_file(record):
    destination = inside_root(record["path"])
    if destination.exists():
        if sha256(destination) != record["sha256"]:
            raise FileExistsError(f"Refusing to overwrite changed file: {destination}")
        return False
    blob = inside_root(record["blob"])
    if sha256(blob) != record["compressed_sha256"]:
        raise ValueError(f"Compressed file hash mismatch: {blob}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent,
                                         prefix=".restoring-", delete=False) as stream:
            temporary = Path(stream.name)
            digest, size = read_blob(blob, stream)
            if digest != record["sha256"] or size != record["bytes"]:
                raise ValueError(f"Original content hash/size mismatch: {blob}")
            stream.flush()
            os.fsync(stream.fileno())
        temporary.chmod(stat.S_IMODE(record["mode"]))
        os.utime(temporary, ns=(record["atime_ns"], record["mtime_ns"]))
        # A link creates a separate restored inode, and cannot overwrite a file
        # that another process created while we were decompressing.
        os.link(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("list", "verify", "restore"))
    parser.add_argument("paths", nargs="*", help="Original paths, relative to the project")
    parser.add_argument("--all", action="store_true", help="Select the complete archive")
    args = parser.parse_args()
    index = json.loads(MANIFEST.read_text())
    rows = index["files"]
    if args.paths:
        prefixes = [str(inside_root(p).relative_to(ROOT)).rstrip("/") for p in args.paths]
        rows = [r for r in rows if any(r["path"] == p or r["path"].startswith(p + "/")
                                      for p in prefixes)]
        if not rows:
            parser.error("No archived files match the requested paths")
    elif args.action == "restore" and not args.all:
        parser.error("Give an original file/directory, or --all")
    if args.action == "list":
        for record in rows:
            state = "restored" if inside_root(record["path"]).exists() else "archived"
            print(f"{state:8} {record['bytes'] / 2**20:8.2f} MiB {record['path']}")
    elif args.action == "verify":
        checked = set()
        for record in rows:
            if record["blob"] not in checked:
                verify_blob(record)
                checked.add(record["blob"])
        print(f"Verified {len(checked)} unique blobs for {len(rows)} original files")
    else:
        restored = 0
        for record in rows:
            if restore_file(record):
                restored += 1
                print(f"Restored {record['path']}", flush=True)
        print(f"Restored {restored} files; compressed copies remain available")


if __name__ == "__main__":
    main()
