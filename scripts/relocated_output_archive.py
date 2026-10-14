"""Run frozen output archives through the verified outputs relocation.

Examples::

    python scripts/relocated_output_archive.py 2026-10-06 verify
    python scripts/relocated_output_archive.py 2026-10-07-redistribution-launch restore outputs/example

The original managers and manifests remain unchanged. Only the approved outputs
root alias and verified device/inode changes are accepted. Packing is excluded.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
APPROVED_TARGET = Path("/mnt/data/Documents/Python/embryogenesis-simulation/outputs")
DEFAULT_MIGRATION = PROJECT_ROOT / "archive/output-relocation/2026-10-07/manifest.json"
MANAGERS = ("2026-10-06", "2026-10-07-redistribution-launch")
STAT_FIELDS = {
    "device": ("device", "dev", "st_dev"),
    "inode": ("inode", "ino", "st_ino"),
    "nlink": ("nlink", "st_nlink"),
    "bytes": ("bytes", "size", "st_size"),
    "mtime_ns": ("mtime_ns", "st_mtime_ns"),
    "mode": ("mode", "st_mode"),
    "uid": ("uid", "st_uid"),
    "gid": ("gid", "st_gid"),
}


def normalized_stat(value):
    result = {}
    for field, aliases in STAT_FIELDS.items():
        present = [value[name] for name in aliases if name in value]
        if not present or any(type(item) is not int for item in present):
            raise ValueError(f"Missing or invalid relocation stat: {field}")
        if any(item != present[0] for item in present):
            raise ValueError(f"Conflicting relocation stat: {field}")
        result[field] = present[0]
    return result


def path_stat(path):
    value = Path(path).lstat()
    return normalized_stat({name: getattr(value, "st_" + name)
                            for name in ("dev", "ino", "nlink", "size", "mtime_ns", "mode", "uid", "gid")})


def file_sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


class RelocationMapping:
    """Validate the one approved root alias and project historical inode stats."""

    def __init__(self, project_root, migration, approved_target=APPROVED_TARGET):
        self.root = Path(project_root).absolute()
        self.source = self.root / "outputs"
        self.target = Path(approved_target).absolute()
        self.migration = Path(migration)
        record = json.loads(self.migration.read_text())
        if record.get("passed") is not True:
            raise ValueError("Output relocation is not verified")
        if (record.get("source_root") != str(self.source) or
                record.get("target_root") != str(self.target)):
            raise ValueError("Relocation roots do not match the approved outputs alias")
        self.files = {}
        for name, entry in record["files"].items():
            raw = Path(name)
            if not name or raw.is_absolute() or ".." in raw.parts or raw.as_posix() != name:
                raise ValueError(f"Invalid relocation file path: {name}")
            if not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256", "")):
                raise ValueError(f"Invalid relocation digest: {name}")
            source = normalized_stat(entry["source_stat"])
            destination = normalized_stat(entry["destination_stat"])
            if not stat.S_ISREG(source["mode"]) or not stat.S_ISREG(destination["mode"]):
                raise ValueError(f"Relocation entry is not a regular file: {name}")
            if any(source[field] != destination[field] for field in STAT_FIELDS
                   if field not in ("device", "inode")):
                raise ValueError(f"Relocation changed file metadata: {name}")
            self.files[name] = dict(sha256=entry["sha256"], source=source, destination=destination)
        self._verified = set()
        self.check_layout()

    def check_layout(self):
        if not self.source.is_symlink() or os.readlink(self.source) != str(self.target):
            raise ValueError("Outputs is not the exact approved relocation symlink")
        for root in (self.root, self.target):
            for part in (root, *root.parents):
                if part.is_symlink():
                    raise ValueError(f"Unexpected root symlink: {part}")
        if not self.target.is_dir():
            raise ValueError("Relocated outputs directory is unavailable")

    def safe(self, name, output=False):
        raw = Path(name)
        if ".." in raw.parts:
            raise ValueError(f"Path traversal: {name}")
        path = raw if raw.is_absolute() else self.root / raw
        relative = path.relative_to(self.root)
        is_output = bool(relative.parts and relative.parts[0] == "outputs")
        if output and not is_output:
            raise ValueError(f"Not an output path: {name}")
        self.check_layout()
        for part in (path, *path.parents):
            if part == self.root.parent:
                break
            if part.is_symlink() and part != self.source:
                raise ValueError(f"Symbolic link: {part}")
        path.resolve().relative_to(self.target if is_output else self.root.resolve())
        return path

    def inside_root(self, name):
        raw = Path(name)
        if raw.is_absolute() or ".." in raw.parts or not raw.parts or raw.parts[0] != "outputs":
            raise ValueError(f"Invalid archive path: {name}")
        return self.safe(name, output=True)

    def metadata(self, path, original):
        """Return historical inode identity only for unchanged verified copies."""
        path = Path(path)
        current = original(path)
        try:
            relative = path.relative_to(self.source).as_posix()
        except ValueError:
            return current
        entry = self.files.get(relative)
        if entry is None:
            return current
        self.safe(path, output=True)
        before = path_stat(path)
        if before != entry["destination"]:
            raise ValueError(f"Relocated file stat changed: {path}")
        if current != {field: before[field] for field in current}:
            raise ValueError(f"Relocated file changed while checking stat: {path}")
        if relative not in self._verified:
            if file_sha256(path) != entry["sha256"]:
                raise ValueError(f"Relocated file hash changed: {path}")
            if path_stat(path) != before:
                raise ValueError(f"Relocated file changed while hashing: {path}")
            self._verified.add(relative)
        # Every other field retains its real current value; only cross-device
        # device/inode identity is projected to the migration's source identity.
        return dict(current, device=entry["source"]["device"], inode=entry["source"]["inode"])

    def install(self, manager):
        if Path(manager.ROOT).absolute() != self.root:
            raise ValueError("Archive manager belongs to a different project")
        if hasattr(manager, "safe"):
            original_metadata = manager.metadata
            manager.safe = self.safe
            manager.metadata = lambda path: self.metadata(path, original_metadata)
        elif hasattr(manager, "inside_root"):
            manager.inside_root = self.inside_root
        else:
            raise ValueError("Unsupported archive manager")


def load_manager(path):
    spec = importlib.util.spec_from_file_location("relocated_archive_manager", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--migration", type=Path, default=DEFAULT_MIGRATION)
    parser.add_argument("manager", help="Archive date or its original manage_archive.py path")
    parser.add_argument("archive_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    chosen = next((name for name in MANAGERS if args.manager == name or
                   Path(args.manager).absolute() == PROJECT_ROOT / "archive/output-compression" / name / "manage_archive.py"), None)
    if chosen is None:
        parser.error("Select 2026-10-06 or 2026-10-07-redistribution-launch")
    allowed = ("list", "verify", "restore") if chosen == "2026-10-06" else ("verify", "restore")
    if not args.archive_args or args.archive_args[0] not in allowed:
        parser.error(f"This wrapper accepts only: {', '.join(allowed)}")
    manager_path = PROJECT_ROOT / "archive/output-compression" / chosen / "manage_archive.py"
    manager = load_manager(manager_path)
    RelocationMapping(PROJECT_ROOT, args.migration).install(manager)
    previous = sys.argv
    try:
        sys.argv = [str(manager_path), *args.archive_args]
        manager.main()
    finally:
        sys.argv = previous


if __name__ == "__main__":
    main()
