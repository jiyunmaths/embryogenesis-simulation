"""Private relocation/archive fixtures; never touches scientific outputs."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

PROJECT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("relocated_archive_wrapper", PROJECT / "scripts/relocated_output_archive.py")
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


def archived_stat(path):
    value = w.path_stat(path)
    return {key: value[key] for key in ("device", "inode", "nlink", "bytes", "mtime_ns", "mode")}


def cli(manager, *args):
    previous = sys.argv
    try:
        sys.argv = ["manager", *args]
        manager.main()
    finally:
        sys.argv = previous


def rejects(call, exception=ValueError):
    try:
        call()
    except exception:
        return
    raise AssertionError("Unsafe or changed fixture was accepted")


def run():
    base = Path(tempfile.mkdtemp(prefix="relocated-archive-check-", dir="/tmp"))
    results = []
    for date in w.MANAGERS:
        root = base / date / "project"
        out = root / "outputs"
        study = out / "complete"
        study.mkdir(parents=True)
        record = root / "archive/output-compression" / date
        record.mkdir(parents=True)
        manager_path = PROJECT / "archive/output-compression" / date / "manage_archive.py"
        manager = w.load_manager(manager_path)
        manager.ROOT, manager.RECORD = root, record
        if hasattr(manager, "MANIFEST"):
            manager.MANIFEST = record / "manifest.json"
        status = study / "status.json"
        status.write_text('{"state":"completed"}\n')
        keeper = study / "keeper.npz"
        keeper.write_bytes(b"immutable scientific checkpoint" * 100)
        none = study / "pinned-none"
        none.write_bytes(b"pin without archived content hash")
        source = root / "source.py"
        source.write_text("frozen_source = True\n")
        raw = study / "trajectory.json"
        raw.write_bytes(b'{"trajectory":"exact"}\n' * 1000)
        raw.chmod(0o640)
        os.utime(raw, ns=(123400000000, 123500000000))
        blob = out / ".compressed-archive" / date / "blob.zst"
        blob.parent.mkdir(parents=True)
        with blob.open("wb") as stream:
            subprocess.run(["zstd", "-q", "-c", str(raw)], stdout=stream, check=True)
        row = dict(path=str(raw.relative_to(root)), blob=str(blob.relative_to(root)),
                   sha256=w.file_sha256(raw), bytes=raw.stat().st_size,
                   compressed_sha256=w.file_sha256(blob), mode=raw.stat().st_mode,
                   atime_ns=raw.stat().st_atime_ns, mtime_ns=raw.stat().st_mtime_ns, codec="zstd")
        raw_bytes = raw.read_bytes()
        raw.unlink()
        rows = [row]
        if date.endswith("launch"):
            duplicate = study / "duplicate.npz"
            rows.append(dict(path=str(duplicate.relative_to(root)), codec="reference",
                             reference=str(keeper.relative_to(root)), reference_guard=archived_stat(keeper),
                             sha256=w.file_sha256(keeper), bytes=keeper.stat().st_size,
                             mode=keeper.stat().st_mode, mtime_ns=keeper.stat().st_mtime_ns,
                             atime_ns=keeper.stat().st_atime_ns))
            audit = dict(protected_paths_sha256={str(keeper): w.file_sha256(keeper), str(none): None,
                                                str(source): w.file_sha256(source)},
                         guard_files_sha256={str(status): w.file_sha256(status)},
                         completed_roots={str(study): dict(status=str(status), state="completed", sha256=w.file_sha256(status))})
            audit_path = root / "audit.json"
            audit_path.write_text(json.dumps(audit))
            index = dict(files=rows, audit="audit.json", audit_sha256=w.file_sha256(audit_path),
                         protected_stat={str(path): archived_stat(path) for path in (keeper, none, source, status)})
        else:
            index = dict(files=rows)
        (record / "manifest.json").write_text(json.dumps(index))
        frozen_manifest = (record / "manifest.json").read_bytes()
        source_stats = {str(path.relative_to(out)): w.path_stat(path) for path in out.rglob("*") if path.is_file()}
        hashes = {str(path.relative_to(out)): w.file_sha256(path) for path in out.rglob("*") if path.is_file()}
        target = base / date / "external-outputs"
        shutil.copytree(out, target)
        shutil.rmtree(out)
        out.symlink_to(target, target_is_directory=True)
        migration = root / "migration.json"
        migration.write_text(json.dumps(dict(passed=True, source_root=str(out), target_root=str(target),
                                            files={name: dict(sha256=hashes[name], source_stat=value,
                                                              destination_stat=w.path_stat(target / name))
                                                   for name, value in source_stats.items()})))
        mapping = w.RelocationMapping(root, migration, approved_target=target)
        mapping.install(manager)
        cli(manager, "verify")
        cli(manager, "restore", "outputs/complete")
        cli(manager, "restore", "outputs/complete")
        assert raw.read_bytes() == raw_bytes
        assert (raw.stat().st_mode, raw.stat().st_mtime_ns) == (row["mode"], row["mtime_ns"])
        assert (record / "manifest.json").read_bytes() == frozen_manifest
        if date.endswith("launch"):
            assert duplicate.read_bytes() == keeper.read_bytes()
            assert duplicate.stat().st_ino != keeper.stat().st_ino
            assert manager.metadata(keeper) == archived_stat(keeper) | {
                "device": source_stats["complete/keeper.npz"]["device"],
                "inode": source_stats["complete/keeper.npz"]["inode"]}
        rejects(lambda: mapping.safe("outputs/../source.py", True))
        nested = out / "nested"
        nested.symlink_to(study, target_is_directory=True)
        rejects(lambda: mapping.safe("outputs/nested/keeper.npz", True))
        nested.unlink()
        raw.write_bytes(b"changed restored original")
        rejects(lambda: cli(manager, "restore", "outputs/complete/trajectory.json"), FileExistsError)
        if date.endswith("launch"):
            with none.open("r+b") as stream:
                stream.write(b"X")
            rejects(lambda: cli(manager, "verify"))
        out.unlink()
        out.symlink_to(base / "unapproved", target_is_directory=True)
        rejects(lambda: mapping.safe("outputs/complete", True))
        results.append(dict(manager=date, passed=True, checks=[
            "blob verify/restore", "mode/mtime preservation", "idempotence", "unchanged original manifest",
            "path traversal rejected", "nested symlink rejected", "changed root alias rejected",
            "changed original overwrite rejected", *(["independent reference restore", "historical inode projection",
                                                       "mutated None-valued pin rejected"] if date.endswith("launch") else [])]))
    return dict(passed=True, fixture_root=str(base), wrapper_sha256=w.file_sha256(PROJECT / "scripts/relocated_output_archive.py"), results=results)


if __name__ == "__main__":
    result = run()
    Path(__file__).with_name("archive-fixture-smoke.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
