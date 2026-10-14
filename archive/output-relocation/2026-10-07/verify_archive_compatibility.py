"""Read-only archive checks against the completed outputs relocation."""
from __future__ import annotations

from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import sys
import time

RECORD = Path(__file__).resolve().parent
PROJECT = RECORD.parents[2]
spec = importlib.util.spec_from_file_location("relocated_archive_wrapper", PROJECT / "scripts/relocated_output_archive.py")
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


def cli(manager, *args):
    previous = sys.argv
    try:
        sys.argv = [str(manager.__file__), *args]
        manager.main()
    finally:
        sys.argv = previous


def run():
    started = time.monotonic()
    manifest = RECORD / "manifest.json"
    expected_sha = "edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02"
    if w.file_sha256(manifest) != expected_sha:
        raise ValueError("Unexpected relocation manifest")
    mapping = w.RelocationMapping(PROJECT, manifest)
    managers = {}
    before = {}
    for date in w.MANAGERS:
        path = PROJECT / "archive/output-compression" / date / "manage_archive.py"
        before[str(path.relative_to(PROJECT))] = w.file_sha256(path)
        before[str(path.with_name("manifest.json").relative_to(PROJECT))] = w.file_sha256(path.with_name("manifest.json"))
        manager = w.load_manager(path)
        mapping.install(manager)
        managers[date] = manager
    validation = PROJECT / "archive/output-compression/2026-10-07-redistribution-launch/restore-validation.json"
    restored = json.loads(validation.read_text())["restored_files"]
    original_stats = {row["path"]: w.path_stat(PROJECT / row["path"]) for row in restored}
    new = managers["2026-10-07-redistribution-launch"]
    new_index = json.loads((new.RECORD / "manifest.json").read_text())
    reference_path = restored[0]["path"]
    reference = next(row for row in new_index["files"] if row["path"] == reference_path)
    blob_path = restored[-1]["path"]
    blob = next(row for row in new_index["files"] if row["path"] == blob_path)
    # Keep the frozen manager's real guard checks, including both full passes.
    original_sha = new.sha256
    progress = dict(files=0, bytes=0)

    def progressing_sha(path):
        digest = original_sha(path)
        progress["files"] += 1
        progress["bytes"] += Path(path).stat().st_size
        if progress["files"] % 128 == 0:
            print(f"Guard hashes: {progress['files']} files, {progress['bytes'] / 2**30:.3f} GiB; elapsed {time.monotonic()-started:.1f}s", flush=True)
        return digest

    new.sha256 = progressing_sha
    print(f"Checking all {len(new_index['protected_stat'])} protected guards and idempotent restore of {reference_path}", flush=True)
    cli(new, "restore", reference_path)
    current_guards = {name: new.metadata(new.safe(name)) for name in new_index["protected_stat"]}
    if current_guards != new_index["protected_stat"]:
        raise ValueError("Historical protected inode identities were not preserved")
    new.stream_source(blob)
    new.stream_source(reference)
    old = managers["2026-10-06"]
    old_index = json.loads(old.MANIFEST.read_text())
    small = min(old_index["files"], key=lambda row: row["bytes"])
    cli(old, "verify", small["path"])
    for row in restored:
        path = PROJECT / row["path"]
        if w.path_stat(path) != original_stats[row["path"]] or w.file_sha256(path) != row["sha256"]:
            raise ValueError(f"Previously restored file changed: {path}")
    after = {name: w.file_sha256(PROJECT / name) for name in before}
    if after != before or w.file_sha256(manifest) != expected_sha:
        raise ValueError("Frozen archive manager, index or relocation manifest changed")
    outputs_pins = sum(Path(name).is_relative_to(PROJECT / "outputs") for name in current_guards)
    result = dict(passed=True, verified_utc=datetime.now(timezone.utc).isoformat(),
                  manifest_sha256=expected_sha, wrapper_sha256=w.file_sha256(PROJECT / "scripts/relocated_output_archive.py"),
                  verified_protected_pins=len(current_guards), relocated_protected_pins=outputs_pins,
                  all_protected_historical_inode_stats_match=True,
                  actual_idempotent_restore=dict(manager="2026-10-07-redistribution-launch", path=reference_path,
                                                  files_written=0, stat_and_sha256_unchanged=True),
                  actual_recoveries_verified=[dict(manager="2026-10-07-redistribution-launch", codec="zstd", path=blob_path),
                                              dict(manager="2026-10-07-redistribution-launch", codec="reference", path=reference_path),
                                              dict(manager="2026-10-06", codec="zstd", path=small["path"])],
                  previous_roundtrip_files_unchanged=[row["path"] for row in restored],
                  frozen_archive_files_sha256=before, frozen_archive_files_unchanged=True,
                  fixture_script="archive/output-relocation/2026-10-07/archive_fixture_smoke.py",
                  fixture_result="archive/output-relocation/2026-10-07/archive-fixture-smoke.json",
                  fixture_log="archive/output-relocation/2026-10-07/archive-fixture-smoke.log",
                  elapsed_seconds=time.monotonic()-started,
                  scope="Read-only real archive guard/recovery checks and idempotent restore of an existing verified file; no output data or original archive records changed.")
    (RECORD / "archive-compatibility.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    run()
