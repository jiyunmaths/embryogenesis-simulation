"""Verify real-file restoration and write the dated report and folder notices."""
import collections
import datetime as dt
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

from manage_archive import MANIFEST, RECORD, ROOT, sha256, verify_blob

index = json.loads(MANIFEST.read_text())
assert index["state"] == "complete" and index["verified_pinned_files"] == 939
files = index["files"]
assert len(files) == 296 and all(r["original_removed"] for r in files)
assert all(not (ROOT / r["path"]).exists() for r in files)
assert all((ROOT / r["blob"]).is_file() for r in files)

examples = ["outputs/development-refinement/time-0.015/history.json",
            "outputs/domain-conservative-extended/grid-40/final_state.npz",
            "outputs/domain-conservative-extended/grid-40/progress_state.npz"]
subprocess.run([sys.executable, str(RECORD / "manage_archive.py"), "restore", *examples], check=True)
restored = []
for name in examples:
    record = next(r for r in files if r["path"] == name)
    path = ROOT / name
    info = path.stat()
    assert sha256(path) == record["sha256"]
    assert info.st_size == record["bytes"] and info.st_mtime_ns == record["mtime_ns"]
    assert stat.S_IMODE(info.st_mode) == stat.S_IMODE(record["mode"])
    assert info.st_nlink == 1
    restored.append(dict(path=name, sha256=record["sha256"], inode=info.st_ino,
                         exact_original_hash=True, mode_and_mtime_preserved=True, independent_inode=True))
assert restored[1]["inode"] != restored[2]["inode"]
# Keep the storage savings after the smoke check. Remove only these restored
# copies, after rechecking their bytes and the retained compressed replacement.
for name in examples:
    record = next(r for r in files if r["path"] == name)
    verify_blob(record)
    assert sha256(ROOT / name) == record["sha256"]
    (ROOT / name).unlink()

groups = collections.defaultdict(list)
for record in files:
    groups["/".join(record["path"].split("/")[:2])].append(record)
for folder, rows in sorted(groups.items()):
    note = ROOT / folder / "ARCHIVED_FILES.md"
    if note.exists():
        raise FileExistsError(note)
    note.write_text("# Some large historical files are compressed\n\n"
                    "These originals are preserved byte-for-byte in the lossless archive. "
                    "Restore this study before opening archived viewers, analysing its raw data, "
                    "or verifying its historical evidence hashes. From the project directory:\n\n"
                    "```bash\n"
                    f"python archive/output-compression/2026-10-06/manage_archive.py restore {folder}\n"
                    "```\n\n"
                    "See `docs/output_compression.md` and "
                    "`archive/output-compression/2026-10-06/manifest.json`.\n\n"
                    "Archived originals:\n\n"
                    + "\n".join(f"- `{r['path'][len(folder)+1:]}`" for r in sorted(rows, key=lambda r:r["path"]))
                    + "\n")

gib = lambda size: size / 2**30
report = ["# Lossless output compression — 6 October 2026", "",
          f"Saved **{gib(index['saved_allocated_bytes']):.2f} GiB** by archiving "
          f"**{len(files)} files from {len(groups)} completed historical studies**. "
          "The live delayed-response run continued throughout the operation.", "",
          f"Original data occupied {gib(index['original_bytes']):.3f} GiB; "
          f"{index['unique_blobs']} unique compressed blobs occupy "
          f"{gib(index['compressed_bytes']):.3f} GiB. The original {len([r for r in files if r['path'].endswith('.npz')])} "
          "duplicate NPZ paths share 59 archive blobs. Large text histories, surface meshes, "
          "and viewers account for the remaining 162 paths. Restored files remain independent.", "",
          f"Free filesystem space changed from {gib(index['disk_free_before']):.2f} to "
          f"{gib(index['disk_free_after']):.2f} GiB during the operation. This is a filesystem-wide "
          "snapshot while the simulation continued writing; per-file savings are measured "
          "separately. Small audit files and folder notices add a little overhead.", "",
          "## Verification", "",
          "- All 212 archive blobs passed complete decompression and original SHA-256 verification "
          "before originals were unlinked; the complete archive subsequently passed the restore tool's `verify` command.",
          "- All **939 live protocol source/input hashes** matched afterward, and their recorded file stats were unchanged. "
          "The protocol hash remained `" + index["active_protocol_sha256"] + "`.",
          "- The recursive scan read 1,477 JSON files and excluded all 36 dependency folders. "
          "Compiled kernels, live outputs, summaries, protocols and analysis reports were retained in place.",
          "- Restoration passed tests for exact binary/text recovery, preserved permissions/timestamps, "
          "separate inodes, refusal of changed originals and corrupt archives, and unsafe-path rejection.",
          "- The command-line restore was also exercised on an actual archived history and two identical "
          "checkpoints. Their hashes matched and their restored inodes were distinct. These three smoke-test "
          "copies were then returned to archive-only storage after re-verification.", "",
          "## Restore before historical analysis", "",
          "From the project directory:", "", "```bash",
          "python archive/output-compression/2026-10-06/manage_archive.py restore outputs/domain-conservative-extended",
          "```", "",
          "Replace the last argument with an archived file or study folder. Use `restore --all` "
          f"for every original, requiring approximately **{gib(index['original_bytes']):.2f} GiB** "
          "of additional free space; compressed copies remain available. "
          "See the [full guide](../../../docs/output_compression.md).", "",
          "## Archived studies", "",
          "| Study | Original files | Original GiB |", "| --- | ---: | ---: |"]
for folder, rows in sorted(groups.items()):
    report.append(f"| `{folder}` | {len(rows)} | {gib(sum(r['bytes'] for r in rows)):.3f} |")
report += ["", "The [manifest](manifest.json) records original and compressed hashes and exact file metadata. "
           "The [audit](audit.json) records dependency exclusions and duplicate selections. "
           "Payloads remain in `outputs/.compressed-archive/2026-10-06/blobs/`; keep them with the manifest.", ""]
(RECORD / "report.md").write_text("\n".join(report))
protocol = ROOT / index["active"] / "protocol.json"
assert sha256(protocol) == index["active_protocol_sha256"]
check = dict(captured_utc=dt.datetime.now(dt.timezone.utc).isoformat(), passed=True,
             real_file_cli_restore=restored, duplicate_restores_have_distinct_inodes=True,
             smoke_test_copies_removed_after_reverification=True,
             all_originals_archived=all(not (ROOT / r["path"]).exists() for r in files),
             active_status=json.loads((protocol.parent / "status.json").read_text()),
             archive_manifest_sha256=sha256(MANIFEST),
             scripts_sha256={p.name:sha256(p) for p in RECORD.glob("*.py")})
(RECORD / "postcheck.json").write_text(json.dumps(check, indent=2) + "\n")
print(json.dumps({k:check[k] for k in ("passed", "all_originals_archived", "active_status")}, indent=2))
