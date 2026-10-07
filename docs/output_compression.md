# Historical output compression — 6 October 2026

Large files from completed historical studies are held in a lossless archive to
reduce local storage. This is a change in storage, not a change to the model or
its scientific evidence. File selection and measured savings are recorded in the
[operation report](../archive/output-compression/2026-10-06/report.md).

The archive contains large JSON histories, surface meshes, HTML viewers, and
byte-identical checkpoint copies from 16 older studies. Summaries, protocols,
status records, and analysis reports remain at their original paths. All 36
output folders found through the live delayed-response experiment's recursive
JSON dependency chain were excluded. Compiled kernels and simulation source
files were also retained in place.

Most NPZ files were already compressed. Only verified duplicate NPZ files were
selected: their exact original ZIP bytes share a compressed blob, rather than
sharing a writable checkpoint file. No arrays were rounded, converted to a lower
precision, downsampled, or regenerated. Each unique blob was fully decompressed
and checked against the original SHA-256 before its originals were removed.

## Restore data before using an archived study

Run these commands from the project directory. The installed `zstd` executable
must be available on `PATH`.

```bash
# See the original file paths and which are currently archived.
python archive/output-compression/2026-10-06/manage_archive.py list

# Restore just one study, including its archived viewers and checkpoints.
python archive/output-compression/2026-10-06/manage_archive.py restore outputs/domain-conservative-extended

# Restore one original file instead.
python archive/output-compression/2026-10-06/manage_archive.py restore outputs/domain-conservative-extended/grid-40/trajectory.json

# Restore every archived original, or verify the compressed archive.
python archive/output-compression/2026-10-06/manage_archive.py restore --all
python archive/output-compression/2026-10-06/manage_archive.py verify
```

Restore a study before rerunning its analysis, opening an archived viewer, or
checking its historical evidence hashes. The scientific code still expects its
original filenames; it does not automatically read compressed blobs.

Restoration checks both compressed and original hashes, writes each file
independently, preserves its permissions and modification time, and refuses to
overwrite a changed existing file. Already restored matching files are skipped.
The compressed copies remain available after restoration, so restoring all data
requires space for the full original files in addition to the archive.

## Archive and audit records

Payloads are stored under `outputs/.compressed-archive/2026-10-06/blobs/`, which
is ignored by Git along with other simulation outputs. The
[manifest](../archive/output-compression/2026-10-06/manifest.json) maps every
original path to a blob and records hashes, sizes, permissions, and timestamps.
Keep the payloads and manifest together when backing up the project. The blob
filenames alone do not record the original directory structure.

The [read-only audit](../archive/output-compression/2026-10-06/audit.json),
[selection script](../archive/output-compression/2026-10-06/audit_script.py),
[compression script](../archive/output-compression/2026-10-06/compress_script.py),
and [restore validation](../archive/output-compression/2026-10-06/restore-validation.json)
document this operation. The dated audit and compression scripts deliberately
refuse an unchanged second invocation. Historical result ledgers and scientific
hash records were not rewritten; restored originals retain their recorded
content hashes.
