# Storage recovered for the redistribution-response test — 7 October 2026

Saved **4.946 GiB** by losslessly archiving **443 older output files**: 409 raw trajectory JSON files and 34 exact duplicate NPZ files. Identical JSONs share 301 zstd blobs. Duplicate NPZs refer to 30 byte-identical protected keepers; restoration creates independent files, never writable hard links to those keepers.

Every payload was fully verified before removing its original path. Both accepted parent review evidence sets, their source/input registries, the completed dosage study and both review folders remained in place. The final checksum/stat pass confirms **2482 protected files and status guards** unchanged. Failed and unresolved scientific decisions were retained. Sources, kernels, ledger and manuscript were not edited.

The [manifest](manifest.json) records exact original hashes, permissions, timestamps and recovery locations. The [read-only audit](audit.json) records candidate selection and measured compression. Keep `outputs/.compressed-archive/2026-10-07-redistribution-launch/blobs/`, the manifest, the restore tool and its protected NPZ reference keepers together.

## Restore before reading an archived historical trajectory

From the project directory:

```bash
python archive/output-compression/2026-10-07-redistribution-launch/manage_archive.py restore outputs/moving-initiation-controls
```

Replace the final argument with a recorded original file or study prefix. `verify` checks all recoverable contents without restoring; `restore --all` recreates all originals with their exact bytes and original mode/mtime, requiring about **5.512 GiB** additional free space. Recovery blobs remain available. Existing changed originals are never overwritten. Keep protected reference files at their recorded locations.

This archive is separate from the [earlier compression](../2026-10-06/report.md); use each stage's own manifest and restore tool. Storage recovery enables launch, not scientific acceptance of the new response experiment.

## Actual recovery smoke test

A real trajectory and six ancestor source checkpoints were restored byte-for-byte, with their original mode/mtime and independent inodes. These copies remain in place for earlier reuse helpers. Net savings after this validation are **4.639 GiB. See [restore validation](restore-validation.json).
