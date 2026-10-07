# Outputs cleanup — 6 October 2026

Removed **six redundant files (32.24 MiB)** and four empty directories. The active network-context experiment, its protocol inputs, completed scientific results, historical failures, compiled kernels, and offline viewers were retained.

The inventory contained approximately **46.63 GiB across 6,055 files in 140 directories**. These are inventory-time figures: the active simulation continues to write files. The audit scanned 4,372 text files across the repository, including ignored outputs, manuscript files, archived documents, manifests, and protocols. It found no read errors. Absence of a literal path reference is not sufficient evidence that simulation data are unused: runners also construct filenames dynamically.

## Removed files

| Files | Size | Why removal was safe |
| --- | ---: | --- |
| `outputs/cell-exchange-moving-incomplete-preparation/{source,initial_states,frozen_reference}.npz` | 22.26 MiB | All three were verified byte-for-byte duplicates of the corresponding files in the retained `outputs/cell-exchange-moving/` run. The abandoned folder had no protocol, status, or results. No individual file was referenced by the scan or active protocol. |
| `outputs/moving-initiation-controls-interim/2026-10-05T023407Z/partial-seed-7_{coarse,fine}_fixed-conductances.json` | 9.96 MiB | Orphan working copies with no assessment or incoming file references. Their 401 and 1 rows exactly matched prefixes of the retained completed coarse and fine histories. The separate assessed interim snapshot was preserved. |
| `outputs/cytokinesis-baseline/__pycache__/model.cpython-312.pyc` | 0.026 MiB | Regenerable Python bytecode; the historical source file was retained and verified. |

The two emptied preparation/snapshot directories and the bytecode directory were removed using `rmdir`. The already empty `outputs/live-transport-validation/` directory was also removed. The completed conservative-transport validation remains intact under `outputs/live-transport-validation-complete/`.

## Retained files and larger storage opportunities

About **33.34 GiB** was stored in directories containing paths referenced by the active protocol, including its live output directory. This is a directory-level total, not a claim that every byte in those folders is individually pinned. The largest inherited dependency, `outputs/phase-carry-exchange-response/`, occupied about **13.30 GiB** and cannot be removed just because its run has finished.

- **Scientific records:** histories, checkpoints, sources, full and pilot comparisons, failed checks, and audit bundles remain available. Completed or superseded runs can still be required to reproduce conclusions. The retention policy in [backends.md](backends.md) explicitly keeps historical benchmark outputs.
- **Compiled caches:** retained. The live experiment and recorded launch checks depend on specific CUDA libraries; native CPU runners also construct library paths dynamically. A directory named “cache” is not automatically disposable.
- **Generated playback:** 37 viewers with no explicit incoming file path references have paired `trajectory.json` files, totaling **0.93 GiB**. They are possible candidates for a later policy that regenerates viewers on demand. They were preserved because their historical templates and labels also matter; a trajectory plus the current renderer does not establish an identical replacement. Viewers explicitly linked from documents were preserved as well.
- **Possible duplicate data:** equal-size files suggest a theoretical maximum of roughly **6.6 GiB** of duplicate storage. This is a candidate screen, not a completed content-hash comparison. Required checkpoint paths must remain valid, and sharing mutable checkpoint files through hard links could corrupt independent runs.
- **Interim assessments:** retained when they have their own evidence hashes and saved assessment scripts. A final review does not by itself make the earlier audit bundle disposable.

The dated results ledger and its archived copies were not rewritten. Their directory inventories remain historical snapshots; the cleanup manifest records the subsequent removals.

## Verification and audit files

The cleanup checked each candidate's size, modification time, inode, references, and SHA-256 before removal. It then verified the retained replacement files, all **772 existing output paths extracted from the active protocol**, and the unchanged protocol hash. No simulation source or pinned document was edited. The worker was still running after cleanup.

- [Removal manifest and retained replacement hashes](../archive/output-cleanup/2026-10-06/manifest.json)
- [Directory decisions and larger storage candidates](../archive/output-cleanup/2026-10-06/candidates.json)
- [Compressed inventory and reference map before cleanup](../archive/output-cleanup/2026-10-06/inventory-before.json.gz)
- [Read-only audit script](../archive/output-cleanup/2026-10-06/audit_script.py)
- [Exact cleanup script used](../archive/output-cleanup/2026-10-06/cleanup_script.py)

The cleanup script is a record of this operation, not a general-purpose deletion tool. Its dated archive destination and candidate checks prevent an unchanged second invocation from proceeding.
