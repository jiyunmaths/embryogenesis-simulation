# Outputs relocation — 7 October 2026

Moved `outputs` to `/mnt/data/Documents/Python/embryogenesis-simulation/outputs` and kept the original repository path as a symbolic link. The original disk recovered **80.037 GiB**, with **89.236 GiB free** at cleanup. Code, Git history, manuscript and the separate `output/presentations` directory remain in the original repository.

The destination was newly created and empty. `rsync -aH` copied the tree; a final stopped-tree synchronization changed no files. Every one of **8,253 regular files**, totaling **80.010 GiB**, passed independent source/destination SHA-256 and stable metadata comparisons. All **14 symlinks** retained their targets and metadata, and all **1,757 directories** matched. Permissions, ownership and modification times were preserved. Original and destination device/inode identities are recorded separately in the immutable migration manifest.

Before removing the original copy, CPU-only recomputation reproduced the existing history-7 pilot and full numerical reports exactly, including their original evidence path keys. All **53 scientific source hashes**, 85 prepared/checkpoint/report hashes and 461 current-study file metadata records remained unchanged. Both historical archive recovery paths passed; all **2,482 protected guards** matched, including 2,416 relocated files. The idempotent restore wrote zero files. All seven retained prior recovery copies remain exact. **52 manuscript/ledger document guards** remained unchanged.

The previous worker PID 229767 was already absent when inspected at 21:33 EDT. Its log last changed at 21:00 EDT; its exact exit time and cause are unknown. No stop signal was sent. Its committed history-8 checkpoint at model **t=513, elapsed=3**, step 273600, was preserved exactly. The unchanged scientific protocol was relaunched at **2026-10-07T21:48:34.314532-04:00**, PID **10839**, on physical GPU 1, using the external path-compatibility launcher and the original numerical thread settings. Independent observation confirmed the process on the qualified GPU and the checkpoint advanced to **t=516, elapsed=6**, step 275200. All 21 original history rows remain exact, and cell IDs, job definition and protocol hash are unchanged. The execution evidence is recorded in `running-verification.json`.

The launcher preserves logical paths only under the relocated outputs root, while normal file operations follow the link. The archive wrapper accepts only the approved root alias and verifies contents/current metadata before translating the historical device/inode identities. No frozen scientific source, equation, parameter, protocol, archive manager or original archive manifest was edited. Numerical and scientific outcomes of the continuing assay remain pending.

## Evidence

- [Verified migration manifest](manifest.json), SHA-256 `edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02`.
- [Initial transfer](initial-copy.json), [transfer log](initial-copy.log), [final sync](final-sync.json), [verification log](verification.log), and [exact verification/cutover source](verify_and_cutover_script.py).
- [Scientific compatibility](scientific-compatibility.json), [archive compatibility](archive-compatibility.json), [archive fixture checks](archive-fixture-smoke.json), and [combined compatibility](compatibility-verification.json).
- [Cutover](cutover.json), [original-copy retirement](retirement.json), [strengthened retirement source](retirement_script.py), and [interruption/checkpoint record](interruption.json).
- [Worker launch](resume.json), [resume source](resume_worker.py), and [independent running observation](running-verification.json).
- [Path-alias tests](path-alias-tests.log) and [storage usage instructions](../../../docs/output_storage.md).
