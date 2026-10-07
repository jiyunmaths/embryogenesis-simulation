# Lossless output compression — 6 October 2026

Saved **2.91 GiB** by archiving **296 files from 16 completed historical studies**. The live delayed-response run continued throughout the operation.

Original data occupied 4.564 GiB; 212 unique compressed blobs occupy 1.651 GiB. The original 134 duplicate NPZ paths share 59 archive blobs. Large text histories, surface meshes, and viewers account for the remaining 162 paths. Restored files remain independent.

Free filesystem space changed from 17.22 to 20.12 GiB during the operation. This is a filesystem-wide snapshot while the simulation continued writing; per-file savings are measured separately. Small audit files and folder notices add a little overhead.

## Verification

- All 212 archive blobs passed complete decompression and original SHA-256 verification before originals were unlinked; the complete archive subsequently passed the restore tool's `verify` command.
- All **939 live protocol source/input hashes** matched afterward, and their recorded file stats were unchanged. The protocol hash remained `d488e3268d154795d8229746ffd4a8adf0c9f689633eff51e68b00a14199442d`.
- The recursive scan read 1,477 JSON files and excluded all 36 dependency folders. Compiled kernels, live outputs, summaries, protocols and analysis reports were retained in place.
- Restoration passed tests for exact binary/text recovery, preserved permissions/timestamps, separate inodes, refusal of changed originals and corrupt archives, and unsafe-path rejection.
- The command-line restore was also exercised on an actual archived history and two identical checkpoints. Their hashes matched and their restored inodes were distinct. These three smoke-test copies were then returned to archive-only storage after re-verification.

## Restore before historical analysis

From the project directory:

```bash
python archive/output-compression/2026-10-06/manage_archive.py restore outputs/domain-conservative-extended
```

Replace the last argument with an archived file or study folder. Use `restore --all` for every original, requiring approximately **4.56 GiB** of additional free space; compressed copies remain available. See the [full guide](../../../docs/output_compression.md).

## Archived studies

| Study | Original files | Original GiB |
| --- | ---: | ---: |
| `outputs/attribute-development` | 40 | 0.180 |
| `outputs/attribute-development-parallel-smoke` | 4 | 0.005 |
| `outputs/attribute-development-smoke` | 4 | 0.005 |
| `outputs/causal-signaling` | 10 | 0.035 |
| `outputs/causal-signaling-confirmation` | 8 | 0.028 |
| `outputs/coupled-resolution` | 10 | 0.022 |
| `outputs/cutoff-dynamics` | 12 | 0.070 |
| `outputs/development-refinement` | 25 | 0.398 |
| `outputs/domain-conservative-extended` | 15 | 1.260 |
| `outputs/exchange-response-histories` | 58 | 0.631 |
| `outputs/fast-mechanics-validation` | 2 | 0.007 |
| `outputs/feedback-long` | 10 | 0.021 |
| `outputs/feedback-polarity-ablation` | 2 | 0.005 |
| `outputs/moving-causal` | 82 | 1.707 |
| `outputs/moving-causal-smoke` | 12 | 0.181 |
| `outputs/refinement-confirmation` | 2 | 0.010 |

The [manifest](manifest.json) records original and compressed hashes and exact file metadata. The [audit](audit.json) records dependency exclusions and duplicate selections. Payloads remain in `outputs/.compressed-archive/2026-10-06/blobs/`; keep them with the manifest.
