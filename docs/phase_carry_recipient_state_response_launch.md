# Recipient-state matching — launched 8 October 2026

The [scientific protocol](phase_carry_recipient_state_response.md) matches both recipient species to either measured reference pair, while preserving each surrounding chemical vector on either shared physical context. Each pulse is compared with its own freely evolving matched control. The target is not clamped after the initial reset.

The study has 45 distinct backgrounds and 93 distinct starting problems: 45 controls and 48 pulses. There are 186 primary plus 93 diagnostic GPU paths, 279 CPU verification paths and 186 independent references. The matched target's 10% fractional and common-amount pulses are identical within each physical context; one challenge is integrated. All paths are new continuations within **three existing mature histories**, with zero new histories.

Thirteen tests qualify the matching design and complete small concurrent solver pipeline. Three further tests qualify the [grouped operator layout](phase_carry_recipient_state_response_grouped.md), agreement with the separately blocked trajectories at all timesteps, and exact checkpoint resume. An initial missing-import error in a layout test was fixed before preparation; its log is retained and no tolerance changed.

The separate literal [input audit](../archive/study-launches/2026-10-08-recipient-state-response-grouped/input-audit.json) confirms recipient matching, unchanged surrounding chemistry, operator/capacity pairing and amount records. Every actual start also passes the short [CUDA/CPU check](../archive/study-launches/2026-10-08-recipient-state-response-grouped/gpu-prefix.json) at all three timesteps. Maximum raw log discrepancy is approximately `1.11e-15`, and maximum normalized discrepancy `1.27e-14`. Full-window numerical checks remain mandatory.

Nine measured operator blocks carry independent chemical trajectory columns instead of 45 repeated matrices. Grouping preserves all scientific states and creates no exchange between columns. The short prefix is roughly four times faster; full timing is not inferred as an accepted benchmark. The earlier separately blocked preparation remains unlaunched at `output/studies/phase-carry-recipient-state-response`, with its original pins and qualification evidence intact.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/home/sji/miniforge3/envs/deeplearning/bin/python -m embryo.phase_carry_recipient_state_response_grouped run \
  --output output/studies/phase-carry-recipient-state-response-grouped --workers 4
```

Physical GPU 1 is logical `cuda:0`. After dual-reference qualification, four single-thread CPU workers verify trajectories concurrently with float64 PyTorch production. Scientific integration, reactions, transport and response thresholds remain unchanged. No custom CUDA mechanics or geometry/polarity kernel is advanced or modified here.

The [launch verification](phase_carry_recipient_state_response_launch_verification.json) pins sources, inputs, test/prefix evidence, layout counts and unchanged ledger/manuscript guards. The [pre-launch snapshot](../archive/study-document-snapshots/2026-10-08-before-recipient-state-response-launch/manifest.json) preserves the accepted parent and reporting state. Raw outputs stay under the small writable workspace study path; parent arrays are read in place and no 3D fields are copied.

Results remain pending full numerical qualification and a separate raw-array reconstruction. This is an immediate-reset, whole-surrounding-network test; it does not match global amounts, isolate nearest neighbors, establish autonomous identity or add fresh carry-corrected development. Launch promotes no ledger/manuscript claim and starts no subsequent experiment.
