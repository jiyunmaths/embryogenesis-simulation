# Shared physical-context response assay — launched 8 October 2026

The [prespecified reciprocal control](phase_carry_shared_context_response.md) crosses the actual t=510 chemical concentration vectors with each measured transport operator and its matching volumes. Whole chemical states are preserved, amount changes are recorded, and pulses occur immediately on reassignment. Both shared physical contexts are tested; no trajectory supplies a reused scientific result.

There are 21 backgrounds, 69 logical cases and 57 byte-counted distinct initial-value problems per timestep. The study runs 138 primary plus 69 diagnostic GPU paths, 207 unchanged CPU verification paths and 138 independent reference trajectories. These are nested interventions in three existing histories, with zero new histories.

All 17 focused tests passed, including the complete small factorial solver pipeline, actual dose normalization, signed interaction algebra, own-control response analysis, amount/capacity bookkeeping and evidence tamper rejection. Every actual start also passed a 0.6-unit CUDA/CPU comparison at all three timesteps. Maximum raw log discrepancy was `8.89e-16`, and maximum pulse-normalized discrepancy was `7.76e-14`. See the [test evidence](../archive/study-launches/2026-10-08-shared-context-response/tests.json) and [GPU-prefix evidence](../archive/study-launches/2026-10-08-shared-context-response/gpu-prefix.json). Prefix checks do not qualify the full scientific horizon.

The [launch verification](phase_carry_shared_context_response_launch_verification.json) confirms 58 scientific source pins, 153 consumed input pins, unchanged ledger/manuscript guards, runtime versions and the [pre-launch reporting snapshot](../archive/study-document-snapshots/2026-10-08-before-shared-context-response-launch/manifest.json). Independent literal preparation and design reviews are retained in the [launch archive](../archive/study-launches/2026-10-08-shared-context-response/).

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/home/sji/miniforge3/envs/deeplearning/bin/python -m embryo.phase_carry_shared_context_response run \
  --output output/studies/phase-carry-shared-context-response --workers 4
```

Physical GPU 1 is logical `cuda:0`. Float64 PyTorch arrays and matrix operations use the unchanged chemical step; the fixed-context assay advances no mechanics or geometry/polarity kernels. CPU references and verification use four workers with one numerical-library thread each. Full-window dual-reference, timestep, normalized-response and convergence gates remain mandatory before independent assessment.

Small chemical outputs remain in the writable repository at `output/studies/phase-carry-shared-context-response`; parent data are read in place and no 3D fields are copied. Launch does not promote ledger or manuscript results. Recipient-state matching remains a separate later intervention.
