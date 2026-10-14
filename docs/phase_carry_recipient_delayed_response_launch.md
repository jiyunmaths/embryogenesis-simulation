# Delayed recipient response — launch qualification, 8 October 2026

The [frozen protocol](phase_carry_recipient_delayed_response.md) applies challenges after a fixed 60-unit wait following local state matching. The accepted unpulsed fine t=570 states define every new initial value. Targets are not rematched at challenge time. Retain both original local references, both shared physical contexts, all six selected recipients and fractional/common-amount depletion.

All **eight design and numerical-pipeline tests pass**, including the full small 141-start pipeline. An initial duplicate `parent_context` keyword failure was fixed before preparation and retained in the launch archive. No scientific threshold changed.

The separate [literal input audit](../archive/study-launches/2026-10-08-recipient-delayed-response/input-audit.json) verifies exact unpulsed sources, operator/volume/ID preservation, both pulse designs and amount records. There are **45 backgrounds, 141 logical starts and 117 distinct initial-value problems**; identical dose preparations retain labels without adding independent evidence.

Every source passes the numerical donor gate. Maximum fine/other-method endpoint log discrepancy is **2.93e-13**. All 45 backgrounds pass the prespecified descriptive settling screen: maximum preceding 24-unit raw log drift is **6.71e-11**, and maximum challenge-time absolute log-concentration rate is **5.73e-14** per model time unit. These measured checks do not establish mathematically exact equilibrium or broader attractor uniqueness.

The [actual CUDA prefixes](../archive/study-launches/2026-10-08-recipient-delayed-response/gpu-prefix.json) qualify every start at all three timesteps. Identical physical operators are grouped into nine independent blocks and 23 chemical columns, with 66 positive-control padding columns excluded from reporting. The equations and shared kernels remain unchanged. Full-window reference, CPU/GPU, refinement, signed-contrast and recovery checks remain mandatory.

The full study computes **423 GPU paths, 423 concurrent CPU checks and 282 reference trajectories** within three existing mature histories, with zero new histories. Dual references qualify before full production. Use four single-thread CPU workers alongside physical GPU 1 as logical `cuda:0`:

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/home/sji/miniforge3/envs/deeplearning/bin/python \
  -m embryo.phase_carry_recipient_delayed_response run \
  --output output/studies/phase-carry-recipient-delayed-response --workers 4
```

The [launch verification](phase_carry_recipient_delayed_response_launch_verification.json) verifies **64 source pins and 369 input pins**, source/settling checks, GPU prefixes, tests and unchanged parent/ledger/manuscript guards. The [pre-launch snapshot](../archive/study-document-snapshots/2026-10-08-before-recipient-delayed-response-launch/manifest.json) preserves reporting and accepted parent documents. Outputs remain small chemical arrays in the writable workspace; no 3D fields are copied.

Delayed scientific results require a separate saved-array reconstruction before interpretation. Target chemistry and surrounding global amounts remain unmatched at challenge time, and physical contexts remain frozen. No ledger/manuscript promotion or subsequent experiment is part of this launch.
