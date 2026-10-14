# Frozen challenge-time response control — launched 8 October 2026

The [prespecified control](phase_carry_frozen_response.md) is launched from the exact t=510 starts of the [accepted moving redistribution-response study](phase_carry_redistribution_response_assessment.md). Each of nine measured transport operators and volume vectors stays fixed while its chemistry continues to t=570. Every original recipient, background, pulse factor and dose is retained. Outcomes remain pending independent assessment.

There are 66 primary production paths at two timesteps, plus 33 finest-timestep diagnostic paths. The 33 logical starting cases contain 27 unique initial-value problems; six fractional/common-amount preparations coincide exactly, and both dose records remain included. There are three existing histories and zero new histories. No moving trajectory supplies a frozen result.

Before launch, 113 relevant existing regression tests passed, with two existing optional moving-GPU tests skipped. The new runner passes 21 additional tests, including independent block transport, actual-pulse normalization, complete CPU/reference qualification and exact committed-state resume. The numerical reviewer found no remaining scientific blocker before freezing.

All 33 actual starting cases were additionally checked on the dedicated GTX 1080 Ti at all three timesteps through 0.6 units. Maximum GPU/CPU raw log discrepancy is `7.78e-16`; maximum pulse-normalized response discrepancy is `7.76e-14`. The [actual-prefix report](../archive/study-launches/2026-10-08-frozen-response/gpu-prefix.json) qualifies these prefixes only. The full reference, production and CPU checks remain mandatory.

An independent weak-dose diagnostic initially failed the normalized reference screen with unrestricted DOP853/Radau steps at ordinary tight tolerances. Using `rtol=1e-13`, `atol=1e-15` and maximum step `0.15` reduced disagreement below `1e-11`, without changing the `1e-8` acceptance screen. Both diagnostics are retained in the [qualification archive](../archive/study-launches/2026-10-08-frozen-response/qualification-evidence.json).

Production uses unchanged float64 PyTorch chemistry with independent, uncoupled context blocks. References and NumPy/SciPy verification run on four CPU workers, each with one numerical-library thread. No 3D fields are copied or advanced, and no custom CUDA mechanics or geometry kernel is modified. The small study is saved under `output/studies/phase-carry-frozen-response`, inside the writable repository.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/home/sji/miniforge3/envs/deeplearning/bin/python -m embryo.phase_carry_frozen_response run \
  --output output/studies/phase-carry-frozen-response --workers 4
```

Physical GPU 1 is logical `cuda:0` in this environment. Host NVIDIA driver access is required; the default execution sandbox cannot initialize CUDA. The initial sandbox attempt stopped at device availability, before any GPU stepping. The qualified host execution then passed the unchanged checks.

The [launch verification](phase_carry_frozen_response_launch_verification.json) records the protocol, 56 scientific source pins, 115 consumed input pins, runtime versions and preserved ledger/manuscript guards. The [pre-launch reporting snapshot](../archive/study-document-snapshots/2026-10-08-before-frozen-response-launch/manifest.json) preserves the accepted parent and prior status. No ledger or manuscript numerical claim is promoted by launch.

This tests whether physical evolution during the challenge window is necessary for a detected response effect given the inherited chemistry and graph. A shared operator and target-state matching remain separate controls. Runtime completion still requires an independent reconstruction before scientific acceptance.
