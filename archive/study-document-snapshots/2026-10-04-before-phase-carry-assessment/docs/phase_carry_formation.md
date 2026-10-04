# Phase-carry confirmation through formation

Two matched moving GPU continuations are running from the original history-9 near-uniform mature restart through **240 elapsed model time units** (physical t=150–390). This extends the promising [six-unit precision diagnosis](geometry_precision.md) across nonlinear formation, including the original largest mismatch near elapsed 118. Results remain pending; the old quantitative failures are preserved.

## What is held fixed

The same 16-cell geometry, polarity, lineage, random streams and chemical start are used at dt=0.001875 and 0.0009375. Keep the 72³ grid, original interface width, conservative transport, beta=2, D_a=0.02, D_b=0.55, zero directional polarity-tension contrast, and activity-dependent tension/adhesion coefficients 0.25/0.35. Contacts still accumulate in float32. The active model has no prescribed A/B labels or downstream fate switch.

The only numerical intervention is the tested float64 carry of phase increments lost to float32 storage. Visible fields, occupancy, spatial geometry and force arithmetic retain their previous conventions. PyTorch owns resident GPU arrays and matrix products; custom CUDA performs mechanics and spatial geometry/polarity. This is not a full float64 backend.

Reuse the exact elapsed 0–6 histories and physical checkpoints from the completed short phase-carry pair. Re-encode only experiment metadata for the new protocol; preserve every physical array, including the accumulated residual. The first new step is from elapsed 6, **without resetting chemistry or residual**. Both transferred states reproduce two subsequent original-checkpoint steps bit for bit on the actual GPU, including their residuals.

There are **two numerical trajectories within one existing developmental history (9)** and zero new histories. These are mature-state formation tests, not fresh zygote-to-identity experiments. Jobs run sequentially on the scientific GTX 1080 Ti.

## Recorded decisions

- Preserve the original maximum absolute log chemical discrepancy limit **0.01** over every cell, both chemicals and all aligned observations. Do not time-align the paths to accept them.
- Retain the inherited limits for polarity, axis ratio, volume, transport, instantaneous growth, onset timing, growth-band crossing and late contrast agreement. Quantitative agreement and qualitative pattern survival are separate outcomes.
- Observe every 0.15 model time units; save carry-preserving atomic checkpoints every 3. Check volume error below 5%, equivalent radius at least four grid spacings, no clipping, finite positive state, dilution amount error at most 2e-14, and sampled boundary occupancy below 0.01.
- Compare the two carry runs over the first 60 and full 240 units. A pilot discrepancy is retained; these diagnostic runs continue through the authorized full window unless a physical/numerical quality screen fails.
- Compare each carry path with its immutable no-carry baseline at the same timestep. An intervention difference is not itself a numerical failure; only the matched within-method refinement tests agreement.
- Apply the existing independent frozen-endpoint chemical assay to each completed endpoint, retaining local basin/stability evidence separately from moving persistence.

The inherited late window is 24 units, with sustained chemical contrast defined by SD(log activator) above 0.1 throughout that window. All tolerances are copied from the previous long protocol. The original no-carry fine/finer raw discrepancy of 0.03750 remains an unresolved historical result even if this new method passes.

## How to interpret the outcome

Passing full-window agreement would support this numerical accumulation method in the tested history and regime. It would not establish full timestep/spatial convergence, full-horizon native/GPU equivalence of the new method, GPU cleavage, arbitrary-parameter acceptance, or autonomous cell identity. A changed onset, pattern or endpoint basin is reported rather than discarded.

The [separate prescribed geometry–chemistry replay](geometry_chemistry_coupling.md) tests time coupling against amount-based ODE references. Its geometry is fixed as a function of time and cannot respond to replayed chemistry; the moving confirmation retains reciprocal coupling. Both experiments are needed to distinguish mechanisms from numerical sensitivity.

## Execution and provenance

The [runner](../embryo/phase_carry_formation.py) creates a new protocol and result directory; no pinned historical scientific source is rewritten. Source/input hashes, accepted and experimental CUDA binary hashes, GPU/PyTorch versions, transfer checks and exact-state restart gates are recorded in `outputs/phase-carry-formation/`. An exclusive lock prevents duplicate coordinators. Interrupted jobs resume with their carry from the last complete checkpoint.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.phase_carry_formation run
```

Root and child `status.json` files show live progress. `summary.json` is written only after both full trajectories and endpoint assays finish. Numerical failures use `completed_with_unresolved_checks`; they are not promoted into the consolidated ledger or manuscript. The user-authorized GPU launch is recorded in `launch.json`.
