# Signaling on replayed geometry

This diagnostic asks why the frozen-geometry causal screen sustains signal contrast while the mature moving full-feedback pilot does not reach its late-contrast criterion. It uses the completed full branch of `outputs/moving-causal`, starting from exactly its saved regulator deviations, fate values, and cell IDs at t=18. It does not draw a new perturbation or change the source mechanics.

## Available geometry and scope

The source stores complete phase fields every six time units through t=78, not at every mechanical step. We reconstruct conductances and measured cell volumes from these eleven checkpoints using the same contact adapter as the live model. Between snapshots, each conductance and volume is interpolated linearly. Symmetric nonnegative conductances and positive capacities are thereby retained, and the operator is assembled as

$$
\Delta(t)=-M(t)^{-1}K(t).
$$

This is an **approximate prescribed-history replay**, not exact replay of the mechanical trajectory. It smooths unresolved contact changes and volume transients. We deliberately test fidelity before attributing a live effect to either of those changes. The protocol and source hashes are written before extracting and integrating trajectories.

## Four diagnostic arms

| Arm | Reaction/exchange operator | Mechanical dilution |
|---|---|---|
| `frozen` | Initial conductances and initial capacities | None; geometry is fixed |
| `replay` | Interpolated conductances and capacities | Apply recorded/interpolated volume ratios |
| `no_dilution` | Same changing operator as replay | Omitted intentionally |
| `frozen_transport` | Initial operator throughout | Apply recorded/interpolated volume ratios |

Every arm uses the full Gierer–Meinhardt reactions and the verified joint signal/fate SSP-RK2 stages. Reaction/exchange uses pre-step geometry; replay dilution then preserves each species' amount during the volume update. This remains first-order splitting of the moving problem. The same mature cell IDs are retained without cleavage.

The last two arms are diagnostic interventions, not complete alternative physical models. Omitting dilution changes amounts as volumes change. A fixed initial operator applied to concentrations with changing volumes need not conserve the current volume-weighted amount either. These interventions can isolate terms algebraically, but their contributions may interact nonlinearly and cannot simply be added. Reactions also change regulator amounts in all arms.

## Predeclared checks

All four arms run at chemical timesteps 0.0075 and 0.00375, and with six-unit geometry snapshots versus twelve-unit decimation. Observations match the original full history every 0.6 time units. Checks require maximum signal RMS difference below 0.01, maximum continuous-fate difference below 0.05, and identical cell-wise final labels for:

- Chemical timestep halving in every arm.
- Snapshot decimation in every arm at the smaller timestep.
- The six-unit replay at the smaller timestep versus the recorded full trajectory.

Only if all prerequisites pass is the replay considered ready for provisional mechanism attribution. Decimation agreement alone is insufficient: two undersampled trajectories could agree with each other while both disagreeing with the source. If fidelity or sampling fails, capture denser geometry before interpreting the ablations.

The output records cell-to-cell activator standard deviation, late contrast over t=63–78, fate vectors and counts, and final regulator amounts. Frozen-equilibrium spectra are computed at each geometry snapshot. The trapezoidal integral of the instantaneous maximum spatial growth rate is reported only as a descriptive spectral quantity: it is not actual amplification, because it omits mode rotation, volume forcing, and nonlinear departures from equilibrium.

## Reproduction and outputs

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.geometry_replay \
  --source outputs/moving-causal --output outputs/geometry-replay-repeat
```

The source must have completed all moving controls with passing numerical-quality checks. Outputs include `protocol.json`, `geometry.npz`, all sixteen chemical trajectories, `spectra.json`, `comparison.json`, `RESULTS.md`, and a contrast comparison plot. Regression checks cover capacity-weighted flux conservation under interpolation, exact reduction to the verified joint solver on a constant history, and a single amount-preserving dilution after reaction.

## Completed sparse-snapshot pilot

The pilot completed all sixteen trajectories. Every chemical timestep-halving check passes, but snapshot decimation and fidelity to the recorded full trajectory fail. Replay versus recorded maximum signal RMS discrepancy is **0.01441** (limit 0.01); maximum fate discrepancy is **2.0160** (limit 0.05), and final cell-wise labels differ.

| Arm | Late activator SD range | Final A / B |
|---|---:|---:|
| Frozen | 0.41836–0.52504 | 7 / 9 |
| Approximate replay | 0.00646–0.01329 | 7 / 9 |
| No dilution | 0.01822–0.02263 | 8 / 8 |
| Frozen transport with volume forcing | 0.40364–0.53030 | 5 / 11 |

These approximate patterns motivate examining contact changes, but **do not establish their causal contribution in the recorded live run**. Source fidelity fails, and the frozen-transport/volume-forcing arm is itself strongly sensitive to snapshot spacing. Refining the chemical timestep does not recover missing mechanical history.

## Every-step capture and automatic replay

`embryo.dense_geometry_replay` captured the same full branch from its exact t=18 checkpoint. It saves symmetric conductances and individual measured volumes at every mechanical step (0.0075 time units), while leaving the equations and source files unchanged. At every original six-unit checkpoint, phase fields, joint chemical/fate state, and polarity must be bitwise identical to the source. Any mismatch stops the experiment.

After capture reaches t=78, it automatically runs the four arms at both chemical timesteps, using every-step geometry and a two-step decimation. It requires signal/fate/label fidelity to the recorded full trajectory at both chemical timesteps, plus the same per-arm refinement tolerances as the sparse pilot. The output separates exact source-continuation verification from approximate chemistry-only replay checks; interpolated substeps can differ from the original split evolution.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dense_geometry_replay \
  --source outputs/moving-causal --output outputs/dense-geometry-replay-repeat
```

The completed run is `outputs/dense-geometry-replay`; `status.json` reports capture/replay/completion, and `capture-audit.json` records checkpoint equality. Final decisions are written automatically to `comparison.json` and `RESULTS.md`. Dense capture requires another full mechanical continuation, so it is slower than the graph-only pilot. No source trajectories or acceptance tolerances are overwritten. This experiment still does not refine the underlying mechanical trajectory or establish multi-seed robustness.

## Dense replay completed

The every-step capture and all four replay controls are complete. All eleven source checkpoints match bitwise, and all timestep, snapshot-decimation, and source-fidelity checks pass. At the original step, replay differs from recorded chemistry by at most 1.97e-14 RMS and fate by 7.87e-13; the half-step discrepancies are 8.84e-6 and 0.01008, respectively, with identical final labels.

| Arm | Late activator standard deviation range |
|---|---:|
| Frozen geometry | 0.41836–0.52504 |
| Full geometry replay | 0.01350–0.02437 |
| Replay without dilution | 0.01763–0.02199 |
| Frozen transport with recorded volume forcing | 0.41239–0.49356 |

In this recorded history, changing transport is sufficient to keep contrast small even when explicit dilution is omitted. Recorded volume forcing with frozen transport does not reproduce that suppression. This supports changing transport as the principal tested contributor to the lost contrast, rather than dilution alone. It does not decompose conductance changes versus changing volume capacities, and the diagnostic interventions have the conservation limitations stated above.

This result belongs to the older fate-based mature geometry/seed. It cannot be automatically transferred to the new fresh-zygote attribute-development experiment. The underlying mechanical trajectory has not been independently refined by these chemistry-only replay checks.
