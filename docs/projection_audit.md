# Audit of the failed spatial field comparison

The completed coupled resolution screen failed one declared check: the 72³/88³ final occupancy fields differ by 1.4306% after linear interpolation to a 56³ probe grid, exceeding the 1% criterion. The original report remains unchanged. This follow-up tests the **measurement**, using saved states; it does not run a different solver or change physical parameters.

## Why calibrate interpolation first?

Every initial grid samples the same analytic collection of diffuse spheres. Let $S_h u_0$ denote those samples and $P_h$ interpolation onto common probe points. Even though the continuous initial state is identical,

$$
P_h S_h u_0 - P_H S_H u_0
$$

need not vanish. It measures sampling and reconstruction differences, without any evolution error. For the original diagnostic, this initial discrepancy is **1.4943%**, already exceeding the final-state acceptance threshold before a single simulation step.

At later times, the projected difference contains both solver and reconstruction errors. Subtracting the two scalar error norms is invalid: errors are fields that can reinforce or cancel. The audit therefore reports their initial/final alignment and an evolution-increment comparison separately, without treating either as a rigorous solver-error correction.

## Reproducible audit

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.projection_audit \
  --source outputs/coupled-resolution --output outputs/projection-audit-repeat
```

Compare linear, cubic-spline, and quintic-spline interpolation on probe grids 56³, 88³, and 112³. For each case, reconstruct the analytic initial state directly on the probe grid as an independent reference. Compare both initial and archived final occupancy fields. Same-grid samples are used directly without interpolation. Spline prefiltering is enabled for orders three and five. Nearest-value extension is used outside the native voxel-center range; these outer strips contain negligible cell tails in this controlled geometry.

Interpolation orders here describe the **diagnostic reconstruction**, not the solver's spatial or temporal accuracy. No interpolation is fed back into dynamics. No clipping hides spline undershoot or overshoot; initial/final value ranges are exported. The original comparison and every source checkpoint are identified by SHA-256 hashes, and the source report is verified unchanged afterward.

## Completed results

| Probe grid | Reconstruction | Initial finest-pair difference | Final finest-pair difference |
|---:|---|---:|---:|
| 56³ | Linear | 1.4943% | 1.4306% |
| 56³ | Cubic | 0.1142% | 0.3152% |
| 56³ | Quintic | 0.0666% | 0.3184% |
| 88³ | Linear | 3.1010% | 2.9286% |
| 88³ | Cubic | 0.1544% | 0.2366% |
| 88³ | Quintic | 0.0729% | 0.2653% |
| 112³ | Linear | 1.4599% | 1.3823% |
| 112³ | Cubic | 0.1354% | 0.2629% |
| 112³ | Quintic | 0.0701% | 0.2663% |

All six higher-order final comparisons are below 1%, spanning **0.2366%–0.3184%**. Their largest single-grid initial reconstruction error against the analytic reference is **0.1833%**. The original linear metric is therefore substantially contaminated by reconstruction error; it cannot be read as a clean estimate of evolution error.

Higher-order splines are not positivity-preserving. Among the finer cases, the smallest reconstructed occupancy is about **−1.79e-4** for quintic interpolation (cubic undershoot is approximately −7.62e-6). These diagnostic values were recorded without clipping and were never fed into the simulation. This is a reason to retain calibration and cross-method checks, rather than assume higher order is automatically correct.

Results and the comparison figure are in `outputs/projection-audit/{RESULTS.md,analysis.json,projection.png}`. The original `outputs/coupled-resolution/comparison.json` hash is unchanged. Four targeted tests verify constant preservation, the analytic-reference calculation, and preservation of the archived result.

## Interpretation and next decision

Use the measured reconstruction error against the analytic initial field and agreement across probe grids/orders to assess whether the original failed metric meaningfully constrains solver error. Calibration on the initial geometry cannot rigorously bound error on the evolved shapes.

This audit is exploratory. It does not retroactively turn the original failed test into a pass. Before accepting a revised diagnostic for production validation, specify an independent finer-grid confirmation (for example 112³ with the same domain, interface width, duration, and time step), keep the 1% physical-field tolerance, and declare reconstruction-calibration and probe-sensitivity checks before seeing the new evolved state. Interface-width, cleavage, and long-time pattern convergence remain unresolved.
