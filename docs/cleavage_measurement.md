# Correcting the cleavage comparison measurement

The original cleavage screen passed 24 of 25 checks. Its failed check was the **initial diagnostic reconstruction**, not conservation, daughter connectivity, event timing, or final shape. Cubic interpolation of occupancy from 56³ samples gave 0.477% relative L2 error against the exact analytic sphere, exceeding the unchanged 0.25% tolerance.

## Reconstruction change

The solver's native field is the phase field, phi. Occupancy is a nonlinear derived quantity, h(phi)=phi²(3−2phi). The original diagnostic interpolates samples of h(phi). The revised diagnostic first interpolates samples of phi and then evaluates h at the common probe points:

$$
\widehat h(\mathbf x)=h\!\left(P_h\phi_h(\mathbf x)\right).
$$

Interpolation and nonlinear transformation do not commute. For these diffuse profiles, reconstructing the smoother native phase field reduces the measured error. On the original 56³ sphere and 88³ probe, cubic reconstruction lowers the error from 0.477% to 0.141%; quintic lowers it to 0.058%. This is a change in diagnostic reconstruction, not in the simulated dynamics or the definition of occupancy.

No solver field is modified. No clipping conceals spline undershoot or overshoot. Native-grid comparisons use native samples. Phase-field ranges in the analytic calibration and occupancy ranges in the evolved comparison are exported. Both interpolation orders use nearest-value boundary extension and spline prefiltering; the cells are far from the boundary.

## Separate validation

After the exploratory original-sphere comparison, a separate protocol fixes cubic reconstruction as the primary diagnostic, quintic as a cross-check, and probe grids 88³, 112³, and 128³. Before evaluating these comparisons, it specifies:

- Every analytic reconstruction must be below the original **0.25%** tolerance, on every source grid and probe, for both orders.
- Analytic geometries include the original sphere, a translated radius-0.36 sphere, and a translated/rotated radius-0.44 ellipsoid with axis factors 0.9, 1.05, and 1.1. The latter two are holdouts, not the geometry used to select the method.
- All finest-pair evolved daughter-field comparisons must remain below the original **1%** tolerance and improve over the coarser pair.
- The spread among interpolation-order/probe results must be below **0.1 percentage point** for each division orientation.
- All 24 other original checks must remain passing, and archived inputs must remain byte-for-byte unchanged.

The analytic holdouts test reconstruction on known continuous shapes. They do not constitute independent evolved simulations or rigorously bound interpolation error for arbitrary evolved cells. Cross-method agreement is supporting evidence, not an exact solution. No scalar initial error is subtracted from an evolved-field discrepancy.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cleavage_measurement prepare \
  --source outputs/cleavage-resolution --output outputs/cleavage-measurement-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cleavage_measurement run \
  --output outputs/cleavage-measurement-repeat
```

The separate validation is in `outputs/cleavage-measurement`; its protocol records SHA-256 hashes of the original protocol, comparison report, and six spatial checkpoints. The original `outputs/cleavage-resolution/comparison.json` and its failed decision are never overwritten. The running developmental-refinement model and its code hashes are unaffected by this new diagnostic module.

## Analytic calibration results

| Geometry | Worst cubic error | Worst quintic error |
|---|---:|---:|
| original sphere | 0.1407% | 0.0575% |
| translated smaller holdout | 0.1503% | 0.0654% |
| rotated ellipsoid holdout | 0.1442% | 0.0694% |

The largest error across all geometries, source grids, and probes is 0.1503% for cubic reconstruction and 0.0694% for quintic. Every calibration is below 0.25%. The smallest reconstructed phase value is approximately −6.35e-7, recorded without clipping. This interpolation is diagnostic only and does not claim to conserve integrated volume.

## Completed revised validation

All **10 revised measurement checks pass**. Across both orders and all three probes, finest-pair evolved occupancy differences are **0.5410–0.5577%** for axial cleavage and **0.3854–0.4051%** for oblique cleavage. Every variant stays below 1% and improves on the coarser pair. Diagnostic spreads are below 0.1 percentage point. The other 24 original checks still pass, and all archived-input hashes are unchanged.

The measurement issue is resolved for this controlled cleavage screen by a separately calibrated native-phase diagnostic. The original failed report remains an accurate historical record. This does not establish full developmental or sharp-interface convergence. Eight regression tests pass, including reproduction of the original failure and corrected calibration.
