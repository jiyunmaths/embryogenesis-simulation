# Separate space and time refinement of the coupled model

The larger domain adds clearance but retains voxel spacing 0.08. The mature developmental runs still have minimum cell radii near 3.97 grid spacings. The next numerical test refines the cells at fixed physical domain and interface width, while checking time discretization independently.

This is a **manufactured 16-cell test**, not a new embryogenesis trajectory. Sixteen analytic diffuse spheres of radius 0.30 occupy a prescribed 4×2×2 arrangement, with centers separated by 0.58. Every grid samples the same logistic sphere profiles directly. The target volume is calculated by radial quadrature of the continuous occupancy profile, rather than inherited from any grid. No coarse checkpoint is interpolated into a finer initial state.

Activator starts from the same weak analytic perturbation evaluated at the prescribed cell centers. Inhibitor starts at one; fate and polarity start at zero. Conservative signaling, fate response, polarity, and mechanical feedback all evolve. Cell division is disabled by starting at the 16-cell cap. This separates local coupled integration from cleavage-time and lineage divergence.

## Protocol

Keep half-width 2.24, interface width 0.085, diffusivities 0.02 and 0.4, and all other constitutive parameters fixed. Evolve for 0.6 model-time units.

| Sweep | Grid | Time step |
|---|---:|---:|
| Spatial, coarse | 56³ | 0.00375 |
| Spatial, intermediate | 72³ | 0.00375 |
| Spatial, fine | 88³ | 0.00375 |
| Temporal, coarse | 72³ | 0.015 |
| Temporal, intermediate | 72³ | 0.0075 |
| Temporal, fine | 72³ | 0.00375 |

The shared 72³ / 0.00375 case runs only once, so there are **five simulations**. The spatial sweep uses the smallest tested step throughout. A temporal check at 72³ is useful evidence of time-error sensitivity, but does not establish that time error is negligible at every spatial resolution.

```bash
python -m embryo.resolution prepare --output outputs/coupled-resolution-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.resolution run \
  --output outputs/coupled-resolution-repeat \
  --after outputs/domain-conservative-extended
```

Preparation immediately evaluates the frozen geometry on all three grids. The `run --after` process waits for the specified domain study to finish. It starts the five cases sequentially only if the largest-domain boundary screen and largest-pair shape screen pass. If those screens fail, its status becomes `blocked`; an upstream execution failure produces `failed`. No physical integration is started while the prerequisite study is running. Omit `--after` for a deliberately independent run.

The completed study is `outputs/coupled-resolution`. Its `status.json` records the current phase and completed cases. The final domain comparison hash is retained when the gate opens. Execution failures are recorded rather than reported as completed science.

## Measurements and declared screens

Each case records metrics every solver step, a final checkpoint, and initial/final surface snapshots. The viewer contains **two endpoint snapshots**, not a continuous animation of this short test.

Compare each adjacent pair using cell occupancy fields, concentrations, axis ratio, and individual volumes. Temporal comparisons use the same voxel coordinates directly. Spatial comparisons use trilinear diagnostic interpolation onto the fixed 56³ probe grid; this projection is never used by the solver. Consequently the spatial L2 measurement includes projection error and can obscure smaller-scale discrepancies. Native-grid volume and shape diagnostics supplement it.

For each sweep, require that the finest pair differs by less than 1% in occupancy-field L2, activator/inhibitor L2, and axis ratio, and that its field difference is smaller than the coarser pair's. Also require the finest spatial case to resolve minimum equivalent radius with at least four grid spacings, every case to remain below 0.01 boundary occupancy and 5% cell-volume error, and no phase-field clipping. All thresholds are stored in the protocol before evolution. Failure remains a result; parameters are not silently adjusted to pass.

These are finite-difference agreement screens, not a proof of convergence. Unequal mesh ratios and possible cancellation preclude a naive log2 order estimate. The full coupled solver uses first-order splitting; the signaling integrator alone being RK2 does not imply second-order coupled time accuracy.

## Completed frozen-geometry preflight

The same analytic cells yield:

| Grid | Voxel spacing | Minimum equivalent radius / spacing | Maximum relative volume-quadrature error | Unstable modes |
|---:|---:|---:|---:|---:|
| 56³ | 0.080000 | 3.9349 | 2.604e-6 | 14 |
| 72³ | 0.062222 | 5.0591 | 3.744e-8 | 14 |
| 88³ | 0.050909 | 6.1834 | 5.214e-9 | 14 |

The conductance-matrix relative difference falls from **7.100e-5** (56→72) to **2.010e-6** (72→88). The maximum positive-eigenvalue relative difference falls from **7.165e-5** to **1.927e-6**, or approximately **0.000193%** for the finer pair. These observations concern quadrature and spectra of fixed phase fields. They do not establish continuum accuracy of the overlap-to-area closure for curved contacts, or convergence of the evolving geometry. The manufactured geometry's 14 unstable modes are not a prediction for the developing embryo.

## What remains after this screen

First read the actual coupled space/time comparison and retain any failed criteria. Then test longer windows, division and lineage sensitivity, geometric rotations, contact cutoff, and interface-width variation independently. A clean short manufactured-state result cannot establish convergence of long-time spontaneous pattern selection or biological cleavage. The imposed shape and chemical perturbation are controls, not evidence of emergent order.

## Completed coupled screen and targeted follow-up

All five cases completed. Eleven of twelve declared checks pass. The spatial occupancy-field check **fails**: the 72³/88³ pair differs by **1.4306%**, exceeding the predeclared 1% criterion, although the 56³/72³ discrepancy was larger (3.7284%). Finest-pair temporal field discrepancy is **0.01320%**, and signal/shape, volume, boundary, resolution, and clipping screens pass.

The failure is retained in the original report. The next test is a [projection-error audit](projection_audit.md): even identical analytic initial geometry produces a nonzero inter-grid discrepancy after diagnostic interpolation. Calibrating that measurement is necessary before interpreting the final discrepancy as solver error or choosing another expensive refinement. The audit does not alter the solver, archived states, threshold, or original pass/fail result.

The projection audit is now complete. The original linear diagnostic gives 1.4943% discrepancy even at the identical analytic initial state; cubic/quintic final comparisons across three probe grids range from 0.2366% to 0.3184%. This identifies substantial reconstruction contamination, but does not retroactively change the archived failure. An independently specified finer-grid confirmation is still needed before adopting a revised validation result.
