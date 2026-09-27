# Controlled small-cell cleavage refinement

The dynamic cutoff experiment passed its sensitivity checks but retained the coarse-cell-resolution failure. The next test targets a missing part of numerical validation: progressive cleavage and abscission. The previous successful 112³ confirmation contained sixteen prescribed cells and no division.

## Controlled geometry and protocol

A single spherical cell with analytic radius **0.4** is sampled independently on grids **56³, 72³, and 88³**, all in the same domain of half-width **2.24**, with interface width **0.085**. Its continuous diffuse occupancy volume is integrated by radial quadrature and used as the common volume target. No coarse checkpoint is interpolated onto a finer grid.

Division begins at time zero with either axis (1, 0, 0) or normalized (1, 2, 3). The latter tests an oblique cut through the Cartesian lattice. The progressive ring dynamics, neck and overlap criteria, and volume projection are unchanged. The cell cap is two, so there is exactly one potential division. Standard coupling remains enabled, but partition noise is zero; two cells remain below the four-cell competence threshold for fate differentiation. This is a controlled mechanical/division test, not a symmetry-breaking or identity-differentiation experiment.

Both orientations use time step **0.00375** for the spatial comparison. Separate oblique 72³ runs use **0.015** and **0.0075**, reusing the 72³/0.00375 case as the finest time step. Eight cases run through **t = 3**. Failure to finish division by then is a reported failure; no unresolved neck is forcibly cut.

## Measurements and declared criteria

The protocol is saved before running. Spatial comparisons use cubic reconstruction on 88³ and 112³ probes. Initial analytic reconstruction must agree within 0.25%. Both probe grids must show a finest-pair daughter occupancy discrepancy below 1%, decreasing from the coarser pair. Native samples are used for temporal comparisons, which must also be below 1% and decrease with refinement. Final comparisons require matching daughter IDs and two completed daughters.

Finest spatial and temporal pairs must differ in abscission time by less than 0.05 model-time units. Finest spatial pairs must agree in final axis ratio within 1%. On the finest grid, axial versus oblique cuts must agree in abscission time within 0.1 and in axis ratio within 1%. No direct rotated voxel-field comparison is made; these are rotation-invariant scalar checks.

An experiment-only subclass records volume, activator amount, inhibitor amount, and aggregate occupancy **immediately before and after the actual abscission method**, separating relabeling conservation from reaction and mechanical changes elsewhere in the time step. Relative jumps must be below 1e-6. Each daughter must have exactly one connected component at the phi=0.5 threshold immediately after abscission and at subsequent 0.15-time observations, using six-neighbor voxel connectivity. A narrowing mother's disconnected lobes are not classified as daughter fragmentation.

Every-step cell-volume error must stay below 5%, with no clipping. Boundary occupancy must remain below 1% at observations. The finest-grid minimum equivalent cell radius must remain at least four grid spacings. All gates remain explicit, including failed ones; these are practical numerical tolerances, not biological calibration.

## Reproduction and outputs

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cleavage_resolution prepare \
  --output outputs/cleavage-resolution-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cleavage_resolution run \
  --output outputs/cleavage-resolution-repeat
```

The production study is in `outputs/cleavage-resolution`. `protocol.json` records configurations and the model hash. `status.json` records the active case and completed cases; each case writes progress to `history.json` and the run log. Final outputs include an abscission audit, full sampled metrics, a final checkpoint, and an explicitly two-endpoint surface viewer. `comparison.json` and `RESULTS.md` are generated automatically after all cases complete.

This test does not validate the entire zygote-to-sixteen-cell trajectory. It isolates one division at a smaller cell scale, with imposed orientation and fixed interface width. Even a pass would leave repeated cleavage, spontaneously chosen axes, long-time refinement, and seed replication unresolved. A failure should identify which feature needs attention before launching an expensive refined developmental trajectory.

## Completed outcome

All eight cases finished; **24/25 checks passed**. The only failure was initial reconstruction calibration: the 56³ initial field differs by about 0.477% on the common probes, above the 0.25% criterion. This failure remains recorded. At the user’s request, [full developmental refinement](development_refinement.md) is proceeding as a separate observable-based study, without treating this calibration failure as resolved.

The subsequent [native-phase measurement validation](cleavage_measurement.md) resolves the diagnostic issue separately: all 10 new checks pass, including translated/deformed holdouts under the unchanged 0.25% calibration tolerance and evolved-field comparisons under 1%. The original 24/25 report remains unchanged.
