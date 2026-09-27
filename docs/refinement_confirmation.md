# Independent finer-grid confirmation

This follow-up adds a genuinely new 112³ evolved state to the archived 72³ and 88³ cases. Its protocol and source hashes were written before starting the new evolution. The original coupled-resolution result and exploratory projection audit remain unchanged.

## Fixed experiment

The manufactured 16-cell aggregate, domain half-width 2.24, interface width 0.085, conservative transport, reaction and mechanical parameters, time step 0.00375, and duration 0.6 are inherited exactly from the archived spatial cases. Only the voxel grid changes. The new initial state is sampled directly from the analytic diffuse spheres; it is not interpolated from a coarse checkpoint. This is a local numerical confirmation with imposed geometry, not a spontaneous developmental experiment.

## Prospective acceptance criteria

Cubic and quintic spline reconstruction are each evaluated on probe grids 88³, 112³, and 128³. Native-grid samples are used directly. Reconstruction never modifies the solver state, and negative spline values are not clipped. All six diagnostic variants must satisfy the following checks:

- The final 88³/112³ occupancy relative L2 difference is below **1%**.
- That difference is smaller than the 72³/88³ difference measured with the same diagnostic.
- Every initial reconstruction of all three cases differs from the analytic reference by less than **0.25%**.
- The spread between the largest and smallest of the six final 88³/112³ differences is below **0.1 percentage point**.

The last two tolerances were specified after the exploratory audit but before seeing the new evolved state. They allocate part of the 1% field-error budget to diagnostic uncertainty; they are practical screening criteria, not rigorous bounds on evolved-field error.

The finest pair must also agree within 1% for activator, inhibitor, and axis ratio. All three cases must stay below 1% boundary occupancy and 5% cell-volume error, with no phase-field clipping. The new grid must keep the smallest cell radius at least four grid spacings. Archived-source hashes must remain unchanged. Every condition must pass; there is no averaging away a failed reconstruction variant.

## Reproduction and status

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.refinement_confirmation prepare \
  --source outputs/coupled-resolution --output outputs/refinement-confirmation-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.refinement_confirmation run \
  --output outputs/refinement-confirmation-repeat
```

The completed study is in `outputs/refinement-confirmation`. `protocol.json` is immutable during execution; `status.json` records preparation, evolution, comparison, completion, or failure. The final report is `RESULTS.md`, detailed measurements are `comparison.json`, and the new case contains its checkpoint, every-step metrics, and two-endpoint surface viewer. Source hashes cover the old protocol, comparison, checkpoints, and analysis files.

A pass would support short-time coupled spatial consistency for this manufactured state at fixed interface width. It would not establish an asymptotic convergence order, remove temporal error at the fixed time step, or validate cleavage, long-time pattern selection, interface-width sensitivity, biological mechanisms, or the overlap-to-contact-area closure.

## Completed outcome

All eleven declared checks passed. The six 88³/112³ final field comparisons range from **0.1619% to 0.2068%**, and each decreases from its corresponding 72³/88³ comparison. Archived inputs are unchanged. The original linear-projection screen remains failed; this is a separately specified confirmation. The next [contact-cutoff sensitivity screen](contact_sensitivity.md) is also complete.
