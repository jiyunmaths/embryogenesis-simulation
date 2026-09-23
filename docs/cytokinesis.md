# Shape-oriented, progressive cytokinesis

The numerical results below predate explicit graph signaling and apical–basal polarity. The cytokinesis mechanism remains in use; see [the current model](graph_signaling.md) for subsequent changes.

The default now uses the longest cell axis for spindle orientation and retains the mother cell throughout a mechanically evolving furrow. The cleavage plane is perpendicular to the spindle. Nearly equal longest axes are handled by sampling within their eigenspace; an isotropic rule remains available as a control.

A prescribed contracting equatorial annulus drives ingression over multiple steps. An interface-local volume correction holds the mother's measured onset volume. Abscission requires a sufficiently thin neck and low prospective daughter overlap; it does not occur merely because the contraction timer expires. The final occupancy partition is conservative, and daughter target volumes reflect actual lobe fractions to avoid imposing new pressure differences.

Read [the equations](model.md#cleavage-and-progressive-cytokinesis) for the precise force, volume constraint, and thresholds. This is an exploratory ring surrogate, not a resolved actomyosin/fluid model or a tensile-stress-based spindle rule.

## Run and inspect

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo --output outputs/new-cytokinesis-run
OPENBLAS_NUM_THREADS=1 python -m embryo.cytokinesis_check --output outputs/new-cleavage-check
```

The already generated 16-cell example is `outputs/progressive/viewer.html`, with a static overview at `outputs/progressive/summary.png`. Playback shows the mother deforming before two IDs appear; the sidebar reports how many cells are undergoing cytokinesis. Lineage records include onset, abscission, spindle axis, and neck/overlap values.

For the old/new comparison recorded here, the previous source was archived as an output artifact before editing. Reproduce the comparison locally with:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.cytokinesis_check \
  --baseline outputs/cytokinesis-baseline/model.py \
  --output outputs/new-old-versus-new-check
```

The optional baseline argument executes that local Python module. It is unnecessary for checking the current solver alone. Existing output directories are not overwritten.

## Matched single-cleavage experiment

One isolated zygote, identical explicit spindle direction $(1,2,3)/\sqrt{14}$, differentiation and fate feedback disabled, duration $3$. Metrics are recorded every integration step, including before and after initiation. The cell cap is two. Thus this comparison isolates cytokinesis; it does not compare orientation rules.

| Model | Grid | $\Delta t$ | Abscission time | Largest cell volume error | Largest contact-proxy increment |
|---|---|---:|---:|---:|---:|
| instantaneous | $40^3$ | 0.015 | 0.0000 | 1.052% | 0.028809 |
| progressive | $40^3$ | 0.015 | 1.0950 | 0.626% | 0.004943 |
| progressive_half_dt | $40^3$ | 0.0075 | 1.0425 | 0.622% | 0.005162 |
| progressive_grid48 | $48^3$ | 0.015 | 1.1250 | 0.618% | 0.004768 |

The largest individual volume error decreased by approximately 41%, and the largest contact-proxy increment by approximately 83%, in the matched default-grid comparison. The old model introduces its largest contact jump at initiation; the new model leaves contacts unchanged at initiation, with a smaller graph change at abscission. The proxy is the existing diffuse interface-shell overlap, **not a physical contact area or measured pressure**.

This small benchmark did not exhibit severe collapse in the old solver. It establishes reduced numerical transients for this case, not universal stability or biological correctness. Volume is held to numerical tolerance while the mother constricts; the remaining reported error develops during daughter relaxation under the original soft volume penalty.

Raw time series and plots: `outputs/cytokinesis-check/comparison.json` and `comparison.png`.

## Refinement limits

Halving $\Delta t$ changes the abscission time from $1.095$ to $1.0425$ (about 4.8%). Increasing the grid from $40^3$ to $48^3$ changes it to $1.125$ (about 2.7%). Final axis ratios change by about 0.4% in each comparison. These are limited refinement checks, not a convergence certification. Abscission uses discrete neck and overlap thresholds, and remains sensitive to resolution and time step. Thin interfaces and the smallest cells require further refinement before interpreting cytokinesis timing quantitatively.

## Full default run

Seed 7, $40^3$ grid, $\Delta t=0.015$, duration 15, cap 16:

- 16 cells at the end; 13 A-like, 2 B-like, and 1 uncommitted activity.
- Maximum sampled aggregate volume error: 1.389%; final aggregate error: −1.221%.
- Final maximum individual volume error: 1.379%.
- Final principal axis ratio: 1.2705.
- Every final cell has one connected interior at $\phi=0.5$.
- No overdue divisions in recorded frames; all divisions completed at the final time.

The embryo still has coarse small-cell resolution (minimum sampled equivalent radius about 3.96 grid spacings). These results do not establish robust fate commitment or a developmental axis. Original outputs in `outputs/demo` and [the earlier results](results.md) describe the old implementation and are retained as historical comparisons.

## Regression checks and compatibility

The test suite covers rotated elongated cells, unbiased orientation fallback for a spherical cell, unchanged fields/contacts at division onset, intermediate furrow states, volume conservation through abscission, post-division volume transients, connected daughter interiors, delayed unresolved cuts, cell-cap reservations, and exact continuation from a checkpoint during cytokinesis. Existing fate and mechanical checks remain.

Checkpoint schema 2 stores active division geometry, onset time and volume, and both RNG states. Older instantaneous-division checkpoints raise a clear error instead of silently continuing with new dynamics. Independent RNG streams separate regulatory noise from mechanical sampling, but geometry feedback can now change cleavage axes and completion times; same-seed feedback ablations are no longer guaranteed to have identical division schedules.
