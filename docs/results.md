# Historical results: instantaneous cleavage

These results describe the previous instantaneous, isotropic cleavage implementation. They are preserved for comparison and do not characterize the current default model. See [the progressive cytokinesis results](cytokinesis.md).

These are exploratory numerical results, not biological validation or evidence of a developmental axis.

## Default demonstration

Configuration: seed 7, $40^3$ grid, interface width 0.085, $\Delta t=0.015$, dimensionless duration 15, 16-cell cap, full feedback. Output: `outputs/demo`.

| Quantity | Result |
|---|---:|
| Final cells | 16 |
| A-like / B-like / uncommitted | 10 / 6 / 0 |
| Final total volume error | −1.261% |
| Maximum sampled absolute total volume error | 1.730% |
| Final maximum individual volume error | 1.512% |
| Final principal axis ratio | 1.25238 |
| Final asphericity | 0.01608 |
| Maximum sampled boundary occupancy | $9.86 \times 10^{-7}$ |
| Minimum sampled equivalent cell radius | 3.969 grid spacings |

The sample has both regulatory tendencies and changing geometry. It is coarse: the smallest cells span only about eight grid spacings in diameter. A/B thresholds do not prove stable fate commitment. Volume is exactly partitioned at division but softly constrained during subsequent evolution.

## Deterministic four-cell refinement

Command: `OPENBLAS_NUM_THREADS=1 python -m embryo.convergence --output outputs/convergence.json`.

To avoid confusing stochastic differences with integration error, differentiation is disabled, cycles have no jitter, division interval is 1.5, and duration is 6. Interface width remains fixed at 0.085 as the grid changes.

| Grid | $\Delta t$ | Final principal axis ratio | Total volume error |
|---|---:|---:|---:|
| $32^3$ | 0.015 | 1.20003 | −0.835% |
| $40^3$ | 0.015 | 1.20289 | −0.838% |
| $48^3$ | 0.015 | 1.20434 | −0.840% |
| $40^3$ | 0.0075 | 1.20286 | −0.838% |

The axis ratio changes about 0.12% between $40^3$ and $48^3$ and about 0.003% when halving $\Delta t$ on $40^3$. This supports the basic four-cell mechanics calculation over this interval. It does not establish convergence of the smaller 16-cell shapes, contact geometry, fate dynamics, or a sharp-interface limit.

## Software and visualization checks

Fourteen automated tests pass, covering both isolated fate-switch attractors, cleavage conservation, contact symmetry, single-cell isotropy, exact checkpoint continuation, matched division sampling in controls, a basic time refinement, and invalid configurations. The demonstration's JSON trajectory and embedded HTML data were checked; JavaScript syntax was checked with Node. The summary figure was visually inspected. All 16 final cells have one connected interior at the $\phi=0.5$ threshold (a final-frame check, not a trajectory-wide topology guarantee).

In-browser interaction testing was unavailable because the browser URL policy blocked the local file. The HTML viewer is self-contained and can be opened by the user locally; playback and camera controls have not been verified in-browser in this environment.

## Interpretation and next experiment

For seed 7, the mechanics-only and one-way models both have final axis ratio 1.242, versus 1.252 for full feedback. Thus cleavage and packing already account for most of the observed anisotropy. Do not describe this result as established fate-driven axis formation.

The next mechanism to implement is a dynamic polarity field and directional cortical mechanics, followed by persistence, rotation, and multi-seed controls. First extend spatial refinement to the full 16-cell model and measure cell connectivity and contact geometry.

## Three-seed paired control screen

Configuration: $40^3$ grid, $\Delta t=0.015$, duration 15, 16-cell cap. Each seed shares division timing and sampled directions across modes. The one-way mode senses geometry but does not alter mechanics from fate.

| Seed | Mode | A / B / uncommitted | Axis ratio | Final volume error |
|---|---|---|---:|---:|
| 7 | mechanics | 0 / 0 / 16 | 1.24152 | -1.209% |
| 7 | one_way | 10 / 6 / 0 | 1.24152 | -1.209% |
| 7 | two_way | 10 / 6 / 0 | 1.25238 | -1.261% |
| 11 | mechanics | 0 / 0 / 16 | 1.12575 | -1.201% |
| 11 | one_way | 12 / 4 / 0 | 1.12575 | -1.201% |
| 11 | two_way | 12 / 4 / 0 | 1.13200 | -1.297% |
| 23 | mechanics | 0 / 0 / 16 | 1.17507 | -1.167% |
| 23 | one_way | 13 / 3 / 0 | 1.17507 | -1.167% |
| 23 | two_way | 13 / 3 / 0 | 1.17781 | -1.288% |

Both fate tendencies appear in all three differentiated seeds. Full feedback changes the final axis ratio only slightly relative to matched mechanics controls. Three seeds are a smoke test of variability, not evidence of robustness or statistical significance. Raw tables and histories are in `outputs/screen`.
