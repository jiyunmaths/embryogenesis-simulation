# Feedback-off repeated development: spatial refinement protocol

This prospective study follows the [independently reviewed timestep ladder](development_time_diagnosis_assessment.md). Its finest pair clears the feedback-off time screen at grid 72, with a narrow spectral margin. Spatial accuracy and the other mechanical controls remain unqualified. This experiment tests spatial refinement at the smaller timestep and temporal accuracy on the finest new grid; it does not add developmental histories or change the model.

## Fixed design

Follow the same analytic seed-7 zygote through naturally scheduled progressive cleavage to sixteen cells at t=18. Retain the domain half-width 2.24, interface width 0.085, radius 0.8, common continuous initial volume target, chemistry, conservative transport, polarity dynamics, cell cycles, pressure correction, daughter partitioning and random streams. Mechanical feedback remains disabled. No labels or chemical patches are supplied.

| Case | Grid | Timestep | Steps | Role |
|---|---:|---:|---:|---|
| off72 | 72³ | 0.00046875 | 38,400 | Exact reuse of completed off_quarter; no rerun or copied checkpoints |
| off88 | 88³ | 0.00046875 | 38,400 | New spatial comparison |
| off104 | 104³ | 0.00046875 | 38,400 | New spatial comparison |
| off120 | 120³ | 0.00046875 | 38,400 | New finest spatial comparison |
| off120_half | 120³ | 0.000234375 | 76,800 | Finest-grid temporal check |

Four new serial GPU paths use resident PyTorch arrays/matrices and the unchanged custom CUDA mechanics/geometry kernels on the designated GTX 1080 Ti UUID. The grid-72 case is exposed by a read-only reuse policy through a directory symlink to the preserved completed study; the runner refuses to execute that case. Its original report, protocol and every saved evidence hash are retained. All other cases begin from independently sampled versions of the same analytic physical state, rather than interpolating a developed trajectory. Initial target, chemistry, numerical corrections, lineage, cycle schedule and RNG states must match apart from grid, timestep and related strides. This is **zero new independent scientific histories**.

## Measurements and unchanged limits

Retain original observations every 0.15 (121 per path), native checkpoints at t=0,3,6,9,12,15,18, and dense read-only network/furrow records every 0.015 plus event/active-parent changes. Record every dividing-mother projection. Use the previously checked instrumentation unchanged; its three-window bitwise gate and full-path noninterference evidence are pinned. Retain eight exact restart steps at each 2→4, 4→8 and 8→16 live window, the independent initial eight-step chemical/dilution check, every-step physical/event checks, sixteen final cells, fifteen events, no active furrows and at least three mature time units.

Compare all three spatial adjacent pairs at the same physical times. Use the inherited primary field reconstruction at 112³, cubic order, at all seven checkpoints. Retain the original endpoint variants (112³/128³, orders 3/5) and maximum endpoint error-spread limit 0.001. These are field-comparison grids, not mechanics grids. Before production, calibrate all four native grids on the original analytic sphere, translated small sphere and rotated ellipsoid; all 48 combinations must remain below the inherited 0.0025 relative L2 calibration limit. No post hoc event alignment or spatial registration is allowed.

Spatial limits remain: chemical distribution and log-distribution distances 0.01; contrast 0.005; axis ratio 0.01; total volume and whole species amounts 0.005; mean polarity 0.01; spectrum 0.05; maximum occupancy/chemical field relative L2 and occupancy-normalized polarity field L2 0.01; cell-count difference at most one; at most thirty unmatched original observations; maturity-time difference 0.6. Report every coarser failure. Check the inherited spatial error-decay rule across the three pairs for concentration/log distributions, axis ratio and all four fields, with the original 1e-10 floor; this is descriptive decay, not a proved order.

Compare off120 with off120_half using all inherited time limits and **native grid-120 fields**: spectrum and field errors 0.01, maturity-time difference 0.05, and the remaining original temporal limits. In addition, require the same-count dense temporal spectrum discrepancy to stay below the same 0.01 limit. Report unequal-count dense observations and edge-support differences separately; do not invent spectra for unequal-dimensional graphs. Dense spatial spectra/conductances are additional diagnostics. This extra temporal sampling strengthens the screen without loosening any original limit.

Every inherited physical/event bound remains unchanged: maximum volume error below 5%, minimum equivalent radius four grid spacings, zero mechanical clipping, positive finite chemistry, boundary occupancy below 0.01, connected nondividing cells, resolved neck/overlap and connected daughters, and the original event target/amount/volume/occupancy and dilution conservation tolerances.

## Decisions and safeguards

The scoped spatial/time screen can clear only if every path passes its physical/event checks, the **104→120 spatial pair** passes every spatial limit and endpoint spread check, the original spatial error-decay rule passes, and the **120-grid timestep pair** passes all temporal gates including the additional dense spectrum check. Passing a finer pair does not accept a failed coarser pair. All results remain provisional pending independent saved-data review. This does not qualify feedback/polarity controls, full-horizon live-event restart replay, a full finer CPU ladder, identity inheritance, a many-cell continuum, or fresh developmental histories.

If spatial fields still fail, preserve those results and use their fixed-time error peaks to specify the next refinement or isolated diagnostic. If the finest-grid time check fails, do not interpret the spatial ladder as a time-qualified spatial result. No model/threshold changes or automatic extension to another grid is authorized by this runner.

Pin the earlier reviewed time screen, failed repeated-development study, accepted scientific evidence, ledger/manuscript, solver/kernel sources and runtime binary. Snapshot reporting documents separately. Preparation and launch refuse existing output/result paths. Check source/input and reuse evidence before production and assessment; use a nonblocking coordinator lock. No other GPU or user's processes are touched. Four new paths run serially to avoid GPU memory contention; no full CPU references are launched.

Prepare with `python -m embryo.development_space_refinement prepare`, calibrate with `calibrate`, then launch `run` in a fresh directory. The launch record will include preparation, calibration, software checks, independent initial-state/provenance audit, pinned protocol and GPU process details. Runtime estimates are provisional until the first new path supplies measured throughput.
