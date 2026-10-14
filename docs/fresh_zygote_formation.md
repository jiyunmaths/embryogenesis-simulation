# Inherited fresh-history chemistry: 60-unit formation pilot

Registered 12 October 2026 before execution, following the [accepted mature-clock screen](fresh_zygote_mature_clock_fix_assessment.md). Test whether the small chemical differences produced by actual zygote development become sustained patterns while cells and their chemical contact network continue to move.

## Histories, starts and controls

Use the twelve actual t=18 developmental preparations from histories **1101, 1102 and 1103**, with four matched branches: feedback off, tension/adhesion (0.25/0.35), polarity tension (0.35), and the full combination. Keep baseline mechanics, polarity evolution, conservative transport, dilution, beta=2, D_a=0.02, D_b=0.55, grid120, extent2.24 and all other coefficients unchanged. There are **three existing independent histories**, not 24 independent biological replicates.

Run each branch at **dt=0.00375 and 0.001875**, giving **24 new numerical paths and twelve full-horizon comparisons**. Start again from the immutable, independently audited t=18 checkpoint of the same clock used in the short screen, and run the complete physical interval **t=18–78 (60 elapsed units)**. This repeats the short window as part of a full-horizon trajectory; none of the short trajectories is silently relabeled as a new history.

Preserve all initial array bytes: chemistry and its correction, phase fields and carry, rounding diagnostic, polarity, targets, IDs/parents, cycle due times, lineage/event audits and random states. Only the declared Config step horizon changes from t=24 to78; physical time, dt, observation stride and step_number do not change. There is no near-uniform chemistry reset, imposed pattern, feedback switch, backend handoff or further division. Mature chemistry and geometry already differ across branches because their mechanical contexts acted during development; these comparisons test entire developmental contexts, not interventions applied to a common mature shape.

## Solver and numerical gates

Use the unchanged compensated resident DevelopmentGpu solver: PyTorch GPU arrays/matrices, custom CUDA mechanics/geometry/polarity, float32 visible fields with float64 phase residual, and float64 chemical state/correction. Sources, inputs, accepted review and runtime binary are pinned. Each retimed actual start receives eight independent NumPy-reference/GPU steps and an exact correction-aware GPU replay at the inherited prefix limits.

Check every-step volume error <5%, equivalent radius ≥4 voxels, zero clipping, positive finite chemistry, dilution amount error ≤2e-14 and inherited total target volume. Keep sixteen cells, fifteen completed cleavage events, no active furrow and unchanged ancestry, cycle and RNG metadata. At regular observations require boundary occupancy <0.01 and one connected component for every cell.

Observe every **0.15** units (401 observations per path), record dense networks every **0.015** (4,001 graphs per path), and save complete checkpoints at elapsed **0,12,24,36,48,60**. Require exact eight-step restart replays at elapsed **12 and36**. All checkpoint arrays retain compensation. Native field comparisons use the same occupancy, activator/inhibitor-density and polarity-density definitions as the accepted short screen.

Retain all inherited regular distribution/shape/volume/amount/polarity and dense spectrum gates. Native density error and regular/dense relative spectral error must be ≤1%. Homologous log chemical error ≤0.01, absolute polarity-component error ≤0.01, relative cell-volume error ≤0.005 and relative transport-matrix norm error ≤0.01. Require exact cell/ancestry correspondence. Report every between-clock contact-support mismatch without hiding graph transitions.

Also require the maximum difference in frozen uniform nonzero-mode growth rate ≤0.001. These growth rates hold the current conductances/volumes fixed and linearize uniform chemistry; they are diagnostics, not stability of the whole moving system.

## Prespecified formation measurements

Contrast is across-cell SD of **natural-log activator**, not a cluster label or biological identity. The primary threshold is **0.1**, with **0.05 and0.2** sensitivity thresholds. These inherit the earlier mature formation convention and are screens for chemical contrast, not lineage-marker measurements.

For each threshold report:

1. First sampled threshold crossing, without shifting/alignment of trajectories.
2. First sustained interval: the first observation starting an interval of at least **24 model-time units** during which every saved contrast remains strictly above threshold. Report its start time only once the entire interval has been observed. A crossing too late to establish 24 units is recorded as a crossing, not sustained formation.
3. Late persistence: contrast remains strictly above threshold throughout elapsed **36–60**. Record late minimum, final contrast and all sampled losses/recoveries separately.

Both clocks must agree on sustained-formation and late-persistence classification at all three thresholds. When sustained onset is observed at both clocks, its sampled start times must differ by ≤0.3 units. An outcome that remains below threshold passes scientifically if numerical checks pass; failed/nonorganizing histories are not replaced. Absence of onset by this horizon is finite-window nonformation, not proof that formation is impossible.

Measure cell chemistry/polarity/volumes, tissue axis ratio, summed undirected conductance, contact-support changes, graph spectrum and supported growing-mode counts. Report contrasts and conductance changes relative to each branch's own inherited start. The feedback-off branch retains passive mechanics, so contact remodeling in it cannot be attributed to activity-dependent material feedback. Matching a branch effect does not identify complete mediation; fixed-conductance and chemical/polarity cue interventions remain later controls if needed.

## Execution, stopping and interpretation

Execute serially on the dedicated GTX1080Ti with one CPU numerical thread; compare each coarse/fine pair before advancing. The six-unit screen suggests roughly **15–18 hours** for the full 24-path design, conditional on numerical gates and long-window speed; update this estimate from the first completed pair. Sparse full checkpoints limit storage while dense matrix diagnostics remain available.

Stop on any prefix, physical, exact restart or full-horizon numerical comparison failure. Preserve failed states and unexecuted cases. Do not automatically lower dt, relax thresholds, replace histories, relaunch a recorded case or continue to240 units. Numerical failures require a separate registered diagnosis. Completed execution remains provisional until an independent saved-array review.

If this pilot passes its own numerical gates, the next separately registered stage tests maintenance through **240 total elapsed units from t=18**. Frozen endpoint coexistence and response/context tests are subsequent questions. This pilot alone does not establish autonomous or inherited cell identity, commitment, calibrated embryo biology, universal convergence, a geometric continuum or chemically caused large-scale shape organization. Manuscript and ledger scientific snapshots remain unchanged by launch.
