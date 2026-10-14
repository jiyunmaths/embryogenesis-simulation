# Matched mechanical-control developmental qualification

## Question and scope

Can the current developmental method resolve the three remaining mechanical branches through naturally scheduled cleavage at the finer settings that passed the feedback-off baseline? Can a saved state during live cleavage reproduce the complete remaining development?

This is numerical qualification of **one synthetic seed-7 preparation through t=0–18**, not a new independent history or an identity-formation experiment. The [independent feedback-off assessment](development_space_refinement_assessment.md) accepts its finest spatial and temporal screens; the historical failed comparisons remain preserved. Qualification of a baseline cannot substitute for testing feedback/polarity branches.

The existing resident PyTorch arrays/matrix operations and custom CUDA mechanics/geometry/polarity kernels, pressure-corrected phase carry, compensated chemical stepping, conservative transport, dilution, longest-axis cleavage and progressive ring remain unchanged. No shared source, compiled kernel, physical model or acceptance threshold is modified.

## Matched design

All starts sample the same analytic radius-0.8 zygote, half-domain 2.24 and interface width 0.085 on their own native grids. They share the continuous target volume, activator 1.03, inhibitor 0.98, initially zero polarity, chemistry, schedules, lineage and random streams. Polarity evolves in every branch. There are no fate labels. The only physical differences are the previously tested mechanical couplings:

| Branch | Activity-dependent tension cγ | Activity-dependent adhesion cA | Polarity tension χ |
|---|---:|---:|---:|
| Feedback off, accepted exact reuse | 0 | 0 | 0 |
| Tension/adhesion only | 0.25 | 0.35 | 0 |
| Polarity tension only | 0 | 0 | 0.35 |
| Full combination | 0.25 | 0.35 | 0.35 |

The names `fate_tension` and `fate_adhesion` are historical configuration keys; coefficients act on measured activator through tanh(a−1), not prescribed identities. Directional polarity tension also depends on the evolving polarity vector. This design changes the mechanical couplings while holding chemical feedback active.

Each of the three new branches has four trajectories:

| Level | Native grid | Timestep |
|---|---:|---:|
| Spatial 88 | 88³ | 0.00046875 |
| Spatial 104 | 104³ | 0.00046875 |
| Candidate 120 | 120³ | 0.00046875 |
| Temporal verification | 120³ | 0.000234375 |

The study adds **12 new numerical trajectories, one exact reused feedback-off candidate, and zero independent scientific histories**. Reuse is through a read-only directory link and retains the original job/protocol/report. It is never executed or overwritten. The prelaunch audit independently verifies analytic phase arrays and matched nongeometric arrays, metadata and random states.

Natural division interval is 2 with the inherited jitter, progressive ring duration 0.9 and ring strength 8; maximum count is 16. Original observations remain every 0.15, regular checkpoints every 3 and dense network/furrow diagnostics every 0.015 plus events. No time registration, field clipping or ancestry alignment is added to obtain a pass.

## Prospective numerical gates

Every path must pass the unchanged physical/event requirements, finish 15 resolved divisions with 16 cells and no active furrows, and retain at least three model units after last abscission. Minimum equivalent radius is four spacings; cell-volume discrepancy must be below 5%, boundary occupancy below 0.01, and mechanical clipping zero. Chemical positivity, target/amount inheritance, constant-preserving and volume-conservative transport, daughter connectivity, every-step dilution and short exact live-window restarts remain required.

For **each branch separately**:

1. Compare 88→104 and 104→120 at the unchanged physical times, using the inherited calibrated 112³ cubic primary reconstruction and all 112/128, order-3/5 endpoint variants. The **104→120 pair must pass every original spatial gate**, including maximum field discrepancy ≤1% and spectrum discrepancy ≤5%.
2. Require the inherited decreasing-error rule across both adjacent pairs for concentration/log distribution distances, axis ratio and all four density fields, with the original 1e−10 floor. The coarser pair may fail its absolute limits and remains explicitly failed; it cannot replace the finest pair.
3. At 120³, halve the timestep and compare native fields, chemical distributions/contrast, amounts/volumes, polarity, counts, maturity timing and spectra at the original times. Every inherited temporal gate must pass; both field and spectrum limits are 1%.
4. Require the dense, equal-cell-count temporal spectrum discrepancy to remain ≤1% on the complete 0.015 clock. Report unequal-count gaps and equal-ID support mismatches explicitly. Dense spatial spectra are additional diagnostics, not a new spatial acceptance gate.

The parent study's **48 analytic calibration records** and its manifest are reused and hash-verified: all selected grids, probes and reconstruction orders were already covered under the same geometry/interface. This is inherited calibration evidence, not a new reconstruction computation. The earlier instrumentation gate and identical compiled method remain pinned. No full finer CPU reference ladder is added; earlier coarse CPU/GPU agreement is a separate scoped result.

## Full remaining-horizon restart replays

For the grid-120 candidate at dt=0.00046875 in **all four mechanical branches**, independently restore each saved mixed-count, active-furrow checkpoint from the 2→4, 4→8 and 8→16 windows. Run each restored state to t=18: **12 full replays** in addition to the original eight-step runtime checks. These are repeat executions of the same numerical histories, not new scientific replicates.

Each replay must retain positive finite chemistry and the same physical/event gates. Compare every later original observation exactly, and compare the complete state bitwise at every later regular checkpoint, including the final state. This includes visible phase, phase carry, chemical concentrations/corrections, targets, polarity, derived geometry, rounding counters, IDs/parents, schedules, lineage, live-event state, event audits, clocks and all random-generator states. Save replay finals and reports in separate directories; never write to the original source directories. A full replay checks restart fidelity of the same solver; it is not an independent physical reference or proof of equality at unsaved intermediate steps.

Replays run after production/refinement regardless of a numerical refinement failure, to isolate restart correctness from accuracy. Any physical/event or restart mismatch stops the coordinator and preserves the failure. Successful restarts cannot accept a failed refinement branch.

## Decision and preservation

The prospective scoped screen requires **all three branch decisions and all 12 distinct branch/stage replays**. Production results remain provisional until an independent saved-array assessment reconstructs their decisions. Preserve every failed coarser or finest result. Do not loosen limits, automatically increase a grid, alter the solver, or launch fresh scientific histories within this runner.

All previously consumed scientific files, the latest accepted review/report/figures, solver/kernel sources, ledger and manuscript are pinned by the [preservation manifest](../archive/study-document-snapshots/2026-10-10-before-development-controls/manifest.json). README/plan reporting copies are separate. Preparation requires a fresh study directory, the coordinator uses a nonblocking lock, and source/input/reuse evidence is checked before and after production. Run serially on the designated GTX 1080 Ti; the other GPU and unrelated processes are untouched.

The completed baseline suggests approximately **6–9 hours** for the 12 production paths, field comparisons and 12 remaining-horizon replays. This is a provisional estimate: mechanics, event times, final cell count and checkpoints affect throughput. Update it using measured branch/replay times.

After independent acceptance, launch a prespecified pilot of new independent zygotes with matched four-branch controls and a separately qualified post-cleavage scientific window. Report histories as the units of replication and measure chemical initiation/maintenance, attribute/response differences and shape separately.

Prepare with `python -m embryo.development_control_qualification prepare`, then launch `run` from the same pinned environment. The [runner](../embryo/development_control_qualification.py) reuses the original developmental dynamics and comparisons; [software guards](../tests/test_development_control_qualification.py) test branch matching, refusal to overwrite evidence and mandatory refinement/replay coverage.
