# Fresh-zygote pilot: matched development before long-window identity assays

## Question and experimental unit

Do the organization mechanisms found in mature tissue remain relevant when the tissue develops from one cell with the current numerical method? Start **three new histories, seeds 1101, 1102 and 1103**, without a supplied mature geometry, chemical pattern or prescribed identity. A history is one stochastic developmental realization; its branches, cells and timestep repeats are nested observations. Three histories constitute a pilot, not a reliable population-frequency estimate.

The first stage develops all three histories through cleavage. Later stages will test initiation and maintenance over longer times and then test cell responses. Passing the first stage does not by itself establish persistent chemistry, biological identity or lineage commitment.

## Matched mechanical controls

| Branch | Activity-dependent tension `c_gamma` | Activity-dependent adhesion `c_A` | Polarity tension `chi` |
| --- | ---: | ---: | ---: |
| Feedback off | 0 | 0 | 0 |
| Tension and adhesion | 0.25 | 0.35 | 0 |
| Polarity tension | 0 | 0 | 0.35 |
| Full combination | 0.25 | 0.35 | 0.35 |

All branches retain activator–inhibitor chemistry, conservative volume-weighted transport, dilution, polarity evolution and the same baseline mechanics. Only the indicated mechanical actions differ. Polarity is a measured continuous vector, not a predefined cell type. Tension and adhesion follow activity instantaneously and therefore cannot count as independent behavioral evidence of identity.

Within each history, the full initial arrays, cell-cycle schedule, lineage, and geometry/chemical random states match exactly across all four branches and both timesteps. Subsequent geometry-dependent cleavage can cause branches to diverge in event order and random-number assignment; matched initial streams do not promise identical daughter noise after such divergence. Event draws and lineage are saved for inspection. No cells, branches or repeats are counted as additional independent histories.

## Initial state and symmetry breaking

Each history begins at time zero with one analytic spherical cell of radius 0.8, interface width 0.085 and the same continuously integrated target occupancy volume as the accepted qualification. The domain is unchanged, extent 2.24 and grid 120 cubed. The zygote has uniform activator 1.03, inhibitor 0.98 and zero polarity. Both chemicals are regulatory species; they are not two identity labels. Unused compatibility fate values remain zero, with the historical fate equation disabled.

No localized chemical cue or imposed fate is added. The initial global chemical displacement is identical in every branch and introduces no spatial axis. The qualified chemical partition rule introduces unbiased zero-mean daughter fluctuations with amplitude 0.001 while preserving mother chemical amount at each division. Cell-cycle jitter is 0.12; shape-based division orientation uses a random tie-break only for nearly degenerate axes. These existing stochastic mechanisms, driven by distinct history seeds, supply developmental variation. A single chemical compartment cannot possess a cell-to-cell spatial pattern before division.

Seed reservation is checked against all existing protocol/launch JSON inventoried under `outputs`, `output/studies` and `archive`. The inventory and its source hashes are frozen before preparation. This establishes non-reuse in recorded local studies, not in undocumented external simulations. Failed or nonorganizing histories will not be replaced.

## First executable stage: development from time 0 to 18

The [accepted matched-control qualification](development_control_qualification_assessment.md) supports grid 120 and developmental timestep **0.00046875** for its one analytic seed-7 preparation. New seeds therefore receive their own paired repeats at **0.000234375**, with no relaxation of the inherited limits. This creates **24 trajectories: 3 histories × 4 branches × 2 timesteps**. Each is a new zygote, not a restart from histories 7–9.

The unchanged resident `DevelopmentGpu` solver uses PyTorch for GPU arrays/matrix operations and custom CUDA for mechanics, geometry, polarity and progressive furrows. Visible phase fields are float32 with a float64 phase residual; chemical state and its rounding correction are float64. Progressive, naturally scheduled divisions stop at sixteen cells. No oxygen or nutrient variables are introduced: adequate resources remain an implicit condition of the prescribed cell cycle.

Regular measurements occur every 0.15 time units, dense network records every 0.015 plus event transitions, and complete checkpoints at 0, 3, 6, 9, 12, 15 and 18. Save event geometry, chemical partition draws, conservative transport/conductance matrices, generalized graph spectra, volume/shape/polarity diagnostics, chemical distributions and full correction/RNG state. Checkpoints allow later cell-attribute reconstruction without assigning clusters or identities.

### Actual-start and numerical gates

Each of the 24 starts receives an eight-step independent NumPy full-state prefix and an exact correction-aware GPU restart check before production. Prefix limits are 2e-5 absolute visible/effective phase, 1e-5 maximum log chemical discrepancy, 1e-5 absolute polarity, 1e-5 relative volume and 1e-5 relative transport discrepancy. A prefix is a component/reference check, not a full fine CPU developmental ladder.

Every production step retains positive finite chemistry, zero clipping, volume error below 5%, equivalent radius at least four voxels, dilution conservation error at most 2e-14 and unchanged total inherited target volume. Observations require boundary occupancy below 0.01 and one connected component for each nondividing cell. Every abscission must pass the inherited neck, overlap, connected-daughter and amount/volume/occupancy conservation limits. Require fifteen completed divisions, sixteen final cells, no active furrows and at least three mature time units before time 18.

Live mixed-count 2→4, 4→8 and 8→16 checkpoints each receive an exact eight-step restart comparison, including phase carry, chemical correction, event and random state. The accepted previous full-horizon replays remain inherited numerical evidence; this new pilot does not duplicate them or claim new full-horizon CPU equivalence.

For each history/branch, retain all inherited timestep limits: maximum density-field discrepancy 1%, relative spectrum discrepancy 1%, chemical/log distribution distance 0.01, contrast difference 0.005, shape-axis discrepancy 1%, total-volume/amount discrepancy 0.5%, polarity discrepancy 0.01 and maturity-time difference 0.05. Regular unequal-count observations remain recorded, with the inherited limit of 30. A separate dense equal-count spectrum check must also remain below 1%; unequal-count dense records and edge-support mismatches remain visible. These distribution-based tests do not assert homologous daughter correspondence after divergent cleavage.

Execute serially on the dedicated GTX 1080 Ti, pairing the two timesteps before advancing to the next branch/history. Stop on any numerical, actual-start or physical failure. Preserve completed and failed data and the unexecuted registered design; do not substitute seeds, change thresholds or continue failed cases into science. Results remain provisional pending an independent saved-array review. No automatic grid, timestep or model extension is authorized by the coordinator.

## Separately gated post-cleavage study

The developmental coordinator deliberately ends at time 18. After independent acceptance, continue the actual saved states, preserving chemistry, geometry, polarity, numerical corrections and random state. Do not reset chemistry, supply a pattern, reopen division or silently change backend. The proposed first mature window is **60 elapsed units**, followed by a **240-unit total continuation** only after the pilot passes its own numerical gates. The final late window is 24 units. These future durations are declared now; their timestep and long-window acceptance are not inherited from the cleavage test.

First qualify efficient mature timesteps using candidates 0.00375, 0.001875 and, if required, 0.0009375 on actual endpoints. Each new mechanical context needs its own CPU/GPU prefix, correction-aware restart, physical checks and matched time refinement. Failure requires a separately recorded diagnostic plan, without retroactively changing this developmental experiment. No longer runs are started by the present coordinator.

The follow-up must separate first sustained formation, continued retention and loss/recovery. Report continuous chemical contrast and prespecified sustained-threshold sensitivity, measured polarity and cell-shape attributes, contact conductances/spectra and tissue shape separately. Compare matched controls within each history; a branch remaining uniform is a scientific outcome if numerical checks pass. Calling a difference a candidate identity additionally requires response/context experiments, and biological commitment or inheritance requires further tests. Frozen endpoint spectra and local uniform/patterned coexistence assess available chemical states, not the stability of the whole moving system.

## Provenance and implementation

[Coordinator](../embryo/fresh_zygote_pilot.py), [design and failure-path tests](../tests/test_fresh_zygote_pilot.py). The protocol pins the unchanged accepted solver/kernel sources, accepted assessment/review, seed inventory, new starts and this design. The coordinator refuses existing output directories, uses an exclusive lock, and never overwrites a recorded case or comparison. Large results reside on the data drive through `output/studies`; a SHA-verified relocation wrapper preserves the canonical paths in historical protocols. Previous failed reports, ledger entries and manuscript claims remain unchanged.
