# Completed chemical-distribution versus dosage assay — 7 October 2026

**Changing where the surrounding cells contain activator and inhibitor changes the recipient's maintained chemical state, even when each species' initial amount is preserved.** Both selected recipients exceed the prespecified late-state screen in each of the three existing histories, at both initial dosage levels. Changing dosage alone stays below that screen in every recipient at both distributions. All 42 paths pass the unchanged numerical and physical criteria.

Dosage does produce early responses. The finding concerns **maintained states over the declared late window**, not the absence of any dosage effect. It controls the initial amount confound in the previous surrounding-cell reset; it does not establish autonomous identity, a polarity-specific cause or a unique molecular mechanism.

## Design and controlled quantities

The [unchanged protocol](phase_carry_context_dosage.md) starts from the common fine **t=450** states of the moving network-context experiment and follows chemistry, mechanics, polarity, conservative transport and dilution to **t=510**. The aggregate has 16 cells; there is no new division or new zygote development in this assay. Model time and amounts are dimensionless and are not calibrated to biological hours or molecule numbers.

For each selected recipient, intervene on **all other cells**, using their original concentrations or the same-ID reference concentrations already defined by the parent assay. The recipient's activator and inhibitor are preserved exactly. The initial saved phase fields, volumes, polarity, residual carry and random-state metadata are shared between arms. Chemistry-dependent tension and adhesion can change immediately, so their derived values are not controlled by preserving the saved physical fields. Mechanics and chemistry then co-evolve.

| Arm | Initial chemical distribution in other cells | Initial amount of each species in other cells |
| --- | --- | --- |
| Sham | Original | Original |
| Raw reset | Reference | Reference |
| Amount-preserving redistribution | Reference, uniformly scaled per species | Original |
| Dosage-only rescaling | Original, uniformly scaled per species | Reference |

For a species, let `x_j` be its original concentration, `r_j` its reference concentration, and `V_j` the common initial measured volume of nonrecipient cell `j`. Define `Q = sum(V_j*x_j)` and `R = sum(V_j*r_j)`, summing only over nonrecipient cells. Redistribution uses `r_j*(Q/R)`; dosage-only rescaling uses `x_j*(R/Q)`. The construction is applied separately to activator and inhibitor. It preserves each species' relative concentration ratios across the intervened cells, not every cell's activator/inhibitor ratio. Because recipient chemistry and volume are fixed initially, matching surrounding amounts also matches the initial whole-aggregate amounts.

These are **initial** controls. Reactions can subsequently change chemical amounts; conservative transport and volume-dilution updates preserve aggregate amounts in their respective steps. No reservoir clamps the totals. For example, final aggregate activator after redistribution differs from its sham total by up to about 8.6%, despite exact initial matching. Later amount differences may be part of the response to redistribution rather than an excluded mediator.

There are **24 new GPU paths and 18 exactly verified reused controls**, at dt=0.00375 and 0.001875. Both timesteps use the same fine preparation for each arm. Six recipient comparisons and their four arms are nested within **three existing histories**, with **no new independent developmental histories**. PyTorch arrays and matrix operations are combined with custom CUDA mechanics and geometry/polarity kernels; the phase fields retain inherited float64 residual carry while chemistry is float64.

## Maintained state effects

The primary metric is the time RMS of the recipient's **signed log concentration difference**, combining activator and inhibitor equally, over elapsed **36–60**. A detected effect is **strictly greater than 0.01**, an operational effect screen, not a cell-type boundary, an equivalence margin or a statistical significance test. Coarse and fine decisions must agree; the unchanged full/late metric-discrepancy tolerance is 0.01.

Compare redistribution with sham to test distribution at original dosage, and raw reset with dosage-only rescaling to test distribution at reference dosage. Compare dosage-only rescaling with sham to test dosage at original distribution, and raw reset with redistribution to test dosage at reference distribution. The interaction is formed from the signed log trajectories as `log(reset) - log(bulk) - log(redistributed) + log(sham)` **before** taking the RMS.

Fine results are:

| History | Recipient | Distribution, original dose | Distribution, reference dose | Dosage, original distribution | Dosage, reference distribution | Interaction | Final activator redistribution/sham |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 22 | 0.015494 | 0.015490 | 9.59e-06 | 5.42e-06 | 4.17e-06 | 1.01550 |
| 7 | 23 | 0.025205 | 0.025205 | 7.04e-05 | 7.07e-05 | 5.88e-07 | 1.02733 |
| 8 | 19 | 0.135458 | 0.135464 | 6.58e-06 | 2.45e-07 | 6.36e-06 | 0.87048 |
| 8 | 20 | 1.957196 | 1.957250 | 0.000195 | 9.56e-05 | 0.000136 | 0.07361 |
| 9 | 29 | 0.136768 | 0.136771 | 6.57e-06 | 3.79e-06 | 2.77e-06 | 0.86931 |
| 9 | 30 | 2.123871 | 2.123882 | 0.000116 | 0.000105 | 1.17e-05 | 0.05712 |

Thus the history-level counts are **2/2 recipients with maintained distribution effects in each history**, **0/2 with maintained dosage effects in each history**, and **0/2 with an above-screen late interaction in each history**. The same decisions hold at both timesteps and both conditional levels of the factorial. Recipients and contrasts are not additional independent histories.

The largest late dosage RMS is **0.000195**, more than 51 times below the screen. Amount-preserving redistribution and raw reset differ by at most **0.000105** in late-state RMS. This supports the conclusion that the initial bulk-amount change is not required to produce these maintained shifts at this tested point; it does not prove exact equivalence or chemical-amount independence in other regimes. The largest late signed-interaction RMS is **0.000136**.

Effects vary substantially across histories and recipients. In history 7, final activator shifts are only about **1.5% and 2.7% upward**. In histories 8 and 9, recipients 19 and 29 shift about **13% downward**, while activator in recipients 20 and 30 ends **13.6-fold and 17.5-fold lower** than its own sham, respectively. Species-specific fold changes come from final log differences; the two-species RMS itself is not a fold-change measure.

![Chemical redistribution shifts maintained recipient states](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/context-dosage.png)

Left: fine late RMS for redistribution at original amounts and dosage-only rescaling at original distribution. Right: activator trajectories for the two strongest recipients relative to their own sham, with the declared late window shaded. Amount matching refers to the initial intervention. Timesteps, arms and selected recipients are not treated as independent replicates. [Editable SVG](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/context-dosage.svg).

## Dosage effects are transient rather than absent

Dosage-only rescaling has a peak instantaneous two-species log distance as large as **0.0811** and a full-window RMS as large as **0.011827**. Several dosage contrasts therefore exceed 0.01 if assessed over the full window. Those are descriptive early responses, not the primary maintained-state outcome: every late dosage contrast stays below 0.01. The interaction is also larger during reorganization than during the late window. The assessment does not move the window or change the threshold after observing the results.

All paths retain strong across-cell chemical contrast at every saved observation: minimum SD of ln(activator) is **1.3253**, above 0.1. This is finite-horizon persistence on a live geometry, not proof of an indefinitely stationary attractor.

## Numerical and provenance assessment

The independent CPU audit reconstructs preparations directly from common initial volumes and concentrations rather than invoking the intervention helper. It verifies **51 source hashes, 1,073 input hashes, all 42 preparations and context gates, and 1,492 immutable evidence/source paths before and after review**. Initial amount relative error is at most **3.84e-16**; formula log error is at most **5.83e-16**. Exact recipient preservation, physical-array dtypes/values, retiming, residual carry and random-state metadata are checked.

All three pilot reports and all three full reports pass. The full stage contains **21 matched mechanical pairs, six recipient comparisons and 36 factorial contrasts**. Raw trajectories, final checkpoints, mechanical deviations and full/late state contrasts, including signed interactions, were reconstructed independently and matched the original reports. Maximum full/late RMS disagreement between timesteps is **4.05e-06**. No temporal realignment or scientific rerun was used.

| Coarse/fine quantity | Maximum discrepancy | Unchanged limit |
| --- | ---: | ---: |
| Log concentration, all cells/species/observations | 4.34e-05 | 0.01 |
| Polarity | 3.12e-06 | 0.01 |
| Relative axis ratio | 4.99e-08 | 0.01 |
| Relative cell volume | 9.72e-07 | 0.005 |
| Relative transport matrix | 5.31e-06 | 0.01 |
| Final phase field | 1.73e-06 | 0.02 |

Maximum recorded volume error is **1.139%**, minimum equivalent radius **5.107 grid spacings**, clipping zero, sampled boundary occupancy at most **1.51e-08**, and dilution amount-conversion error at most **4.44e-16**.

The audit reconstructs conservation, volume-weighted symmetry/positive semidefiniteness, Laplacian spectra and uniform-state Gierer–Meinhardt growth diagnostics across **16,842 saved observations**. Graph diagnostic errors are below 1.7e-15. This checks the implemented operator, not the physical accuracy of diffuse-interface conductance closure. Recorded eight-step independent chemistry and four-step exact-restart prefix gates are verified from immutable evidence rather than rerun. Pilot quality records and comparisons are verified; long continuations intentionally supersede new pilot checkpoint files. Saved observations do not independently reconstruct every-step extrema.

## Meaning for the identity question and next test

This factorial strengthens **network-supported context sensitivity of maintained chemical states**. A cell can change its maintained activity because chemical resources are allocated differently across the surrounding network, even when the initial aggregate amount is unchanged. Global chemical abundance alone does not predict the tested recipient outcome.

A plausible explanation comes from conservative contact transport: the input to a recipient depends on contact-weighted concentration differences from specific cells. Equal total amount does not imply equal local input. Nonlinear activator–inhibitor reactions and mechanically evolving conductances can then sustain different trajectories. This is an interpretation consistent with the model; the assay does not uniquely separate direct chemical transport from subsequent geometry, polarity, or total-amount mediation. It changes every nonrecipient cell, so it also does not isolate the nearest neighbors or establish network necessity.

Chemical state and perturbation-response behavior remain different phenotype attributes. The previous delayed-pulse assay used raw resets, which also changed initial amounts. **The next recommended test is delayed negative activator pulses from the new amount-preserving redistribution endpoints**, with each background's own continuing moving control, fractional and feasible common-amount doses, and paired timestep gates. Compare those responses with the sham endpoints and retain below-screen outcomes. This would ask whether controlled chemical redistribution changes behavior as well as maintained levels. No pulse assay was performed by this factorial or launched by this review.

Afterward, frozen-geometry or replay controls can separate direct chemical-context effects from mechanical mediation. Autonomous identity, inheritance through division, fresh carry-corrected zygote development, general spatial/transport convergence and broader parameter robustness remain separate requirements. These mature continuations inherit earlier no-carry development; they do not broadly replace the ledger or manuscript snapshot.

## Evidence

- [Frozen protocol JSON](../outputs/phase-carry-context-dosage/protocol.json), SHA-256 `547d648ebb819e8f892c523328832f57dea7694ccf21b22fcb0848cc3ed02a2a`.
- [Original completed status](../outputs/phase-carry-context-dosage/status.json) and [summary](../outputs/phase-carry-context-dosage/summary.json).
- [Independent full assessment](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/assessment.json), [saved audit script](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/assessment_script.py), [compact assessment](phase_carry_context_dosage_assessment.json), and [audit provenance](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/review-provenance.json).
- [Recipient table](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/recipient-summary.csv), [plot data](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/plot-data.json), and [figure provenance](../outputs/phase-carry-context-dosage-review/2026-10-07T211547Z-final/figure-provenance.json).
- [Pre-review reporting snapshot](../archive/study-document-snapshots/2026-10-07-before-completed-context-dosage-review/manifest.json).

Scientific sources, kernels, tests, protocol documents, preparation records, histories, checkpoints, prior assessments, ledger and manuscript remain unchanged. README and the working plan now point to this completed review. No new simulation was launched.
