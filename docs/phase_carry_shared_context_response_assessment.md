# Reciprocal shared physical-context response test — independent assessment

**Accepted:** different inherited chemical states retain selective response differences when the physical context is shared. Cells 19 and 29 clear the prespecified response screen in both measured physical contexts, at both doses, all three timesteps and both reference solvers. Counts remain **0/2, 1/2 and 1/2** for histories 7, 8 and 9 in each physical context. All numerical gates pass. This supports conditional chemical-state sensitivity within the tested model; it does not establish recipient-autonomous identity.

## What was controlled

The [prespecified protocol](phase_carry_shared_context_response.md) crosses each actual t=510 sham/redistributed whole chemical vector with both measured transport operators and their matching cell volumes. Each crossed arm continues its own chemistry on fixed physics through t=570. Pulses occur immediately at reassignment and are compared with that arm's own unperturbed continuation. Chemistry includes both concentration levels and their allocation across cells. Transfers preserve concentrations, and explicitly record per-cell and total amount changes. Operator and capacity vector stay paired.

The common removed amount is recomputed within each shared physical context. Canonical concentration-based factors are exactly identical across the two physical contexts; amounts rescale with target volume. Those doses generally differ from the older native common-dose convention. Every trajectory was computed anew. Native comparison paths are numerical continuity checks, not scientific reuse.

There are **21 backgrounds, 69 logical starts and 57 distinct initial-value problems per timestep**. Twelve fraction/common preparations coincide exactly, with both labels retained. There are 138 primary plus 69 diagnostic GPU paths, 207 CPU verification paths and 138 independent reference trajectories. The independent units remain **three existing mature histories**, with two selected recipients nested within each and **zero new histories**. No mechanics, geometry/polarity, dilution or division is advanced in these frozen chemical continuations.

## Response results

The full 60-unit effect is the target's two-species time RMS of the difference between signed pulse/own-control log responses, normalized by each pulse's actual absolute initial log amplitude. The primary screen is strictly greater than 0.01. Each comparison below holds its physical context fixed between chemical states. Below-screen is not an equivalence claim.

| History / cell | Fraction: sham physics | Fraction: redistributed physics | Common amount: sham physics | Common amount: redistributed physics | Above screen at both doses / contexts |
| --- | ---: | ---: | ---: | ---: | --- |
| 7 / 22 | 0.009245 | 0.009272 | 0.009472 | 0.009498 | No |
| 7 / 23 | 0.000033 | 0.000030 | 0.000024 | 0.000027 | No |
| 8 / 19 | 0.070197 | 0.072793 | 0.070682 | 0.073269 | Yes |
| 8 / 20 | 0.002320 | 0.001866 | 0.001964 | 0.002891 | No |
| 9 / 29 | 0.072946 | 0.073022 | 0.073446 | 0.073519 | Yes |
| 9 / 30 | 0.005572 | 0.002443 | 0.007159 | 0.003727 | No |

Cells 19/29 have effects approximately **0.0702–0.0735** on either shared physical context. Their differing inherited operators and volume vectors are therefore unnecessary for these *detected response differences* at this horizon, conditional on either tested shared physical context. This does not imply that transport is unnecessary: each arm retains conservative exchange through its fixed network.

The physical contrast is also measurable: its largest RMS is **0.007515**, for cell 20 at redistributed chemistry under fractional depletion. The largest signed interaction RMS is **0.003576**, for cell 30 under the common-amount pulse. These are descriptive continuous measures, with no additional scientific screen. They must not be called absent or converted into a percentage of chemical causation. The exact signed vector decompositions hold to `5.56e-17`; scalar RMS effects do not add because their waveforms can reinforce or cancel.

All **48 fine pulse trajectories recover** toward their corresponding continued controls: target recovery is 2.25–4.35 units and volume-weighted network recovery is 2.25–5.40. Recovery does not erase differences between the chemical backgrounds themselves.

## Baseline states and reassignment limits

Whole-network chemical differences remain substantial during the last 24 observed units, with volume-weighted two-species log RMS **0.6398–0.9962**. These are finite-horizon observations, not an equilibrium or identity-persistence acceptance screen.

| History / cell | Final activator R/S: sham physics | Final activator R/S: redistributed physics | Late network log RMS: sham physics | Late network log RMS: redistributed physics |
| --- | ---: | ---: | ---: | ---: |
| 7 / 22 | 1.016016 | 1.015349 | 0.820945 | 0.820593 |
| 7 / 23 | 1.030038 | 1.029990 | 0.639840 | 0.641828 |
| 8 / 19 | 0.867502 | 0.868067 | 0.930717 | 0.908015 |
| 8 / 20 | 0.055364 | 0.060362 | 0.978543 | 0.972990 |
| 9 / 29 | 0.863418 | 0.868414 | 0.926083 | 0.914642 |
| 9 / 30 | 0.054804 | 0.047922 | 0.978050 | 0.996244 |

Cells 20/30 retain approximately **16.6–20.9-fold lower** activator than their shared-physics sham-chemistry controls, despite below-screen normalized response contrasts. Chemical levels and normalized response behavior remain different aspects of the phenotype.

Crossed controls relax after the instantaneous physical-context transfer. Maximum final volume-weighted network log drift from a reassigned start is **0.05191** (history 8, cell-20 redistributed chemistry on sham physics). The pulses consequently measure immediate-transplant conditional responses, not responses after an assumed equilibration.

Challenge-time chemical amounts and target chemistry remain unmatched. Redistributed/sham total activator amounts on shared physics range approximately 0.9218–1.0860; inhibitor amounts range 0.9316–1.0821. Switching physical contexts at fixed concentrations changes total amounts by approximately −0.0194% to +0.0176%; exact per-cell amounts and volume ratios are retained. The chemical factor cannot be interpreted as allocation alone, a dosage-independent effect, or isolated local memory.

## Numerical qualification and evidence

| Independent check | Result | Prespecified limit |
| --- | ---: | ---: |
| DOP853/Radau raw log discrepancy | `6.42e-12` | `1e-8` |
| DOP853/Radau normalized pulse discrepancy | `5.91e-11` | `1e-8` |
| DOP853/Radau signed contrast discrepancy | `3.22e-11` | `1e-8` |
| Production/reference raw log discrepancy | `1.23e-6` | `0.01` |
| CPU/PyTorch raw log discrepancy | `2.45e-13` | `1e-8` |
| Production normalized pulse discrepancy | `1.15e-5` | `0.01` |
| Production signed contrast discrepancy | `4.41e-6` | `0.01` |
| Maximum RMS contrast discrepancy | `2.31e-7` | `0.01` |
| Maximum relative response-area discrepancy | `7.72e-6` | `0.02` |
| Recovery-time discrepancy | `0` | `0.3` |
| Observed order, all 69 resolved cases | `2.00085–2.00699` | `1.7–2.3` |
| Maximum relative pulse amount error | `4.66e-15` | `2e-13` |

All 138 raw timestep pairs, 432 normalized pulse comparison rows and 540 signed contrast comparison rows qualify. Primary chemical screen decisions agree across all numerical methods and timesteps. The unchanged PyTorch float64 step checks both internal Euler stages; independently saved positivity and checkpoint consistency are reconstructed separately. Unsaved internal stages are runtime/code evidence, not independently recovered observations.

All 58 scientific source pins and 153 input pins verify, and **1,124 distinct consumed artifacts remain unchanged**. An independent script importing no project analysis helpers reconstructs 8,685 recorded quantities to maximum discrepancy `2.22e-16`. The 27 exactly matched native cases reproduce old paths at all three timesteps and both reference methods with zero log discrepancy. A separate [native continuity audit](../archive/study-reviews/2026-10-08-shared-context-response/native_continuity.json) checks 81 case/timestep pairs against 117 saved old curves, all byte-identical; duplicate old dose labels do not add independent cases.

GPU production took **144.14 seconds**; CPU verification used 763.46 summed worker-seconds and references 47.78 summed worker-seconds on four workers. These CPU sums are not elapsed wall time. The [launch report](phase_carry_shared_context_response_launch.md) retains tests, actual GPU prefixes, frozen source/runtime hashes and pre-launch guards.

The [comparison figure](../archive/study-reviews/2026-10-08-shared-context-response/shared-context-response.png), [recipient CSV](../archive/study-reviews/2026-10-08-shared-context-response/recipient-summary.csv), [independent reconstruction](../archive/study-reviews/2026-10-08-shared-context-response/independent_review.json) and [compact assessment](phase_carry_shared_context_response_assessment.json) retain the complete effects and numerical results. Raw chemical paths remain under `output/studies/phase-carry-shared-context-response`; no 3D fields were copied.

## What comes next

**Match the selected recipient's own activator and inhibitor while preserving the distinct surrounding chemical states on the same physical context.** This would test whether response differences remain attributable to surrounding-network chemistry once target chemistry is controlled. It is a separate intervention and has not been launched here.

Retain the limitations of three selected mature histories, historical no-carry developmental inheritance, unmatched challenge-time amounts, unchanged conductance closure and frozen mechanics. These results establish neither autonomous/inherited biological identity nor fresh carry-corrected zygote differentiation, population prevalence or a universal polarity mechanism. The results ledger, manuscript numerical Results, baseline code and shared kernels remain unchanged.
