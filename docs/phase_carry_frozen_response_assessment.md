# Completed frozen challenge-time response control — 8 October 2026

**Stopping physical evolution at the challenge time does not remove the detected delayed response differences.** Cells 19 and 29 remain above the prespecified response screen, at both original pulse doses and all three timesteps. History-level counts remain **0/2, 1/2 and 1/2** for histories 7, 8 and 9. All pulses recover toward their own continuing controls.

The completed study passes independent reconstruction of 99 production paths, 99 CPU verification paths and 66 reference trajectories. The independent units remain **three existing mature histories**. These paths include 66 primary coarse/fine paths and 33 finest-timestep numerical diagnostics; they do not represent new developmental replication.

## What was stopped and what was retained

The [prespecified control](phase_carry_frozen_response.md) starts from the **actual t=510 chemical states, measured volumes and conservative transport matrices** of the [accepted moving redistribution-response experiment](phase_carry_redistribution_response_assessment.md). It holds each background's operator and volumes fixed while its unchanged activator–inhibitor reactions and chemical exchange continue to t=570. Mechanical and polarity updates, changing transport and geometric dilution stop together. No chemical reservoir or continuing clamp is introduced.

Each of the nine backgrounds retains its own chemical continuation as the pulse control. The two original challenge designs remain a 10% recipient activator depletion and a common feasible removed amount, computed from the smaller target amount in the sham/redistributed pair. Each response uses its **actual initial log pulse amplitude**, including the much weaker common-amount sham pulses. Source chemistry, doses, cell order and observation times match the moving experiment exactly.

There are 33 logical starting cases: nine controls and 24 pulses. Six fractional/common-amount pairs coincide exactly, leaving **27 unique initial-value problems**. Both dose records are retained. Their arrays agree exactly at every numerical level and reference method. No moving path is reused as a frozen result.

The intervention originally matched chemical amounts at **t=450**, not at the t=510 challenge. Each frozen background still inherits different chemistry, contacts and volumes. This control asks whether their **continuing physical evolution during the challenge window** is necessary; it does not remove their earlier history or match the physical contexts to one another.

## Response comparison

The response effect is the target's two-species, pulse/own-control normalized log-response difference between redistributed and sham backgrounds, averaged as time RMS over the full 60 units. The operational screen remains **strictly greater than 0.01**. Below-screen is not equivalence or absence of any effect.

| History / cell | Fractional: moving → frozen | Common amount: moving → frozen | Frozen above screen at both doses |
| --- | ---: | ---: | --- |
| 7 / 22 | 0.009352 → 0.009387 | 0.009578 → 0.009613 | No |
| 7 / 23 | 0.000052 → 0.000050 | 0.000101 → 0.000099 | No |
| 8 / 19 | 0.072664 → 0.072162 | 0.073146 → 0.072643 | Yes |
| 8 / 20 | 0.006403 → 0.006365 | 0.007979 → 0.007940 | No |
| 9 / 29 | 0.072965 → 0.072470 | 0.073469 → 0.072974 | Yes |
| 9 / 30 | 0.006645 → 0.006612 | 0.008251 → 0.008217 | No |

The effects in cells 19 and 29 are approximately **0.7% smaller** after freezing. The largest absolute moving/frozen RMS-effect change is 0.0005031, for cell 19 at the common-amount dose. This is a small, numerically resolved change; unchanged classifications do not imply identical waveforms. The audit also preserves the full signed-waveform comparisons. A ratio of response norms is **not a fraction of an effect mechanically mediated**.

![Moving and frozen response effects and chemical levels](../archive/study-reviews/2026-10-08-frozen-response/frozen-response.png)

Panels A/B show the full-window effects under the two doses; panel C shows final unperturbed activator ratios. All bars use the fine timestep, `dt=0.001875`. [SVG figure](../archive/study-reviews/2026-10-08-frozen-response/frozen-response.svg) and [six-recipient CSV](../archive/study-reviews/2026-10-08-frozen-response/recipient-summary.csv).

All 24 fine pulse paths recover, with matching coarse, diagnostic and reference recovery times. Target recovery occurs within **2.25–4.35 units**, and network recovery within **2.25–5.40 units**. The definition requires at least 24 subsequent observed units within the 10%-of-initial-displacement band. Recovery is toward each background's continuing frozen control, rather than restoration of the sham state or exact equality.

Approximately 99.77%–99.98% of the integrated squared response differences occurs in the first six units. This is descriptive; the declared full-window measure remains primary. No early-window substitution or response-screen change was used.

## Chemical levels remain a different attribute

Cells 20 and 30 retain large concentration differences despite below-screen normalized response effects: frozen final activator is **13.6-fold and 17.5-fold lower** than sham, respectively. Their moving counterparts were 16.3-fold and 17.0-fold lower. Freezing therefore retains the large deficits but changes their values, particularly in cell 20. Maintained concentration and response behavior should continue to be reported separately.

Cells 19 and 29 have smaller frozen activator shifts, to approximately 0.871 and 0.870 times sham, while their fractional-pulse activator response areas are about **34% larger** than their respective sham responses. Across-cell log-activator contrast remains strong; its minimum over all fine paths is 1.355. These baseline measures are descriptive companions, not newly added scientific acceptance conditions.

## Numerical and provenance acceptance

The unchanged float64 PyTorch SSP-RK2 chemistry runs on the dedicated GTX 1080 Ti. Nine uncoupled blocks preserve the exact operators; positive duplicate-control padding is excluded from path counts. The coarse/fine/diagnostic timesteps are `0.00375`, `0.001875` and `0.0009375`. Every context uses one internal substep at each level, so batching does not alter the positivity restriction.

DOP853 and Radau independently integrate every logical start with `rtol=1e-13`, `atol=1e-15` and maximum step 0.15. Both trajectories are retained. Before launch, a looser reference diagnostic failed small-pulse normalized agreement; tightening integration qualified the unchanged screen. The [launch record](phase_carry_frozen_response_launch.md) retains both diagnostics, the 113 existing regression passes, 21 new test passes and actual CUDA/CPU prefix checks.

| Independently reconstructed quantity | Maximum / range | Declared limit |
| --- | ---: | ---: |
| Reference raw log disagreement | 4.49e-12 | 1e-8 |
| Reference normalized response disagreement | 4.26e-11 | 1e-8 |
| Production/reference raw log error | 1.20e-6 | 0.01 |
| CPU/GPU raw log disagreement | 2.45e-13 | 1e-8 |
| Coarse/fine raw log disagreement | 8.98e-7 | 0.01 |
| Fine/diagnostic raw log disagreement | 2.24e-7 | 0.01 |
| Resolved convergence orders | 2.00086–2.00658 | 1.7–2.3 |
| Normalized response discrepancy | 1.14e-5 | 0.01 |
| Relative activator-response area error | 7.48e-6 | 0.02 |
| Recovery-time disagreement | 0 | 0.3 |
| Numerical RMS-effect discrepancy | 2.28e-7 | 0.01, with matching decisions |

All 33 reference cases, 99 production accuracy cases, 33 resolved order tests, 66 raw timestep pairs, **216 pulse-response comparisons and 108 effect comparisons** pass. Eight independently generated pre-launch golden reference trajectories match exactly. The audit imports no project functions and performs no integrations; it reconstructs the checks directly from saved arrays, including 401-point clocks, source/dose equality, checkpoint mappings and the 66 matched moving histories.

All **56 scientific source pins, 115 input pins and 510 consumed artifacts** remain unchanged. The 2,415 independently recomputed recorded metrics agree within 1.25e-16. Saved concentrations are positive and finite. Unsaved internal Euler-stage positivity remains a runtime check in the pinned code; the sampled-data audit cannot independently regenerate those stages. No 3D mechanics is advanced here, and this result does not validate the physical accuracy of the inherited diffuse-interface conductance closure.

GPU production took approximately **142 seconds** for all three levels. Four-worker CPU verification consumed 361 summed worker seconds, and references 20 summed worker seconds; these sums are not total elapsed study time. The small chemical study avoids copying 3D phase fields.

## Interpretation and next control

The detected context-sensitive response differences **do not require continuing physical evolution after t=510**, given the inherited chemical state, graph and volumes. This narrows the mechanism: response behavior can persist through chemical dynamics on an already established physical environment. It does not isolate whether the inherited chemical configuration or the inherited geometry is responsible, or whether earlier mechanics was needed to establish either.

The next control should put the different chemical configurations on a **shared transport operator and its matching volume vector**. The operator and capacities must move together to preserve conservative exchange. Specify whether concentrations or amounts are preserved during reassignment; preserving concentrations can change amounts when volumes differ. Recipient-state matching is a separate subsequent intervention for testing surrounding-network sensitivity at equal recipient chemistry. No additional experiment is launched by this assessment.

This is a finite-horizon result in three inherited mature histories. Autonomous/inherited biological identity, irreversible fate, nearest-neighbor specificity, fresh carry-corrected zygote development and population prevalence remain untested. The ledger and manuscript numerical Results retain their previous qualifications.

The [compact assessment](phase_carry_frozen_response_assessment.json), [independent reconstruction and hashes](../archive/study-reviews/2026-10-08-frozen-response/independent_review.json), [reconstruction script](../archive/study-reviews/2026-10-08-frozen-response/independent_review.py) and [review provenance](../archive/study-reviews/2026-10-08-frozen-response/review-provenance.json) preserve the evidence. Scientific parent outputs, prior assessments, baseline kernels, the ledger and manuscript remain unchanged; README and the working plan record completion.
