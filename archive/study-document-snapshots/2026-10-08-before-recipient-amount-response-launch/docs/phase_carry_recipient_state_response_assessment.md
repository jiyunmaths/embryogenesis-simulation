# Recipient-state matching — accepted assessment

**Matching the recipient's initial chemistry does not remove surrounding-network effects.** Cells **19, 20, 29 and 30** clear the prespecified response screen using both measured recipient reference states and both shared physical contexts. The counts are **0/2, 2/2 and 2/2** recipients in histories 7, 8 and 9, respectively, in every reference/context combination. All numerical gates pass. These are immediate post-reset responses in three existing mature histories, not settled response traits or autonomous cell identities.

## What was controlled

The [prespecified protocol](phase_carry_recipient_state_response.md) asks whether surrounding chemistry changes a recipient's response when its own initial activator and inhibitor are equal. For each selected recipient, both surrounding starts use the same measured pair: either its sham pair or its redistributed pair. Every other cell's concentrations remain exactly as measured. The two chemical starts share the same saved conservative transport operator, matching volumes and cell-ID order. Repeat this on both measured physical contexts.

The matching and 10% activator depletion occur together at t=510. The target then evolves freely; each pulse is compared with its own matched, unpulsed continuation through t=570. There is no continuing clamp or pre-pulse equilibration. The target concentration, volume and removed activator amount are equal between surrounding starts within each reference/physical-context comparison. Fractional and common-amount depletion therefore coincide; there is one challenge design, not two independent doses. Global chemical amounts remain unmatched. Preparation records separately preserve physical-volume reassignment and local-matching amount changes for both species.

There are **45 distinct backgrounds and 93 distinct initial-value problems**: 45 controls and 48 target pulses. The sham/sham native control is shared across the two recipients within each history; all starts are calculated anew. Three timesteps produce **279 GPU paths**, **279 CPU checks**, and **186 DOP853/Radau reference trajectories**. Forty-five unchanged cases reproduce accepted parent paths at five numerical settings, providing 225 continuity comparisons. No parent path is counted as a new calculation or new history. Six selected recipients are nested within **three histories**, with **zero new developmental histories**.

Only Gierer–Meinhardt chemistry and conservative exchange advance. Transport and volumes are frozen: mechanics, polarity, geometric dilution and division are not advanced. The [qualified grouped GPU layout](phase_carry_recipient_state_response_grouped.md) uses independent chemical columns on nine measured physical blocks, with 42 duplicated positive-control padding columns excluded from reported paths. It changes no equation or shared scientific kernel.

## Response results

The primary effect is the full 60-unit time RMS of the difference between the target's signed, two-species pulse/own-control log responses, normalized by the actual absolute initial log pulse amplitude. The screen is strictly **greater than 0.01**. Each table entry holds initial target chemistry and physics equal between the two surrounding starts. Below-screen is not equivalence.

| History / cell | Sham match, sham physics | Sham match, redistributed physics | Redistributed match, sham physics | Redistributed match, redistributed physics | Above screen in all four |
| --- | ---: | ---: | ---: | ---: | --- |
| 7 / 22 | 0.006939 | 0.006967 | 0.006710 | 0.006737 | No |
| 7 / 23 | 0.000897 | 0.000894 | 0.000906 | 0.000903 | No |
| 8 / 19 | 0.053137 | 0.055205 | 0.066995 | 0.069693 | Yes |
| 8 / 20 | 0.097906 | 0.092006 | 0.052608 | 0.047396 | Yes |
| 9 / 29 | 0.054889 | 0.054997 | 0.069517 | 0.069662 | Yes |
| 9 / 30 | 0.098431 | 0.100553 | 0.048896 | 0.049131 | Yes |

Cells 19/29 retain detected differences after target matching. Cells 20/30, previously below screen with unmatched target chemistry, now have effects **0.0474–0.1006** under all four matching/context combinations. Holding the recipient's starting concentrations equal therefore does not eliminate the influence of the rest of the chemical network. Matching changes the subsequent trajectory; comparing old and new effects does not yield an additive percentage attributable to local versus surrounding chemistry.

All **48 fine pulse paths recover** toward their own controls. Target and volume-weighted network recovery each range from **0.45 to 5.40 units**, using the prespecified sustained-recovery definition and remaining observation requirement. Pulse recovery does not erase differences between the two chemical backgrounds.

The largest physical-context and signed interaction RMS values are **0.011955** and **0.009881**, respectively. They were prespecified as descriptive continuous measures with no scientific screen; the primary chemical screen must not be applied to them. Signed waveform identities hold to `1.11e-16`, but their scalar RMS magnitudes are not additive.

## Target states reappear, but immediate responses depend on the reset

Every matched control starts with redistributed/sham target ratios exactly **one** for both species. Later, differences reappear. The following fine results use the sham matching reference; the alternative reference returns to the same endpoints to numerical precision in the observed window.

| History / cell | Final activator R/S: sham physics | Final activator R/S: redistributed physics | Late network log RMS: sham physics | Late network log RMS: redistributed physics |
| --- | ---: | ---: | ---: | ---: |
| 7 / 22 | 1.016016 | 1.015349 | 0.820945 | 0.820593 |
| 7 / 23 | 1.030038 | 1.029990 | 0.639840 | 0.641828 |
| 8 / 19 | 0.867502 | 0.868067 | 0.930717 | 0.908015 |
| 8 / 20 | 0.055364 | 0.060362 | 0.978543 | 0.972990 |
| 9 / 29 | 0.863418 | 0.868414 | 0.926083 | 0.914642 |
| 9 / 30 | 0.054804 | 0.047922 | 0.978050 | 0.996244 |

Cells 20/30 recover approximately **16.6–20.9-fold lower activator** under redistributed surrounding chemistry despite initially equal target states. The surrounding network can restore local differences in this model. Whole-network log separation during elapsed 36–60 units remains approximately **0.6398–0.9962**. These are descriptive finite-horizon state differences, not an equilibrium or biological identity acceptance test.

A separate [post-run endpoint diagnostic](../archive/study-reviews/2026-10-08-recipient-state-response/reference-endpoint-diagnostic.json) compares the two matching references at fixed surrounding chemistry and physics. Across 24 unpulsed whole-network pairs, maximum final log discrepancy is `2.61e-13`, and maximum discrepancy in the last 24 units is `6.53e-11`. This diagnostic integrates no new trajectories and introduces no acceptance threshold. It shows agreement for these measured references and this horizon; it establishes neither unique attractors nor general equivalence.

The immediate pulse effects differ between the two matching references even though their unpulsed endpoints agree. Cells 20/30 show their large response differences during early reset relaxation, as the [figure](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.png) illustrates. Thus initial matching reveals surrounding-network sensitivity, but it does not yet establish different response behavior after the controls have settled. Delaying the pulse is the next necessary distinction.

## Numerical qualification

| Saved-array check | Maximum / range | Prespecified limit |
| --- | ---: | ---: |
| DOP853/Radau raw log discrepancy | `6.42e-12` | `1e-8` |
| DOP853/Radau normalized pulse discrepancy | `5.91e-11` | `1e-8` |
| DOP853/Radau signed contrast discrepancy | `2.79e-11` | `1e-8` |
| Production/reference raw log discrepancy | `7.25e-6` | `0.01` |
| CPU/PyTorch raw log discrepancy | `4.96e-14` | `1e-8` |
| Production normalized pulse discrepancy | `1.14e-5` | `0.01` |
| Production signed contrast discrepancy | `4.07e-6` | `0.01` |
| RMS contrast discrepancy | `1.66e-7` | `0.01` |
| Relative activator-response area discrepancy | `7.47e-6` | `0.02` |
| Recovery-time discrepancy | `0` | `0.3` |
| Observed order, all 93 resolved cases | `2.00056–2.00699` | `1.7–2.3` |
| Relative requested/measured pulse amount error | `9.99e-16` | `2e-13` |

All 93 raw dual-reference checks, 279 production accuracy checks, 186 raw timestep comparisons, 432 pulse comparison rows and 540 contrast comparison rows pass. Primary decisions agree at all three timesteps and both reference solvers. Saved states are positive; padding and checkpoint histories agree exactly. Internal Euler-stage checks are unchanged runtime/source evidence, not independently reconstructed unsaved observations.

A separate [raw-array audit](../archive/study-reviews/2026-10-08-recipient-state-response/independent_review.json), importing no project analysis helpers, reconstructs **6,357 recorded quantities** to maximum discrepancy `2.22e-16`. All **62 source pins, 250 input pins and 1,162 consumed artifacts** verify unchanged. The 225 native case/method comparisons have zero log discrepancy. The [launch evidence](phase_carry_recipient_state_response_launch.md) retains 13 design tests, three grouping tests, literal prepared-input checks and all actual CUDA prefixes. The initial grouping-test missing-import failure is retained; it was fixed before preparation, without changing scientific criteria.

GPU production took **165.82 seconds**. Full coordinator elapsed time was **287.80 seconds** (about 4.8 minutes), with four single-thread CPU verification workers overlapping GPU production. The grouped layout was **4.12 times faster in the short prefix benchmark**; no full ungrouped timing experiment was run. Worker-time sums are not elapsed wall time.

The [PNG](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.png), [SVG](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.svg), [24-comparison CSV](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-summary.csv), [compact assessment](phase_carry_recipient_state_response_assessment.json) and review provenance preserve the reported result. Raw chemical trajectories remain in `output/studies/phase-carry-recipient-state-response-grouped`; no 3D fields were copied. The earlier ungrouped preparation remains unlaunched and unchanged.

## Next priority and limits

**Delay the challenge after target matching**, keeping each evolving control and the same crossed chemical/physical contexts. Compare immediate and delayed responses to determine which differences describe reset relaxation and which remain after local chemical states have re-established. This follow-up has not been launched. A subsequent initial-amount-matching control would separate surrounding allocation from total network dosage.

The present intervention changes the chemistry of the other fifteen cells as a whole and leaves their global amounts unmatched. It establishes neither nearest-neighbor specificity nor dosage-independent allocation effects. Retain the limits of three selected mature histories, historical no-carry developmental inheritance, sixteen-cell resolution, frozen mechanics and the inherited conductance approximation. No autonomous or inherited biological identity, fresh carry-corrected zygote differentiation, population prevalence or universal polarity mechanism is established. The results ledger, manuscript, accepted parent, baseline simulator and shared kernels remain unchanged; no subsequent experiment is launched here.
