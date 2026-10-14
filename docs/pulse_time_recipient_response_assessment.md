# Challenge-time recipient matching: independent completed review

The completed study is **accepted at its tested scope**. On the actual t=570 graphs produced by the moving parent, surrounding chemistry can change a recipient's immediate pulse-response trajectory even when its starting activator/inhibitor concentrations and physical context are equal. This response sensitivity occurs in both selected recipients of histories 8 and 9, and stays below the declared screen in history 7. All pulses recover toward their own continuing controls. The evidence supports context-sensitive chemical behavior during a reset; it does not establish autonomous or inherited cell identities.

## What was controlled

The study uses twelve independently accepted fine waiting endpoints from the [moving initial-allocation experiment](moving_recipient_amount_response_assessment.md), comprising six selected recipients nested in **three existing developmental histories**, 7–9. It adds zero histories. Each endpoint's actual conservative transport operator and volumes are frozen through t=570–630.

For each recipient, compare whole-network chemical allocations S and R on the same endpoint physical context. S denotes the sham chemical allocation; R denotes redistribution. Each physical context is the geometry that co-evolved with that allocation in the parent. Both parent allocations began with the original sham physical fields at t=510, so these new endpoint labels are not the older historical sham/redistributed-physics branches.

Three interventions are compared: retain the native recipient state, match both recipient species to its measured S pair, or match both to its measured R pair. Every other concentration is retained exactly. At the reset time, apply a −10% activator pulse and compare it with its **own unpulsed, identically reset/reassigned control**. Both then evolve freely. Matched local concentrations and common recipient volume give equal removed activator amounts across the allocations; fractional and equal-amount pulses therefore coincide and are not duplicated as dose evidence.

The target is equal **at the pulse**, not throughout the response. The controls themselves relax after matching. This is an immediate reset/pulse experiment, not a later settled challenge. Reassignment preserves concentrations rather than amounts when capacities change; matching is a declared external chemical intervention. Both changes are accounted for explicitly. Neither whole-network nor surrounding species amounts are matched at t=570.

## Scientific outcomes

The primary effect E is the 60-unit time RMS difference between the two allocations' signed recipient responses. Each response contains activator and inhibitor log ratios relative to its own unpulsed control, divided by the magnitude of its actual log pulse. E averages over the two species. The strict **E > 0.01** screen is a predefined operational criterion, not a statistical significance test or a percentage difference in concentration. Values below it are not evidence of equivalence.

The following fine-timestep ranges span the two shared endpoint physical contexts:

| History | Recipient | Native, no reset | Matched to S pair | Matched to R pair |
|---|---:|---:|---:|---:|
| 7 | 22 | 0.007434–0.007456 | 0.005763–0.005787 | 0.005618–0.005641 |
| 7 | 23 | 0.0000287–0.0000314 | 0.000858–0.000861 | 0.000866–0.000869 |
| 8 | 19 | 0.063481–0.066559 | 0.049706–0.052208 | 0.061071–0.064249 |
| 8 | 20 | 0.002345–0.005484 | 0.096365–0.098605 | 0.048450–0.048838 |
| 9 | 29 | 0.067313–0.067898 | 0.052076–0.052602 | 0.064670–0.065348 |
| 9 | 30 | 0.002581–0.005598 | 0.098034–0.100137 | 0.048770–0.049067 |

With no local reset, above-screen counts are **0/2, 1/2 and 1/2** recipients in histories 7, 8 and 9. After matching to either reference, the counts are **0/2, 2/2 and 2/2**. Every classification agrees across both physical contexts, all three timesteps and both reference solvers. These fractions describe nested recipients in three histories, not population success probabilities.

Cells 19/29 retain their detected response differences after local matching. Cells 20/30 show larger differences again immediately after matching, while the native arms on those same frozen contexts remain below screen. Their matching-reference choice matters: E is roughly 0.10 with the measured S pair and 0.049 with the R pair. The classification is consistent across references; the magnitude is not reference-independent. This does not overturn the accepted below-screen delayed response measurements for cells 20/30. It shows that a reset places them in a different response regime.

Matched-minus-native effects are computed from the **signed contrast waveforms**, then reduced to RMS. They are not obtained by subtracting scalar E values or interpreting an RMS reduction as an additive causal percentage. The largest such change is 0.10348. Physical-context and interaction contrasts remain descriptive, with maxima 0.004474 and 0.003788; no new primary detector is assigned to them.

All **72 labeled fine pulses**, representing 48 distinct physical pulse problems, recover toward their own continuing controls. Recipient and whole-network recovery times range from 0.45 to 5.40 units. Recovery requires every remaining observation to remain within the specified threshold, with at least 24 units of follow-up. Timesteps and reference methods give identical recovery times on the 0.15-unit observation grid. Maximum late whole-network log RMS from the continuing control is 4.31e-11. These are finite-window recoveries, not proofs of equilibrium or pulse-induced lasting identity.

The local reset does not erase the contextual endpoints. Under each particular chemical allocation and physical context, native, S-matched and R-matched controls return to the same saved final network state to within 6.67e-16 in maximum log concentration. Paired allocations restore different recipient states, including **17.54–20.02-fold activator separation in cell 20** and **18.16–20.61-fold in cell 30**, depending on endpoint geometry. Both species initially agree exactly in every matched allocation pair. Strong between-cell contrast remains throughout every control; the minimum whole-window log-activator SD is 1.3535. This is maintenance and recovery in mature tissue, not differentiation newly demonstrated from a zygote.

## Independent numerical and provenance checks

The archived review reconstructs saved arrays using literal NumPy/SciPy formulas, with no project assessment or integration helpers. A separate archived literal preparation auditor verifies the accepted parent's actual endpoint arrays, target-only resets, conserved graph structure, reassignment accounting and exact pulses. The earlier actual-start prefix/restart audit remains pinned and unchanged. Twelve physical blocks each carry twelve independent chemical columns, with no padding or communication between columns.

| Check | Independently reconstructed evidence |
|---|---|
| Scientific sources and inputs | 77/77 source pins and 331/331 input pins unchanged |
| Consumed pinned artifacts | 1,660 hashes verified |
| Logical versus distinct starts | 144 labeled starts; 96 distinct physical initial-value problems |
| Actual-start GPU/CPU/restart qualification | All 432 cases previously independently audited; evidence reverified |
| GPU paths / NumPy paths / reference trajectories | 432 / 432 / 288 |
| Production/reference/CPU accuracy rows | 432/432 pass |
| Raw timestep comparisons | 288/288 pass |
| Resolved order checks | 144/144 pass; observed orders 2.00056–2.00672 |
| Response comparison reports | 9/9 pass: 648 pulse rows and 1,026 signed contrast rows |
| Duplicate preparation checks | 240/240 pass; trajectories agree bit for bit |
| Maximum dual-reference raw log discrepancy | 4.81e-12; limit 1e-8 |
| Maximum GPU/NumPy raw log discrepancy | 6.66e-16; limit 1e-8 |
| Maximum production/reference raw log discrepancy | 7.21e-6; limit 0.01 |
| Maximum normalized response discrepancy | 1.15e-5; limit 0.01 |
| Maximum signed contrast / E discrepancy | 4.20e-6 / 2.24e-7; both limits 0.01 |
| Maximum relative response-area discrepancy | 7.85e-6; limit 0.02 |
| Maximum recovery-time discrepancy | 0.0; limit 0.3 |
| Recorded scalar comparisons | 12,486; maximum reconstruction difference 1.11e-16 |
| Full coordinator runtime | 301.52 seconds, or 5.03 minutes |

Every saved concentration trajectory is finite, positive and float64. Checkpoints agree with saved histories, final states, clocks, method metadata and committed hashes. Numerical corrections are finite float64 arrays with the required shape and magnitude; independent CPU final corrections also pass. Corrections are numerical rounding bookkeeping, not biological memory. Unsaved intermediate SSP-RK2 stages remain source/runtime evidence rather than independently reconstructed observations. These temporal checks do not establish spatial convergence or a many-cell limit.

The parent's **non-input manuscript PDF export exception** is preserved explicitly. Its broad historical audit found 2,181/2,183 entries unchanged; both differing entries refer to `output/pdf/embryogenesis_o1.pdf`. The original expected hash is `8e90238617f1e3f67d139a699cfc1d5ae404cd1a58a76509019d9a7159b1940b`; the export currently has `0c60e5a0ec48698e69800a883722f6875022153994235e85bb5adb3d0de032af`. This review rechecks both the current export and the original expected bytes preserved in `manuscript/main.pdf`. The export is not a study source or input. No PDF is rebuilt, restored or substituted. This exception prevents claiming all historical artifacts are unchanged; it does not affect the accepted numerical gates.

## Interpretation and next focused control

The supported conclusion is **a surrounding chemical context can alter response behavior despite equal starting recipient chemistry and a common measured transport geometry**. Matching local chemistry does not isolate a permanent intrinsic cell identity: the surrounding cells drive the recipient away from that matched state, and its response depends on the chosen reset state. The no-reset arms preserve selective delayed-style classifications; matching exposes additional transient response sensitivity in cells 20/30. This is consistent with network-supported, state-dependent behavior.

Current dosage remains a limitation. At the challenge, paired global activator amounts have R/S ratios of 0.9211–1.1000 and inhibitor amounts 0.9304–1.0954 across the design. Surrounding-only ratios range from 0.9218–1.1089 and 0.9358–1.1032. Thus these chemical-context differences include both spatial allocation and current total amount. They cannot yet be attributed specifically to allocation at equal current dosage, or to immediate neighbors rather than the whole network.

The next focused experiment should **match both surrounding species amounts at t=570**, while retaining equal recipient chemistry and the common measured operator/volumes. Include the present unmatched-budget and native arms, paired own resetting controls, explicit chemical interventions and unchanged actual-start/numerical qualification. That would distinguish allocation from current total dosage in this immediate assay. Initial budget matching would not maintain equal amounts during reactions or make the response a settled challenge. Moving continuation after the new reset, a settled challenge, nearest-neighbor specificity, fresh carry-corrected development, autonomous/inherited identity and spatial convergence remain separate questions. **No new test is launched by this review.**

Accepted parents, prospective protocol, simulation sources, results ledger and manuscript remain unchanged. The coordinator's historical pending-review flag is preserved; this separate assessment records acceptance.

[Prospective protocol](pulse_time_recipient_response.md), [launch/completion record](pulse_time_recipient_response_launch.md), [machine-readable assessment](pulse_time_recipient_response_assessment.json), [independent review](../archive/study-reviews/2026-10-09-pulse-time-recipient-completed/independent_review.json), [review source](../archive/study-reviews/2026-10-09-pulse-time-recipient-completed/independent_review.py), [36 primary comparisons](../archive/study-reviews/2026-10-09-pulse-time-recipient-completed/response-comparisons.csv), [24 signed matching changes](../archive/study-reviews/2026-10-09-pulse-time-recipient-completed/matching-changes.csv), [response/control/amount arrays](../archive/study-reviews/2026-10-09-pulse-time-recipient-completed/independent_review.npz).
