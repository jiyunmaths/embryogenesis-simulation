# Completed delayed moving-response assay — 7 October 2026

**All 66 paths pass the unchanged numerical and physical criteria across three existing histories.** New pulses delivered 60 model units after the surrounding-cell reset still reveal a response change in one of the two selected recipients in histories 8 and 9. Neither recipient in history 7 exceeds the prespecified full-window response-effect screen. The classifications agree at both timesteps and with both pulse-dose designs. All tested pulses recover toward their own moving controls.

The result supports **selective context sensitivity of later responses**, alongside maintained differences in chemical levels. It narrows the earlier immediate-response finding: five recipient comparisons exceeded the immediate screen, while two exceed the delayed screen. Below-screen responses are not evidence of exact equivalence. This is a finite-horizon assay of continuous chemical phenotypes, not a new biological identity label or an autonomous/inherited fate.

## What was tested

The [unchanged protocol](phase_carry_delayed_response.md) begins at the actual fine **t=510** endpoints of the completed moving reset assay. The original reset occurred at t=450. Chemistry, mechanics, polarity, conservative transport, measured volumes and dilution continue together to **t=570**, with no further reset and no division. Each background has its own unperturbed moving control.

Compare two negative activator challenges in each selected recipient:

- **Fractional pulse:** remove 10% of that background's initial target activator amount.
- **Matched amount:** remove the same feasible amount in sham and reset, equal to 10% of the smaller of their target amounts. Each response is normalized by its own actual initial log pulse, rather than assuming both received a 10% concentration change.

Inhibitor and other cells are initially unchanged by the pulse. The common-amount dose can be small in the activator-rich background: the smallest fractional depletion is about 0.57%. Neither dose design controls the chemical amounts changed by the original surrounding-cell reset.

There are **66 new paths: 18 controls and 48 pulse paths**, at dt=0.00375 and 0.001875. Both timesteps start from the same fine physical preparation within each background. There are six selected recipient comparisons nested in **three existing histories**, with **no new developmental histories**. Model time is dimensionless and has not been calibrated to biological hours.

## Delayed response effects

The signed target response is the log pulse/control concentration ratio for each species, divided by the magnitude of that pulse's initial log amplitude. The reported effect is the time RMS of the reset-minus-sham response difference over the full 60-unit window, including both species equally. The prespecified screen is **strictly greater than 0.01**, an operational effect screen rather than a cell-type boundary or statistical significance test.

| History | Recipient | Fractional-pulse RMS | Matched-amount RMS | Above 0.01 at both doses |
| --- | --- | ---: | ---: | --- |
| 7 | 22 | 0.009363 | 0.009588 | No |
| 7 | 23 | 0.000052 | 0.000101 | No |
| 8 | 19 | 0.072665 | 0.073148 | Yes |
| 8 | 20 | 0.006403 | 0.007979 | No |
| 9 | 29 | 0.072962 | 0.073466 | Yes |
| 9 | 30 | 0.006645 | 0.008251 | No |

Thus the history-level counts are **0/2, 1/2 and 1/2 recipients**, respectively. Cells 19 and 29 retain clear context-dependent delayed responses. Their normalized activator-response areas are about **33% larger** in reset than sham under either dose design. The other four comparisons stay below the full-window RMS screen; they can still have measurable area or recovery differences.

History 7 recipient 22 is especially useful as a retained negative outcome: its six-unit pilot RMS is about 0.030, but its full-window values remain below 0.01. The pilot is a numerical qualification window, not a replacement scientific effect criterion. The smaller effects in cells 20 and 30 likewise cross the screen in the short pilot but not over the declared full window. No time alignment or threshold adjustment was used.

![Maintained levels and selective delayed response shifts](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/delayed-response.png)

Left: final unperturbed reset/sham concentration ratios, showing activator and inhibitor separately. Right: full-window response effects from the earlier immediate assay and this new delayed assay. Fine data are shown without treating timesteps or recipient challenges as independent replicates. The immediate and delayed challenges occur at different model ages and physical contexts. [Editable SVG](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/delayed-response.svg).

## Chemical levels and behavior are different attributes

These state comparisons describe the unperturbed backgrounds during the new assay; they are not a newly added state-acceptance test. Recipient state log RMS combines activator and inhibitor concentration differences between reset and sham. Species-specific ratios provide the fold-change interpretation.

| History | Recipient | Final state log RMS | Final activator reset/sham | Final inhibitor reset/sham | Fractional response-area reset/sham |
| --- | --- | ---: | ---: | ---: | ---: |
| 7 | 22 | 0.01192 | 1.01250 | 1.01145 | 0.951 |
| 7 | 23 | 0.02273 | 1.02506 | 1.02072 | 1.001 |
| 8 | 19 | 0.11896 | 0.88294 | 0.89300 | 1.328 |
| 8 | 20 | 2.06098 | 0.06134 | 0.43208 | 0.870 |
| 9 | 29 | 0.12146 | 0.88064 | 0.89089 | 1.325 |
| 9 | 30 | 2.09223 | 0.05898 | 0.42235 | 0.865 |

Activator in recipients 20 and 30 ends approximately **16.3-fold and 17.0-fold lower** than sham, yet their delayed waveform effects are below 0.01. In recipients 19 and 29, much smaller baseline concentration changes accompany the clearly detected delayed effects. Chemical abundance and perturbation-response behavior therefore should remain separate attributes when defining a model phenotype. The data do not justify assigning an identity from activator concentration alone.

All 24 distinct fine pulse/control comparisons recover; all corresponding coarse paths agree. Target recovery ranges from **2.25 to 4.35** model units and network recovery from **2.25 to 5.40**, with at least 24 subsequent observed units within the declared 10%-of-initial-displacement recovery band. Recovery means return toward each background's own moving control, not restoration of the original sham chemical state or exact equality.

As a descriptive decomposition, **99.75%–99.98%** of the integrated squared delayed response difference lies within the first six units of the *new pulse*. A later background can therefore retain a different response profile while the perturbation itself is transient and recovers. This decomposition does not replace the full-window screen or prove indefinitely maintained stationary states.

## Numerical and provenance review

The independent CPU audit verifies **49 scientific source hashes, 890 input hashes, all 66 preparations and context gates, and 1,811 evidence/source paths before and after review**. It checks physical-array/dtype and residual-carry preservation, exact retiming and random-state metadata, recipient-only pulses, dose accounting, completed histories/checkpoints, pilot/full snapshots, positivity, cell IDs, conservative graph validity and unchanged numerical decisions. State trajectories, normalized two-species waveform distances, response areas and recovery metrics were independently reconstructed from raw histories. No scientific simulation was rerun.

All three per-history pilot reports and all three full reports pass. Each stage contains 33 matched mechanical pairs and 24 matched response comparisons. Full-window maxima are:

| Coarse/fine quantity | Maximum discrepancy | Unchanged limit |
| --- | ---: | ---: |
| Log concentration | 1.69e-06 | 0.01 |
| Polarity | 1.47e-07 | 0.01 |
| Relative axis ratio | 2.54e-08 | 0.01 |
| Relative cell volume | 4.69e-07 | 0.005 |
| Relative transport matrix | 5.32e-07 | 0.01 |
| Final phase field | 1.19e-07 | 0.02 |
| Pulse-normalized response, all cells/species/times | 1.76e-05 | 0.01 |
| Relative activator-response area | 1.42e-05 | 0.02 |
| Recovery time | 0 | 0.3 |

Maximum response-effect RMS disagreement is **4.2e-07**. Both timesteps agree on all effect classifications. All paths retain strong across-cell contrast at every saved observation: minimum SD of ln(activator) is **1.3551**, above 0.1. Maximum recorded target-volume error is **1.102%**, minimum equivalent radius **5.107 grid spacings**, clipping zero, sampled boundary occupancy at most **1.93e-08**, and dilution amount-conversion error at most **4.44e-16**.

A separate direct reconstruction confirms all 24 full-network response-normalization errors exactly and checks conservation, symmetric volume-weighted Laplacians, spectra and uniform-state growth diagnostics across **26,466 saved observations**, with graph errors below 2.1e-15. This verifies the implemented operator on saved data, not the physical accuracy of its diffuse-interface conductance approximation. Recorded eight-step independent chemistry and four-step exact-restart gates were verified from hashed evidence, not rerun. Saved observations do not independently reconstruct every-step extrema.

## Interpretation and next test

The immediate context effect is partly specific to reorganization, but it is **not entirely restricted to that phase**: later challenges retain detected response changes in one recipient in each of two histories. Large maintained concentration differences and a strong delayed response effect need not occur in the same recipient. Perturbation recovery and a different response profile can coexist.

The original reset alters **chemical distribution and total chemical dosage**; the backgrounds also develop different geometry and polarity. Equalizing the later pulse amount rules out different removed pulse amounts as the sole explanation of the two detected effects, but does not isolate spatial organization, network necessity, chemical mediation or a polarity-specific cause.

The next priority is an **amount-preserving surrounding reset together with a bulk-dosage-matched control**. Keep the recipient and initial physical state unchanged; compare changing the distribution at fixed surrounding chemical totals with changing those totals while preserving the original relative distribution. Carry, positivity, actual-volume amount accounting, own controls and paired refinement must remain part of that test.

Autonomous identity, inheritance through division, fresh carry-corrected zygote development, general spatial/transport convergence and broader parameter robustness remain separate requirements. The earlier ledger/manuscript scope is preserved, and no follow-up simulation was launched by this review.

## Evidence

- [Protocol JSON](../outputs/phase-carry-delayed-response/protocol.json), SHA-256 `d488e3268d154795d8229746ffd4a8adf0c9f689633eff51e68b00a14199442d`.
- [Original completed status](../outputs/phase-carry-delayed-response/status.json) and [summary](../outputs/phase-carry-delayed-response/summary.json).
- [Independent full review and saved script](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/assessment.json), [supplemental reconstruction](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/supplement.json), and [compact assessment](phase_carry_delayed_response_assessment.json).
- [Recipient table](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/recipient-summary.csv), [plot data](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/plot-data.json), and [figure provenance](../outputs/phase-carry-delayed-response-review/2026-10-07T105523Z-final/figure-provenance.json).
- [Pre-review reporting snapshot](../archive/study-document-snapshots/2026-10-07-before-completed-delayed-response-review/manifest.json).

Scientific sources, kernels, tests, protocol documents, trajectories, checkpoints, preparation records, prior assessments, ledger and manuscript remain unchanged. README and the working plan now point to this completed assessment.
