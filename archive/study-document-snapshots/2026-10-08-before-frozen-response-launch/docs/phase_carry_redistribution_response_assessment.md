# Completed delayed response after amount-preserving redistribution — 8 October 2026

**Changing the initial distribution of chemicals while preserving their initial amounts can alter a recipient's later response to a challenge. The effect is selective:** one of two selected recipients exceeds the declared response screen in histories 8 and 9, while neither exceeds it in history 7. Both pulse-dose designs and both timesteps agree. All tested pulses recover toward their own continuing controls. All 66 logical paths pass the unchanged numerical, physical and independently verified provenance gates.

These results extend the maintained-state finding to later response behavior in two recipients. They also retain an informative limitation: large chemical-level changes can coexist with below-screen response changes. The independent units remain **three existing developmental histories**, with two selected recipients nested in each.

**What was tested.** The [frozen protocol](phase_carry_redistribution_response.md) starts from actual fine t=510 endpoints of the [completed distribution/dosage study](phase_carry_context_dosage_assessment.md). Its original intervention occurred at **t=450**: the recipient's concentrations were preserved, while concentrations in **every other cell** were redistributed and rescaled separately for activator and inhibitor to preserve the aggregate's initial species amounts. This was a one-time intervention without a reservoir or continuing clamp. It affected all nonrecipients, rather than isolating nearest neighbors.

The backgrounds were allowed to reorganize for 60 units before the new challenge. Their amounts are **not matched at t=510**: redistribution/sham aggregate amount ratios then range from 0.922 to 1.086 for activator and 0.932 to 1.082 for inhibitor. Chemistry, volumes, geometry and polarity can all differ by challenge time.

Continue the actual backgrounds to t=570 with live chemistry, mechanics, polarity, conservative transport and geometric dilution. Each has its own unperturbed moving control. Apply either a 10% recipient activator depletion, or the same feasible removed amount in both backgrounds, equal to 10% of their smaller target amount. Inhibitor and other cells remain unchanged at the pulse. Both timesteps begin from the identical fine preparation within each background.

There are **66 logical paths: 48 new paths and 18 exactly reused sham control/fractional paths**, including 18 controls and 48 pulses. No new zygote history or division is included. The tested diffusivity ratio is D_b/D_a=27.5 and polarity-tension contrast chi=0.35 on a 72-cubed grid. Phase fields retain float64 update-residual carry and chemistry remains float64. Model time and chemical quantities remain dimensionless.

**Primary delayed response.** Each signed two-species response is the pulse/own-control log concentration ratio divided by that pulse's actual initial log amplitude. The effect is the full 60-unit time RMS difference between redistributed and sham target responses, weighting activator and inhibitor equally. The declared screen is **strictly greater than 0.01**. It is an operational comparison, rather than a biological identity boundary, equivalence test or significance test.

| History | Recipient | Fractional-pulse RMS | Common-amount RMS | Above screen at both doses | Final activator redistribution/sham |
| --- | ---: | ---: | ---: | --- | ---: |
| 7 | 22 | 0.009352 | 0.009578 | No | 1.01249 |
| 7 | 23 | 0.000052 | 0.000101 | No | 1.02499 |
| 8 | 19 | 0.072664 | 0.073146 | Yes | 0.88294 |
| 8 | 20 | 0.006403 | 0.007979 | No | 0.06134 |
| 9 | 29 | 0.072965 | 0.073469 | Yes | 0.88064 |
| 9 | 30 | 0.006645 | 0.008251 | No | 0.05898 |

Thus the history-level counts are **0/2, 1/2 and 1/2**, at either dose. Recipient 22 is close to the screen, but its smaller below-screen margin is 0.000422, compared with a coarse/fine effect discrepancy of 1.72e-08 for that dose. Its classification is numerically resolved under this protocol; this does not establish equivalence of its responses.

![Maintained chemical levels and selective delayed response differences](../archive/study-reviews/2026-10-08-redistribution-response/redistribution-response.png)

[Editable figure](../archive/study-reviews/2026-10-08-redistribution-response/redistribution-response.svg), with [raw-derived plot data](../archive/study-reviews/2026-10-08-redistribution-response/plot-data.json).

**Chemical levels and response behavior remain distinct attributes.** Activator in cells 20 and 30 ends approximately **16.3-fold and 17.0-fold lower** than sham, despite their below-screen waveform effects. Their late-window two-species baseline state log RMS values are 2.070 and 2.097. In contrast, recipients 19 and 29 have much smaller maintained concentration shifts but clearly different delayed responses. Their normalized activator-response areas under fractional pulses increase by about 33%. Defining a phenotype from concentration alone would miss this distinction.

All 48 pulse paths recover, with matching coarse/fine recovery classifications and times. Target recovery occurs within 2.25–4.35 units; network recovery within 2.25–5.40 units. Recovery requires at least 24 subsequent observed units inside the declared 10%-of-initial-displacement band. It means return toward each background's continuing control, rather than restoration of the sham background or exact equality.

The six-unit pilots have five above-screen recipient comparisons, whereas the full-window primary assessment has two. This is expected from the **prespecified window dependence of time RMS**, not a numerical failure: the response differences mostly occur early and then decay. Approximately 99.75%–99.98% of their integrated squared differences lies within the first six units. The pilot qualifies numerics; it cannot replace the declared full-window outcome or justify changing the screen.

**Numerical assessment.** Independent calculations reconstruct all 66 raw histories, 33 full mechanical pairs, 24 pulse-response refinement comparisons and the twelve full-window recipient/dose effects at each timestep. All reconstructed gates and saved decisions agree. Maximum full-window discrepancies are:

| Coarse/fine quantity | Maximum | Declared limit |
| --- | ---: | ---: |
| Log concentration | 3.78e-04 | 0.01 |
| Polarity | 1.47e-07 | 0.01 |
| Relative axis ratio | 2.54e-08 | 0.01 |
| Relative volume | 4.69e-07 | 0.005 |
| Relative transport | 5.32e-07 | 0.01 |
| Final phase field | 1.19e-07 | 0.02 |
| Pulse-normalized response, all cells/species/times | 1.76e-05 | 0.01 |
| Relative activator-response area | 1.42e-05 | 0.02 |
| Recovery time | 0 | 0.3 |

Maximum primary-effect disagreement is 4.21e-07. Maximum recorded cell-volume error relative to prescribed target volume is 1.102%, minimum equivalent radius 5.107 grid spacings, clipping zero, sampled boundary occupancy 1.93e-08 and dilution conversion error 4.44e-16. Pulse accounting agrees within 4.22e-15 relative error. Across-cell contrast stays above its 0.1 screen at every observation; the minimum SD of log activator is 1.355.

These checks independently reconstruct sampled quantities and endpoint volumes. **Internal-step extrema remain cumulative runtime audit records**, whose consistency and hashes are checked; saved observations cannot reproduce every unsaved step. Conservative-transport conservation checks validate the implemented operator on saved data, rather than establish the physical accuracy of the diffuse-interface conductance approximation.

**Interpretation and next test.** The same recipients, 19 and 29, exceeded the delayed screen after the earlier raw reset. They still do so after an initial amount-preserving redistribution. Thus the raw reset's change in total initial chemical amount is **not required** for the detected delayed effects at this tested point. This supports selective context sensitivity arising after redistribution. It does not demonstrate universal behavioral changes, autonomous or inherited identity, or independence from challenge-time chemistry and geometry.

Next, freeze each **actual t=510 transport operator and measured volumes** and repeat its own-control pulse comparisons with continuing chemical reactions, both original pulse doses, and independent chemical-solver checks. This tests whether ongoing mechanical evolution during the challenge is necessary for the moving assay's response differences. Each frozen background still inherits its distinct chemistry, graph and volumes, so persistence under freezing would not separate those initial differences. A shared transport operator would further isolate chemical-state differences from inherited physical-context differences; recipient-state matching would separately test surrounding-network sensitivity at matched recipient chemistry. These require distinct future controls. These mature continuations inherit earlier no-carry development and do not establish fresh carry-corrected zygote differentiation, inherited identity, or general spatial/transport convergence. The current review launches no further simulation.

**Evidence.** The [response reconstruction](../archive/study-reviews/2026-10-08-redistribution-response/response_reconstruction.json) and [numerical verification](../archive/study-reviews/2026-10-08-redistribution-response/numerical_verification.json) include saved scripts and provenance. The protocol SHA-256 is `5574abd1b40bb711676685611cb84b511bdb342a1c377345bbef5d752d78ccc1`. The [independent provenance audit](../archive/study-reviews/2026-10-08-redistribution-response/provenance_verification.json) verifies 53 scientific source hashes, 2,412 input hashes, two pinned runtime binaries, all 66 exact preparations/context gates, and 3,091 immutable paths before and after review. The user-authorized experimental-literature manuscript revision after launch is recorded separately from the frozen scientific pins; it does not change this assay's equations or outcomes. No manuscript file is a frozen source/input pin for this assay. Prior ledger and manuscript claims remain unchanged pending separate promotion.


The [compact assessment](phase_carry_redistribution_response_assessment.json), [six-recipient CSV](../archive/study-reviews/2026-10-08-redistribution-response/recipient-summary.csv), and [pre-review reporting snapshot](../archive/study-document-snapshots/2026-10-08-before-completed-redistribution-response-review/manifest.json) preserve the decisions and reporting provenance. Scientific sources, kernels, tests, protocols, preparations, histories, checkpoints, prior assessments, ledger and manuscript remain unchanged by this review. README and the working plan point to the completed assessment; no new simulation is launched.
