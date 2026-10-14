# Moving challenge-time amount control: interim assessment

**Historical snapshot:** the study is now complete and [independently accepted](moving_pulse_time_amount_response_assessment.md). The partial assessment below is preserved as recorded; use the completed review for current results.

Latest snapshot: 2026-10-09 23:36 UTC / 19:36 America/New_York. This is a partial completed-path review, not acceptance of the whole study. [Prospective protocol](moving_pulse_time_amount_response.md), [launch record](moving_pulse_time_amount_response_launch.md), [literal saved-array reconstruction](../archive/study-reviews/2026-10-09-2335-moving-pulse-time-amount-interim/independent_interim.py), [machine-readable assessment](../archive/study-reviews/2026-10-09-2335-moving-pulse-time-amount-interim/interim.json). The [earlier 16:07 assessment](../archive/study-reviews/2026-10-09-2335-moving-pulse-time-amount-interim/prior-interim-assessment.md) is preserved, along with its original script/JSON archive.

## Progress

The coordinator is active on the qualified GTX 1080 Ti. **30/48 full paths have completed**, all 48 have completed their six-unit pilots, and approximately 66.5% of the total model-time work is recorded. The active path is history 8, cell 20, R allocation, fine timestep, unpulsed control. Seven of twelve full background timestep reports pass. Three of six recipients have complete eight-path S/R comparisons at both timesteps: cells 22, 23 and 19. Cell 20's R fine paths and all history-9 full paths remain unfinished.

All 48 native/GPU start checks and the separate prefix reconstruction passed, followed by all 48 pilots, twelve pilot background reports and six signed pilot contrasts. This review verifies all 83 source and 908 input pins, complete-path checkpoint/history hashes, and literal pulse-response and recovery formulas for the seven completed paired backgrounds. The whole study remains subject to a final independent review.

## Full-window findings available now

E is the time-and-species RMS difference between S/R recipient pulse responses, each normalized by the actual log pulse relative to its own moving control. The primary window is 60 units and the operational screen is E>0.01. The table uses the fine timestep and compares the exact frozen M=S/B=S/P=S reference:

| History | Recipient | Moving E | Frozen E | Relative to screen | Descriptive E change |
|---|---:|---:|---:|---|---:|
| 7 | 22 | 0.019185 | 0.019203 | Above | -0.09% |
| 7 | 23 | 0.002258 | 0.002305 | Below | -2.05% |
| 8 | 19 | 0.038318 | 0.037593 | Above | +1.93% |

Both timesteps agree on all three classifications. Maximum signed waveform discrepancy across these completed S/R comparisons is 2.37e-6; maximum RMS discrepancy is 2.37e-7, each against 0.01. **The newly detected cell-22 response difference survives 60 units of live mechanics.** Cell 23 remains below screen and cell 19 remains above. Moving E differs from frozen E by about −0.09%, −2.05% and +1.93%, respectively; these are descriptive changes of RMS, not additive causal fractions attributable to mechanics.

All seven completed fine pulse paths recover toward their own continuing controls under the predeclared sustained threshold, with target/network recovery times spanning 2.70–4.80 units and at least 24 subsequent units verified. Final whole-network log RMS from controls is nonzero (up to 5.03e-5), so operational recovery is not exact equality or proof that a mechanical perturbation leaves no residual effect.

The remaining recipients still have their earlier pilot results only: cell 20 and cells 29/30 are above screen over six units. Pilot RMS values must not be compared directly with 60-unit RMS or treated as completed primary outcomes. History 7 now has a complete 1/2 result; history 8 has one completed above-screen recipient and one pending; history 9 remains pending. These are nested measurements in **three existing histories, zero new histories**. They support immediate allocation-sensitive behavior conditional on mature preparations, not inherited/autonomous biological identity or fresh-zygote organization.

## Numerical health

All available full background reports pass. Maximum raw chemical timestep discrepancy is 4.24e-6 against 0.01; normalized response discrepancy is 1.68e-5 against 0.01; relative response-area discrepancy is 1.37e-5 against 0.02. Recovery times agree on the 0.15-unit observation grid.

Current recorded quality extrema: maximum volume error 1.12% against 5%; minimum cell radius 5.107 grid intervals against four; boundary occupancy 2.39e-8 against 0.01; dilution amount mismatch 4.44e-16 against 2e-14; clipping zero. No failed numerical or physical gate is observed. These quality records are not an independent mechanical replay of every timestep. About 62 GiB remains free, well above the storage reserve.

## Updated time estimate

At this snapshot, 5.46 hours have elapsed since launch. Remaining-time extrapolation uses every unfinished path's pilot speed, calibrated against measured full continuations and crediting all completed pilot work. It gives **2.55 hours remaining**, including a heuristic allowance for initialization, checkpointing and final paired checks. Central predicted finish is **22:09 New York time October 9** (02:09 UTC October 10). A practical planning window is approximately 22:00–22:30 if current throughput holds; further slowdown or a failed gate would change it. This is not a statistical confidence interval or a guarantee. Independent completed-array review follows simulation completion.

No running job, pinned scientific input, solver, threshold, ledger entry or manuscript conclusion was changed. Fresh-zygote development remains the next major priority after division-specific qualification.
