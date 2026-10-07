# Completed moving chemical-network reset with phase carry

**All 42 paths pass the unchanged full-study criteria across three existing histories.** Resetting the chemistry outside a preserved recipient changes its maintained chemical state in both selected cells of every history. Immediate pulse responses change above the prespecified screen in one of two recipients in history 7 and both recipients in histories 8 and 9. Both recipients remain closer to their donor response reference after the reset in every history. All tested pulses recover toward their own moving controls.

These results support **maintained, context-sensitive chemical phenotypes with a resilient response tendency associated with prior chemical preparation**. They do not establish autonomous identity or show that the tested surrounding state is necessary for keeping the donor-nearer response.

## What was changed

The [unchanged protocol](phase_carry_network_context.md) uses fine exchanged endpoints at model time 450 in mature, nondividing 16-cell histories 7, 8 and 9. For each selected recipient, its two chemical concentrations and every cell's initial physical state are preserved exactly. Every other cell is reset once to its same-ID chemistry from the matched untouched endpoint. Geometry, polarity, mechanics, conservative signaling transport and dilution then continue to evolve freely through model time 510.

Each reset has an unperturbed moving control and an immediate -10% activator pulse in the preserved recipient. A sham leaves the exchanged chemistry unchanged and shares its unperturbed control between the two recipient pulses. Recipient pulses match both their initial fraction and absolute removed amount between sham and reset. The reset itself externally adds or removes separately recorded chemical amounts; it is not a conservative chemical transplant.

There are **33 new paths and nine exact reused fine sham controls/pulses**, with dt=0.00375 and 0.001875 starting from the same fine physical preparation. Cell targets, contexts, pulses and timesteps are nested interventions. The replication unit is **three existing histories**, not 42 independent embryos or six independent recipients. No new zygote history is added. Model time has not been calibrated to developmental hours.

## Results by history

| History | Recipients | Late state shifts | Immediate response shifts | Donor-nearer after reset | Fine pulse/control recoveries |
| --- | --- | --- | --- | --- | --- |
| 7 | 22 / 23 | 2/2 | 1/2 | 2/2 | All four |
| 8 | 19 / 20 | 2/2 | 2/2 | 2/2 | All four |
| 9 | 29 / 30 | 2/2 | 2/2 | 2/2 | All four |

The state and response screens are both 0.01 in their respective log-based measures. They are operational effect screens, not biological cell-type definitions. The two timestep levels agree on every classification. All 42 paths retain strong across-cell chemical contrast for the whole 60-unit window; the minimum SD of log activator is **1.3246**, above the 0.1 contrast criterion.

![Maintained state shifts, immediate response shifts, and reference distances](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/network-context.png)

Left: recipient state deviation from its sham control, with the final 24 units shaded. The logarithmic axis omits the exactly zero initial deviation. Middle: full-window reset-versus-sham response differences, each pulse measured relative to its own background. Right: post-reset distances from the same fixed fine donor/destination response references. Smaller distance determines the nearer reference; this is not an equivalence test. Fine data are shown; the coarse/fine checks are separate acceptance evidence. [Editable SVG](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/network-context.svg).

## Maintained state and immediate behavior are separate findings

State deviation is the RMS of the recipient's activator and inhibitor log differences between the reset control and sham control. The prespecified screen uses its **maximum over the final 24 units**. As an additional descriptive check, its minimum over that window also exceeds 0.01 for all six recipients; none of these state effects depends on one isolated late sample. This does not establish an indefinitely sustained plateau or a stationary equilibrium.

For a pulse, the signed response of each species is its log pulse/control concentration ratio divided by the magnitude of the initial log pulse, abs(ln(0.9)). Response change is the full-window time RMS of the difference between reset and sham responses, with both species included equally. The table reports fine-path values.

| History | Recipient | Final state log RMS | Full response-change RMS | Final activator reset/sham | Activator-response area reset/sham |
| --- | --- | ---: | ---: | ---: | ---: |
| 7 | 22 | 0.01477 | 0.011256 | 1.01550 | 1.055 |
| 7 | 23 | 0.02485 | 0.000621 | 1.02741 | 1.004 |
| 8 | 19 | 0.13250 | 0.071299 | 0.87049 | 1.387 |
| 8 | 20 | 1.93987 | 0.084763 | 0.07361 | 2.765 |
| 9 | 29 | 0.13381 | 0.071075 | 0.86931 | 1.377 |
| 9 | 30 | 2.11811 | 0.091500 | 0.05713 | 2.935 |

The largest changes occur in recipients 20 and 30: activator ends at **7.36% and 5.71%** of the sham level, approximately **13.6-fold and 17.5-fold lower**. Their inhibitor levels are 42.81% and 41.37% of sham. Their activator-response areas increase by factors of 2.77 and 2.94. A state log RMS of 1.94 or 2.12 is not a fold change; the species-specific ratios provide that interpretation.

History 7 gives smaller effects: final activator differs by about +1.55% in cell 22 and +2.74% in cell 23. Cell 23's immediate response change remains below the prespecified screen, despite its maintained state shift. This negative result is retained.

All twelve distinct fine pulse/control comparisons—six sham and six reset comparisons—recover under the declared 10%-of-initial-displacement rule, with at least 24 later observed units. Target recovery times range from 2.7 to 5.55; network recovery ranges from 2.7 to 6.45. Recovery means convergence toward that background's own control, not recovery of the original sham chemical state. Pulse recovery and a lasting difference between the backgrounds can coexist.

**The current pulses challenge the initial reorganization period.** A post hoc decomposition places **99.30%–99.98%** of the time-integrated squared reset-minus-sham response difference within the first six units. The backgrounds still have different late chemical states, but these immediate pulse trajectories do not establish how a new pulse would behave after those states develop. This decomposition is descriptive and does not replace the prespecified full-window response criterion.

## Numerical and provenance acceptance

The independent CPU review verified **47 scientific source hashes, 764 input hashes, and 1,292 evidence/source paths before and after review**. It checked all completed histories/checkpoints, source-array and metadata preservation during retiming/reset/pulse preparation, recipient preservation, actual reset and pulse amounts, native/GPU prefix histories and fields, volume-weighted graph conservation/spectra, pilot/full snapshots, and raw timestep differences. State shifts, waveforms, area and recovery metrics were independently recomputed from the raw histories. No scientific simulation was rerun.

All 33 new context gates and nine exact reused parent gates remain valid. Their recorded eight independent chemical steps and four exact carry-restart steps are verified from their original hashed evidence, rather than rerun in this review. Three pilot comparisons and three full comparisons pass; the full comparisons contain 21 mechanical pairs and twelve response pairs.

| Full-window coarse/fine quantity | Maximum discrepancy | Unchanged limit |
| --- | ---: | ---: |
| Log chemistry | 1.86e-05 | 0.01 |
| Polarity | 2.97e-06 | 0.01 |
| Relative axis ratio | 6.55e-08 | 0.01 |
| Relative cell volume | 9.72e-07 | 0.005 |
| Relative transport matrix | 0.0024 | 0.01 |
| Final phase field | 1.67e-06 | 0.02 |
| Pulse-normalized response, all cells/species/times | 2.41e-05 | 0.01 |
| Relative target activator area | 1.44e-05 | 0.02 |
| Recovery-time discrepancy | 0 | 0.3 |

Maximum target-state-shift disagreement is 2.1e-07; maximum response-shift disagreement is 3.53e-07. The smallest positive response-effect margin is approximately 0.00126, well above its numerical discrepancy. Every reset response remains donor-nearer, with donor/destination distance margins of 0.0731–0.1716; nearest-reference signs are not marginal at the tested timesteps.

Recorded every-step physical extrema pass: target-volume error **1.117%** against 5%, minimum equivalent radius **5.107 grid spacings** against four, zero clipping, maximum sampled boundary occupancy **1.51e-08**, and dilution amount-conversion error **4.44e-16** against 2e-14. Sparse saved histories cannot independently reconstruct every-step extrema. Reconstructed weighted-mass conservation error is 8.33e-17; graph-spectrum reconstruction error is 1.78e-15.

No recorded graph supports a growing spatial mode about uniform chemistry during these maintained-state trajectories. This is consistent with the initiation/maintenance distinction; it does not qualify a new frozen-endpoint coexistence assay.

The new moving-path audits total **3.68 sequential GPU hours**, excluding reused parent paths and additional coordinator/context-check overhead. The relaunch preserved eight completed pilots and restarted one exactly source-equivalent time-zero state; [the operational record](phase_carry_network_context_relaunch.md) preserves the interrupted artifacts. This did not change the scientific implementation or thresholds.

## Interpretation and limits

Changing other cells' chemical preparation changes a recipient that was initially left untouched. This supports a causal effect of the **collective chemical-context intervention within the coupled model**. The subsequent route can include diffusion, chemistry-dependent mechanics and changed geometry; the assay does not isolate those mediators.

At the same time, the donor reference remains nearer after every reset. Thus the test shows context sensitivity alongside resilience to this challenge. It does not show that the original surroundings are necessary for the response tendency, that the recipient is autonomous, or that its response equals the donor's. The network remains present and coupled in all arms.

The reset changes both spatial chemical distribution and total amounts: activator changes range from about -7.88% to +11.38% of the sham total, and inhibitor from -4.00% to +10.39%. No matched bulk-dosage control was included. Attribution specifically to the spatial neighborhood rather than dosage requires a separate amount-controlled intervention. Every other cell was reset, not only immediate contact neighbors.

This qualifies maintained chemical-state and immediate negative-pulse behavior at the tested mature point and time window. It does not establish discrete biological types, gene-regulatory commitment, inheritance through division, fresh carry-corrected differentiation from a zygote, indefinitely persistent identity, general parameter/spatial robustness, or a continuum limit of the geometric transport approximation. The original no-carry ledger/manuscript snapshot and older positive-pulse, pre-relaxed, reservoir and neighbor claims are not retroactively accepted.

## Next test

Apply **delayed pulses from the actual t=510 sham/reset endpoints**, with each endpoint's own continuing moving control and the same native/carry and paired timestep gates. This tests behavior after the initial reorganization period, without assuming the moving state is stationary.

Include both a -10% fractional pulse and a **common feasible removed amount**. A suitable shared removal is 10% of the smaller initial target activator amount in the two delayed backgrounds; use their actual measured volumes. Reusing the original absolute dose could make a strongly depleted recipient negative. Normalize each response by its actual log pulse and record doses separately.

Then test an amount-preserving surrounding reset and/or a matched bulk-dosage rescaling to separate chemical distribution from dosage. These are proposed follow-ups, not simulations launched by this review. Preserve the six recipient comparisons within their three histories when reporting replication.

## Evidence

- [Unchanged protocol](../outputs/phase-carry-network-context/protocol.json), SHA-256 `155d8dfd673cb3d68630297e6e9a2fdb531607140f9a9340fc06fb6556466278`.
- [Original completed summary](../outputs/phase-carry-network-context/summary.json) and [status](../outputs/phase-carry-network-context/status.json).
- [Independent review and saved assessment script](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/assessment.json).
- [Compact machine-readable assessment](phase_carry_network_context_assessment.json).
- [Recipient table](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/recipient-summary.csv) and [plot data](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/plot-data.json).
- [Figure provenance](../outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/figure-provenance.json).
- [Pre-review reporting snapshots](../archive/study-document-snapshots/2026-10-06-before-completed-carry-network-context-review/manifest.json).

Original scientific sources, tests, kernels, protocol document, preparation records, trajectories, checkpoints, earlier assessments, ledger and manuscript remain unchanged. README and the working plan now point to this completed review.
