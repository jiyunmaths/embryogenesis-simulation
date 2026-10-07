# Completed carry-corrected moving chemical exchange and response

**Exchanging chemical states transfers response tendencies across initially different geometric contexts in all three tested histories.** Both recipients are closer to the donor response than to their original response after 60 units of coupled chemical and mechanical evolution. All 48 paths pass the original numerical and physical gates.

The concentrations themselves adjust after exchange: all three exchanged pairs are classified as reorganized by the prespecified state-distance rule. Donor-nearer behavior therefore means a transferable response tendency, not preservation of an exact chemical state or conversion into a supplied cell type.

## Design and reporting units

The [unchanged protocol](phase_carry_exchange_response.md) starts from the qualified fine-timestep, chi=0.35 maintenance endpoints at model time 390 in existing histories 7, 8 and 9. These are mature, nondividing 16-cell preparations. They inherit original-method development; this experiment does not generate new carry-corrected zygote trajectories.

Within each history, select the largest log-activator contrast among pairs with at least the median standardized geometric-context distance. Context features are exposure, conservative exit rate and measured cell volume. Both species are exchanged with pairwise volume corrections that conserve each species' initial total amount. Geometry, polarity, cell IDs, lineage, RNGs and float64 phase carry are initially preserved. Reactions subsequently change chemical amounts.

The exchange stage runs from t=390 to 450 in untouched and freshly exchanged backgrounds. From each own completed endpoint, an unperturbed control and separate -10% activator pulses in both selected cells run from t=450 to 510. Pulse removal is an external intervention with recorded amount; it is not an amount-preserving exchange. Each pulse is compared with its own moving-background control and with same-age untouched donor and destination responses.

Both stages use dt=0.00375/0.001875 from identical initial physical states and inherited phase residuals. Resident PyTorch arrays/matrix operations and custom CUDA mechanics and spatial geometry/polarity remain active. Chemistry, mechanics, polarity, measured volumes, conservative conductances and dilution co-evolve. There is no fixed-network or reservoir override.

Parameters are beta=2, D_a=0.02, D_b=0.55, chi=0.35, c_gamma=0.25, c_A=0.35, grid 72 cubed, extent 2.24 and interface width 0.085. Model times and concentrations are nondimensional; response gains and waveform distances below are also dimensionless.

**Report replication as three existing developmental histories.** Twelve exchange paths and 36 response paths, two timesteps, selected cells and pulse interventions are nested measurements. The six recipient comparisons are not six independent histories. No new history or population success probability is added. The runner calls the exchange stage `formation`; here it measures exchange retention/relaxation, not fresh initiation.

## Chemical states persist but reorganize

Contrast S is the across-cell standard deviation of natural-log activator concentration. Every path remains above S=0.1 at every saved observation over its 60-unit window, including the final 24-unit window. Observation spacing is 0.15, so brief unobserved transients are not excluded.

| History | Selected IDs, initially low/high | Untouched final S | Exchanged final S | Destination state ratio | Transferred state ratio |
|---|---|---:|---:|---:|---:|
| 7 | 22 / 23 | 1.499060 | 1.612086 | 1.127193 | 0.170378 |
| 8 | 19 / 20 | 1.493240 | 1.695634 | 1.162363 | 0.212425 |
| 9 | 29 / 30 | 1.603563 | 1.777735 | 1.164199 | 0.212190 |

Fine-timestep t=450 values are shown. State ratios are maximum pairwise two-species volume-weighted log distances over the final 24 units, divided by the initial exchange separation. The destination is the untouched moving control. The transferred reference is the same-time conservative exchange of that control, using fixed t=390 masses for descriptive weighting. Actual current volumes govern transport and dilution.

The exchanged pairs stay closer to the transferred reference than to the destination reference, but neither ratio is below 0.1. All therefore meet the original **reorganized** label. This is a finite-horizon comparison with a moving reference, not a claim of stationary cell identity or a new endpoint bistability assay.

Context contrast is relative to this preparation: selection uses standardized features and an above-median pair distance. Both selected cells are highly exposed (about 98% on the model exposure diagnostic), with exposure differences of only 0.29–0.32 percentage points. Their conservative exit rates differ more clearly, but this is not an inner-versus-outer transplantation test or evidence across arbitrary environments.

| History | Initial exposure, low/high | Initial exit rate, low/high | Initial measured volume, low/high |
|---|---|---|---|
| 7 | 0.981707 / 0.984613 | 1.468677 / 1.247048 | 0.134769 / 0.136313 |
| 8 | 0.982236 / 0.979155 | 1.416039 / 1.611788 | 0.136467 / 0.134393 |
| 9 | 0.982981 / 0.979786 | 1.371306 / 1.553833 | 0.136016 / 0.135007 |

## Responses are donor-nearer in both exchange directions

For each species, the target response is its signed pulse/control log ratio divided by |ln(0.9)|. Waveform distance is the time RMS of the two-species difference over all 60 units. Donor and destination references must themselves differ by more than 0.01; all six comparisons are informative, with reference separation 0.168–0.173.

| History | Recipient receives from donor | Donor distance | Destination distance | Donor/destination ratio |
|---|---|---:|---:|---:|
| 7 | 22 receives from 23 | 0.087447 | 0.191150 | 0.457479 |
| 7 | 23 receives from 22 | 0.001814 | 0.173401 | 0.010461 |
| 8 | 19 receives from 20 | 0.039443 | 0.169647 | 0.232500 |
| 8 | 20 receives from 19 | 0.003380 | 0.167877 | 0.020135 |
| 9 | 29 receives from 30 | 0.041421 | 0.170164 | 0.243421 |
| 9 | 30 receives from 29 | 0.003399 | 0.168616 | 0.020157 |

Both recipients are donor-nearer in each history and at both timesteps. The high-state recipients have larger donor distances than the low-state recipients, especially in history 7. This is evidence against exact response equivalence: the new coupled context still modifies the waveform.

![Carry-corrected moving exchange response waveforms](../outputs/phase-carry-exchange-response-review/2026-10-06T105823Z-final/carry_exchange_responses.png)

Solid curves show activator; dashed curves show inhibitor. Orange is the exchanged recipient, blue the untouched donor, gray the untouched destination. The first ten elapsed units are displayed; the reported distances use the full 60-unit observation. [PDF figure](../outputs/phase-carry-exchange-response-review/2026-10-06T105823Z-final/carry_exchange_responses.pdf).

A -10% activator pulse produces a substantial inhibitor response in recipients that received the high state: peak normalized inhibitor magnitude 1.47–1.64, versus about 0.003 in the original low-state destinations. Recipients that received the low state have peak inhibitor magnitude 0.00079–0.00134, versus 1.62–1.72 in the original high-state destinations. The transfer changes behavior as well as concentration.

Every pulse returns within the declared 10%-of-initial-displacement tolerance relative to its own evolving control, with at least 24 later units observed. Across the twelve fine-level pulse paths, target recovery is sampled at 2.7–5.55 elapsed units and whole-network recovery at 2.7–6.15. Recovery from this small challenge is compatible with persistence of the reorganized background; it is not reversal to the original unexchanged state or proof of permanent commitment.

The uniform-state spectrum has no unstable spatial modes at any saved observation in these 48 paths. Finite-amplitude differences continue despite that absence. This concerns initiation about uniform chemistry, not the stability of the complete moving coupled system.

## Numerical and physical qualification

All 48 context-specific checks pass: strict native/carry prefixes, eight independent NumPy chemical steps, amount conversion and four exact carry-restart steps. All six grouped pilot and six grouped full-stage reports pass, covering 24 coarse/fine path comparisons in each window. Responses receive an additional pulse-minus-control comparison; there is no time alignment or threshold relaxation.

| Full-window diagnostic | Maximum measured error | Original limit |
|---|---:|---:|
| Raw log-concentration discrepancy | 7.1649388e-05 | 0.01 |
| Pulse-normalized response discrepancy | 0.00067328766 | 0.01 |
| Relative activator-response area discrepancy | 2.2455748e-05 | 0.02 |
| Recovery-time difference | 0 | 0.3 |
| Exchange state-ratio difference | 6.0759894e-07 | 0.05 |

Polarity, axis, volume, conductance and final field discrepancies also pass. Across all 48 paths, target-volume error is at most **1.1218%** against 5%, minimum equivalent radius is **5.10680 grid spacings** against four, clipping is zero, sampled boundary occupancy is at most **1.52e-08** against 0.01, and recorded dilution amount-conversion error is at most **4.44e-16** against 2e-14.

The independent CPU review verifies 45 source hashes and 572 pinned input hashes, plus all completed evidence. It restores checkpoints, verifies carry and nonchemical handoff preservation, independently ranks selection candidates, reconstructs conservative spectra and mass balance, checks saved native prefix trajectories/fields, and recomputes state distances, waveforms, areas and recovery. All original evidence hashes remain unchanged.

The original hashed records establish the eight-step independent chemistry and four-step exact GPU restart checks; these were not rerun during review. Every-step physical audit extrema are verified from checkpoint/result records and cannot be reconstructed from sparse observations. Pilot field snapshots are hash-verified; full paired fields are also compared with their completed checkpoints.

The measured moving-path runtime totals 5.69 sequential GPU hours; coordinator elapsed time including response-context checks is 6.75 hours. This review launches no scientific simulation.

## Interpretation and remaining limits

The result supports **persistent, transferable chemical/response differences within a co-evolving network**. Initial geometric position does not uniquely fix the later state or response: changing chemical preparation while preserving the initial physical state changes both. This strengthens the behavioral part of the project's operational identity question.

It does not establish autonomous identity, discrete biological cell types, gene-regulatory commitment, inheritance through division, indefinite persistence or fresh identity formation from a carry-corrected zygote. Donor-nearer responses can arise from nonlinear chemical dynamics supported by the surrounding network. Local context also changes the response. Spatial/developmental convergence, other parameter values and the general geometric transport approximation remain separate questions.

This qualifies only fresh two-species exchange and -10% activator responses in the new t=390–510 window at the tested point. It does not retroactively accept original t=150–270 amplitudes, positive pulses, pre-relaxed exchanges, neighbor/reset or reservoir results, the 262.1-fold moving/frozen ratio or the full original ledger/manuscript snapshot.

## Next scientific priority

Test network support with a carry-corrected moving neighbor-context/reset control. Preserve one recipient's chemical state and the full initial physical state, while resetting other cells' chemistry to the matched untouched context. Compare with an unchanged-neighbor control and challenge each resulting background against its own unperturbed control. Account explicitly for chemical amounts added or removed by the reset. This asks whether donor-like response tendencies survive a changed chemical neighborhood; it distinguishes transferable state from network dependence without claiming autonomy from a single reset.

Use these qualified endpoint preparations and require their own native/carry, restart, physical and timestep gates. Delayed-pulse and additional-amplitude tests can follow; broaden parameters and return to fresh carry development only as separately qualified experiments. No follow-up is launched by this review.

## Evidence

- [Original protocol](../outputs/phase-carry-exchange-response/protocol.json), SHA-256 `c8fcbfc071b094b8cce292fdd06c079cff6c50176750b860b8276a4e6903ab27`.
- [Completed original summary](../outputs/phase-carry-exchange-response/summary.json).
- [Read-only independent review and saved script](../outputs/phase-carry-exchange-response-review/2026-10-06T105823Z-final/assessment.json).
- [Compact machine-readable assessment](phase_carry_exchange_response_assessment.json).
- [Figure provenance](../outputs/phase-carry-exchange-response-review/2026-10-06T105823Z-final/figure-provenance.json).
- [Pre-review reporting snapshots](../archive/study-document-snapshots/2026-10-06-before-completed-carry-exchange-response-review/manifest.json).

The original runner, tests, kernels, pinned protocol document, preparation records, trajectories, checkpoints, earlier assessments, ledger and manuscript remain unchanged.
