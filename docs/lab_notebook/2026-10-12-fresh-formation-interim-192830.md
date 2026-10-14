# Fresh-history formation: interim assessment

Observed **12 October 2026 at 15:28:30 America/New_York (19:28:30 UTC)**. This is a dated snapshot of the [registered 60-unit formation pilot](../fresh_zygote_formation.md), not final numerical or scientific acceptance. [Calculation data and evidence hashes](2026-10-12-fresh-formation-interim-192830.json); [read-only analysis script](../../archive/study-reviews/2026-10-12-fresh-formation-interim/assess.py).

## Execution

**Nine of 24 trajectories and four of twelve paired timestep comparisons have completed; all their recorded gates pass.** History 1101 has completed all four branches at both clocks. History 1102 has completed feedback-off at the coarse clock; its fine counterpart is running at physical t=52.8, elapsed34.8 of60. The other history-1102 contexts and all history-1103 continuations are pending. These are nested paths within **three existing independent developmental histories**, not 24 independent histories; only one history has all its comparisons completed.

Process675358 remains active. The dedicated GTX1080Ti was at96% utilization with about2.9GiB allocated; the process/GPU inventory shows one scientific coordinator and no competing compute process on that GPU. The execution log contains no error output. No solver, frozen protocol, launcher or old status/report was changed for this assessment.

## Chemical formation in the completed history

Contrast is the unweighted standard deviation across cells of the **natural logarithm of activator concentration**, a dimensionless chemical-state measure. The primary threshold is0.1. Sustained onset means the beginning of a subsequently observed, uninterrupted24-model-unit interval above that threshold; an isolated threshold crossing is insufficient. Late persistence requires every saved observation in elapsed36–60 to remain above threshold.

| History1101 inherited mechanical context | Sustained onset, elapsed since t=18 | Physical onset time | Final contrast, finer clock | Total conductance change | Endpoint growing chemical modes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Feedback off | 11.70 | 29.70 | 1.4084 | −42.88% | 5 |
| Tension and adhesion | 11.70 | 29.70 | 1.3615 | −43.80% | 5 |
| Polarity tension | 12.15 | 30.15 | 1.4364 | −49.94% | 5 |
| Full coupling | 12.15 | 30.15 | 1.4548 | −50.85% | 4 |

All four branches form sustained differences and remain above threshold throughout the late window. Both clocks agree on onset at the0.15-unit observation resolution, with no sampled loss/recovery. The0.05 and0.2 sensitivity thresholds also give sustained and late-persistent outcomes in every branch, at both clocks. Late-window minimum contrast is1.247–1.273, comfortably above even the0.2 threshold.

**Polarity-associated contexts delay the primary onset by0.45 model units in this history; they do not prevent formation.** Formation in the feedback-off branch rules out a claim that polarity feedback is necessary for chemical differences in this preparation. Each branch inherits its own developmental chemistry and geometry, so this comparison concerns entire developmental contexts, rather than a feedback switch applied to one common mature state. One completed history cannot establish prevalence across histories or overturn earlier, different mature-reset experiments.

Contacts remodel substantially in every branch. Total undirected conductance, reconstructed as the sum over i<j of V_i Δ_ij, falls43–51%. Nevertheless, the frozen uniform-state calculation retains4–5 growing nonzero chemical modes at the endpoint, with maximum growth0.0955–0.1711. The observed conductance decline therefore does not remove all supported unstable modes in this history. These instantaneous chemical diagnostics do not establish stability of the full changing mechanochemical system.

History1102 feedback-off at the coarse clock also completes with sustained onset at elapsed13.8 and final contrast1.4188. Its active finer counterpart has crossed the primary threshold at the same sampled time and currently has contrast1.3306, but the full24-unit persistence window and late window are still pending. It is not yet a qualified second-history result.

## Numerical evidence and limits

The four completed paired comparisons pass all registered gates. The largest homologous chemical log discrepancy is0.000419 against0.01; the largest density-field discrepancy is0.001561% against1%; the largest dense spectral discrepancy is0.09378% against1%. **One contact-support mismatch** occurs in history1101 feedback-off at physical t=40.35:40 versus41 edges. It remains within the spectral and transport limits; the hard contact cutoff still produces a visible discontinuity and is not silently treated as identical topology.

Across the completed paths, maximum recorded cell-volume error is1.228% against5%, minimum equivalent radius is8.51 voxels against4, and clipping is zero. Every completed path records the required exact eight-step restarts at elapsed12 and36, unchanged ancestry/cycle/random-state bookkeeping and sixteen cells. The active path's recorded physical gates also pass.

This review independently recomputes contrast, sustained/late classifications, conductances, mass-weighted graph spectra and frozen growth rates from3,609 completed saved regular observations plus233 active observations, without importing simulation modules or using the GPU. The recalculated contrast and spectrum agree exactly with their saved values; maximum growth-rate discrepancy is2.22e−15. All63 frozen source pins match. Checkpoint density-field errors, every-step audits and dense paired maxima above are the coordinator's recorded results; full independent saved-array/dense-graph reconstruction remains a final-review requirement.

## Time estimate

Measured full60-unit production times average **17.4minutes at dt=0.00375** and **29.4minutes at dt=0.001875**. Actual-start prefixes average1.42minutes; each paired comparison takes about13seconds after its final trajectory report. Fourteen unstarted paths remain, seven at each clock, plus about13minutes of the active fine run.

Adding those measured costs and a two-minute final-integrity allowance gives **6.07hours remaining**, a central finish of **21:33 America/New_York on12October** (01:33UTC on13October). Plan for **6–7hours**, approximately **21:30–22:30 local time**, if all registered gates pass and GPU/I/O availability remains comparable. A failed gate stops the study. Final independent review takes additional time. The revised total execution estimate is about10hours; the original15–18-hour launch estimate extrapolated short runs that had disproportionately high checkpoint overhead.

The cohort remains provisional. Biological lineage identity, irreversible commitment and the separately gated240-unit maintenance extension are not established or started. README science claims, results ledger and manuscript are unchanged by this interim assessment.
