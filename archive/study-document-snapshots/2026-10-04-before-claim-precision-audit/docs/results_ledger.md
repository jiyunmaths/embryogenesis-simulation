# Results ledger

Evidence snapshot: 2026-10-04T00:06:48.974327+00:00.

**Replication unit: developmental history.** The core attribute-model evidence represents **three distinct histories (7, 8, 9)**, reused across experiments. Cells, graphs, branches, pulses, perturbation seeds, and timestep/backend repeats are nested measurements. Do not sum histories across rows or report intervention fractions as population success probabilities.

Initiation, maintenance, and context sensitivity are separate claims. Numerical acceptance is separate from a scientific outcome; failed gates and unsettled states remain evidence. The history column counts registered attribute-model histories. The historical fate-model developmental-refinement row uses one separate history with five nested numerical repeats; analytic reference problems have no developmental history.

## Initiation

| Study | Attribute histories | Evidence and current conclusion | Validation |
|---|---:|---|---|
| `initiation` | 1 | No-feedback development produces strong chemical contrast; full coupling suppresses it in the paired zygote comparison. [attribute-development/comparison.json](../outputs/attribute-development/comparison.json) | quality_pass |
| `polarity-components` | 1 | Polar mechanics alone suppresses initiation; removing its mechanical contribution from tension+adhesion restores formation on the matched geometry. [feedback-long/comparison.json](../outputs/feedback-long/comparison.json), [feedback-polarity-ablation/comparison.json](../outputs/feedback-polarity-ablation/comparison.json) | quality_pass |
| `homogeneous-spectrum` | 1 | The direct branch loses graph-supported homogeneous instability during development. [feedback-spectrum/spectra.json](../outputs/feedback-spectrum/spectra.json), [feedback-spectrum/counterfactuals.json](../outputs/feedback-spectrum/counterfactuals.json) | descriptive_spectral_analysis |

- **initiation:** One paired history; formation contrast is not a discrete identity or shape-causality result.
- **polarity-components:** One mature geometry, two prepared chemical families, twelve continuations; chemical modulation of polarity remains active.
- **homogeneous-spectrum:** Instantaneous frozen-graph linearization; not a nonlinear-maintenance or moving-system stability proof.

## Maintenance

| Study | Attribute histories | Evidence and current conclusion | Validation |
|---|---:|---|---|
| `moving-maintenance` | 3 | All three histories retain developed chemical contrast and cell association after switching full coupling on; paired keep-off controls also retain it. [feedback-survival-validation/summary.json](../outputs/feedback-survival-validation/summary.json), [feedback-survival/comparison.json](../outputs/feedback-survival/comparison.json), [survival/comparison.json](../outputs/feedback-survival-validation/seed-8/survival/comparison.json), [survival/comparison.json](../outputs/feedback-survival-validation/seed-9/survival/comparison.json) | quality_pass |
| `frozen-coexistence` | 3 | Both endpoint graphs in each of the three histories support locally stable uniform and patterned chemical states. [feedback-endpoint-bistability/results.json](../outputs/feedback-endpoint-bistability/results.json), [feedback-endpoint-bistability-seed-8/results.json](../outputs/feedback-endpoint-bistability-seed-8/results.json), [feedback-endpoint-bistability-seed-9/results.json](../outputs/feedback-endpoint-bistability-seed-9/results.json) | chemical_acceptance_pass |
| `original-frozen-recovery` | 1 | Small perturbations recover on both original graphs; only the no-feedback graph retains strong chemical contrast. [attribute-persistence/results.json](../outputs/attribute-persistence/results.json) | tolerance_and_recovery_pass |

- **moving-maintenance:** Finite horizon t=90-150 after the 16-cell cap; not direct-feedback formation, continued division, or inherited identity.
- **frozen-coexistence:** Six graphs and 252 nested starts are not independent histories; no full moving-system stability or bifurcation classification.
- **original-frozen-recovery:** Two graphs from one paired history; recovery of uniform chemistry is not differentiation.

## Context

| Study | Attribute histories | Evidence and current conclusion | Validation |
|---|---:|---|---|
| `all-pair-exchange` | 1 | Among 96 informative chemical transplants, half return and half reorganize; none meets transferred-concentration likeness. [attribute-exchange/results.json](../outputs/attribute-exchange/results.json), [attribute-exchange/assessment.json](../outputs/attribute-exchange/assessment.json) | tolerance_pass |
| `common-environment` | 1 | Isolation erases differences in the tested releases; common reservoir exchange can support multiple stable chemical states. [attribute-common-environment/results.json](../outputs/attribute-common-environment/results.json) | independent_checks_pass |
| `finite-reservoir` | 1 | An evolving finite reservoir can support differences; formation and release depend on exchange regime. [attribute-finite-reservoir/results.json](../outputs/attribute-finite-reservoir/results.json), [attribute-finite-reservoir/summary.json](../outputs/attribute-finite-reservoir/summary.json) | qualified_with_retained_failures |
| `frozen-exchange-response` | 3 | Exchanged-cell responses favor the donor reference on both endpoint graphs in each history. [cell-response-exchange/results.json](../outputs/cell-response-exchange/results.json), [cell-response-exchange/assessment.json](../outputs/cell-response-exchange/assessment.json) | independent_solver_and_settling_pass |
| `moving-exchange-response` | 3 | Each history has eight donor-nearer moving response comparisons; exchanged chemistry reorganizes rather than retaining transferred values unchanged. [cell-exchange-response-moving/comparison.json](../outputs/cell-exchange-response-moving/comparison.json), [exchange-response-histories/comparison.json](../outputs/exchange-response-histories/comparison.json) | completed_quality_and_backend_checks |
| `neighbor-context` | 3 | Conservative neighbor averaging shifts the fresh low-state recipient in all three histories; broader challenges alter late states and pulse responses. [neighbor-context/results.json](../outputs/neighbor-context/results.json), [neighbor-context/summary.json](../outputs/neighbor-context/summary.json) | independent_solver_and_sham_pass |
| `delayed-pulses` | 1 | Both selected seed-8 immediate-pulse nonrecoveries recover when challenged after settling, under percentage and original-amount controls; alternative endpoints also resist these pulses. [neighbor-context-delayed/results.json](../outputs/neighbor-context-delayed/results.json) | independent_solver_and_control_pass |

- **all-pair-exchange:** 120 pairs within one fixed graph; post-run endpoint stability is exploratory, not independent replication.
- **common-environment:** 2256 sampled cell states descend from one history; fixed reservoirs supply external support. An initial strength-4 tolerance failure is retained; further refinement passes, with one selected independent Radau check per arm.
- **finite-reservoir:** Original strong-release numerical failure, unstable symmetric endpoints, and an unsettled follow-up remain exceptions; not new developmental histories.
- **frozen-exchange-response:** 48 reference comparisons nested in three histories; donor-nearer is not unchanged transfer or autonomous identity.
- **moving-exchange-response:** Prepared mature patterned basins; 24 comparisons reuse the same three histories, not new formation or population-frequency evidence.
- **neighbor-context:** 36 challenges and 72 comparisons nested in three histories; immediate pulses occur during reorganization; geometry is frozen.
- **delayed-pulses:** Post-selected follow-up of two cases in one reused history; 24 pulses are not a replication rate, a measured delay window, or moving-context evidence.

## Numerical

| Study | Attribute histories | Evidence and current conclusion | Validation |
|---|---:|---|---|
| `gpu-validation` | 1 | Four full accepted mature CPU/GPU control/pulse replays pass; additional contexts are checked before scientific use. [gpu-backend-validation/comparison.json](../outputs/gpu-backend-validation/comparison.json) | backend_pass |
| `mature-timestep` | 2 | Both existing histories pass mature exchange/retention and selected same-state response timestep halving; classifications and sampled recoveries agree. [exchange-response-histories-refined/refinement.json](../outputs/exchange-response-histories-refined/refinement.json) | refinement_pass |
| `survival-timestep` | 1 | The paired t=90-150 seed-7 maintenance continuation passes timestep halving. [feedback-survival-validation/refinement.json](../outputs/feedback-survival-validation/refinement.json) | pass |
| `conservative-reference` | 0 | Conservative volume-weighted transport converges on the controlled regular reference; random-walk scaling does not give fixed physical diffusivity. [live-transport-validation-complete/report.json](../outputs/live-transport-validation-complete/report.json) | pass |
| `geometric-closure` | 0 | Flat orthogonal reference checks pass; all three declared general-geometry closure checks fail. [geometry-transport-validation/comparison.json](../outputs/geometry-transport-validation/comparison.json) | fail |
| `positive-skew-prototype` | 0 | The separate conservative-positive flux prototype passes its specified patches. [positive-skew-flux/comparison.json](../outputs/positive-skew-flux/comparison.json) | pass |
| `boundary-skew-prototype` | 0 | The separate corrected-boundary prototype passes its prescribed patches. [skew-boundary-correction/comparison.json](../outputs/skew-boundary-correction/comparison.json) | pass |
| `historical-development-refinement` | 0 | The historical full-development refinement completes five cases but fails the combined predeclared gate. [development-refinement/comparison.json](../outputs/development-refinement/comparison.json) | fail |
| `cleavage-measurement` | 0 | A separate calibrated diagnostic resolves the cleavage measurement failure without rewriting the original failed screen. [cleavage-measurement/comparison.json](../outputs/cleavage-measurement/comparison.json) | pass |

- **gpu-validation:** Backend agreement in the tested mature regime, not an extra history, developmental convergence, or GPU division validation.
- **mature-timestep:** 18 repeated continuations, four full backend replays, ten context checks; no new history, spatial convergence, or refined development.
- **survival-timestep:** Same physical starts; no validation of preceding cleavage.
- **conservative-reference:** Regular reference does not validate diffuse-contact conductance closure.
- **geometric-closure:** Curvature, diffuse gaps, and nonorthogonal flux remain physical approximation limits.
- **positive-skew-prototype:** Not coupled to live mechanics; does not repair the current conductance approximation.
- **boundary-skew-prototype:** Not live-integrated; no general moving-tissue acceptance.
- **historical-development-refinement:** Historical fate model; not three new attribute-model histories, and no full-development acceptance transferred.
- **cleavage-measurement:** Controlled analytic cleavage, not full developmental or sharp-interface convergence.

## History-level context assessment

| History | State shifts / challenges | Response shifts / comparisons | Fresh recipient conservative-mixing ratio |
|---|---:|---:|---:|
| 7 | 4 / 12 | 17 / 24 | 0.342 |
| 8 | 4 / 12 | 15 / 24 | 0.321 |
| 9 | 3 / 12 | 14 / 24 | 0.583 |

The counts describe nested interventions within each history. The conservative recipient shift repeats across all three histories. Both immediate-pulse nonrecoveries are in history 8; its delayed follow-up is selected and does not add a fourth history.

## Unsupported claims and decisions

- Chemical differences and donor-nearer responses do not establish autonomous cell types, a prescribed type count, or biological function.
- Frozen chemical local stability is not full moving-system stability. A finite-horizon recovery is not permanent memory or inheritance.
- Mature continuation timestep/backend passes do not replace the failed historical developmental gate or establish attribute-development spatial convergence.
- Conservative exchange does not resolve the failed general geometric conductance closure. Passed flux prototypes are not silently substituted into live results.
- Polarity-dependent initiation suppression has one matched geometry; necessity is conditional on the tested coefficient set. The chemical multiplier remains active.

## Reproduce the consolidation

Run `python -m embryo.results_ledger` to deliberately refresh the snapshot and `python -m embryo.results_ledger --verify` to check its current evidence/source hashes. Refreshing records report values; it does not rerun assays or strengthen their inference. [Machine-readable ledger](results_ledger.json) contains study IDs, explicit history IDs, details, SHA-256 evidence hashes, and the complete output-directory inventory. Uncurated folders and caches are not promoted to accepted studies by their names or completion status.

No new scientific experiment was launched during consolidation. The next scientific decision remains moving-geometry confirmation of the context interventions, after reviewing this ledger and manuscript.
