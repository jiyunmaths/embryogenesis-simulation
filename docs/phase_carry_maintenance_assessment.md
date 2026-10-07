# Completed carry-corrected moving maintenance

**Developed chemical patterns survive evolving polarity-dependent mechanics in all three tested histories, at both directional-tension settings.** All twelve new moving paths retain contrast at every saved observation from model time 150 to 390. Their timestep and physical checks pass, and all twelve new frozen endpoint graphs support local uniform/patterned coexistence.

This qualifies the maintenance side of the initiation-versus-maintenance explanation. Near-uniform chemistry still fails to initiate at chi=0.35 in these same three starting contexts, whereas developed chemistry persists. The original polarity-specific initiation suppression remains **1/3 histories, history 9**. Histories 7/8 fail to initiate at both original contrasts, so their shared moving-context restriction remains a separate result.

## Design and reporting units

The [prespecified protocol](phase_carry_maintenance.md) uses the supplied mature 16-cell geometry, polarity and chemical preparations from existing histories 7, 8 and 9. New developed-start paths compare chi=0/0.35 at dt=0.00375/0.001875 for 240 elapsed units. Twelve qualified near-uniform initiation paths are reused. Chemistry, polarity, mechanics, measured volumes, conservative live conductances and dilution evolve throughout; there is no fixed-network or reservoir override. Setting chi=0 removes directional tension action while retaining polarity dynamics and activity-dependent tension/adhesion.

Parameters remain beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35, grid 72 cubed, extent 2.24 and interface parameter 0.085. Resident PyTorch arrays/matrix operations and custom CUDA carry mechanics and spatial geometry/polarity remain unchanged.

**Developmental histories are the replication units.** Twelve new paths, twelve reused controls, two timesteps, 24 endpoint graphs and 120 chemical preparations are nested within three existing histories. This adds no new history or population success estimate.

## Maintenance and matched initiation outcomes

Contrast S is the across-cell standard deviation of natural-log activator activity. The primary rule requires S>0.1 at every saved observation over elapsed 216–240. Continuous retention separately requires S>0.1 at every saved observation over elapsed 0–240. All twelve developed-start paths satisfy both; none has an observed first downward crossing or a loss-and-recovery episode. Observation spacing is 0.15, so unobserved brief transients are not excluded.

| History | Developed, chi=0 | Developed, chi=0.35 | Near-uniform, chi=0 | Near-uniform, chi=0.35 |
|---|---:|---:|---:|---:|
| 7 | 1.304081 | 1.459202 | 1.88427e-06 | 8.06451e-09 |
| 8 | 1.320426 | 1.456939 | 8.10362e-07 | 8.66772e-09 |
| 9 | 1.423889 | 1.564542 | 1.58318 | 1.36317e-08 |

Fine-timestep final contrasts are shown. Both timestep levels agree on every declared outcome. Developed-start whole-window minima equal their initial contrasts: 1.181858, 1.203618 and 1.289623 for histories 7, 8 and 9. The smallest developed-start late-window contrast across both levels is 1.296700.

![Matched initiation and maintenance trajectories in three histories](../outputs/phase-carry-maintenance-review/2026-10-06T021641Z-final/carry_maintenance.png)

Solid curves start developed; dashed curves start near-uniform. Top: chemical contrast, with dotted threshold 0.1. Bottom: actual chemical-graph snapshot growth about uniform chemistry, with dotted zero. Shading marks elapsed 216–240. Fine timestep shown. [PDF figure](../outputs/phase-carry-maintenance-review/2026-10-06T021641Z-final/carry_maintenance.pdf).

The positive directional-tension branches end with more contrast than their matched zero-contrast developed branches in all three histories. This descriptive comparison contradicts a claim that polarity tension universally suppresses chemical differences. It does not establish universal amplification or equivalence between branches.

## Linear initiation and nonlinear maintenance

Every developed-start graph eventually loses its unstable spatial modes, but the finite-amplitude pattern remains. The following fine-level crossing times refer to the instantaneous uniform-state growth diagnostic on each branch's actual evolving chemical operator:

| History | Developed chi=0 crossing | Developed chi=0.35 crossing | Final developed chi=0.35 contrast |
|---|---:|---:|---:|
| 7 | 44.216 | 17.125 | 1.459202 |
| 8 | 44.514 | 17.580 | 1.456939 |
| 9 | 145.191 | 50.094 | 1.564542 |

Uniform-state linear stability concerns whether small deviations grow. It does not rule out a separate stable finite-amplitude chemical state. The nonlinear feedback can maintain the supplied patterned preparation after the initial growth opportunity closes. Different chemical starts can remodel subsequent geometry differently, so the moving start comparison alone is not a demonstration of two stable attractors on one identical moving trajectory.

The frozen endpoints provide the separate same-operator test: each of the twelve new actual endpoint graphs supports stable uniform chemistry and a locally stable pattern. The developed state and two small patterned perturbations return locally to the pattern; two near-uniform preparations settle to uniform chemistry. All five preparations settle and pass the independent-solver checks. Across all 24 new/reused endpoints, 14 support coexistence and ten return to uniformity from their sampled starts. Those ten do not exclude an untested finite-amplitude basin.

These are **local chemical coexistence tests on frozen geometry**, rather than global bistability, indefinite moving maintenance or stability of the complete coupled system. The prior [conductance intervention](polarity_conductance_controls_assessment.md) supplies the complementary formation evidence: preserving G(0) restores initiation at both contrasts in all three histories, while mechanics and dilution remain active.

## Numerical and physical qualification

All **six new pilot and six new full-window timestep comparisons pass**. Including reused initiation, all twelve pilot and twelve full comparisons pass. Comparisons use raw matched times with no time alignment or relaxed thresholds. Continuous-retention and loss-presence decisions also agree.

| New maintenance context | Maximum log-concentration discrepancy | Original limit |
|---|---:|---:|
| History 7, chi=0 | 1.70957e-06 | 0.01 |
| History 7, chi=0.35 | 4.10077e-06 | 0.01 |
| History 8, chi=0 | 1.44429e-05 | 0.01 |
| History 8, chi=0.35 | 3.51091e-05 | 0.01 |
| History 9, chi=0 | 1.78427e-06 | 0.01 |
| History 9, chi=0.35 | 3.97816e-06 | 0.01 |

The maximum new maintenance discrepancy is **3.51091e-05**, in history 8 at chi=0.35. Across all twelve full pairs, the maximum is **0.000409631**, from the reused history-9 zero-contrast formation transient. The new maximum uniform-growth crossing-time discrepancy is 0.000651182, below 0.3. Polarity, axis, volume, transport and growth discrepancies pass their unchanged limits.

Across the twelve new paths: target-volume error is at most **1.1601%** against 5%; minimum equivalent radius is **5.10614 grid spacings** against four; clipping is zero; sampled boundary occupancy is at most **8.34e-09** against 0.01; amount-conversion error is at most **4.44e-16** against 2e-14.

All **24 new/reused frozen endpoint assays pass**. Maximum recorded DOP853/Radau log disagreement is **2.45006e-09**. The CPU review verifies all 43 scientific source and 462 input hashes, restores all completed checkpoints, checks float64 carry residuals, independently reconstructs weighted spectra/conservation and recomputes paired metrics, endpoint RHS residuals, Jacobian stability, return distances and classifications. The secondary solver paths were not saved: their original hashed agreement records are verified without rerunning the solves.

Twelve actual developed-context checks and the three independent float64 preparation comparisons also remain verified. Measured new moving runtime totals 5.94 sequential GPU hours. No new simulation was launched by this review.

## What this result does and does not establish

The accepted explanation is conditional chemical organization: mechanics can close an initiation opportunity while a developed pattern survives. History 9 demonstrates polarity-specific initiation suppression alongside maintenance at both contrasts. Histories 7/8 demonstrate maintenance despite failed initiation in both moving controls. This strengthens the distinction between starting a pattern and maintaining one.

The starts still inherit original-method development and supplied mature chemical preparations. This does not qualify fresh carry development from a zygote, autonomous or inherited biological identities, broad parameter robustness, general geometric transport closure, full spatial/developmental convergence or a new global shape axis. Aggregate axis ratios remain close to their already anisotropic starts. The original no-carry failures and 22-study ledger snapshot remain unchanged; this qualified continuation is a separate result.

## Next scientific priority

Use carry-qualified developed states for matched chemical-state exchange and response tests across different geometric contexts. Compare donor and recipient responses with unexchanged controls, account for imposed chemical amounts, and require context-specific numerical/refinement checks. Follow with moving neighbor-context and delayed-pulse comparisons. Older exchange/response amplitudes, reservoir results, other polarity contrasts and the 262.1-fold moving/frozen ratio still require their own carry revalidation.

## Evidence

- [Original protocol](../outputs/phase-carry-maintenance/protocol.json), SHA-256 `3a293494c46c6fbe12a540c9ce2109c3c1f488494206af0fbbe0c836975db573`.
- [Completed original summary](../outputs/phase-carry-maintenance/summary.json).
- [Read-only final verification and saved script](../outputs/phase-carry-maintenance-review/2026-10-06T021641Z-final/assessment.json).
- [Compact machine-readable assessment](phase_carry_maintenance_assessment.json).
- [Figure provenance](../outputs/phase-carry-maintenance-review/2026-10-06T021641Z-final/figure-provenance.json).
- [Pre-review reporting snapshots](../archive/study-document-snapshots/2026-10-06-before-completed-carry-maintenance-review/manifest.json).

The original maintenance protocol document, source files, kernels, inputs, trajectories, checkpoints, endpoint records, parent assessments, ledger and manuscript remain unchanged.
