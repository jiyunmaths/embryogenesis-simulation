# Completed carry polarity initiation study

**Directional polarity mechanics suppresses initiation in one of three tested developmental histories (history 9). Histories 7 and 8 fail to initiate even with directional polarity tension disabled, although their matched frozen references form patterns.** All three history comparisons pass the original numerical checks. The physical cause of the moving-versus-frozen restriction in histories 7 and 8 remains unresolved.

This is the completed assessment of the [matched carry experiment](phase_carry_polarity.md), applying the [reporting rule recorded before final assessment](moving_geometry_initiation_controls.md). Twelve numerical paths are nested within three existing histories: ten new moving runs and two reused history-9 zero-contrast runs. No new developmental histories were generated. A fraction of these three histories is not a population probability.

## What was tested

Each continuation starts from the same mature 16-cell state at physical time 150 within its history, with the same near-uniform chemistry. Directional tension contrast is chi=0 or 0.35; activity-dependent tension and adhesion remain 0.25 and 0.35. Polarity continues to evolve in both branches. Thus chi=0 disables its directional tension contribution; it does not remove polarity or its earlier influence on the starting geometry.

Both dt=0.00375 and 0.001875 run for 240 elapsed units, ending at physical time 390. Gierer–Meinhardt parameters are beta=2, D_a=0.02 and D_b=0.55. The grid is 72 cubed, extent 2.24, with interface parameter 0.085. PyTorch chemistry and matrices remain resident on the GPU; custom CUDA computes mechanics, geometry and polarity. Visible fields and contact accumulation remain float32, while phase-update residuals are retained in float64.

Formation requires SD(log activator) above 0.1 throughout elapsed 216–240. It measures sustained chemical heterogeneity; it does not identify biological cell types.

## History-level results

Values below use the fine timestep; the coarse timestep gives the same classifications and passes full-window agreement.

| History | chi=0: final SD(log activator) | chi=0.35: final SD(log activator) | Matched initial frozen reference: final SD | Interpretation |
|---|---:|---:|---:|---|
| 7 | 1.8843e-6; no formation | 8.0645e-9; no formation | 1.18455; sustained formation | Moving restriction persists without directional polarity tension. |
| 8 | 8.1036e-7; no formation | 8.6677e-9; no formation | 1.19023; sustained formation | Moving restriction persists without directional polarity tension. |
| 9 | 1.58318; sustained formation | 1.3632e-8; no formation | 1.47042; sustained formation | Polarity-dependent suppression supported. |

History 9's zero-contrast branch crosses the contrast threshold at elapsed 75.4163 and remains above it throughout the late window, whose minimum SD is 1.57521. None of the other five branches crosses the threshold. Tiny residual spreads near uniformity are reported as numerical measurements, not distinct phenotypes.

The frozen references were recomputed on the **actual exported initial GPU operator and volumes** for each history. All four initial moving paths within a history have identical chemistry, cell IDs, volumes and transport matrices. The earlier native operators differ slightly from these exported operators, so old references were not assumed to be exact matches. DOP853 and Radau agree to at most 1.84e-10 in log concentration; all three references satisfy the same late-window formation criterion. This verifies the moving-versus-frozen comparison conditional on these starts.

## Numerical qualification

All six 60-unit pilot comparisons and all six complete 240-unit timestep comparisons pass. The comparisons use raw concentration trajectories at matched times, without temporal alignment or relaxed tolerances.

| History | chi | Maximum absolute log-concentration discrepancy | Original limit | Result |
|---|---:|---:|---:|---|
| 7 | 0 | 7.5592e-7 | 0.01 | Pass |
| 7 | 0.35 | 7.5412e-7 | 0.01 | Pass |
| 8 | 0 | 7.8496e-7 | 0.01 | Pass |
| 8 | 0.35 | 7.8093e-7 | 0.01 | Pass |
| 9 | 0 | 4.0963e-4 | 0.01 | Pass |
| 9 | 0.35 | 7.7509e-7 | 0.01 | Pass |

All original polarity, volume, axis-ratio, transport, growth-rate, crossing-time and onset-time checks also pass. The ten new context gates pass, including exact carry checkpoint/restart checks. The three short positive-contrast mechanics comparisons against the independent float64 reference pass; their 0.15-unit scope does not establish long-horizon full-float64 equivalence.

Across all twelve moving paths, maximum relative volume error is **1.143%** against a 5% limit, minimum equivalent cell radius is **5.107 grid spacings** against a four-spacing minimum, maximum sampled boundary occupancy is **7.61e-9** against 0.01, and no clipping occurs. Maximum dilution amount error is **4.44e-16**, below 2e-14. The review verifies the pinned 37 source hashes and 200 input hashes, completed checkpoint/history consistency, and immutable comparison records.

## What the spectra suggest

The following times are newly derived from the carry trajectories. They refer to when the largest instantaneous spatial growth rate about uniform chemistry crosses zero. These are frozen-snapshot diagnostics, not stability results for the complete moving system.

| History | Frozen initial geometry: contrast onset | Moving chi=0: spatial-growth crossing | Moving chi=0.35: spatial-growth crossing |
|---|---:|---:|---:|
| 7 | 177.612 | 63.160 | 19.624 |
| 8 | 156.662 | 59.520 | 18.324 |
| 9 | 58.024 | 136.968 | 46.206 |

Histories 7 and 8 need much longer to form contrast on their fixed initial graphs than the interval over which their moving snapshots support positive growth about uniform chemistry. History 9 forms before that crossing without directional tension, whereas positive contrast brings the crossing earlier and formation fails. This is consistent with a restricted opportunity for initiation. It does **not** separate contact-conductance remodeling from changing volumes, dilution, inherited geometry, or their combined effects. The integral of positive instantaneous growth is descriptive and does not prove amplification in a changing network.

## Actual frozen endpoints

All twelve actual endpoint assays pass independent-solver and settling checks. Maximum DOP853/Radau log disagreement is 7.49e-10. The two history-9 zero-contrast endpoints support local coexistence of stable uniform and patterned chemical states. At the remaining ten endpoints, sampled preparations settle to uniform chemistry.

Those ten assays start from near-uniform moving endpoints and their prescribed perturbations; they do not establish that no separately prepared finite-amplitude pattern or alternative basin exists. Frozen chemical stability is conditional on the endpoint graph and is not full moving-system stability. These are twelve nested endpoint measurements, not twelve developmental histories.

## What this resolves and what it does not

The carry correction resolves the tested raw timestep failure while retaining the earlier initiation classification in this mature-state regime. The result supports polarity-dependent suppression in history 9 and confirms a moving-versus-frozen restriction beyond directional polarity tension in histories 7 and 8. It does not yet identify the physical cause of that restriction or establish the same cause in both histories.

The starts still inherit no-carry development. This study does not validate fresh carry development from a zygote, repeated division, inheritance, autonomous cell identity, general conductance closure, spatial convergence, or a full float64 backend. Initial aggregate axis ratios already range from 1.349 to 1.433, and subsequent whole-aggregate shape changes are modest; this is not evidence of a newly generated large geometric axis.

This assessment supplements the [original results ledger and precision audit](results_ledger.md#phase-update-precision-audit). Its 22 original claims, decisions, evidence hashes and carry dispositions are preserved. The tested near-uniform chi=0/0.35 continuation subset is now qualified; developed-start maintenance, other contrasts, exchange, neighbor and reservoir assays still need their own carry revalidation. The older 262.1-fold ratio at another contrast/window remains an old-method measurement.

## Next discriminating test

Histories **7 and 8 both qualify** for the already specified [dilution/contact-conductance experiment](moving_geometry_initiation_controls.md#first-follow-up-dilution-and-contact-remodeling): cross dilution on/off with evolving/fixed initial symmetric contact conductances while mechanics and chemistry co-evolve at chi=0. Reuse four accepted baseline paths and run twelve new paths across the three additional arms and two timesteps.

Fixed conductances must retain actual changing compartment volumes: Delta(t)=-M(t)^(-1)K(0). Holding Delta(0) fixed against changing volumes would generally lose amount conservation. Dilution-off arms must account for their explicit chemical-amount source. Keep the original 60-unit pilot, 240-unit horizon and raw refinement limits. Different rescue patterns in the two histories imply different mechanisms; the same qualified rescue pattern would support a common contribution in these tested contexts. No follow-up simulations were launched during this review.

## Evidence and preservation

- [Original completed protocol](../outputs/phase-carry-polarity/protocol.json), SHA-256 `7d31576dcde9ddc809aaefac438dbaea7841124bcc0177f7b33132330b31b20c`.
- [Original completed summary](../outputs/phase-carry-polarity/summary.json).
- [Independent final verification, actual-input frozen references and review script](../outputs/phase-carry-polarity-review/2026-10-05T005029Z/assessment.json).
- [Compact machine-readable assessment](phase_carry_polarity_assessment.json).
- [Pre-review document snapshots and hashes](../archive/study-document-snapshots/2026-10-04-before-completed-carry-polarity-review/manifest.json).

The original simulation protocol, kernels, histories, checkpoints, comparison reports, endpoint assays and ledger snapshot were not rewritten by this review. Raw evidence under outputs remains local and ignored by Git; the compact assessment records its provenance.
