# Mature-timestep qualification: interim assessment

Observed **12 October 2026 at 10:35:10 America/New_York** (14:35:10 UTC). This is a dated snapshot, not live acceptance. [Data snapshot](2026-10-12-mature-timestep-interim-143510.json), [launch](../fresh_zygote_mature_clock_fix_launch.md), [coordinator status](../../output/studies/fresh-zygote-mature-clock-fix/status.json).

**22 of 24 continuations and 11 of 12 timestep comparisons are complete. All recorded production gates pass.** Both levels of feedback off, tension/adhesion, polarity and full coupling are complete for histories1101 and1102. History1103 has completed off, tension/adhesion and polarity; its full-feedback coarse run is active at physical t=21, with the fine run still pending. Twenty-three actual-start prefixes are recorded, including the declared reused prefix. These remain branches and numerical repeats within three existing independent histories, with zero new histories.

The host process656696 was alive at14:33:10 UTC. GPU utilization was0% during CPU/reference/checkpoint work at that sample, which does not mean the coordinator stopped. Production subsequently advanced to t=21. The runtime is the pinned compensated DevelopmentGpu binary on the dedicated GTX1080Ti, with PyTorch transport and custom CUDA mechanics/geometry/polarity. All58 source pins were checked again and remain unchanged.

## Numerical agreement

| History | Mechanical branch | Maximum dense spectral difference | Maximum checkpoint-density difference |
|---|---|---:|---:|
| 1101 | full | 0.000783% | 0.000647% |
| 1101 | off | 0.000873% | 0.000725% |
| 1101 | polarity | 0.000784% | 0.000647% |
| 1101 | tension/adhesion | 0.000873% | 0.000725% |
| 1102 | full | 0.000788% | 0.000654% |
| 1102 | off | 0.000886% | 0.000733% |
| 1102 | polarity | 0.000792% | 0.000654% |
| 1102 | tension/adhesion | 0.000889% | 0.000733% |
| 1103 | off | 0.000981% | 0.000738% |
| 1103 | polarity | 0.000829% | 0.000635% |
| 1103 | tension/adhesion | 0.000980% | 0.000737% |

Both displayed criteria have an unchanged **1% limit**. The largest dense spectral discrepancy is **0.000981%**, and largest density-field discrepancy is **0.000738%**. Maximum per-cell log-concentration difference is **1.05e-06**, versus its0.01 limit. Maximum relative transport-matrix discrepancy is **0.001347%**, versus1%. There are no between-timestep contact-topology mismatches in the eleven completed pairs; this does not mean their transport weights are constant over time.

All22 completed paths preserve lineage/cell-cycle/random states and pass exact eight-step checkpoint replay at elapsed3. Each finishes at physical t=24 with sixteen cells. Maximum volume error is **1.161%** against the5% limit; minimum resolved radius is **8.51 voxels** against four. Clipping is zero, and dilution-stage amount error is at most4.44e-16. This is conservation for the dilution operation, not a claim that reaction kinetics conserve total chemical amount.

These are production reports and read-only source checks, not the independent final saved-array audit. Native field discrepancies are linked from production comparison records and have not been independently reconstructed in this interim assessment.

## Early behavior and its limits

Across the eleven completed fine-clock endpoints, chemical contrast—defined here as the standard deviation of natural-log activator concentration over the sixteen cells—increases **4.98–5.78 times** from t=18 to24. Its endpoint values are about0.0107–0.0209. Growth occurs in every completed branch, including passive mechanics with chemical/mechanical feedback disabled. Small inherited differences are being amplified in this window; neither persistent cell identity nor polarity-specific suppression is established.

The sum of pair conductances falls **9.6–12.9%** over the same window. This quantity is calculated from off-diagonal entries of MΔ, with M=diag(cell volumes), and expressed relative to its t=18 value; counting each undirected edge twice cancels in the ratio. Tissue axis ratios also change. Thus chemistry is growing on a remodeling contact network, without a chemistry reset. The short qualification cannot tell whether growth persists, saturates or reverses over60/240 units.

## Completion estimate

At the observation, measured mean production times were **3.12 minutes coarse** and **4.33 minutes fine**, plus about **86.3 seconds** per non-reused CPU/GPU prefix and **5.1 seconds** per paired comparison. Subtracting current-stage progress and allowing for final pin verification gives about **8.3 minutes remaining**. Allow **8–12 minutes from10:35**, approximately **10:43–10:47 America/New_York**, if gates pass. Final independent review is separate and is not included in that compute estimate.

No thresholds were relaxed, no finer fallback or new job was launched, and no README, results-ledger or manuscript claim is promoted by this snapshot. After completion, independently review saved graphs, observables, checkpoint fields and compensation-aware replay evidence before accepting the candidate mature timestep for the separately registered longer experiments.
