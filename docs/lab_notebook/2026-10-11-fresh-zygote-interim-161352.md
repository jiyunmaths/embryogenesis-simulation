# Fresh-zygote pilot: interim assessment, 11 October 2026, 16:13:52 EDT

At this observation time, **15/24 developmental trajectories and 7/12 paired timestep comparisons are complete**. All recorded path and pair gates pass. These are nested controls and numerical repeats within **three registered independent histories**, not 24 independent samples. History 1101 has all four pairs complete; history 1102 has three pairs complete and its full-coupling fine trajectory is running; history 1103 has not started.

The active path is `h1102_full_development_half` at **t=12**, with **12 cells**, 11 completed abscissions and 4 live furrows. Coordinator PID 527626 was confirmed alive at 16:10 EDT; its dedicated GPU was at 94% utilization. Approximately 301 GiB remained free on the data filesystem at 16:11 EDT.

[Data snapshot](2026-10-11-fresh-zygote-interim-161352.json), [frozen design](../fresh_zygote_pilot.md), [launch](../fresh_zygote_pilot_launch.md), [earlier 11:47 assessment](../fresh_zygote_pilot_interim.md). This entry does not update the accepted ledger or manuscript claims.

## Numerical status

| History | Mechanical branch | Maximum density-field error | Dense equal-count spectral error | Pair gate |
| --- | --- | ---: | ---: | --- |
| 1101 | off | 0.3699% | 0.9196% | Pass |
| 1101 | tension_adhesion | 0.3373% | 0.6722% | Pass |
| 1101 | polarity | 0.3579% | 0.9461% | Pass |
| 1101 | full | 0.3271% | 0.8829% | Pass |
| 1102 | off | 0.4368% | 0.7512% | Pass |
| 1102 | tension_adhesion | 0.4369% | 0.7508% | Pass |
| 1102 | polarity | 0.4888% | 0.9244% | Pass |

The density-field and dense spectral limits are both 1%. The tightest spectral comparison remains history 1101, polarity-only: **0.9461%**, just **0.0539 percentage points** below the limit. Its peak occurs at t=8.61 with four cells and five versus six contact edges. History 1102 polarity-only peaks at 0.9244%, also at a contact-edge difference. Passing these registered limits is provisional numerical evidence, not an independent final acceptance.

All regular observation cell counts match for these seven pairs. Dense traces retain occasional unequal-count observations near abscission: t=11.715 in history 1101 off/tension–adhesion, and t=8.895 in history 1102 tension–adhesion/polarity. These observations are reported and excluded from equal-dimension spectral comparisons; they are not silently imputed.

Every completed path ends at t=18 with 16 cells, 15 abscissions and no live furrows. The largest recorded volume error is **1.662%** (5% limit), minimum effective cell radius **8.497 voxels** (four-voxel limit), clipping zero, and dilution amount error at most 4.44×10⁻¹⁶ (2×10⁻¹⁴ limit). Boundary occupancy remains far below its 0.01 limit. All 16 recorded actual-start CPU/GPU prefixes pass; all 15 completed paths have exact eight-step live-cleavage restart checks at the three registered stages. These short checks do not replace the pending independent saved-array review.

## What the early endpoints show

| History | Branch | SD of log activator at t=18 | Shape axis ratio | Unstable modes on the frozen endpoint |
| --- | --- | ---: | ---: | ---: |
| 1101 | off | 0.003589 | 1.32446 | 9 |
| 1101 | tension_adhesion | 0.003611 | 1.32442 | 9 |
| 1101 | polarity | 0.003541 | 1.32703 | 8 |
| 1101 | full | 0.003564 | 1.32700 | 8 |
| 1102 | off | 0.001942 | 1.32065 | 9 |
| 1102 | tension_adhesion | 0.001956 | 1.32061 | 9 |
| 1102 | polarity | 0.001923 | 1.32292 | 8 |

Chemistry remains close to uniform: fine endpoints have log-activator SD approximately 0.0019–0.0036. The coarse full branch of history 1102 is also near uniform, but its fine partner and pair gate are unfinished. These short endpoints do **not** establish persistent chemical states or cell identities.

Shape anisotropy is already present in feedback-off as well as coupled controls, with axis ratios about 1.32. Polarity produces a small matched shift in shape and transport spectrum; it has not yet established a large shape effect or a causal identity difference in these fresh histories. At the completed fine endpoints, off/tension–adhesion have nine frozen positive-growth modes, whereas polarity/full have eight. This is an instantaneous calculation about the homogeneous equilibrium on a fixed endpoint, with the same chemistry parameters. It excludes continued mechanics, mode rotation, polarity dynamics and dilution, and does not predict persistence in the moving system.

## Runtime forecast

The observed mean runtime is **20.09 minutes** at dt=0.00046875 and **38.14 minutes** at dt=0.000234375. The estimate uses these completed full-path timings, subtracts observed elapsed runtime on the active fine path, and adds measured coordinator overhead. It does not assume that equal intervals of developmental time cost equally much; later stages contain more cells.

There are 9 unfinished paths: the active fine full-coupling path plus eight paths in history 1103. Estimated remaining computation is **4.27 hours**, with a practical range of **4–5 hours**. Central completion is **about 20:30 EDT on 11 October**, or **00:30 UTC on 12 October**. This covers t=0–18 development and automatic timestep comparisons only; independent review and the separately gated 60/240-unit continuations are excluded. A numerical stop or interruption invalidates the completion forecast.

## Historical-document check requiring reconciliation

All **38 source pins and 382 input pins** still match. A separate read-only check found **1,769 of 1,774 historical preservation entries unchanged**; the five differences are the manuscript README, sections Markdown/LaTeX, introduction-and-methods Markdown and bibliography modified by the authorized introduction revisions. Their original bytes are preserved in the manuscript archive and were verified against the old hashes. These documents are not inputs to the active simulation.

[Authorized revision record](../../archive/study-document-snapshots/2026-10-11-authorized-restart-introduction-edit/record.json). The unchanged coordinator checks the original document locations at completion, so **its final historical preservation gate is expected to flag these edits**, even if every numerical pair passes. No frozen source, protocol, preservation manifest or running process was changed during this assessment. Final review must explicitly reconcile the authorized editorial revisions and retain the distinction between a historical-document mismatch and a numerical failure. Scientific acceptance remains pending.
