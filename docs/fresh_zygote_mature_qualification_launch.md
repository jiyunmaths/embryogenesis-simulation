# Fresh-history mature-timestep qualification launch

**Subsequent status: stopped at elapsed zero on a checkpoint-scheduling error.** The CPU/GPU prefix passed; the six-unit trajectory did not advance. See the [stop assessment](fresh_zygote_mature_qualification_assessment.md) and [separate scheduler-fix design](fresh_zygote_mature_clock_fix.md). The launch observations below are retained as history.

Launched **12 October 2026 at 08:19:43 America/New_York** (12:19:43 UTC), process **654752**, on the dedicated GTX 1080 Ti. [Process record](../archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/process-launch.json), [run log](../archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/run.log).

This is the next numerical gate after [the accepted fresh-development control completion](fresh_zygote_control_completion_assessment.md). It continues the actual t=18 endpoints of histories 1101, 1102 and 1103, retaining chemistry, geometry, polarity, compensation, lineage and random states. Only integer clocks and dt change. See the [registered design](fresh_zygote_mature_qualification.md).

| Item | Registered scope |
|---|---|
| Physical time window | t=18 to24; six elapsed model-time units |
| Candidate timesteps | 0.00375 and 0.001875 |
| Mechanical contexts | off, tension/adhesion, polarity, full |
| Numerical paths / paired comparisons | 24 / 12 |
| Independent developmental histories | Three existing histories; zero new histories |
| Solver | Unchanged compensated DevelopmentGpu, PyTorch transport and custom CUDA mechanics/geometry/polarity |
| GPU scheduling | Serial paths on the dedicated GTX1080Ti; one CPU numerical thread |

All **58 software checks passed**, followed by an [independent pre-launch NumPy/hash audit](../archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/prelaunch-audit.json). The audit verified **53 source pins, 486 input pins and 432 inherited arrays across all 24 starts**. State arrays match the accepted developmental endpoints byte for byte. All four branches retain their original couplings and all inherited thresholds remain unchanged. The original failed pilot remains failed. About 285.6 GiB was free on the data drive at launch; inherited starts occupy 5.54 GiB.

Frozen protocol SHA-256: `f6c099718c3d368646912660fb91fc2ce1280ce11130fb105cf07790f9811131`.

## Initial status

At 12:19:59 UTC the host process was alive and using one CPU core for the launch-time integrity check. The coordinator status still read `prepared`, because GPU prefixes begin only after that check. This launch observation does not report a completed prefix or mature trajectory. A completion forecast will be recorded after production timing is available.

## Monitoring and decision

- [Coordinator status](../output/studies/fresh-zygote-mature-qualification/status.json).
- [Frozen protocol](../output/studies/fresh-zygote-mature-qualification/protocol.json).
- [First actual-start prefix](../output/studies/fresh-zygote-mature-qualification/h1101_off_coarse/actual-start-prefix.json) and [first trajectory progress](../output/studies/fresh-zygote-mature-qualification/h1101_off_coarse/gpu-progress.json), written when those stages reach them.
- [Software test record](../archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/test-verification.json).

The coordinator stops at any prefix, physical, restart or paired-comparison failure. There is no automatic finer-timestep substitution or long-window extension. Passing execution remains provisional pending independent saved-array review. The 60- and 240-unit scientific experiments have not been launched. The README claims, results ledger and manuscript evidence are unchanged by this launch.
