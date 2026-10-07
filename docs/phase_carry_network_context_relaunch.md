# Network-context study relaunch — 6 October 2026

The interrupted network-context experiment was relaunched at **08:11 EDT** on the dedicated GTX 1080 Ti (physical GPU 1), using the existing protocol and scientific source files. The new worker PID is **9027**. Its current record is [current-worker.json](../outputs/phase-carry-network-context/current-worker.json); the [original launch record](../outputs/phase-carry-network-context/launch.json) is preserved unchanged.

Before relaunch, the worker recorded in the original launch was absent, the coordinator lock was available, and GPU 1 had no compute process. The relaunch verified all **47 source hashes and 764 input hashes**, the recorded runtime and CUDA libraries, available context-gate evidence, and nine committed checkpoints.

Eight new paths had completed six-unit pilots; their states and histories remain available. One interrupted path, `seed-7_coarse_reset-cell-23_cell-23_neg10`, had committed only its initial state. Its attempted three-unit temporary save was an incomplete ZIP archive and was not used as a checkpoint.

The frozen runner's history validator requires a positive elapsed duration and therefore cannot directly resume a time-zero checkpoint. Every prepared-source array in that checkpoint—including chemistry, geometry, metadata, polarity, the float64 phase carry, and rounding counters—was verified exactly identical to the pinned source. The initial checkpoint, history, status, and incomplete temporary file were archived with their hashes. This path restarts from its unchanged prepared source. The scientific implementation and protocol were not modified.

An earlier relaunch inspection stopped at the zero-duration validation before launching any worker or moving files. Its [attempt record](../outputs/phase-carry-network-context/relaunch/2026-10-06T120835Z/attempt.json) distinguishes this operational interruption from a failed scientific acceptance gate.

The experiment still comprises three existing histories, 33 new paths and nine exact reused paths. Every required context gate and paired pilot comparison remains in effect before long continuation. Relaunch verification does not qualify unfinished scientific outcomes.

- [Successful relaunch record](../outputs/phase-carry-network-context/relaunch/2026-10-06T121057Z/launch.json)
- [Verification and archived artifact hashes](../outputs/phase-carry-network-context/relaunch/2026-10-06T121057Z/verification.json)
- [Resumed worker log](../outputs/phase-carry-network-context/relaunch/2026-10-06T121057Z/run.log)
- [Study status](../outputs/phase-carry-network-context/status.json)
