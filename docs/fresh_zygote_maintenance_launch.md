# Fresh-history long maintenance: launch

Launched **12 October2026 at22:45:30 America/New_York** (13October02:45:30UTC), process **740991**, on the dedicated GTX1080Ti. [Registered design](fresh_zygote_maintenance.md), [process record](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/process-launch.json), [confirmed launch](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/launch-verification.json), [execution log](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/run.log).

## Question and experimental units

Continue patterns that arose during actual development, while chemistry, mechanics, polarity and conservative contact transport continue to evolve. Start from the **24 independently accepted t=78 endpoints** of the formation study, without resetting chemistry or changing mechanical couplings. The preceding60-unit formation stage remains frozen; it is not repeated or relabeled as a new history.

Run another **180 model-time units**, physical **t=78–258**, reaching **240 total elapsed units since t=18**. Histories **1101–1103** retain four nested contexts: feedback off, tension/adhesion, polarity tension and the full combination. Run every context at **dt=0.00375 and0.001875**, grid120 and extent2.24. There are **three existing independent histories**, **24 new numerical trajectories**, **twelve timestep comparisons**, and **zero new independent histories**.

Measure continuous retention, all sampled losses/recoveries, complete24-unit sustained losses/recoveries and persistence during the final24 units. Use SD(natural-log activator)>0.1, with0.05/0.2 sensitivity thresholds. Both clocks must agree on each registered outcome and corresponding first loss/sustained-loss/recovery times within0.3 units. Reproducible loss is a valid scientific result; execution does not require the pattern to survive.

## Completed launch checks

All **100 software checks passed**, including the actual180-unit coordinator loops at both clocks in a minimal CPU fixture:1,201 regular observations,12,001 dense graphs, six checkpoints, two exact eight-step restarts and correct physical/maintenance/total clocks. Tests cover strict thresholds, full loss/recovery windows, paired outcome sensitivity, stop-on-failure and provisional completion. These software fixtures do not substitute for scientific GPU validation.

The [independent NumPy/hash audit](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/prelaunch-audit.json) verifies **68 source pins,602 input pins and432 inherited arrays**. It checks original reviewed endpoint hashes, complete correction/carry inheritance, unchanged model/gates and retained prior failures. Only metadata `config.steps` changes. The unchanged solver uses **PyTorch arrays/matrices with custom CUDA mechanics, geometry and polarity**, plus compensated phase and chemical state. Large outputs stay on the existing data-drive output link. The drive had183.8GiB free before launch, and the scientific GPU was idle with8MiB allocated.

Protocol SHA-256: `78f7a98bdb2a3d0ecdaa825f94898d9111b664ba197654014d2ecc02ed656bfd`.

## Confirmed first production observation

At **22:50:12 America/New_York (02:50:12UTC)**, the history1101 feedback-off coarse path had reached **physical t=84.45, maintenance elapsed6.45**, with sixteen cells. Its newly executed eight-step NumPy/GPU prefix passed, including exact correction-aware GPU replay. Maximum log-chemical discrepancy was3.476e-8 and relative transport discrepancy5.997e-7. The accepted runtime binary remains `878c308c87e3b29a369dc1f94cf7de71b0e302aca24a980a08c71fcdbc504693`.

Recorded physical checks pass: maximum volume error1.0647%, minimum radius8.5143 voxels, zero clipping and dilution amount error3.33e-16. The initial full checkpoint is saved. Host process and GPU use are confirmed in the launch verification. **No complete path or timestep comparison is yet finished** at this observation, and no long-maintenance scientific claim is accepted.

## Monitoring and decision

- [Coordinator status](../output/studies/fresh-zygote-maintenance/status.json) and [frozen protocol](../output/studies/fresh-zygote-maintenance/protocol.json).
- [First trajectory progress](../output/studies/fresh-zygote-maintenance/h1101_off_coarse/gpu-progress.json) and [actual-start prefix](../output/studies/fresh-zygote-maintenance/h1101_off_coarse/actual-start-prefix.json).
- [Software verification](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/test-verification.json) and [independent audit log](../archive/study-launches/2026-10-13-fresh-zygote-maintenance/audit.log).

Allow provisionally **28–32hours from launch**, approximately the early morning of14October America/New_York, conditional on unchanged speed and passing gates. The measured preceding run projects about28.7hours when trajectory duration is tripled and prefix overhead is retained. Re-estimate after the first completed coarse/fine pair; early steps alone do not establish whole-run speed.

Execute serial matched pairs with one CPU numerical thread. Every path receives a new actual-start prefix, and each full pair is compared before another context starts. Stop on any prefix, physical, exact restart or numerical/outcome comparison failure. No automatic fallback, tolerance relaxation, replacement history, rerun or downstream endpoint/identity assay. Completed execution remains provisional pending an independent saved-array review.

README science claims, manuscript and results ledger remain unchanged by launch. This test concerns chemical-pattern retention; it cannot establish autonomous memory, committed biological identity or universal numerical convergence. Dated status belongs in the lab notebook.
