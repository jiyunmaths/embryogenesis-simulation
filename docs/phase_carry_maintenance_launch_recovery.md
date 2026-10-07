# Maintenance launch review and recovery

The preparation completed successfully, but the original supervisor failed **before any of the twelve scientific maintenance paths started**. There were no new maintenance results to assess. The launch handoff has now been repaired and the first pilot is advancing; scientific maintenance acceptance remains pending.

## What failed

At 2026-10-05T18:14:27.924036+00:00, the supervisor read `state=prepared` before `prelaunch-verification.json` existed. It accepted only `preparing` while waiting, so it raised `Unexpected preparation state: prepared`. The verification file was written 3.462 seconds later and records a successful preparation. Preparation completion therefore did not mean the scientific experiment had run.

The corrected supervisor waits in either `preparing` or `prepared` until both prerequisites are present. Preparation failures and unexpected states still stop launch. The original launcher still verifies the frozen source/input hashes, fresh prepared status and successful qualification; the scientific coordinator retains its exclusive lock. This changes launch orchestration only. Scientific equations, kernels, inputs, timesteps and acceptance thresholds are unchanged.

The failed supervisor status, traceback, original script, launch metadata and prior documentation snapshot are preserved in `/home/sji/Documents/Programming/Python/embryogenesis-simulation/outputs/phase-carry-maintenance/launch-recovery/2026-10-05T200032Z`. The corrected supervisor and exact regression/verification script are archived alongside them. The [machine-readable recovery record](phase_carry_maintenance_launch_recovery.json) identifies their hashes and the new launch.

## Qualification reviewed

- All twelve actual developed-state context checks passed, including chemical stepping, amount accounting, native/carry coupling and exact checkpoint restart.
- All three prespecified independent float64 mechanics checks passed. These are component checks, not full moving-trajectory acceptance.
- All six reused initiation pilot comparisons and six reused full comparisons remain qualified.
- Read-only verification confirmed 43 scientific source hashes, 462 input hashes and 192 parent evidence hashes, with the protected prior documents unchanged.
- Nine handoff regression cases passed, including the exact `prepared`-before-verification sequence and rejection of failed or unexpected states.

Maximum preparation chemical log error and amount-accounting error were both 2.22e-16. Maximum short native/carry transport discrepancy was 3.55e-6, within its original 1e-5 limit. These checks qualify launch; they do not establish long-term pattern maintenance or broadly revalidate the old ledger.

The scientific protocol hash remains `3a293494c46c6fbe12a540c9ce2109c3c1f488494206af0fbbe0c836975db573`. The [prespecified maintenance protocol](phase_carry_maintenance.md) is unchanged.

## Recovered execution

Actual scientific launch: **2026-10-05T20:02:02.326510+00:00**, coordinator PID **2859520**, physical GPU 1 (GTX 1080 Ti), with two CPU endpoint workers. The original supervisor has finished the handoff; the detached scientific coordinator continues independently.

At the documentation snapshot (2026-10-05T20:03:25.852667+00:00), `seed-9_coarse_chi-0_pattern` had reached elapsed **11.4** toward its 60-unit pilot. Zero of twelve new full paths were complete. Each paired pilot must pass before its 240-unit continuation. Check `outputs/phase-carry-maintenance/status.json`, the current job status and `run.log` for live progress; this document is a dated snapshot.

The planning estimate is about **7–9 hours from the actual restart**, approximately 2026-10-06 03:02 UTC to 2026-10-06 05:02 UTC, conditional on passing pilot gates. The earlier measured estimate was 7.42 hours of moving computation; endpoint work, qualification and runtime variation affect the finish time. The short restarted pilot is insufficient for a new calibrated ETA.

## What can be concluded now

Preparation is qualified and scientific execution has begun. There is **no new maintenance conclusion**. The latest completed scientific result remains the [polarity/conductance initiation study](polarity_conductance_controls_assessment.md): preserving initial chemical conductances restored initiation in the three existing histories, while original polarity-specific suppression was demonstrated only in history 9.

The restarted experiment asks whether developed chemical differences persist under evolving conductances and chi=0/0.35 mechanics. It will report continuous retention separately from loss followed by late recovery and will compare those outcomes with the qualified near-uniform starts. Twelve new maintenance paths and twelve reused initiation paths remain nested within **three developmental histories**. They do not establish biological cell types, autonomous identity or new developmental replicas.
