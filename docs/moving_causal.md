# Moving-geometry feedback controls

The production pilot is running in `outputs/moving-causal`. No completed shape-causality result is available yet. It uses the independently validated [joint signal/fate integrator](joint_fate.md); preparation rejects a failed validation or a different source checkpoint.

All branches start from the same resolved 72³, 16-cell checkpoint at developmental time 18. The domain extent is 2.24 and time step 0.0075. Geometry, target volumes, and lineage are retained. Activator/inhibitor receive identical seed-7, volume-weighted zero-mean perturbations of RMS 0.001 around equilibrium; fate and polarity are reset to zero. There is no further cleavage. Each branch runs to time 78, with observations every 0.6 and checkpoints and surface frames every 6 time units.

| Arm | Intervention | Comparison purpose |
|---|---|---|
| Full | Full signaling, fate, polarity, and mechanical feedback | Reference |
| No feedback | Disable fate/polarity-dependent mechanical properties; retain regulatory dynamics | Measure the contribution of mechanical feedback |
| No self-activation | Replace activator reaction by 1/h − a; retain coupling | Test the role of autocatalysis |
| No signal-to-fate | Set activator forcing of fate to zero; retain signaling and polarity | Separate the fate-mediated pathway |

Contacts and conservative transport are reconstructed as geometry moves. Signaling and fate share RK stages; existing polarity and mechanics follow. Regulator amounts are preserved across each mechanical volume change by applying dilution exactly once. Intervention-specific graph Jacobians are reported. Checkpoints preserve the deviation variables and arm; resume with `JointMovingSimulation.restore` to retain the experiment's equations rather than the base solver.

The protocol declares a late window of time 63–78. It compares the full branch's mean principal axis ratio with each control. A feedback-specific or self-activation-specific elongation requires an excess of at least 0.05 over the respective control. Full-loop persistent activator contrast requires weighted standard deviation at least 0.1 throughout the late window. The no-signal-to-fate difference is also reported, without making it an additional success gate. This axis-ratio test addresses elongation; it does not exhaust all possible forms of spatial organization.

Numerical quality is reported separately: maximum relative volume error below 0.05, minimum effective cell radius at least four grid spacings, zero phase-field clipping, and sampled boundary occupancy below 0.01. Failure to produce excess elongation is a scientific outcome, distinct from numerical failure. Every arm retains signal and fate vectors, shape histories, a final checkpoint, and a surface viewer. The runner writes `comparison.json` and `RESULTS.md` automatically after all branches finish.

This is one chemical seed on one already anisotropic mature geometry. It tests additional deformation caused by the interventions, not spontaneous axis formation from a zygote. Resetting regulators does not erase the geometry's developmental history. Moving contacts and dilution can generate signaling differences even in controls that homogenize on a frozen graph. The geometry/operator split still needs its own time refinement, and seed replication and the remaining transport-closure checks are required before broader claims.

The end-to-end two-step smoke test completed all four arms with passing quality checks. Forty-one relevant tests passed, including dilution amount balance, single fate advancement, intervention spectra, and exact checkpoint replay.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.moving_causal prepare \
  --output outputs/moving-causal-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.moving_causal run \
  --output outputs/moving-causal-repeat
```

Preparation freezes source and input hashes. Keep those files unchanged while a run is in progress. This experiment leaves the core live solver and the existing developmental refinement study unchanged.

## Concurrent completion and assessment

The full branch was already in progress when the user requested completion and assessment. `benchmarks/finish_moving_causal.py` runs the three remaining controls in separate processes using the unchanged `moving_causal.run` function and frozen source hashes. Each worker selects one original arm; it defers between-arm comparisons until collection. The original full branch is retained through t=78. Its redundant serial continuation is then stopped, with any partial no-feedback output archived rather than mixed into the worker trajectory. The original scientific protocol and criteria are unchanged. Execution records and worker logs are under `outputs/moving-causal/parallel-workers`.

A separate watcher calls `benchmarks/assess_moving_causal.py` after all branches finish. The assessment verifies complete observation times and matched initial shape, regulator, and fate values. It reports final and late-window shape, signal contrast, fate counts, numerical-quality extrema, and the original acceptance decisions. It writes `ASSESSMENT.md`, `assessment.json`, `comparison.png`, and `comparison.pdf` in the output directory. `assessment-status.json` distinguishes waiting, assessing, completed, and failed states. Both collection and assessment run automatically; their production results remain pending until the status files report completion.

The assessment script has been exercised against all four completed smoke-test branches, including the case of undefined fate separation before commitment. Smoke-test plots are workflow checks, not scientific results.
