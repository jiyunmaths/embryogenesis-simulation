# Development GPU single-cleavage qualification: completed assessment

The independent saved-array review accepts all **16 CPU/GPU cases and eight timestep pairs** under the unchanged prospective limits. The isolated developmental extension is qualified for the prescribed single-cleavage design at grid 72 through t=3. This is **zero new scientific developmental histories**. Naturally scheduled repeated cleavage and developmental spatial/time refinement remain unqualified.

[Prospective protocol](development_gpu_qualification.md), [launch record](development_gpu_qualification_launch.md), [literal independent review](../archive/study-reviews/2026-10-10-development-gpu-cleavage/independent_review.json), [review code](../archive/study-reviews/2026-10-10-development-gpu-cleavage/review_saved_arrays.py), [machine-readable assessment](development_gpu_qualification_assessment.json).

## Method and controls

Resident PyTorch arrays and chemical/transport matrix operations are retained, with a separate custom CUDA contractile-ring mechanics kernel. The accepted mature sources are unchanged. The new backend preserves phase carry through pressure correction and daughter partition, conserves chemical amounts and volume targets at abscission, resizes cell buffers, and retains polarity, lineage, cycle/partition RNG states and numerical corrections on restart. Missing developmental schema/correction fields are rejected.

An analytic radius-0.4 diffuse cell starts with uniform activator 1.03 and inhibitor 0.98. Two prescribed cleavage directions, four mechanical branches and dt=0.00375/0.001875 produce sixteen synthetic cases. All start one progressive division at t=0, are capped at two cells and are followed through t=3. The branches are coupling off, tension/adhesion only, directional polarity tension only and the full combination; chemistry remains active throughout. The oblique cases have a specified initial polarity vector. These axes and polarity vectors exercise numerical mechanics and are not spontaneous-axis evidence.

The independent NumPy reference implements the developmental equations without calling the GPU solver. A separate 0.06-unit historical native prefix checks the numerical extension against prior arithmetic. The mother stays a well-mixed chemical compartment until abscission; gradients between connected lobes are omitted. Beginning-geometry splitting remains first order. Pressure correction and occupancy splitting preserve the represented field, extending historical float32 cleavage rather than claiming identical old-method trajectories.

## Numerical results

Each entry gives maximum observed discrepancy / declared limit. Cellwise comparisons use identical IDs; ancestor-averaged chemistry and whole amounts cover all observations. **No unmatched observations occurred** in the sixteen backend or eight timestep comparisons.

| Metric | NumPy/GPU | Coarse/fine GPU |
|---|---|---|
| Cellwise log concentration | 3.16e-09 / 1e-05 | 0.000184 / 0.01 |
| Polarity component | 2.64e-09 / 1e-05 | 0.000321 / 0.01 |
| Relative tissue axis ratio | 1.38e-09 / 0.0001 | 0.000881 / 0.01 |
| Relative cell volume | 3.16e-09 / 1e-05 | 0.000136 / 0.005 |
| Relative transport-operator norm | 1.05e-06 / 1e-05 | 0.00174 / 0.01 |
| Final phase field | 5.96e-08 / 2e-05 | 0.000778 / 0.02 |
| Ancestor-averaged log concentration | 1.79e-09 / 1e-05 | 0.000177 / 0.01 |
| Relative whole species amount | 2.25e-11 / 1e-05 | 4.62e-05 / 0.005 |

The largest native-prefix phase discrepancy is 4.17e-07 against 2e-5, and its largest chemical log discrepancy is 8.25e-09 against 1e-5. All sixteen independent eight-step chemical/dilution checks pass the 1e-11 limit. The existing mature branch and the new no-cleavage branch agree bit-for-bit for eight steps, including phase carry, chemical correction and geometry.

Every case completes one resolved division: coarse abscission at t=0.90375 and fine at t=0.901875. CPU/GPU event times agree; the coarse/fine shift 0.001875 is below 0.0075. Every daughter has one connected component at the event boundary and saved nondividing observations. The pinching mother is excluded from that connectivity requirement.

Maximum visible-volume jump at abscission is 1.34e-08 against 1e-6; recorded chemical-amount jump is 4.44e-16 against 2e-14. Independent long-double reconstruction of the GPU represented chemical pairs gives maximum amount jump 1.81e-16. Target-volume jump is 0; represented occupancy discrepancy is 3.33e-16 against 1e-14. Maximum every-step dilution closure error is 4.44e-16. Maximum volume-target error is 1.237% against 5%, minimum equivalent radius is 5.225 grid intervals against four, and no mechanical clipping occurs. Every saved boundary screen passes. All saved transport operators preserve constants and volume-weighted amounts, with symmetric volume-weighted conductances and nonnegative exchange coefficients.

## Restart and preservation evidence

All sixteen paths pass eight exact restored/uninterrupted steps during the live furrow and eight immediately after abscission. A separate full replay from t=0.45 through t=3 reproduces the final tensors, corrections, geometry, cell count, IDs, lineage, events, clocks and RNG states bit-for-bit in **all sixteen cases**. [Replay record](../archive/study-reviews/2026-10-10-development-gpu-cleavage/restart-replay.json). The initial replay script compared the decimal clock to exact 0.45; correcting that audit-only comparison to the existing 1e-12 clock tolerance resolved the script assertion. The saved-array review initially required every due time to be finite, overlooking the model's deliberate +infinity marker for an already dividing mother. Its corrected check permits that marker only for active mothers. Both first attempts and scripts are archived; solver code, data and acceptance thresholds were unchanged.

Two additional tests on the GPU pass: an ellipsoid divides along its long axis without changing fields at division initiation; and naturally scheduled development reaches four cells with an exact restart from the three-cell stage while another mother is dividing. These are software checks, not repeated-development accuracy or new histories. Twenty-eight unit tests and five explicitly enabled existing CUDA backend tests also pass.

All 23 source pins, 19 input pins and 1052 previous scientific/ledger/manuscript pins verify. The review records 1285 consumed hashes and the runtime binary hash. Original production reports retain their independent-review-pending flag; this separate review supplies the acceptance decision. The preliminary unsuffixed prepared directory was never run, and its earlier source version is preserved. The qualification took 21.93 minutes with four parallel CPU references; this is a qualification runtime, not a full-development speed benchmark.

## Scope and next gate

Immediate GPU event amounts, targets and occupancy volumes are reconstructed from saved arrays. The pointwise represented-occupancy jump, intermediate component counts, every-step extrema and CPU event jumps are recorded runtime evidence rather than raw-field reconstruction of every step. Final CPU/GPU and immediate GPU daughter fields are independently relabeled. Restart replay uses the candidate solver; the independent physical comparison is the separate NumPy trajectory.

Next qualify **naturally scheduled repeated cleavage from one cell toward sixteen**, retaining live-event restarts, matched controls, temporal refinement, cell/interface resolution and boundary screens. Then declare unused seeds and a separate protocol for three new independent developmental histories. This single-cleavage pass does not establish spatial convergence, a tissue continuum, chemical pattern initiation, inherited identity or biological lineage commitment. Existing scientific claims and manuscript remain unchanged.
