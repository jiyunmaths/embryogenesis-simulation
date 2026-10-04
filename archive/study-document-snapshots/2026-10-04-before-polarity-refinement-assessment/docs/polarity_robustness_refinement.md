# Targeted refinement of zero-contrast pattern initiation

The [completed 240-unit polarity study](polarity_robustness_assessment.md) reproduces history 9's initiation at both tested timesteps, but its zero-directional-tension formation transient fails the declared maximum chemical discrepancy criterion: **0.0257403 against 0.01**. Close final states do not resolve that transient failure. This test adds one further timestep halving, with the same physical start and unchanged tolerances.

The implementation is [polarity_robustness_refinement.py](../embryo/polarity_robustness_refinement.py). Its output directory is `outputs/polarity-robustness-refined/`; the previous study, protocol, sources, histories, checkpoints, and unresolved status remain unchanged.

**Launch snapshot:** the persistent runner started on **2026-10-04 at 07:56 EDT** on the free GTX 1080 Ti. The exact-context finer native/GPU check passes, and the 60-unit timestep pilot is running. Long-horizon refinement and endpoint outcomes are pending. [Launch verification](polarity_robustness_refinement_verification.json) records this pending status separately from the completed parent assessment.

## Design and independent units

| Quantity | Setting |
|---|---|
| Developmental history | Existing history 9; one selected numerical check, zero new histories |
| Chemical start | Exact original near-uniform perturbation, not a new random realization |
| Geometry, polarity, lineage and random streams | Original mature physical state at $t=150$ |
| Chemistry and material laws | $\beta=2$, $D_a=0.02$, $D_b=0.55$, $\chi=0$; activity-tension/adhesion unchanged |
| Existing fine reference | Completed GPU trajectory at $\Delta t=0.001875$, reused read-only |
| New finer trajectory | $\Delta t=0.0009375$; 256,000 new coupled steps |
| Full horizon | 240 elapsed model-time units, ending at $t=390$ |
| Observation / checkpoint interval | 0.15 / 3 model-time units |
| Pilot horizon | 60 elapsed units, before the original formation transient |
| Sustained contrast | Across-cell SD of $\log a$ above 0.1 throughout the final 24 units |

The new run starts at the original $t=150$ state. It does not switch timesteps at elapsed 60, 75, or 117.75, and it does not start from the developed endpoint. Retiming changes only the discrete step count, timestep and physically equivalent observation/end bookkeeping. Exact array and random-stream comparisons verify the physical start after applying the identical chemical and mechanical intervention. The initial chemical file is copied without regeneration.

This is an outcome-selected numerical diagnosis, not a new independent replication of initiation. Report it as **one additional trajectory within history 9**, nested in the existing three-history polarity study.

## Backend and pilot gates

The computation uses resident **PyTorch arrays and matrix operations**, with unchanged custom CUDA mechanics, geometry and polarity calculations. The completed four-run full native/GPU baseline validation at timestep 0.001875 is reverified, including source, trajectory, hardware, software and binary evidence.

The finer timestep has its own **0.6-unit native C++/GPU comparison**, using four native CPU threads, the exact history-9 starting geometry and chemistry, and the original strict backend tolerances. This is explicitly short-horizon agreement at 0.0009375; it does not claim a full native/GPU replay at the new timestep.

Only after this context gate passes does the finer trajectory advance to elapsed 60. Its observations are compared against the corresponding immutable portion of the completed fine reference. The long continuation starts from that same evolved finer checkpoint only if every pilot criterion passes. Pilot evidence is saved separately and hashed before appending the long history. A failed backend or pilot gate stops continuation automatically.

## Unchanged timestep acceptance criteria

| Maximum fine/finer discrepancy | Tolerance |
|---|---:|
| Raw absolute chemical log error over both species, cells and sampled times | 0.01 |
| Absolute polarity component error | 0.01 |
| Relative aggregate axis-ratio error | 0.01 |
| Relative cell-volume error | 0.005 |
| Relative transport-matrix Frobenius error | 0.01 |
| Absolute largest frozen modal growth-rate error | 0.001 |
| Spectral zero-crossing time error | 0.30 |
| Contrast-onset time error | 0.30 |

Crossing/onset presence and late sustained-contrast classification must agree. The pilot uses its original 12-unit late window; the full comparison uses 24 units. Individual volume/radius, clipping, positive state, boundary and dilution screens remain active at their existing thresholds. No tolerance is increased if the finer run fails.

At completion, five dual-solver frozen chemical trials assay the finer run's own endpoint, using the existing stationarity, local-return and Jacobian criteria. This is a separate chemical endpoint check, not stability of the coevolving system.

## Decision and scope

The final report retains both timestep pairs: the previous coarse/fine failure and the new fine/finer result. It reports the new maximum chemical discrepancy, its timing, late discrepancy, onset and spectral-crossing errors, outcome agreement, and the discrepancy relative to the earlier pair. A smaller discrepancy alone is insufficient: the new pair must meet every original criterion. Agreement of these two finer trajectories does not establish an asymptotic convergence order, spatial convergence, general contact-geometry accuracy, fresh zygote differentiation, or autonomous/inherited identity.

If the full comparison passes, the finer pair provides a numerically accepted version of this selected initiation trajectory under the declared sampled tolerances. The earlier coarse/fine failure remains in the evidence. If it fails, retain `completed_with_unresolved_checks` and diagnose the discrepancy before extending the parameter study.

## Execution and monitoring

Before preparation, **28 focused software checks passed**, including stage gating, exact-state retiming, restart behavior, immutable pilot evidence and independently recomputed numerical decisions.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.polarity_robustness_refinement prepare

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.polarity_robustness_refinement run

# Read-only root status; the finer job and its prefix also have status.json files.
python -m embryo.polarity_robustness_refinement status

# Reassess only after completion, between coordinator writes.
OPENBLAS_NUM_THREADS=1 python -m embryo.polarity_robustness_refinement assess
```

`protocol.json` pins sources, inputs, reference evidence and tolerances. `seed-9/physical-start-verification.json` records the retiming checks. `pilot-refinement.json` and `long-refinement.json` preserve separate numerical decisions. `summary.json` combines those decisions with endpoint results without altering the parent. Interrupted runs resume from the new trajectory's own physical checkpoint.

The existing fine trajectory took approximately **24.3 minutes** for its moving integration. Doubling its step count estimates about **49 minutes** of finer integration, plus roughly 2–3 minutes for the native/GPU prefix and endpoint/launch overhead. Allow **50–60 minutes overall** on the free GTX 1080 Ti; the 60-unit pilot should finish roughly 15 minutes after launch. These are measured-workload estimates, not completion guarantees.
