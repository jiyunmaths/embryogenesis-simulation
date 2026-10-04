# Completed assessment of the third initiation timestep

The history-9, zero-directional-tension repeat completes all 240 elapsed units. Its backend context check, 60-unit pilot, moving-quality screens and five frozen endpoint trials pass. **The full-horizon timestep criterion fails again:** the fine/finer maximum chemical log discrepancy is **0.0375031**, greater than both the original **0.01** limit and the previous coarse/fine discrepancy **0.0257403**. Retain `completed_with_unresolved_checks`.

The scientific outcome remains consistent across the three sampled timesteps: the exact near-uniform start develops sustained chemical contrast. This strengthens agreement in that selected qualitative outcome; it does not establish quantitative convergence of formation, add an independent developmental history, or establish autonomous cell identity.

## What was tested

The [protocol](polarity_robustness_refinement.md) preserves the original history-9 physical state at $t=150$, near-uniform chemical perturbation, polarity, cell order, lineage, random streams and material laws. Chemistry uses $\beta=2$, $D_a=0.02$, $D_b=0.55$. Directional polarity-tension contrast is $\chi=0$; polarity dynamics and activity-dependent tension/adhesion remain active.

One new trajectory uses timestep 0.0009375. Completed trajectories at 0.00375 and 0.001875 are reused read-only. All end at $t=390$ with observations every 0.15 units. These are **three numerical trajectories within one existing history**, nested in the three-history polarity experiment, not three new embryos.

| Timestep | First contrast crossing | Final SD of log activator | Minimum SD over final 24 units | First frozen-growth zero crossing |
|---|---:|---:|---:|---:|
| 0.00375 | 75.35797 | 1.583051 | 1.575091 | 137.11458 |
| 0.001875 | 75.32093 | 1.582970 | 1.575018 | 137.20748 |
| 0.0009375 | 75.26496 | 1.582849 | 1.574908 | 137.34952 |

Times are elapsed model-time units from the restart. Every trajectory remains above the original sustained-contrast threshold of 0.1 throughout the final 24 units. The zero crossings describe uniform chemistry on instantaneous frozen graphs, not stability of the complete moving system.

![Three-timestep outcomes, transient errors and descriptive timing diagnostic](images/polarity-refinement-completed.png)

## Quantitative acceptance remains unresolved

| Fine/finer measure | Discrepancy | Original limit | Decision |
|---|---:|---:|---|
| Maximum raw chemical log error | 0.0375031 | 0.01 | **Fail** |
| Maximum polarity component error | 0.00014538 | 0.01 | Pass |
| Maximum relative axis-ratio error | 0.00001255 | 0.01 | Pass |
| Maximum relative volume error | 0.00001332 | 0.005 | Pass |
| Maximum relative transport-matrix error | 0.00026651 | 0.01 | Pass |
| Maximum frozen modal growth-rate error | 0.00014675 | 0.001 | Pass |
| Contrast-onset time error | 0.055972 | 0.30 | Pass |
| Frozen-growth crossing time error | 0.142035 | 0.30 | Pass |
| Sustained-contrast classification | Both true | Same classification | Pass |

The new maximum discrepancy is **1.457 times** the previous adjacent-pair discrepancy, an increase of about **45.7%**. Comparing the coarsest and finest trajectories gives an even larger maximum, **0.0632434**. There is no demonstrated decrease in raw trajectory error over these three timestep levels; an asymptotic convergence order cannot be inferred.

All three pair maxima occur at elapsed **117.75**, in cell **27**'s activator. Its activities are 0.155195, 0.151251 and 0.145684 as the timestep decreases. The new adjacent-pair discrepancy corresponds to a concentration ratio difference of approximately **3.82%**. It exceeds the limit at **233/1,601 sampled times**, in two observed intervals: 82.8–88.8 and 93.0–121.65. In the final 24 units, its maximum chemical log discrepancy is much smaller, **0.00038704**.

The 60-unit pilot passes with maximum chemical discrepancy **0.00019427**, but it ends before the nonlinear formation transient. Its success therefore does not contradict the full-horizon failure. No threshold, late window or acceptance rule is changed after seeing the result.

## Timing explains much of the mismatch, without accepting it

A supplementary diagnostic fits one temporal offset over elapsed 75–125, using all cells and both species equally in log space. It compares the finer states at time $t$ with linearly interpolated fine log states at $t+s$, searching $s\in[-0.3,0.3]$.

The fitted offset is **0.082326 units**. It reduces log-state RMS discrepancy from **0.0045933 to 0.0006520**, and maximum discrepancy within that fitted window from **0.0375031 to 0.0065084**. The fine/finer paths thus resemble the same transition with a small timing displacement; cell 27's falling activator makes that displacement conspicuous in raw relative concentrations.

This fit is **post hoc and descriptive**. It does not replace the required comparison at identical physical times, establish the source of the timing error, or prove convergence of the chemical path after reparameterizing time. The raw criterion still fails. Residual state differences remain, and the fitted shift increases rather than decreases across successive adjacent timestep pairs: 0.055637 for coarse/fine versus 0.082326 for fine/finer.

## Endpoint basins and numerical quality

The finer endpoint's conservative graph has no supported growing uniform mode; its largest uniform growth rate is **−0.045617**. Nevertheless, its developed pattern and two 1% log perturbations return to a stable patterned equilibrium, while two 0.1% near-uniform perturbations settle near the uniform equilibrium. All five trials pass the original solver, stationarity, Jacobian and local-return checks.

Maximum independent-solver log discrepancy is **$3.06\times10^{-10}$**. The patterned endpoint Jacobian has largest real part approximately **−0.70031**; the near-uniform trials have largest real parts approximately **−0.04562**. This supports **local uniform/patterned coexistence on the frozen finer endpoint**, not global bistability or stability of the complete moving embryo.

The finer moving audit records maximum volume error **1.143%**, minimum equivalent radius **5.110 voxels**, zero clipping, boundary occupancy below **$7.17\times10^{-11}$**, and dilution amount error **$4.45\times10^{-16}$** or less. The exact-context native/GPU prefix has maximum chemical log discrepancy **$1.13\times10^{-8}$**, transport discrepancy **$2.31\times10^{-7}$**, and phase-field discrepancy **$5.97\times10^{-8}$**. Backend agreement is short-horizon at the new timestep, distinct from full-horizon timestep convergence.

The moving integration took **43.7 minutes**; the native/GPU prefix took **2.4 minutes**. The completed job is a numerical diagnosis within history 9, not another successful history in a statistical denominator.

## What to test next

The most useful next step is to separate numerical mechanisms before another full timestep halving:

1. Hold a measured contact graph and volumes fixed, and compare the actual chemical stepping rule at the three timesteps against independent ODE references through pattern formation. Follow with an identical prescribed geometry/volume path if necessary. A replay built from 0.15-spaced observations needs its own interpolation/refinement check; it does not recover all original per-step geometry.
2. Run a short, matched precision control for mechanics and geometry at fixed timestep. The current kernel computes forces and chemistry in float64 but stores each updated phase field in float32; contact-shell products also originate in float32. Test higher-precision field updates and contact accumulation separately, retaining the PyTorch/custom-CUDA architecture and original backend unchanged.
3. If those controls identify a source, verify a separately implemented correction against the existing references, then repeat this full formation trajectory with the unchanged criteria before expanding the parameter map.

Accumulated rounding, contact quadrature sensitivity, split-coupling error and nonlinear sensitivity are **hypotheses**, not established explanations. Agreement with a native prefix does not isolate errors shared by both implementations. No additional scientific simulation is launched by this assessment.

The broader [polarity finding](polarity_robustness_assessment.md) remains conditional: stronger directional mechanics restricts initiation at the tested mature starts, while developed states survive. This repeat supports qualitative persistence of the successful zero-contrast branch but leaves its quantitative formation trajectory unresolved.

## Verification and reproducibility

Assessment verifies pinned source/input hashes, exact-state retiming, all three physical checkpoint/history pairs, context gates, endpoint paths, and both pilot/full decisions. It recomputes **4,803 stored spectra and contrast observations**, with exact spectral agreement at saved precision, and recomputes the five endpoint classifications from saved chemical paths. The original tolerances and failed evidence remain preserved.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.polarity_refinement_assessment
```

The [assessment module](../embryo/polarity_refinement_assessment.py) generates `outputs/polarity-robustness-refined/assessment.json` and the figure above. [Completed verification](polarity_refinement_completed_verification.json) records this completed failure separately from the historical launch snapshot. The consolidated ledger and manuscript remain their earlier evidence snapshots. Planning documents are updated; their prior bytes are retained under `archive/study-document-snapshots/2026-10-04-before-polarity-refinement-assessment/`.
