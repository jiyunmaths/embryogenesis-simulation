# Longer fixed-ratio polarity assessment

All **21 moving jobs, 21 exact-context native/GPU checks, and 21 frozen endpoint assays complete**. Every moving-quality and endpoint-chemistry check passes. The study remains **`completed_with_unresolved_checks`** because the full-horizon history-9, zero-directional-tension timestep comparison fails its declared maximum chemical error tolerance. The original tolerances and status are preserved.

The main observed result is conditional: **removing polarity's directional tension permits initiation in history 9, while established patterns survive all tested polarity contrasts in all three histories**. Both timesteps give the history-9 formation outcome. This is qualitative agreement; the failed transient error criterion prevents full quantitative timestep acceptance.

## What was compared

The [protocol](polarity_robustness.md) fixes $\beta=2$, $D_a=0.02$, $D_b=0.55$, and activity-dependent tension/adhesion. Only directional polarity-tension contrast changes: $\chi=0,0.35,0.7$. Within each existing history, geometry, polarity, lineage and each chemical start are identical at $t=150$. The runs extend for 240 model-time units to $t=390$. Six unchanged coarse $\chi=0.7$ first-60-unit prefixes are reused without resetting their evolved states.

There are **three developmental histories**, 7, 8 and 9. The 18 coarse interventions and three fine history-9 near-uniform repeats are nested measurements, not 21 independent embryos. These are mature-state formation/maintenance tests, not new zygote-to-identity trajectories. Setting $\chi=0$ retains polarity dynamics and activity-dependent tension/adhesion; it removes only polarity's directional mechanical action.

## Formation and maintenance diverge

The contrast criterion requires across-cell SD of $\log a$ above 0.1 throughout the final 24 units.

| Directional contrast $\chi$ | Near-uniform starts with sustained contrast | Developed starts retaining contrast | Final developed contrast range |
|---|---:|---:|---:|
| 0 | 1/3 histories: history 9 | 3/3 histories | 1.304–1.424 |
| 0.35 | 0/3 histories | 3/3 histories | 1.457–1.564 |
| 0.70 | 0/3 histories | 3/3 histories | 1.573–1.687 |

History 9's zero-contrast near-uniform run crosses 0.1 at elapsed **75.36**, ending at contrast **1.58305** with late minimum **1.57509**. Its timestep-halved counterpart crosses at **75.32**, ending at **1.58297** with late minimum **1.57502**. Both produce sustained differences. The fine repeats at 0.35 and 0.70 remain weak, as do their coarse counterparts.

The nine developed coarse starts all retain strong differences. In this sampled study, stronger directional tension accompanies larger late developed contrast, even though it restricts formation from the small perturbations. This does not establish a monotonic law beyond these three points or permanent stability of the coevolving system.

![Matched formation references, formation-window timings, maintenance and timestep sensitivity](images/polarity-robustness-completed.png)

The top panels share the exact near-uniform starts with their frozen references. The lower left reports the first zero crossing of the largest uniform-state growth rate on each instantaneous frozen graph; it is not a full moving-system stability calculation. The lower center reports developed-state late minima. The lower right retains the quantitative refinement failure rather than hiding it behind matching final classifications. Lines connecting the three sampled contrasts are visual guides, not an interpolated phase diagram.

## The longer references rule out a short frozen observation window

All three matched near-uniform frozen references now meet the sustained-contrast criterion. DOP853 and Radau agree with maximum chemical log error below $1.75\times10^{-10}$ on these paths.

| History | Frozen first contrast crossing | Frozen contrast at 240 | Moving $\chi=0$ contrast at 240 |
|---|---:|---:|---:|
| 7 | 177.61 | 1.18455 | $1.95\times10^{-6}$ |
| 8 | 156.66 | 1.19023 | $8.37\times10^{-7}$ |
| 9 | 58.02 | 1.47042 | 1.58305 |

The earlier 60-unit frozen references were too short to show strong formation in histories 7 and 8. That concern is resolved for fixed geometry at 240. Their moving branches still lose contrast even without directional tension. Thus **polarity is not the only contributor to suppression**, and removing it from mature, previously polarity-conditioned geometry does not undo that geometry's history.

For history 9, motion without directional tension delays the measured onset relative to its frozen reference but permits formation. The matched nonzero-contrast interventions do not. This supports a polarity-dependent change in formation opportunity at fixed chemical parameters on this history, with the numerical qualification below.

## Spectral opportunity closes before, or after, formation

All near-uniform starts initially have positive largest uniform-state growth rates. The first positive-to-nonpositive crossings are:

| History | $\chi=0$ | $\chi=0.35$ | $\chi=0.70$ |
|---|---:|---:|---:|
| 7 | 63.31 | 19.64 | 11.65 |
| 8 | 59.66 | 18.34 | 10.90 |
| 9 | 137.11 | 46.25 | 27.14 |

Times are elapsed from the restart, interpolated between observations. History 9 forms a strong pattern before its zero-contrast graph loses uniform-state instability. Both nonzero contrasts close that opportunity earlier, while its chemistry is still weak. Histories 7 and 8 have smaller initial growth rates and lose the opportunity even at zero contrast. The recorded chemical trajectories, rather than crossing time alone, establish these outcomes.

At the endpoint, all near-uniform-derived coarse graphs have their largest transport eigenvalue below the lower instability-band edge **4.32503**. The successfully formed history-9 pattern nevertheless persists. These spectra concern uniform chemistry on a frozen snapshot; changing eigenvectors, nonlinear reactions, dilution and chemical feedback on geometry are not included in that linear diagnostic. The successful branch's late geometry has itself developed under nonuniform chemistry.

Measured total undirected conductance declines from the common initial graph:

| History | $\chi=0$ decline | $\chi=0.35$ decline | $\chi=0.70$ decline |
|---|---:|---:|---:|
| 7 | 14.5% | 43.6% | 58.7% |
| 8 | 14.5% | 44.3% | 59.5% |
| 9 | 20.4% | 45.5% | 60.4% |

This is consistent with directional mechanics weakening the contact transport and shifting the spectrum. Late edge-weight variation decreases in these sampled branches, but endpoint homogenization alone does not establish the cause or timing of the earlier spectral change. The contact-area proxy $g_{ij}\ell_{ij}$ is defined by the existing conductance closure; it is not an independent physical-area measurement. General curved, separated and nonorthogonal contact accuracy remains unresolved.

The first sampled nonpositive-growth point gives a more specific check of the homogenization hypothesis. In history 9 at $\chi=0.70$, total conductance is already **16.2% lower**, while the coefficient of variation of positive edge weights has increased slightly, from **0.54618 to 0.55035**. History 7 at the same contrast also has a small increase in variation at its crossing. Thus decreasing edge-weight variation is not a prerequisite for losing linear instability in every branch. Contact weakening and the full weighted spectrum are more informative than uniformity alone; these diagnostics do not independently isolate a single geometric cause. The contact measurements use the first saved nonpositive-growth observation, whereas the crossing times above are interpolated.

## Endpoint local attraction

All nine developed coarse endpoint graphs support local uniform/patterned chemical coexistence under the declared perturbation, stationarity and Jacobian checks. The successfully formed history-9 zero-contrast endpoint supports the same local coexistence at both timesteps. These are whole-network frozen chemical results.

The remaining near-uniform-derived endpoints return to uniform chemistry from their sampled starts. Their assays do not introduce a separately prepared patterned chemical state, so their negative bistability flags do not exclude every patterned basin. Prepared and near-uniform moving branches can develop different endpoint geometries.

## Numerical acceptance and the specific failure

All three 60-unit coarse/fine pilots pass. Over the full 240-unit window, the 0.35 and 0.70 comparisons still pass every declared criterion. Only the zero-contrast history-9 comparison fails:

| History-9 coarse/fine measure | Measured discrepancy | Declared limit | Decision |
|---|---:|---:|---|
| Maximum raw chemical log error, $\chi=0$ | 0.0257403 | 0.01 | **Fail** |
| Contrast-onset time, $\chi=0$ | 0.03704 | 0.30 | Pass |
| Uniform-growth crossing time, $\chi=0$ | 0.09290 | 0.30 | Pass |
| Late sustained-contrast classification, $\chi=0$ | Both true | Same classification | Pass |
| Maximum raw chemical log error, $\chi=0.35$ | $1.9133\times10^{-5}$ | 0.01 | Pass |
| Maximum raw chemical log error, $\chi=0.70$ | $4.9095\times10^{-6}$ | 0.01 | Pass |

The failed maximum occurs at elapsed **117.75** in cell **27**'s activator: coarse 0.155195 versus fine 0.151251, a concentration ratio difference of about **2.61%**. The error exceeds 0.01 at 191 of 1,601 sampled times, in four intervals between elapsed 84.15 and 121.05. In the final 24 units, maximum chemical log error is only **0.00024618**. All its polarity, shape, volume, transport and growth-rate comparison criteria pass.

Therefore, the failure is concentrated in the nonlinear formation transient, with close late chemical agreement and matching outcomes. This does not justify relaxing the original tolerance or calling the whole trajectory quantitatively resolved. Another timestep halving is required to determine whether the transient chemical error decreases adequately.

Across all 21 moving jobs, maximum volume error is **1.161%**, minimum equivalent radius is **5.106 voxels**, clipping is zero, sampled boundary occupancy is below $2.53\times10^{-8}$, and dilution amount error is below $4.45\times10^{-16}$. Exact-context native/GPU prefixes have maximum chemical log error $1.86\times10^{-7}$ and phase-field error $5.97\times10^{-8}$. All **105 frozen endpoint trials** pass, with maximum independent-solver log discrepancy $2.69\times10^{-9}$ and endpoint RHS magnitude $3.81\times10^{-9}$.

Assessment independently recomputes **33,621** stored graph spectra and contrast observations, with exact agreement at saved precision. It verifies protocol/source/input hashes, final histories/checkpoints, native/GPU evidence, frozen paths, endpoint reports and both refinement decisions. Short-prefix backend agreement remains distinct from full-horizon backend equivalence at every new parameter.

## Scope and next decision

The result supports different requirements for **initiating** and **maintaining** chemical organization, and a history-dependent effect of directional mechanics. It does not establish autonomous or inherited cell identity. No new lineage program or population classifier was introduced. Aggregate axis ratios change by at most **3.77%** from their mature starts; this study does not demonstrate a new large chemical-specific global shape axis.

The immediate priority is a targeted further timestep halving for history 9's near-uniform $\chi=0$ start, retaining exact initial state, physical times, backend checks and the original tolerances. Compare the existing fine run with the new finer run before expanding the parameter map. The failing branch should remain explicitly unresolved if the transient does not pass; no additional scientific simulation is launched by this assessment.

Reproduce the derived assessment and figure with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/embryo-mpl python -m embryo.polarity_robustness_assessment
```

Full values, error traces and evidence hashes are in `outputs/polarity-robustness/assessment.json`, generated by [the assessment module](../embryo/polarity_robustness_assessment.py). The [completed verification](polarity_robustness_completed_verification.json) records this assessment's status separately from the original launch snapshot. The consolidated ledger and manuscript retain their earlier evidence snapshot; this study is not silently promoted to fully accepted results.
