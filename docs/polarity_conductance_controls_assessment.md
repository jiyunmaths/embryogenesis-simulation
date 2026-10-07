# Completed polarity and conductance controls

**Preserving the initial chemical contact conductances restores sustained patterns at directional polarity tension chi=0.35 in all three tested histories.** Matched fixed-conductance controls also form at chi=0. In history 9, this removes the previously demonstrated polarity-associated failure to initiate. In histories 7 and 8, it rescues chemical organization from two contexts where neither original branch formed.

The original reporting remains **polarity-specific suppression in 1/3 histories, history 9**. Successful interventions in histories 7/8 do not change that fraction. This completed [prespecified experiment](polarity_conductance_controls.md) comprises eight new GPU paths and sixteen qualified reused paths within **three existing developmental histories**, not 24 developmental replicas or a population success estimate. All required numerical, physical and frozen-endpoint checks pass.

## Formation outcomes

All paths start with matched mature 16-cell geometry, polarity and near-uniform chemistry at physical time 150. The observation lasts 240 model units, ending at time 390. Chemistry uses beta=2, D_a=0.02 and D_b=0.55; activity-dependent tension/adhesion coefficients remain 0.25/0.35. The only crossed interventions are directional tension chi=0/0.35 and evolving/preserved chemical conductances. chi=0 removes polarity's directional tension action while retaining polarity dynamics. Both timesteps, 0.00375 and 0.001875, give the same formation classifications.

Chemical contrast means the across-cell standard deviation of natural-log activator activity. Formation requires contrast above 0.1 throughout elapsed 216–240. The table shows final contrasts at the fine timestep; all conductance-preserving branches satisfy the full late-window rule.

| History | Evolving, chi=0 | Evolving, chi=0.35 | Fixed, chi=0 | Fixed, chi=0.35 |
|---|---:|---:|---:|---:|
| 7 | 1.88427e-06 | 8.06451e-09 | 1.374710 | 1.373882 |
| 8 | 8.10362e-07 | 8.66772e-09 | 1.397166 | 1.397173 |
| 9 | 1.58318 | 1.36317e-08 | 1.471130 | 1.471137 |

The positive-contrast fixed branches cross the threshold at nearly the same times as their matched zero-contrast controls:

| History | Fixed chi=0 onset | Fixed chi=0.35 onset | Absolute difference | chi=0.35 late minimum contrast |
|---|---:|---:|---:|---:|
| 7 | 120.529209 | 120.548968 | 0.019759 | 1.181235 |
| 8 | 130.623773 | 130.642693 | 0.018920 | 1.199722 |
| 9 | 57.095659 | 57.098181 | 0.002522 | 1.471136 |

These are interpolated threshold crossings in dimensionless model time. Adding directional tension changes the fixed-conductance onset by less than 0.02 units. The largest descriptive final state difference across the two fixed contrasts is 0.00765403 in absolute log concentration. Across the full window the largest such difference is 0.0207964; the trajectories are not identical. Similarity is descriptive: no equivalence criterion was specified for this cross-intervention comparison. The patterns in histories 7 and 8 are still evolving at the endpoint; a successful late window is not an indefinitely stationary moving state.

![Chemical contrast and actual chemical-operator spectra in all three histories](../outputs/polarity-conductance-controls-review/2026-10-05T170941Z/polarity_conductance_controls.png)

Top: contrast; the dotted line marks the formation threshold. Bottom: largest instantaneous growth rate about uniform chemistry on the actual chemical operator. Shading marks elapsed 216–240. Fine-timestep curves are shown; both timestep levels pass. Near-overlap of the two fixed-conductance curves is an observed result. [PDF figure](../outputs/polarity-conductance-controls-review/2026-10-05T170941Z/polarity_conductance_controls.pdf).

## What was held fixed

Let G(t) be the symmetric matrix of exchange conductances between cells. Define graph stiffness K(t) and cell-volume matrix M(t) by

$$
K(t)=\operatorname{diag}(G(t)\mathbf{1})-G(t),
\qquad M(t)=\operatorname{diag}(V_1(t),\ldots,V_n(t)).
$$

The evolving and fixed chemical operators are

$$
\Delta_{\mathrm{evolving}}(t)=-M(t)^{-1}K(t),
\qquad
\Delta_{\mathrm{fixed}}(t)=-M(t)^{-1}K(0).
$$

The intervention preserves **exchange capacities G(0)**, while per-volume transport rates still change with measured volumes. It retains mechanical concentration conversion by V_old/V_new, with no artificial volume amount source. Geometry, polarity, activity-dependent tension/adhesion and chemical feedback all continue to evolve; the geometric contact network continues to govern mechanics and polarity. Diffusive amount exchange is conservative with the actual volumes. Reactions can still produce and remove chemicals.

The chemical conductances can therefore differ from the live geometric contacts in the diagnostic branch. This deliberate intervention tests the consequence of retaining the initial communication network; it is not a proposed physical closure for a developing embryo. Frozen spectra and endpoint assays use the actual **chemical** operator, rather than substituting the geometric one.

## A supported route from mechanics to signaling

For these reaction and diffusion parameters, the linear instability band is approximately **4.3250 < lambda < 42.0386**, where lambda is an eigenvalue of -Delta. A finite graph supports only its discrete eigenvalues; an unstable continuous band alone is insufficient.

Every evolving-conductance branch eventually loses its unstable spatial modes. Directional polarity tension makes the instantaneous growth rate cross zero earlier. Preserving conductances instead retains **two unstable spatial modes at every saved observation**, at both contrasts and in all three histories.

| History | Evolving chi=0 growth crossing | Evolving chi=0.35 growth crossing | Final largest lambda, evolving chi=0.35 | Final largest lambda, fixed chi=0.35 | Final growth, fixed chi=0.35 |
|---|---:|---:|---:|---:|---:|
| 7 | 63.160 | 19.624 | 2.552958 | 4.646150 | 0.035852 |
| 8 | 59.520 | 18.324 | 2.529660 | 4.631785 | 0.034351 |
| 9 | 136.968 | 46.206 | 2.792667 | 5.126936 | 0.081325 |

History 9 illustrates the distinction between initiation and maintenance. With evolving conductances and chi=0, contrast forms at elapsed 75.416, before the growth rate crosses zero at 136.968. The finite-amplitude pattern survives after uniform chemistry becomes linearly stable. At chi=0.35 the crossing occurs at 46.206, and contrast never reaches the formation threshold. With fixed conductances, both contrasts form at approximately elapsed 57.10 and retain a positive snapshot growth rate. Histories 7/8 lose their baseline growth opportunity without formation even at chi=0; fixed conductances restore formation at both contrasts.

The evolving branches also show a larger initial-to-final reduction in total geometric conductance with directional tension:

| History | Conductance decline, evolving chi=0 | Conductance decline, evolving chi=0.35 | Maximum individual volume change, fixed chi=0.35 |
|---|---:|---:|---:|
| 7 | 14.503% | 43.613% | 0.2893% |
| 8 | 14.538% | 44.321% | 0.2890% |
| 9 | 20.412% | 45.539% | 0.2875% |

Conductance totals count each undirected contact twice; their fractional changes are unchanged by counting once. Volume changes in this table compare initial and final values, rather than the maximum error relative to target volume over a trajectory. The larger conductance decline is descriptive. The intervention and restored formation support a tested conductance-mediated contribution; they do not identify network homogenization or prove that all polarity effects act through transport. Mode growth and integrated positive growth from frozen snapshots are not a full stability analysis of the time-dependent coupled system.

## Numerical qualification and endpoints

All **12 pilot and 12 full-window timestep comparisons pass**, without time alignment or relaxed limits. The four newly computed contexts have these full-window raw errors:

| History | Directional tension | Maximum log-concentration discrepancy | Original 0.01 gate |
|---|---|---:|---|
| 7 | chi-0.35 | 4.44935e-06 | Pass |
| 8 | chi-0.35 | 5.17091e-06 | Pass |
| 9 | chi-0.35 | 5.92462e-06 | Pass |
| 9 | chi-0 | 5.94475e-06 | Pass |

The maximum discrepancy across **all twelve full pairs**, including reused evidence, is **0.000409631**. Its largest value comes from the reused history-9 evolving chi=0 formation transient. The new contexts have maximum **5.94475e-06**. Polarity, axis ratio, volumes, chemical and geometric transport, growth, onset and crossing-time checks all pass their original gates.

Across all 24 paths, maximum target-volume error is **1.1521%** against 5%, minimum equivalent radius is **5.10658 grid spacings** against four, clipping is zero, maximum sampled boundary occupancy is **7.61e-09**, and maximum amount-conversion accounting error is **4.441e-16** against 2e-14.

All **24 frozen endpoint assays** pass settling and independent-solver checks. Their 120 chemical preparations are nested measurements. The twelve fixed-conductance endpoints have unstable uniform chemistry and stable patterned states from the sampled starts: **initiation and maintenance**, rather than uniform/patterned bistability. The two history-9 evolving chi=0 endpoints support local uniform/patterned coexistence. The other ten evolving endpoints settle to uniformity from their sampled starts; this does not rule out an untested finite-amplitude patterned basin. Different endpoint preparations can select different patterns.

Maximum recorded independent DOP853/Radau log disagreement is **7.49236e-10**. The final CPU review verifies all 41 source and 369 input hashes, restores all completed checkpoint records, and recomputes chemical contrasts, every saved chemical graph spectrum, all paired metrics, endpoint residuals/Jacobian stability, pattern-return distances and classifications. Secondary solver trajectories were not saved; their agreement is verified from the original hashed records without rerunning those solves. The compute path retains resident PyTorch chemistry/matrices with unchanged custom CUDA carry mechanics and geometry/polarity.

## Interpretation and next priority

The model can organize near-uniform chemical states, but **mechanically evolving exchange networks can close the opportunity for initiation**. Preserving those networks removes the observed polarity-associated loss in history 9 and the shared moving-context restriction in histories 7 and 8. This is stronger mechanistic evidence than a polarity-versus-no-polarity observation alone. It establishes a contribution under the tested interventions, not unique or universal causal mediation.

The assessment concerns mature continuations on geometries inherited from no-carry development. It does not establish biological cell identities, autonomous memory, inheritance, fresh carry zygote development, full spatial/developmental convergence or general geometric transport closure. Aggregate axis ratios remain close to their already anisotropic starts; rescued chemical organization does not establish a newly generated large shape asymmetry. The broader old-method ledger remains pending carry revalidation.

The next priority is to **revalidate maintenance of already developed patterns under carry mechanics**, with matched chi=0/0.35 moving controls in the same three histories and their own numerical gates. This would qualify the maintenance side of the initiation-versus-maintenance narrative before revisiting exchange/neighbor-context behavior or expanding the parameter map. No follow-up simulations are launched by this review.

## Evidence

- [Original protocol](../outputs/polarity-conductance-controls/protocol.json), SHA-256 `bc6b748601c96c357c47d12a6a77912fd0613b4b7eb3f448505bb7d8e1899934`.
- [Completed original summary](../outputs/polarity-conductance-controls/summary.json).
- [Full CPU verification, exact metrics and saved script](../outputs/polarity-conductance-controls-review/2026-10-05T170941Z/assessment.json).
- [Compact machine-readable assessment](polarity_conductance_controls_assessment.json).
- [Pre-review documentation snapshots](../archive/study-document-snapshots/2026-10-05-before-completed-polarity-conductance-review/manifest.json).

The original protocol document, launch verification, kernels, trajectories, checkpoints, endpoint assays, parent assessments, reporting rules and ledger remain unchanged. The eight new paths used **4.94 hours** of measured sequential GPU runtime. Raw outputs and figures are local and ignored by Git.
