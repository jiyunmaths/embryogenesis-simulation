# Completed moving-initiation controls

**Preserving the initial chemical contact conductances restores sustained patterns in both tested histories, 7 and 8. Removing mechanical dilution alone does not restore them.** This identifies changing exchange conductances as a shared contribution to the failed initiation in these two mature-state contexts, while cell mechanics, polarity, volumes and chemical feedback continue to evolve.

This is a completed assessment of the [factorial experiment](moving_initiation_controls.md), under the [decision rules specified before execution](moving_geometry_initiation_controls.md). The replication units are **two existing developmental histories**. Twelve new numerical paths and four reused baseline paths are nested within them; the interventions, two timesteps and endpoint trials do not add developmental replicas. Both histories pass all required checks. This is not a population success estimate.

## Results

Each path begins from matched mature 16-cell geometry and near-uniform chemistry at physical time 150 and ends at 390, after 240 elapsed units. Parameters are beta=2, D_a=0.02, D_b=0.55 and chi=0; activity-dependent tension/adhesion remain 0.25/0.35. Directional polarity tension is disabled, but polarity dynamics, mechanics and their prior influence on the inherited starts are retained. Both dt=0.00375 and 0.001875 give the same formation decisions.

Values below use the fine timestep. Chemical contrast is the across-cell standard deviation of natural-log activator activity. Formation requires contrast greater than 0.1 throughout elapsed 216–240, rather than a large final value alone.

| Intervention | History 7: final contrast | History 8: final contrast | Late-window outcome |
|---|---:|---:|---|
| Baseline | 1.88427e-06 | 8.10362e-07 | No formation in either history |
| Dilution removed | 2.59299e-07 | 3.44127e-07 | No formation in either history |
| Fixed chemical conductances | 1.37471 | 1.39717 | Formation in both histories |
| Fixed conductances + dilution removed | 1.18544 | 1.19884 | Formation in both histories |

The late-window minimum contrasts are **1.18123 and 1.20013** with fixed conductances, and **1.18544 and 1.18588** with the combined intervention. These comfortably exceed the original 0.1 threshold. The baseline and dilution-only branches never cross it.

With fixed conductances and dilution retained, the threshold is crossed at elapsed **120.529 and 130.624** in histories 7 and 8. With dilution also removed, crossings occur at **176.329 and 155.712**, respectively. Thus removing dilution does not accelerate initiation here; it delays it in the conductance-preserving branches. These are interpolated threshold-crossing times in model units, not biological hours. Their mechanism is not separately established: omitting dilution changes chemical amount accounting and can change both transient perturbations and later feedback.

![Fine-timestep chemical contrast and actual chemical-operator growth rates](../outputs/moving-initiation-controls-review/2026-10-05T113534Z/initiation_controls.png)

Top: contrast trajectories, with the formation threshold marked. Bottom: instantaneous growth about uniform chemistry on each actual chemical operator. Shading marks the required final 24-unit window. Both timestep levels pass; the fine trajectories are shown. Spectral curves diagnose frozen snapshots, not stability of the full moving system. [PDF figure](../outputs/moving-initiation-controls-review/2026-10-05T113534Z/initiation_controls.pdf).

## Why the conductance intervention matters

The symmetric conductance matrix G describes exchange capacities between cells. K=diag(G 1)-G and M=diag(V_i) describe the graph stiffness and measured cell volumes. The intervention freezes **G(0)** and K(0), while retaining the evolving volumes:

$$
\Delta_{\mathrm{fixed}}(t)=-M(t)^{-1}K(0).
$$

Holding conductances fixed therefore does not freeze cell geometry or per-volume exchange rates. The live geometric network still governs mechanics and polarity. The chemical network conserves diffusive amount exchange with the current volumes. Dilution remains active in the fixed-conductance single intervention, so its rescue does not depend on the artificial amount source introduced by dilution removal.

In the baseline branches, the sum of the symmetric conductances declines by **14.503% and 14.538%** from the initial to final snapshot in histories 7 and 8. Here the sum counts each undirected contact twice; the fractional change is the same when counted once. Measured individual volumes change by at most about **0.192%** between those snapshots. This comparison concerns initial-to-final changes, not the maximum volume-target error across the trajectory.

The frozen-snapshot spectrum supplies a consistent route from contact changes to failure. The continuous linear instability band for these fixed chemical parameters is approximately **4.3250 < lambda < 42.0386**, where lambda is an eigenvalue of -Delta. The largest supported eigenvalue starts at **4.64586 and 4.63212**, inside the band. Baseline final values fall to **3.75550 and 3.74463**, below its lower boundary, leaving no unstable spatial mode. Instantaneous spatial growth about uniform chemistry crosses zero at elapsed **63.160 and 59.520**. Removing dilution produces nearly the same spectral evolution and fails to restore formation.

The fixed-conductance branches instead retain two unstable spatial modes through the observed window. Their final largest eigenvalues are **4.64622 and 4.63186**, and their final largest growth rates are **0.035859 and 0.034359**. The combined branches retain the same opportunity. Geometry and capacities still move, but the chemical exchange network no longer loses its initial conductance structure.

The intervention supports a causal contribution of conductance remodeling in these contexts. It does not identify network homogenization, a single contact-area change, or a unique universal cause. Integrated positive instantaneous growth does not prove amplification in a time-dependent graph; changing modes and nonlinear trajectories require separate analysis. The evolving baseline and the conductance-preserving branches also develop different chemical feedback after their states diverge.

## Numerical qualification

All **eight 60-unit pilot comparisons and eight full 240-unit timestep comparisons pass**. Comparisons use raw matched-time trajectories, without temporal alignment or relaxed limits.

| History | Intervention | Maximum log-concentration discrepancy | Original 0.01 gate |
|---|---|---:|---|
| 7 | Baseline | 7.5592e-07 | Pass |
| 7 | Fixed conductances + dilution removed | 2.1945e-06 | Pass |
| 7 | Dilution removed | 2.2835e-08 | Pass |
| 7 | Fixed chemical conductances | 4.3268e-06 | Pass |
| 8 | Baseline | 7.8496e-07 | Pass |
| 8 | Fixed conductances + dilution removed | 3.4241e-06 | Pass |
| 8 | Dilution removed | 4.9266e-08 | Pass |
| 8 | Fixed chemical conductances | 5.0894e-06 | Pass |

The maximum full-window discrepancy is **5.08937e-06**, far below 0.01. All polarity, volume, axis-ratio, chemical/geometric transport, growth and timing checks also pass. Across all sixteen paths, the maximum target-volume error is **1.1491%** against 5%; the minimum equivalent radius is **5.10658 grid spacings** against four. There is no clipping, maximum sampled boundary occupancy is **1.64e-09** against 0.01, and maximum compartment/total conversion-accounting error is **4.441e-16** against 2e-14.

The final review verifies **39 pinned source hashes and 260 input hashes**, all sixteen completed checkpoint/history records, immutable pair snapshots, initial conductances, explicit dilution-off sources, and the actual chemical operator used at each endpoint. Chemical contrasts and all saved graph spectra were recomputed from arrays. The earlier implementation gates remain separately recorded in the [launch verification](moving_initiation_controls_verification.json).

## Frozen endpoint behavior

All **sixteen endpoint assays** pass their original independent-solver and settling checks; their eighty chemical preparations remain nested in two histories. Maximum recorded DOP853/Radau log disagreement is **5.39144e-10**. The review recomputes endpoint spectra, residuals, chemical Jacobian stability, pattern-return distances and classifications from saved paths. Secondary solver trajectories were not stored; their agreement is verified from the original hashed records rather than rerunning the solves.

The eight baseline/dilution-only endpoints have stable uniform chemistry and settle to uniformity from all sampled starts. Those starts are near-uniform moving endpoints and their prescribed perturbations; they do not exclude an untested finite-amplitude patterned basin.

The eight conductance-preserving endpoints have unstable uniform chemistry. All sampled preparations settle into stable patterned chemical states, supporting initiation and maintenance on the actual frozen chemical operator. They do **not** support stable uniform/patterned coexistence: uniform chemistry is unstable. Different preparations can select different spatial patterns, so this is not a unique-pattern result. Frozen chemical stability remains distinct from stability of the complete moving system.

## Scientific conclusion and limits

The prior [carry polarity assessment](phase_carry_polarity_assessment.md) showed polarity-specific suppression in **one of three histories, history 9**, plus moving-versus-frozen failures without directional tension in histories 7 and 8. This new factorial study identifies a common tested contribution in those latter two: **evolving contact conductances restrict initiation, and preserving them restores it even with dilution and moving mechanics retained.** The polarity-specific reporting fraction remains one of three. This study uses chi=0 and does not establish that the chi=0.35 suppression operates through the same route.

This improves the explanation of when the activator–inhibitor equations can organize initially similar chemical states: pattern initiation depends on whether the changing exchange network preserves a growth opportunity long enough. It does not demonstrate biological cell identity, autonomous state memory, inheritance, or development from a fresh carry zygote. The mature starts inherit no-carry development. Subsequent axis ratios remain close to their already anisotropic starts; rescued chemical organization is not a newly generated large shape asymmetry. Spatial/developmental convergence, geometric conductance closure and the broader ledger carry audit remain separate unresolved work.

The natural next mechanistic control is to cross directional tension (chi=0/0.35) with evolving/fixed conductances on the existing histories, keeping dilution and changing capacities. Reuse qualified paths and include matched chi=0 fixed-conductance controls. Rescue at chi=0.35 would support conductance remodeling as a route for the additional polarity-associated suppression; nonrescue would retain volume dependence and other moving interactions as candidates. This is a proposed follow-up, **not a newly launched experiment**. Broader maintenance, response and context claims still need carry revalidation before expansion.

## Evidence

- [Original protocol](../outputs/moving-initiation-controls/protocol.json), SHA-256 `f25296cd2b58acbcbb5b06a0ec90ab9300c462e9e5a587266f22edad5fe60ddc`.
- [Original completed summary](../outputs/moving-initiation-controls/summary.json).
- [Read-only full final verification and script](../outputs/moving-initiation-controls-review/2026-10-05T113534Z/assessment.json).
- [Compact assessment and precise metrics](moving_initiation_controls_assessment.json).
- [Pre-review documentation snapshots](../archive/study-document-snapshots/2026-10-05-before-completed-initiation-controls-review/manifest.json).

Original protocols, kernels, histories, checkpoints, comparisons, endpoint assays, ledger decisions and reporting rules remain unchanged. The new moving paths took **8.33 hours** of measured sequential GPU runtime. Raw evidence and figures under outputs remain local and ignored by Git.
