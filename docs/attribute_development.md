# Development with continuous attributes and no prescribed identities

The first production experiment is running in `outputs/attribute-development`. It asks whether a developing regulatory–mechanical system produces distinct, persistent **phenotypes without a supplied two-state fate switch**. There is no requested number of types, no clustering, and no A/B labels in its observations or surface exports. Results are pending.

This is an experimental branch implemented in `embryo/attribute_development.py`. The historical core/dashboard model remains available for comparison, and its existing fate equation is still documented elsewhere. Neither that code nor the ongoing dense-geometry replay is changed by this experiment.

## What is removed, and what remains

There is no integration of the downstream drift $f-f^3$, no fate noise, no fate partition perturbation, and no fate thresholds used to classify observations. Activator and inhibitor remain continuous chemical variables; they are not cell types. A zero-valued `fate` array is retained solely for compatibility with base simulation bookkeeping/checkpoint layout. It cannot influence the new mechanical law and is omitted from phenotype records. A/B counts, uncommitted counts, and fate separation are removed from exported metrics.

The supplied ingredients still include chemical feedback, cell mechanics, polarity, division rules, and a new explicit constitutive coupling. Removing a bistable fate switch does not make every assumption emergent. In particular, the number of chemical species is fixed, cells are well mixed chemically, and divisions follow the existing shape/timer rules.

## Two fresh developmental branches

Both branches start from the same analytic radius-0.8 zygote, with initial activities equal to one, zero polarity, identical initial random streams, and seed 7. They use a 72³ grid, domain half-width 2.24, interface width 0.085, timestep 0.0075, and a 16-cell cap. They develop through time 90 without growth between divisions. Signals receive the existing small amount-balanced partition perturbations at cleavage. There is no inherited mature geometry, preassigned axis, chemical patch, or inherited fate pattern.

| Branch | Chemical and polarity dynamics | Mechanical coupling |
|---|---|---|
| `direct` | Full activator–inhibitor dynamics, conservative contact exchange, dilution, and polarity | Instantaneous continuous activity-dependent tension/attraction and existing polarity-dependent tension |
| `no_feedback` | Same equations and initial random streams | Constant baseline tension/attraction; polarity cannot change tension |

Different mechanics can change division timing, orientation, and the assignment of later random draws to descendants. Matching the starting seed does not ensure identical developmental histories. This comparison tests the combined regulatory mechanical feedback, not individual adhesion, tension, or polarity contributions.

## Direct continuous mechanical response

Define an instantaneous bounded activity response

$$
r_i=\tanh(a_i-1).
$$

It has no independent dynamics or stored state. Its zero is the homogeneous chemical reference, not an identity threshold. Mechanical coefficients are

$$
\gamma_i^0=\gamma_0(1+c_\gamma r_i),
\qquad
A_{ij}=A_0(1+c_A r_i r_j),
$$

with $\gamma_0=1$, $A_0=4$, $c_\gamma=0.25$, and $c_A=0.35$. Existing polarity-dependent directional tension is then applied in the direct branch. The contrast values are retained from the prior mechanics for an initial controlled comparison; the old configuration names `fate_tension` and `fate_adhesion` store these numbers, but no fate variable enters their calculation here.

This law assumes that activity changes material properties promptly and that responses of the same sign favor attraction. It can encourage activity-dependent rearrangement and must eventually be compared with alternative constitutive laws. It does not specify two stable identities or discretize activity into categories. Any apparent populations must be assessed independently of this assumed coupling. No-feedback mechanics is identical to the historical no-feedback kernel; for the direct branch, substitution of the same coefficient values reproduces the core mechanical update exactly in regression tests.

## What we measure

The primary per-cell attribute vector contains:

1. Log activator activity.
2. Log inhibitor activity.
3. Polarity magnitude.
4. Cell covariance axis ratio minus one.
5. Cell covariance asphericity.

Logs express chemical deviations relative to the unit reference scale. Polarity magnitude and shape descriptors are independent of absolute orientation; a global rotation alone should not create a new identity. These features are reported separately, without combining their different scales into a clustering distance.

Position, lineage, division state, exposure, measured/target volume, and raw contact weight are recorded as **context**. Tension and attraction are recorded as derived material properties, but excluded from the primary vector because they are already determined by activity in this constitutive law. Counting those derived properties as independent evidence would exaggerate chemical differentiation.

Observations are saved every 0.6 time units. Full checkpoints and cell-surface meshes are saved every six units. The final comparison reports per-feature cell-to-cell dispersion over the last 15 time units and the median within-cell attribute range over that window, provided cell IDs remain stable. It also reports aggregate shape and development completion. Stable feature differences may be explained by stable exposure or geometry, so these descriptive measurements cannot alone establish identities.

A valid outcome may be near-uniform cells, a continuous range of phenotypes, or candidate separated populations. No result is forced into two clusters. The next analysis can examine population structure and its relationship to context, followed by perturbation/relocation tests of persistence. Establishing cell identity will require reproducible responses and memory under those tests, plus independent seeds and numerical refinement.

## Numerical and workflow checks

The protocol freezes configuration, constitutive laws, feature definitions, source hashes, and quality limits before execution. Quality limits are maximum per-cell relative volume error below 5%, minimum equivalent cell radius at least four grid spacings, zero clipping, and sampled boundary occupancy below 1%. Completion requires 16 cells with no active division. These are acceptance checks for a pilot, not a certification of full developmental convergence; the historical trajectory study still has unresolved refinement discrepancies.

Thirty-four relevant tests pass, covering core development plus new material substitution, absence of fate drift/labels, exact checkpoint continuation in both branches, and progressive cleavage without daughter fate assignment. The parallel two-step 72³ smoke test passes volume, resolution, boundary, and clipping checks and completes export/comparison. It intentionally cannot pass the 16-cell development-completion check; that short run is not reported as a full biological experiment.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.attribute_development prepare \
  --output outputs/attribute-development-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.attribute_development run \
  --output outputs/attribute-development-repeat
```

The two branches run in separate worker processes. Per-branch `status.json` reports time. `history.json` contains attributes and context; `surface-*.json` contains meshes with continuous chemical values; `lineage.json` records divisions. Final `comparison.json`, `RESULTS.md`, and `comparison.png` are generated automatically. Resume an experimental checkpoint with `AttributeSimulation.restore`, not `Simulation.restore`, to retain the new equations. The old dashboard has not been converted to this experiment's attribute-based schema.
