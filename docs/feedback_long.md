# Long matched feedback-component continuations

This study follows the short feedback-mechanism screen with longer moving-geometry controls. It asks separately which mechanical terms affect the emergence of chemical differences and which affect their persistence. The production study has completed all ten continuations in `outputs/feedback-long`, with all declared numerical-quality checks passing. The completed assessment is summarized below.

## Matched geometry and two chemical starting states

All runs start from the same no-feedback attribute-development checkpoint at t=18: sixteen resolved cells on a 72³ grid, domain half-width 2.24, interface width 0.085, and dt=0.0075. There are no active divisions and no further cell divisions. Geometry, measured volumes, lineage, and polarity are identical across all starts. This source has a positive maximum homogeneous-state spatial growth rate of approximately 0.20160, giving an initiation test a meaningful opportunity to grow before the geometry evolves.

The two starting conditions are:

- **Formation:** activator and inhibitor are each perturbed around the homogeneous equilibrium of one, using seed 7 and volume-weighted RMS amplitude 0.001. Each perturbation has zero volume-weighted mean.
- **Persistence:** those same initial perturbations evolve for 240 units on the frozen source graph. Their resulting patterned concentrations are transplanted back onto the common t=18 mechanical state. The chemical preparation does not advance the mechanical clock or change polarity.

The persistence preparation has final log-activator SD 0.95379, maximum chemical derivative 1.72e-11, and maximum concentration log discrepancy 2.71e-6 between standard and tightened DOP853 integration. It passes the predeclared contrast, stationarity, and accuracy gates. This is a controlled preconditioning procedure, not a pattern that developed during the moving trajectory.

The initial chemical means and amounts can differ between formation and persistence because the preconditioning reactions create and consume chemicals. Consequently, the two conditions compare different chemical states, not spatial contrast alone at identical chemical means. Within each condition, all mechanical arms have identical starting chemicals and amounts.

## Five mechanical arms

The original linear activity-dependent material law is retained; no exponential extension or larger coefficient is used.

| Arm | Activity–tension contrast | Activity–adhesion contrast | Polarity–tension contrast |
|---|---:|---:|---:|
| Baseline | 0 | 0 | 0 |
| Tension only | 0.25 | 0 | 0 |
| Adhesion only | 0 | 0.35 | 0 |
| Polarity only | 0 | 0 | 0.35 |
| Full coupling | 0.25 | 0.35 | 0.35 |

Polarity dynamics remain active in every arm. The polarity-only intervention switches its mechanical action; it does not introduce a different initial polarity. Every run uses the same mechanics code path with the indicated coefficients, conservative chemical transport, and amount-balanced dilution under volume changes. There is no downstream fate switch, fate noise, or A/B classification.

Five arms under two chemical starting conditions give ten matched continuations. Each runs for 60 units, from t=18 to t=78, with observations every 0.6 units and checkpoints every six units. This is one geometry and one chemical perturbation seed, not a full zygote-to-embryo ablation ensemble.

## Measurements and declared checks

Each observation saves the five continuous per-cell attributes, their context, aggregate shape, and the instantaneous transport spectrum. The primary endpoint is minimum across-cell log-activator SD over the final fifteen units (t=63–78). A value above 0.1 constitutes sustained chemical contrast for this screen. For formation it measures successful emergence and maintenance; for persistence it measures retention of contrast. It does not establish cell identity or a number of cell types.

The comparison also reports late mean contrast and aggregate axis ratio. Full time series of spectral growth indicate whether changes in chemical contrast accompany entry into or exit from the homogeneous-state instability band. Spectral stability of uniform chemistry is not equivalent to stability of an established nonlinear pattern. Likewise, persistent contrast can occur while individual cell states change; attribute histories must be inspected before inferring cellular memory.

Numerical quality is checked separately: maximum individual volume error below 5%, minimum equivalent radius at least four grid spacings, zero clipping, and sampled boundary occupancy below 1%. Volume, radius, and clipping are audited every step. An arm stops if it exceeds a numerical quality limit; the result is reported as failed rather than treated as biological loss of pattern. Boundary occupancy is checked at observations.

The protocol freezes code/input hashes, starting-state preparation, coefficients, observation times, and thresholds before production. Timestep and spatial refinement of this long study are subsequent validation tasks; the selected short-screen timestep check does not certify these longer trajectories.

## Execution and restart

Two workers execute independent arms. Every arm writes `status.json` and `history.json`; atomic `latest_state.npz` checkpoints contain the simulation state together with the matching audit history, arm identity, and protocol hash. A restart restores that internally consistent bundle and recomputes only work after its last saved checkpoint. Completed arms are reused. Code or input changes are rejected by the frozen hashes.

When all ten arms finish successfully, the runner automatically writes `comparison.json`, `RESULTS.md`, and `comparison.png`. Failed arms are listed in the top-level status file, and the study is not reported as complete with passing results. The mechanism still needs scientific assessment even if every numerical check passes.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.feedback_long prepare \
  --output outputs/feedback-long-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.feedback_long run \
  --output outputs/feedback-long-repeat --workers 2
```

The same `run` command resumes a stopped study, provided its original code and inputs remain unchanged. Do not start two drivers against the same output directory simultaneously. Ten relevant tests pass, including matching of initial geometry/polarity across arms and bitwise-exact checkpoint continuation with rejection of mismatched protocol or arm metadata.

## Interpretation boundaries

This experiment can attribute differences among matched mature-geometry trajectories to the selected mechanical terms. It cannot establish that the same term dominates embryogenesis from the zygote, where cleavage and geometry histories may diverge. Independent seeds, full developmental ablations, and refinement remain necessary. The live diffuse contact-conductance approximation is also unchanged.


The two-step 72³ smoke run completed all ten arms and passed every numerical quality check, including final comparison and checkpoint export. Its tiny time horizon verifies workflow only and is not evidence of successful formation or persistence. Production runs use the full predeclared sixty-unit horizon.


The [t=90 feedback-switch experiment](feedback_survival.md) separately addresses survival of the actual developmentally produced pattern, without chemical preconditioning. It runs alongside this component study and directly tests whether turning feedback on destroys an existing pattern while geometry moves.

## Completed assessment

All five preconditioned persistence arms retain contrast over time 63–78. Formation passes with baseline, tension-only, and adhesion-only mechanics, but fails with polarity-only and full coupling. Polarity-only therefore reproduces the qualitative formation suppression on this matched initial geometry; the [tension-plus-adhesion control](feedback_polarity_ablation.md) is now running to test whether polarity mechanics is necessary in the full combination. These are one-geometry intervention results, not fresh developmental replicates.

| Mechanical arm | Formation | Persistence |
|---|---:|---:|
| baseline | 0.12322 | 1.14875 |
| tension | 0.12593 | 1.14786 |
| adhesion | 0.12289 | 1.16040 |
| polarity | 0.02212 | 1.20381 |
| full | 0.02243 | 1.21569 |

See the [28 September completed-component and survival reassessment](feedback_reassessment_2026-09-28.md) for numerical values, interpretation, independent seed-8 survival, and remaining limitations.
