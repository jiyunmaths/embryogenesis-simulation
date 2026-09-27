# Dynamic signaling-cutoff interventions

The frozen contact-cutoff screen retained instability counts but failed its quantitative growth-rate tolerance at four developmental times. This follow-up tests whether that sensitivity changes the actual coupled trajectory.

## Isolated intervention

All three branches restore the same conservative, boundary-cleared 56³ checkpoint at time 30, including cell fields, regulator concentrations, fate states, polarity, lineage, and random-generator states. Signaling cutoffs are 0.01, 0.02, and 0.04. Polarity neighbor filtering is explicitly fixed at 0.02 in every branch. Diffusivities, geometry parameters, domain half-width 2.24, interface width 0.085, and time step 0.015 remain unchanged. Each branch runs through time 90, covering the baseline's signal-growth and subsequent decay phases.

The checkpoint already contains 12 A and 4 B cells. This is an identity-retention and signal/shape-sensitivity experiment, not a test of initial differentiation. Sixteen cells already reach the configured cap, and no division is pending. No perturbation or signal clamp is added.

`Config.polarity_contact_cutoff` defaults to **−1**, meaning inherit `graph_contact_cutoff`. Existing checkpoints and default trajectories therefore keep the historical shared-cutoff behavior. A value in [0, 1) explicitly sets the polarity cutoff. The experiment uses 0.02. This separates the direct filtering intervention; polarity can still change downstream in response to changed signals and geometry, as intended by the coupled model.

## Declared decisions

Observations every 0.6 model-time units record individual activator/inhibitor concentrations, continuous fate, identity labels, polarity vectors, and geometric/quality metrics. Relative shape and identity comparisons use the same cell IDs at matched times. For each perturbed branch versus the 0.02 baseline, all sampled comparisons must satisfy:

- Activator and inhibitor per-cell RMS differences below **0.01**, using homogeneous equilibrium concentration 1 as the reference scale.
- Absolute differences in cell-to-cell concentration standard deviations below **0.005**.
- Relative axis-ratio difference below **1%**.
- Maximum per-cell continuous-fate difference below **0.05**.
- No identity-label disagreements at any sampled time.

Numerical quality is a separate set of decisions: boundary occupancy below 1% at observations; cell-volume error below 5%, no phase-field clipping, and minimum equivalent radius of at least four grid spacings checked every step. The source already narrowly fails the radius criterion (approximately 3.969 spacings), so this study cannot pass the full quality screen. It can still measure paired sensitivity, but cannot establish spatially converged developmental robustness. No threshold is relaxed to hide this limitation.

These are prospective practical tolerances, not biologically fitted ones. Agreement over this finite window and one seed does not establish ensemble robustness or long-time persistence.

## Reproduction and artifacts

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cutoff_dynamics prepare \
  --checkpoint outputs/domain-conservative/grid-56/final_state.npz \
  --until 90 --output outputs/cutoff-dynamics-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cutoff_dynamics run \
  --output outputs/cutoff-dynamics-repeat
```

The completed study directory is `outputs/cutoff-dynamics`. The protocol records the source checkpoint and model hashes before starting. Branches run sequentially, beginning with the baseline. `status.json` and each branch's `history.json` update every observation; progress checkpoints are saved every six model-time units. Completed branches include final checkpoints, full sampled histories, and viewers explicitly showing two endpoint surfaces. After all branches finish, the runner writes `comparison.json`, `RESULTS.md`, and a four-panel trajectory plot, separating sensitivity outcomes from quality failures. Errors set the status to `failed`.

The production run is complete. All five sensitivity checks pass; the only failed quality check is the previously identified minimum cell radius. Tests establish that an explicit baseline polarity cutoff reproduces inherited-cutoff fields bit for bit, that cutoff routing is independent, and that checkpoint continuation and a short three-branch workflow work.

## Follow-up

Inspect the predeclared sensitivity outcomes without suppressing the known resolution failure. If trajectories differ, identify when signals diverge and whether fate or geometry follows. Any claim about developmental robustness requires refined developmental controls and additional seeds; this fixed-grid intervention does not replace them. Testing a shared signaling/polarity cutoff is a distinct experiment and must be labeled accordingly.

## Completed results

Across both perturbed branches, maximum activator RMS difference is **0.002027**, maximum activator-contrast difference is **0.001961**, maximum relative axis-ratio difference is **1.186e-6** (0.000119%), and maximum continuous-fate difference is **0.002352**. All sampled cell identity labels match the baseline. Boundaries, volumes, and clipping pass; minimum cell radius does not. The next [controlled cleavage refinement](cleavage_resolution.md) tests division numerics before full developmental refinement.
