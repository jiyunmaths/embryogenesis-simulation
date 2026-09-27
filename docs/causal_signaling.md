# Does activator–inhibitor feedback cause organization?

This first causal experiment separates **persistent signaling contrast**, **cell-identity differentiation**, and **shape**. It tests the first two on a frozen, resolved embryo graph. Shape causality is not tested here and remains a separate moving-geometry experiment.

## Controlled starting point

The source is the 72³ developmental checkpoint at time 18, after all sixteen cells have formed and no cytokinesis is pending. Its minimum equivalent cell radius is approximately 5.104 grid spacings, and it clears the boundary screen. The live conservative contact graph is frozen. This graph comes from an already organized geometry; freezing it does not establish that the geometry arose from signaling.

For each of twenty perturbation seeds, activator and inhibitor are reset to equilibrium 1 plus independently sampled Gaussian cell perturbations. Each species has zero volume-weighted perturbation mean and volume-weighted RMS 0.001. Every intervention receives the same perturbations for a given seed. Fate is reset to zero, so established labels are not inherited. No spatial axis is imposed, but the fixed graph itself is not rotationally symmetric.

These are twenty chemical perturbation replicates on **one geometry**, not twenty independent developmental embryos. The known overlap-to-conductance closure limitations remain. The experiment tests causality within this particular mathematical model.

## Interventions

All reaction changes retain the homogeneous equilibrium (a,h)=(1,1), avoiding a trivial comparison caused by changing basal equilibrium. With beta=2, the controls are:

| Arm | Activator reaction | Inhibitor reaction | Other change |
|---|---|---|---|
| Full loop | a²/h − a | beta(a² − h) | None |
| No self-activation | 1/h − a | beta(a² − h) | Remove activator dependence from its production |
| No induced inhibitor | a²/h − a | beta(1 − h) | Remove activator induction of inhibitor |
| No inhibitor action | a² − a | beta(a² − h) | Remove inhibition of activator production |
| No transport | Full | Full | Both diffusivities zero |
| Equal diffusion | Full | Full | Both diffusivities 0.02 |
| No signal-to-fate coupling | Full | Full | Fate forcing gain zero |

Frozen-mode growth is computed separately for each modified Jacobian and diffusivity pair. Removing induced inhibitor or inhibitor action makes the equilibrium locally unstable (maximum local growth +1). Those are dysregulation controls; their growing signals cannot be labeled diffusion-driven pattern formation.

Signaling uses positive SSP-RK2 stages with the same production/loss bound as the existing conservative integrator. Downstream fate obeys the existing deterministic drift

$$
\dot z_i=\rho\left[z_i-z_i^3+g(a_i-1)\right],
$$

with no fate noise, neighbor inhibition, or exposure bias in this screen. It uses SSP-RK2 with the updated activity held fixed for each fate step; this sequential coupling requires its own time-refinement check. It is not a claim of exact replay of the moving live model's Euler fate update.

## Prospective observations and criteria

Every arm runs for sixty model-time units, recording complete cell vectors every 0.5. Persistent signal contrast requires volume-weighted activator standard deviation at least 0.1 throughout the final twenty units, in at least sixteen of twenty trials. Suppression requires contrast below 0.01 throughout that window. Both-fate outcomes require at least one cell above +0.55 and one below −0.55 at the end. These separate criteria do not equate having two labels with spatial organization.

A regulator reaching 10 ends that replicate as a **finite dysregulation exit**, never as a clipped or successful pattern. Its exit time and value are recorded; later held values are not treated as evolving observations. A concentration ceiling is a diagnostic stopping rule, not part of the model's dynamics.

The initial dt=0.01 study includes dt=0.005 comparisons for the full loop, no self-activation, and equal diffusion. Limits are 0.01 maximum per-cell signal RMS discrepancy, 0.05 maximum continuous-fate discrepancy, and unchanged final A/B counts. All original decisions are retained.

## Completed first screen

| Intervention | Persistent signal contrast / 20 | Both fate labels / 20 | Completed / 20 |
|---|---:|---:|---:|
| Full loop | 20 | 20 | 20 |
| No self-activation | 0 | 20 | 20 |
| No induced inhibitor | 0 | 0 | 0 |
| No inhibitor action | 0 | 0 | 0 |
| No transport | 0 | 20 | 20 |
| Equal diffusion | 0 | 20 | 20 |
| No signal-to-fate coupling | 20 | 0 | 20 |

Full-loop final contrast ranges from approximately 0.465 to 0.537. The stable signaling ablations decay to contrasts of order 1e-15. Both inhibition knockouts reach the diagnostic ceiling between about times 6.26 and 7.28; their incomplete fate results are not evidence that inhibition is necessary for fate differentiation.

**The loop is necessary for persistent signal contrast among these tested controls, but it is not necessary for producing two fate labels in the present fate model.** The drift z−z³ makes zero fate unstable. Tiny transient activity differences can seed eventual commitment even after signals homogenize. Removing signal-to-fate forcing leaves exactly zero fate unchanged under these noiseless initial conditions, while the signal pattern itself remains. This supports a causal link from transient signaling to the switch, not a requirement for a sustained Turing pattern to generate two identities.

## Temporal limitation and separate confirmation

Signal trajectories agree under the initial step halving, but continuous-fate differences in the no-self-activation and equal-diffusion controls reach approximately 0.225 and 0.204, exceeding 0.05. Their final fate counts match, yet commitment timing/state accuracy is not established by that initial comparison. The original combined screen therefore reports **False**, rather than ignoring this failure.

A separate finer-time confirmation is specified before execution: dt=0.00125 versus 0.000625 for full, no self-activation, no transport, and equal diffusion, using the same twenty seeds, geometry, duration, perturbations, and unchanged tolerances. It preserves the original source report and checkpoint hashes. Results are stored separately in `outputs/causal-signaling-confirmation`.

## Reproduction and scope

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.causal_signaling prepare \
  --checkpoint outputs/development-refinement/space-72/state-18.npz \
  --output outputs/causal-signaling-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.causal_signaling run \
  --output outputs/causal-signaling-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.causal_confirmation \
  --source outputs/causal-signaling-repeat --output outputs/causal-confirmation-repeat
```

The initial production screen is `outputs/causal-signaling`, with per-arm trajectories, exact intervention Jacobians, spectral growth rates, per-seed outcomes, and explicit exits. No embryo model code or ongoing developmental run was changed.

Next, matched moving-geometry controls should test whether sustained signaling generates shape changes beyond mechanics and division history. Separately, test fate-switch sensitivity to transient pulses, noise, and the z−z³ bistability. Do not change that switch merely to force the activator–inhibitor hypothesis to appear necessary. This screen does not establish biological causality, an ensemble developmental result, or shape emergence.

## Completed finer-time confirmation

The independent confirmation also reports **False**, retaining its numerical limitation. All signals agree closely and all final fate counts match across the compared steps. However, maximum continuous-fate discrepancies are 0.00270 (full), 0.01620 (no self-activation), **0.23440 (no transport)**, and **0.08388 (equal diffusion)**; the latter two exceed 0.05. Full-loop persistent contrast and both final fate labels occur in 20/20 trials at both finer steps. Stable ablations retain both final labels in 20/20 while persistent signal contrast is absent.

These outcomes support the signal-pattern mechanism and reproduce the categorical fate result, but do not certify commitment trajectories or timing in all controls. The next numerical check for fate is a joint signal/fate integrator and an independent high-accuracy ODE reference, separating splitting error from sensitivity near the unstable uncommitted state. Preserve both failed temporal reports; do not relax the tolerance. Moving-geometry shape ablations are still needed before attributing geometry to the loop. Sixty-seven relevant regression tests pass.

## Resolution and moving controls

The subsequent [joint signal/fate validation](joint_fate.md) resolves the integration sensitivity against independent tightened DOP853 references. Maximum continuous-fate error at dt=0.0075 is 6.63e-5, and all five tested arms and twenty seeds pass the unchanged signal/fate tolerances with matching cell-wise final labels. The original failed reports above remain intact. Corrected trajectories retain persistent full-loop signaling, suppression in the stable signaling ablations, and both fate labels even in those ablations. Four matched [moving-geometry controls](moving_causal.md) are now running through t=78; shape-causality results are pending.
