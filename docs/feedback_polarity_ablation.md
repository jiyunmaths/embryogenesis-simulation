# Tension and adhesion together without polarity mechanics

Status: **completed; removing polarity mechanics restores formation**. Both arms pass all declared numerical-quality screens. This completes the missing combined control in the [matched mechanical-component study](feedback_long.md).

The completed study found sustained formation with baseline, tension-only, and adhesion-only mechanics, but not with polarity-only or full coupling. Polarity mechanics was sufficient to reproduce suppression on that initial geometry. It remained possible that tension and adhesion together would also suppress formation without polarity mechanics.

## Intervention

Restore the exact same no-feedback developmental time-18 geometry and the saved formation/persistence chemical arrays from `outputs/feedback-long/`. Retain direct activity-dependent tension and adhesion at their original contrasts, 0.25 and 0.35. Set only polarity-tension contrast from 0.35 to **zero**.

Polarity still evolves, including its chemical modulation. This intervention removes its mechanical action; it does not clamp chemistry or remove the polarity state. Geometry, transport, volume dilution, random streams, lineage, and all other coefficients match the corresponding full-coupling initialization. No further cleavage occurs at the existing sixteen-cell cap.

The two arms run from time 18 to 78 on the original 72³ grid, domain half-width 2.24, interface parameter 0.085, and timestep 0.0075:

- **Formation:** the identical near-uniform perturbation used in all original formation arms.
- **Persistence:** the identical pattern prepared on the frozen time-18 graph and used in all original persistence arms.

This is a matched mature-geometry intervention, not a new zygote history. The persistence start was chemically preconditioned and differs from the actual time-90 developed-state switch assay.

## Assessment fixed before results

The primary criterion is unchanged: minimum across-cell SD of log activator must exceed 0.1 at every observation in time 63–78. The persistence arm supplies the corresponding maintenance comparison.

- If removing polarity mechanics restores formation while full coupling suppresses it and baseline permits it, the result supports **necessity for suppression by this particular full combination on this starting state and horizon**.
- If formation remains suppressed, tension and adhesion together can suppress it without polarity mechanics; polarity is not necessary in this setting.
- Numerical failure or inconsistent reference outcomes make the causal comparison inconclusive.

Neither outcome establishes general necessity across parameters, developmental histories, or geometries. The test also does not resolve whether chemical modulation of polarity is needed; that requires a separate intervention on the polarity activation factor.

## Execution and checks

```bash
python -m embryo.feedback_polarity_ablation prepare
python -m embryo.feedback_polarity_ablation run
```

The driver runs the two new arms in parallel, checks frozen source and input hashes, and automatically compares results with all five completed reference arms. It reuses the original stepping, observation, numerical-audit, and checkpoint/restart code. The extra arm is registered only inside its worker process; existing simulation source files and running studies are unchanged.

Output directory: `outputs/feedback-polarity-ablation/`. Per-arm histories/statuses update every 0.6 model-time units, with restart checkpoints every six units. The completed driver will write `comparison.json`, `RESULTS.md`, and `comparison.png`. Re-running `run` resumes a stopped driver; do not run two drivers against the same directory.

Quality criteria remain maximum per-cell volume error below 5%, radius at least four grid spacings, zero clipping, and sampled boundary occupancy below 0.01. These are quality screens, not spatial or temporal convergence certification of this new control.

Four relevant tests pass, including exact single-coefficient ablation, preservation of starting fields and random streams, real worker execution/restart, unchanged original arm registration, and alternative outcome classification.

## Completed outcome

Formation has late minimum log-activator SD 0.12558 (threshold 0.1), compared with 0.02243 under full coupling. Prepared-pattern persistence also passes, with minimum 1.15970. Removing polar mechanics therefore rescues formation while retaining tension and adhesion. Together with the polarity-only control, this supports necessity and sufficiency for suppression in this single tested setting, not general necessity or independence from chemical modulation. See the [completed assessment](feedback_completed_assessment.md).
