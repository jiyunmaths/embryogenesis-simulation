# Moving-geometry survival: timestep and developmental replication

Status: **completed; both timestep comparisons pass, and both additional histories retain patterns in their paired moving continuations**. The completed seed-7 experiment retained chemical contrast and initial-cell association in both moving branches from model time 90 to 150. This study tests whether that result depends on the timestep or the developmental history. See [the original survival experiment](feedback_survival.md).

## Fixed experimental design

| Study | Starting state | Timestep | Interval | Branches |
|---|---|---:|---|---|
| Existing reference | Actual seed-7 developmental checkpoint at time 90 | 0.0075 | 90–150 | Switch feedback on; keep feedback off |
| Temporal refinement | Exactly the same physical checkpoint | 0.00375 | 90–150 | Same paired branches |
| History replication | Fresh zygotes, independently developed with seeds 8 and 9 | 0.0075 | 0–90 development, then 90–150 | No-feedback development, then paired switch/control |

All studies use the original conservative signaling operator, 72³ grid, domain extent 2.24, interface width 0.085, sixteen-cell cap, and original mechanical coupling coefficients. No chemical reseeding, frozen preconditioning, or prescribed cell identities are introduced. Random seeds vary development itself, including division histories. Each paired continuation starts from the same checkpoint within its own history.

The refinement changes the integrator timestep and rescales the integer step counter so that the next step advances from physical time 90, rather than inadvertently restarting the clock at 45. Fields, chemical concentrations, polarity, lineage, division schedules, and random-stream states are preserved. Observation times remain separated by 0.6 model units.

## Assessment fixed before results

For the refinement, match observations by physical time and cell ID. Let the normalized initial cell volumes be the weights:

$$
w_i = \frac{V_i(90)}{\sum_j V_j(90)}.
$$

The chemical discrepancy combines activator and inhibitor on a logarithmic scale:

$$
E_{\mathrm{chem}}(t) =
\sqrt{\frac{1}{2}\sum_i w_i\left[
\left(\log\frac{a_i^{\mathrm{coarse}}(t)}{a_i^{\mathrm{fine}}(t)}\right)^2+
\left(\log\frac{b_i^{\mathrm{coarse}}(t)}{b_i^{\mathrm{fine}}(t)}\right)^2
\right]}.
$$

Here (a_i) and (b_i) are activator and inhibitor concentrations. The logarithm measures proportional discrepancies; weighting prevents small cells from disproportionately influencing the comparison. The practical acceptance limit is a maximum error of **0.02 across all sampled times** (roughly a 2% typical multiplicative concentration discrepancy), separately for each branch. This is an operational tolerance, not a derived error bound.

Aggregate shape is compared using its principal-axis ratio (R): the maximum of |R(coarse)/R(fine) − 1| must be at most **1%**. This checks overall elongation, not every local surface feature. Both chemical-contrast survival and initial-cell association classifications must also agree between timesteps. Two timesteps test sensitivity at these tolerances; they cannot establish an asymptotic convergence order.

The original survival criteria remain unchanged: over time 135–150, log-activator standard deviation must remain above 0.1, and initial-cell log-activator correlation must remain at least 0.8. Both arms must pass numerical quality checks: per-cell volume error below 5%, equivalent cell radius at least four grid spacings, zero phase-field clipping, and sampled boundary occupancy below 0.01.

For developmental replication, every planned seed is reported. An embryo without a mature sixteen-cell state or sufficient initial contrast at time 90 is recorded explicitly; it is not replaced by another seed. Survival is conditional on having a pattern to challenge, so report both the number of patterned starting histories and the number retaining their patterns. Numerical failures remain a separate category. Seeds 7, 8, and 9 form a small pilot cohort, insufficient for precise population-level probabilities.

## Execution and outputs

```bash
python -m embryo.feedback_survival_validation prepare
python -m embryo.feedback_survival_validation run
```

The prepared experiment is in `outputs/feedback-survival-validation/`. The driver runs refinement and the developmental cohort concurrently, using at most four new mechanics workers. The cohort develops one seed at a time and automatically runs its paired survival continuations before advancing to the next seed. Existing independent experiments continue unchanged.

- `protocol.json`: seeds, tolerances, input and source hashes.
- `refined/`: fine-timestep paired histories, checkpoints, and comparison.
- `seed-8/` and `seed-9/`: development and paired survival outputs.
- `refinement.json`: matched-time continuous errors and acceptance gates.
- `cohort.json`: all attempted developmental histories, including unsuccessful cases.
- `summary.json`: final combined results, including the original seed-7 comparison.
- `status.json` and `run.log`: progress and failures.

Checkpoints are written every six model-time units, and completed work can be resumed with the same `run` command after the previous driver has stopped. Hash verification rejects changed source or inputs. A completed process does not itself mean the scientific acceptance gates passed; inspect `refinement.json` and every cohort outcome.

This experiment does not resolve spatial-grid convergence, the geometric transport closure, indefinite moving-system stability, or biological cell-type identity. It asks a narrower question: is finite-horizon survival of an established chemical pattern robust to timestep halving and a small set of independently developed geometries?

## Completed assessment

Both seed-7 timestep comparisons pass all declared gates. Seeds 8 and 9 complete development and both survival branches, retaining contrast and initial-cell association. All three planned histories therefore reproduce finite-horizon maintenance. See the [completed assessment](feedback_completed_assessment.md) for values, denominators, and scope. The quality checks do not resolve spatial convergence, the live conductance closure, or biological identity.
