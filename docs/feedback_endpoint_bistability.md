# Chemical attractor coexistence after moving-geometry survival

The t=90–150 switch experiment retains the developed chemical pattern with both feedback on and feedback off. Both final graphs nevertheless stabilize homogeneous chemistry. This follow-up tests whether each final geometry supports **at least two locally stable chemical attractors**: a uniform state and a patterned state.

## Controlled test

The experiment restores each completed t=150 checkpoint, verifies that the moving study completed with passing quality checks, and freezes its measured geometry and conservative operator. The original Gierer–Meinhardt reactions and diffusivities are retained. There is no change of parameters between the two initial-state families on a given graph.

Each endpoint graph has 42 trajectories:

- The actual developed chemical state, continued without perturbation.
- The exact homogeneous reference (1,1).
- Twenty positive perturbations of the developed state, with log-noise scale 0.01.
- Twenty positive perturbations of the homogeneous state, with log-noise scale 0.001.

Perturbations are separately rescaled to preserve each species' initial total amount. The uniform and patterned families themselves may have different amounts; reaction dynamics do not conserve these totals. For each developmental seed, these are repeated interventions on two endpoints of one history, not independent embryos. The assay has now been completed for developmental seeds 7, 8, and 9. No labels or population count are supplied.

Every trajectory runs for 240 chemical time units. The declared criteria require maximum chemical derivative below 1e-6 at the endpoint, negative largest real eigenvalue of the complete two-species chemical Jacobian, and numerical agreement. Patterned trajectories must retain log-activator SD above 0.1 and return within volume-weighted log RMS 1e-4 of the unperturbed patterned continuation. Uniform starts must return within maximum absolute log deviation 1e-4 of (1,1).

The patterned-state Jacobian is evaluated directly, rather than applying the homogeneous Turing dispersion relation to a nonuniform state:

$$
J_{\mathrm{chem}}=
\begin{pmatrix}
\operatorname{diag}(2a_i/b_i-1)+D_a\Delta_V &
-\operatorname{diag}(a_i^2/b_i^2)\\
2\beta\operatorname{diag}(a_i)&-\beta I+D_b\Delta_V
\end{pmatrix}.
$$

This matrix contains all cell-to-cell chemical coupling on the frozen graph. It excludes derivatives with respect to phase fields, geometry, and polarity, so it is not the stability matrix of the complete moving simulation.

## Results

| Developmental seed | Endpoint graph | Patterned control: final log-activator SD | Patterned-state largest real eigenvalue | Uniform-state largest real eigenvalue | Perturbed patterned starts returning | Perturbed uniform starts returning |
|---|---|---:|---:|---:|---:|---:|
| 7 | Feedback switched on | 1.18186 | -0.30409 | -0.18052 | 20/20 | 20/20 |
| 7 | Feedback kept off | 1.10734 | -0.16378 | -0.08755 | 20/20 | 20/20 |
| 8 | Feedback switched on | 1.20362 | -0.37154 | -0.18244 | 20/20 | 20/20 |
| 8 | Feedback kept off | 1.14039 | -0.21548 | -0.09654 | 20/20 | 20/20 |
| 9 | Feedback switched on | 1.28962 | -0.58994 | -0.12015 | 20/20 | 20/20 |
| 9 | Feedback kept off | 1.22481 | -0.41868 | -0.04669 | 20/20 | 20/20 |

Seeds 8 and 9 were tested with exactly the seed-7 protocol: unchanged perturbation scales, integration horizon, numerical tolerances, and acceptance thresholds. All **168 new trajectories** pass, giving **252 passing trajectories on six endpoint graphs from three developmental histories**. The perturbed starts are basin probes, not additional developmental replicates. No seeds or trajectories were excluded.

The small derivative residuals and negative chemical Jacobian eigenvalues support locally stable uniform and patterned equilibria on every endpoint graph. This is numerical evidence for at least two attractors per graph; it does not show that exactly two exist or determine their complete attraction basins.

![Seed 7: uniform and patterned starts approach distinct chemical equilibria on the same endpoint graphs.](images/feedback-endpoint-bistability.png)

DOP853 at relative/absolute tolerances 1e-8/1e-10 is compared with 1e-11/1e-13 for every trajectory. Across all three histories, maximum absolute log discrepancy is 2.79e-6, below the unchanged 1e-5 threshold. An independent Radau integration checks the most tolerance-sensitive initial state on each graph; maximum discrepancy is 2.05e-9. Maximum endpoint chemical derivative is 7.55e-10, below the 1e-6 limit. All stored trajectories were checked for positivity and finiteness, their endpoint contrast was independently recomputed, and recorded source/checkpoint hashes were verified. A finite-difference unit test previously validated the patterned Jacobian and checked its homogeneous limit against the dispersion relation.

## Mechanistic interpretation

The results distinguish **loss of a linear pattern-initiation mechanism** from **loss of an existing nonlinear patterned state**. When the uniform state becomes stable, small perturbations can no longer grow by that linear instability, but another stable chemical state can remain available. Initial chemical history then determines which state the system approaches on the same fixed graph.

This explains why spectral stabilization can coexist with the observed moving-geometry survival. It also narrows the earlier feedback-suppression conclusion: the tested feedback suppresses formation in one developmental history, but does not universally erase established organization. A cyclic parameter sweep would be required to demonstrate hysteresis explicitly; coexistence alone is not such a sweep.

The moving experiment supplies finite-horizon persistence evidence, while this frozen-endpoint experiment identifies a chemical mechanism consistent with it. Neither proves stability of the full moving system or autonomous single-cell identity. Replication now supports attractor coexistence across the three tested developmental histories, but does not establish generality across model parameters or arbitrary geometries. The moving switch study has also completed seed-7 timestep refinement and seeds 8–9 replication; see the [completed feedback assessment](feedback_completed_assessment.md).

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.feedback_endpoint_bistability \
  --output outputs/feedback-endpoint-bistability-repeat
```

The runner requires completed, quality-passing switch-study outputs and refuses to overwrite an existing directory. It records source/input hashes and thresholds before evolution. `results.json` contains every trajectory's checks; the per-graph NPZ files retain starting states, complete sampled chemical trajectories, operator, volumes, and cell IDs. The production moving trajectories remain unchanged.

To repeat either added developmental history, choose a fresh output directory (example for seed 8; substitute 9 for seed 9):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.feedback_endpoint_bistability \
  --source outputs/feedback-survival-validation/seed-8/survival \
  --output outputs/feedback-endpoint-bistability-seed-8-repeat
```

The new results are stored in `outputs/feedback-endpoint-bistability-seed-8/` and `outputs/feedback-endpoint-bistability-seed-9/`, each with a trajectory figure and complete per-trial records. `outputs/feedback-endpoint-bistability-cohort-summary.json` aggregates all six endpoints and hashes their input reports and trajectory files. Simulation outputs are local artifacts excluded from Git.
