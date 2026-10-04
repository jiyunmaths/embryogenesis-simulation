# Three-timestep phase-carry check and independent mechanics reference

The [completed formation pair](phase_carry_formation.md#completed-formation-assessment) passes its original quantitative limits, but two nearby timesteps cannot establish a convergence trend. This follow-up adds **dt = 0.00375** to the completed **0.001875 and 0.0009375** paths. It also checks the mechanics against a separate float64 CPU calculation. **The study is complete and passes its declared checks.** Its scope remains one mature-state history at zero directional-tension contrast.

## Shared state and units of replication

All three moving paths start from the same mature 16-cell geometry at physical time 150 in history 9, with the same polarity, lineage, random streams and near-uniform chemical reset. They run for 240 elapsed units, ending at physical time 390. The grid is 72³, $D_a=0.02$, $D_b=0.55$, $\beta=2$ and directional polarity-tension contrast $\chi=0$. Activity-dependent tension and adhesion remain 0.25 and 0.35. Polarity itself continues to evolve in the moving runs.

There is **one existing developmental history**, three nested numerical paths, and **no new histories**. The two completed paths are reused read-only. Only the coarse carry path is new. Its residual starts at zero at the shared initial state; later checkpoints retain it. The earlier no-carry failures are unchanged.

## Independent float64 mechanics check

Before launching the moving trajectory, the runner compares the existing CUDA mechanics with [a separate NumPy/SciPy reference](../embryo/mechanics_float64_reference.py). The reference stores phase fields and evaluates volumes, centroids, overlaps, coefficients, forces and updates in float64. It uses an independently evaluated no-flux face operator. Tests check both face conservation and the negative gradient of the discrete energy, holding the directional tension field fixed during differentiation as required by this explicitly coupled model.

Two actual geometries are tested: the shared initial geometry and the finest completed patterned endpoint. Each receives a 0.15-unit **mechanics-only** evolution at all three timesteps. Concentrations and polarity are held fixed; cell shapes, overlap forces and volume penalties evolve. Baseline and carry CUDA paths are compared with the CPU field after every step. Component-assay residuals start at zero on the same visible fields.

The predeclared carry limits are maximum phase error $5\times10^{-6}$, relative field norm error $2\times10^{-6}$, relative volume error $2\times10^{-6}$, and zero clipping. The adjacent float64 reference field discrepancies must decrease by at least a factor of 1.5. These six short comparisons do **not** supply a high-precision reference for the full co-evolving chemical/polarity system.

## Moving-run decision

Both adjacent moving pairs must pass every inherited limit, including maximum absolute log chemical discrepancy 0.01, transport discrepancy 0.01, volume discrepancy 0.005, and matching persistent-contrast classifications and crossing times. The ratio of coarse/fine to fine/finer maximum chemical discrepancies must be at least 1.5. This ratio describes the observed three-point trend; it is not an absolute solution-error estimate or proof of asymptotic order.

The new endpoint also receives the original dual-solver frozen-graph assay: the developed state, two small patterned perturbations and two near-uniform perturbations. Stationarity and numerical checks remain required. Local bistability is assessed rather than imposed as a passing outcome.

## Evidence and follow-up

The frozen design is `outputs/phase-carry-convergence/protocol.json`; root/child `status.json` files report progress. Sources, inputs, software versions and CUDA binaries are verified before execution. The coordinator uses an exclusive lock and resumable carry-preserving checkpoints. Reference field files and results have separate hash verification.

The exact-context GPU implementation gate has passed: 16 baseline steps match the accepted kernel bit for bit, the first zero-carry step matches it, and four checkpoint/restart steps preserve the state and residual exactly. CPU tests of the independent operator, energy gradient, reference timestep trend and retiming also pass. All six independent short mechanics comparisons and the complete three-timestep formation comparison now pass.

## Completed results

| Moving pair | Maximum absolute log chemical discrepancy | Original limit | Decision |
|---|---:|---:|---|
| 0.00375 versus 0.001875 | 0.000409631 | 0.01 | Pass |
| 0.001875 versus 0.0009375 | 0.000205004 | 0.01 | Pass |

The adjacent-error ratio is **1.99816**, giving an observed pair order of 0.99867. This is consistent with first-order temporal behavior over the tested steps and window; it is not proof of asymptotic convergence. Both pairs also pass the original transport, polarity, volume, shape, growth and crossing-time limits. All three paths form sustained chemical contrast, with final log-activator spread about 1.58318.

The independent float64 mechanics comparisons have maximum carry field error **5.32e-8**, below the declared 5e-6 limit. Their reference field discrepancies decrease by factors 2.01772 and 2.00276 under halving, respectively, on the initial and patterned geometries. These remain 0.15-unit mechanics-only checks with chemistry and polarity fixed. The new coarse frozen endpoint passes the dual-solver and stationarity checks and supports local uniform/patterned coexistence.

The completed comparisons, checkpoints, sources, inputs, reference fields and endpoint provenance were independently rechecked before the next launch. The prior no-carry numerical failures remain failures.

The [matched polarity initiation study](phase_carry_polarity.md) is now running across histories 7, 8 and 9, with zero and positive directional-tension contrasts at two timesteps. Full-system high-precision equivalence, spatial convergence, initiation from fresh zygotes, and autonomous cell identity remain separate questions. This qualification does not automatically promote the backend or expand the parameter map.
