# Does polarity still suppress initiation with phase-update carry?

The [completed three-timestep qualification](phase_carry_convergence.md#completed-results) passes its declared numerical checks. The next experiment asks whether the earlier suppression of chemical pattern formation by directional polarity mechanics persists after retaining small phase updates. It uses matched controls on three existing developmental histories, rather than adding new identities or nutrient assumptions.

## Design and starting conditions

Each history starts at physical time 150 with the exact same mature 16-cell geometry, polarity, lineage, random streams and near-uniform chemical preparation used in the earlier polarity study. Within each history, the only physical intervention is directional tension contrast $\chi=0$ versus $\chi=0.35$. Polarity dynamics remain active in both branches; zero contrast removes polarity's directional modulation of tension, not polarity itself. Activity-dependent tension and adhesion remain 0.25 and 0.35.

Both branches are checked at dt=0.00375 and 0.001875. The intended horizon is 240 elapsed units, ending at physical time 390. The grid remains 72³ with extent 2.24 and interface parameter 0.085; $\beta=2$, $D_a=0.02$ and $D_b=0.55$ are fixed. Visible fields and contact accumulation remain float32, while phase-update residuals are retained in float64.

| Design factor | Values |
|---|---|
| Existing histories | 7, 8, 9 |
| Directional tension contrast | 0, 0.35 |
| Moving timesteps | 0.00375, 0.001875 |
| Chemical start | Same near-uniform preparation within each history |
| Numerical paths | 12; two completed history-9 zero-contrast paths reused read-only |
| New moving jobs | 10 |

The replication unit is the **developmental history**: three existing histories, no new histories. The paths, cell measurements, chemical perturbations and numerical repeats are nested. This is a mature-state formation-opportunity assay, not fresh zygote development or a maintenance experiment.

## Gates before interpretation

1. Recheck the completed carry qualification, full accepted backend evidence, software/binary hashes and matched physical starts.
2. Test positive-contrast mechanics against the independent float64 CPU reference on the history-9 initial state. Three 0.15-unit mechanics-only comparisons use the previously declared field/volume/clipping limits and the reference timestep trend. This short check does not validate full-system long-horizon equivalence.
3. Require exact-context native/GPU baseline agreement. Six unchanged historical prefix gates can be reused after checking their original evidence and identical starting parameters; four new fine-step contexts receive their own 0.6-unit comparisons. All ten new jobs receive accepted/experimental baseline checks and exact carry restart checks.
4. Check each coarse/fine pair through 60 elapsed units. Only pairs passing the original raw refinement limits receive new long continuations. A failed pilot remains explicit and blocks that pair; other prespecified contexts can still be tested.
5. Assess complete 240-unit agreement and run the original dual-solver frozen assay on each actual endpoint. Endpoint classifications are measured outcomes; bistability is not imposed as a passing requirement.

Pilot observations are saved as immutable comparison files, separate from histories that later grow. Checkpoints preserve the numerical residual and co-evolving state. Old protocols, kernels and failed results are unchanged.

## Scientific decision

Sustained formation means log-activator spread exceeds 0.1 throughout the final 24 units. Record onset, shape, transport spectra, uniform-state growth and frozen endpoint behavior separately. The integral of positive frozen modal growth is a descriptive measure of the available instability window, not an amplification theorem for a changing network.

A history supports polarity-specific suppression only when the zero-contrast control forms sustained contrast and the positive-contrast branch does not, with passing numerical checks. If both branches remain uniform, that history lacks a successful initiating control and does not establish suppression. If both form patterns, or only the positive branch does, report that outcome directly. Supporting the hypothesis is not required for a numerical pass.

The [history-level reporting addendum and moving-geometry controls](moving_geometry_initiation_controls.md) fix the aggregate wording before the final three-history assessment. If the old outcome repeats with qualified carry comparisons, report polarity-dependent suppression in **1/3 tested histories**, and moving-versus-frozen initiation failure persisting without directional polarity tension in the other two. Their common physical cause remains unresolved. Do not treat the two unsuccessful zero-contrast controls as two extra demonstrations of polarity-specific suppression. The first proposed follow-up crosses dilution on/off with evolving/fixed conservative contact conductances in the eligible histories; it is specified but not launched. The active protocol and runner remain unchanged.

## Execution and scope

The study is launched in `outputs/phase-carry-polarity/`. Its `protocol.json` fixes sources, inputs, parameters and decisions before execution; root and child `status.json` files report progress. Moving jobs run sequentially on the qualified GTX 1080 Ti using PyTorch arrays/matrix operations and custom CUDA mechanics/geometry/polarity. Two separate CPU processes run frozen endpoint assays alongside moving work. The initial estimate for new moving computation is about 6.2 hours; gate and endpoint costs and machine load add uncertainty.

The runner is [phase_carry_polarity.py](../embryo/phase_carry_polarity.py). Eleven orchestration/reference tests pass, covering balanced design, proper interpretation of unsuccessful controls, failed-pilot blocking, immutable comparison records and carry-dependent checkpoint recovery. All three positive-contrast float64 component checks now pass, with maximum carry field error 4.56e-8. The first new history-9 positive-contrast context gate also passes, and its coarse 60-unit moving pilot has started. Moving outcomes remain pending.

Full developmental spatial convergence, physical transport closure, inheritance, autonomy, biological cell-type identification and nutrient-dependent growth remain outside this test. The ledger and manuscript retain their earlier reviewed evidence scope until these new results are assessed.
