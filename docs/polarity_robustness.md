# Fixed-ratio polarity test with longer observation

The [completed parameter study](parameter_robustness_assessment.md) found that developed chemical patterns survive moving geometry, while near-uniform starts remain weak. At the combined point $(D_b/D_a,\chi)=(27.5,0.7)$, the moving contact spectra leave the initial linear instability band. That comparison changed diffusivity and polarity contrast together. This follow-up holds diffusivity fixed and varies polarity contrast alone, with four times the observation window and targeted timestep checks.

This is a **mature-state initiation/maintenance experiment**, nested in **three existing developmental histories** (7, 8, 9). It does not add independent zygote histories or claim new biological cell identities. The code is [polarity_robustness.py](../embryo/polarity_robustness.py); the declared protocol and live status are under `outputs/polarity-robustness/`.

**Completed assessment:** all 21 moving jobs, context checks and endpoint assays finish. Developed patterns survive every contrast across all three histories; only history 9 initiates from its near-uniform start at $\chi=0$, with matching outcomes at both timesteps. The full study remains `completed_with_unresolved_checks`: its nonlinear formation transient exceeds the declared chemical timestep-error tolerance. [Results, figure and numerical qualification](polarity_robustness_assessment.md).

**Historical launch snapshot:** 40 targeted software tests pass and one opt-in test is skipped. All six longer frozen references complete with verified trajectory hashes. Near-uniform starts reach sustained contrast in all three frozen histories: first 0.1 crossings occur at approximately 177.61, 156.66 and 58.02 elapsed units for histories 7, 8 and 9, respectively. Maximum dual-solver chemical log error for those three paths is $1.75\times10^{-10}$. [Preparation and frozen-reference verification](polarity_robustness_verification.json) preserves the original pending moving/timestep status.

## Matched experimental design

| Quantity | Declared setting |
|---|---|
| Chemistry | Gierer–Meinhardt; $\beta=2$, $D_a=0.02$, $D_b=0.55$, ratio 27.5 |
| Polarity-tension contrast $\chi$ | 0, 0.35, 0.7 |
| Activity-dependent tension and adhesion | Existing laws and coefficients held fixed |
| Initial geometry and polarity | Exact mature $t=150$ direct-feedback state within each history |
| Initial chemistry | Exact parent near-uniform perturbation or prepared developed pattern, paired across $\chi$ and timestep |
| Moving and matched frozen observation | 240 model-time units, ending at $t=390$ |
| Coarse moving timestep | 0.00375; 18 jobs: three histories × three contrasts × two starts |
| Targeted fine moving timestep | 0.001875; three history-9 near-uniform jobs, one per contrast |
| Observation / checkpoint spacing | 0.15 / 3 model-time units |
| Sustained moving contrast | Across-cell SD of $\log a$ above 0.1 throughout the last 24 units |
| Pilot gate | All three history-9 coarse/fine near-uniform pairs must agree through 60 units before any continuation beyond 60 |

All times and coefficients retain the model's nondimensional conventions. Setting $\chi=0$ removes directional cortical tension while retaining polarity evolution and activity-dependent tension/adhesion. It does not remove all mechanical feedback. The mature initial geometry already carries its developmental mechanical history; this intervention does not undo that history.

The six unchanged coarse $\chi=0.7$ first-60-unit runs are reused from the parent study. Their parent protocol, histories, physical checkpoints, chemical starts, lineage, parameters and numerical checks are verified. The evolved state at $t=210$ is preserved; chemistry, geometry and polarity are not reset. New checkpoint metadata records the new protocol and transfer provenance.

The finer timestep is installed by retiming the discrete clock at fixed physical $t=150$. Spatial fields, chemistry, polarity, target volumes, cell order, lineage and random streams remain identical. History 9 is selected before new outcomes because it had the strongest initial uniform growth and largest frozen near-uniform contrast in the parent study. This targeted check covers near-uniform starts in that history, not patterned starts or all three histories.

## Frozen references and measured spectra

Six unique frozen references, one per history and chemical start, run to 240 with both DOP853 and Radau. They use the identical initial measured conservative operator and chemistry. Since frozen chemistry contains neither $\chi$ nor a mechanics timestep, these references are shared across those interventions rather than counted as separate evidence.

Each moving observation retains the actual conservative operator $\Delta=-M^{-1}K$ and its mass-symmetric spectrum. The analysis records the largest uniform-state modal growth rate, its first positive-to-nonpositive crossing, and the first upward crossing of chemical contrast 0.1. Crossing times are interpolated between recorded observations. These are descriptive timings, not proofs of coupled moving-system stability. A developed start is already patterned; its lack of a later upward crossing does not imply failed initiation.

Longer frozen references distinguish slow growth from missing pattern formation. Comparing moving branches at the same ratio isolates the intervention on directional tension. Comparing moving and frozen chemistry tests evolving geometry, transport and volume dilution together; it does not isolate dilution from transport.

## Numerical gates and restart behavior

The resident **PyTorch/custom CUDA** backend is retained. Both the accepted coarse and fine full-horizon baseline GPU validations must pass their original source, evidence, hardware and binary checks. Before every new job, a separate 0.6-unit native C++/GPU comparison tests its exact parameters and chemical start with the unchanged strict backend tolerances. These short checks do not assert full-horizon native/GPU equivalence at new parameters.

The 60-unit pilot and final 240-unit timestep comparisons use the same declared criteria:

| Maximum coarse/fine discrepancy | Tolerance |
|---|---:|
| Absolute chemical log error, over both species and all cells/times | 0.01 |
| Absolute polarity component error | 0.01 |
| Relative aggregate axis-ratio error | 0.01 |
| Relative volume error | 0.005 |
| Relative transport-matrix Frobenius error | 0.01 |
| Absolute largest modal growth-rate error | 0.001 |
| Spectral zero-crossing time error | 0.30 |
| Contrast-onset time error | 0.30 |

Both timesteps must agree on the presence of a spectral crossing and contrast onset, and on the late sustained-contrast classification. If both lack a crossing, its timing discrepancy is zero. A crossing in only one trajectory fails. The pilot late window is 12 units; the final late window is 24. Raw chemical log error is used rather than dividing contrast error by a near-zero contrast.

If any pilot comparison fails, the queue stops before longer scientific continuations. No tolerance or horizon is automatically relaxed. If full-horizon refinement disagrees, the completed study is marked as having unresolved checks. Fine/coarse disagreement is a numerical limitation, not a mechanism result.

Existing volume/radius, clipping, positive chemistry, boundary and dilution screens remain active. Histories must preserve exact initial chemistry, IDs and physical clocks. Interrupted jobs resume from their saved physical states with no duplicated observations. Completed results are checked and reused without rewriting their evidence. Pilot observations are saved separately and hashed, so appending a long history cannot change the evidence underlying the pilot gate.

## Endpoint basins and interpretation

Every completed moving job freezes its own endpoint graph for five chemical trials: the current state, two 1% log perturbations of that state and two 0.1% near-uniform perturbations. The original dual-solver, stationarity, Jacobian and local-return criteria are reused. The declared chemical horizon is 240, with extension to 960 only for unresolved endpoint classification, as in the parent assay.

We assess three distinct questions: whether small perturbations develop sustained contrast; whether an established pattern persists; and whether the final frozen graph supports locally stable chemical states. A near-uniform endpoint's assay does not supply a separately prepared patterned start, so a negative bistability flag does not exclude every patterned basin. Conversely, frozen local bistability is not stability of the complete moving system.

If weak polarity permits initiation while stronger polarity closes the formation window, the result would support directional mechanics limiting formation opportunity at fixed chemical rates. If all branches remain weak despite long frozen growth, other shared moving effects or inherited starting geometry remain candidates. If all moving branches initiate with longer observation, the earlier 60-unit result was horizon-dependent. These alternatives are declared before new moving results.

This study retains the existing geometric conductance approximation and does not resolve its general curved/gapped/nonorthogonal closure failures. It also does not establish spatial convergence, later division inheritance, autonomous identity, or a new persistent global shape axis.

## Execution and outputs

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.polarity_robustness prepare

# Run on the accepted GTX 1080 Ti; the coordinator prevents duplicate runners.
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.polarity_robustness run

OPENBLAS_NUM_THREADS=1 python -m embryo.polarity_robustness assess
```

`assess` verifies completed evidence; it does not require a GPU. Run it between coordinator writes to avoid racing the same summary file. Root `status.json` reports the active stage; each job's status reports current elapsed time and numerical audit. `pilot-refinement.json`, immutable `pilot-history.json` files, `long-refinement.json`, frozen paths, endpoint assays and the final summary preserve separate evidence. Source and input hashes are fixed before launch.

The projected total runtime is approximately **5–6 hours** on the GTX 1080 Ti, based on the parent's measured 3.31 seconds per coarse model-time unit, reused prefixes, doubled fine step counts and validation overhead. The first timestep gate should take approximately 30 minutes. These are estimates, not completion promises. The consolidated results ledger and manuscript remain earlier snapshots until this new evidence is assessed and accepted.
