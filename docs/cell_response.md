# Cell response, recovery, and chemical-state exchange

The methodology progresses from frozen endpoint responses to chemical exchange and matched moving responses. All cells share the same regulatory equations, and no identity classifier is supplied. The frozen assays, seed-7 moving exchange/response study, and targeted response timestep refinement are complete. The [current GPU history replication](cell_response_moving.md#replication-across-developmental-histories) adds seeds 8 and 9; complete response results remain pending. The frozen protocol below is retained with its original horizon and acceptance rules.

## Question and scope

Do developed chemical differences predict a cell's response to the same fractional stimulus, beyond its geometric context? This is a test of continuous response heterogeneity, not a cell-type classifier. Every cell has the same reaction equations. Network-level coexistence of uniform and patterned equilibria does not imply two identities available to each isolated cell.

## Fixed protocol

The experiment uses both time-150 endpoint geometries from developmental seeds 7, 8, and 9. Source/checkpoint hashes and successful bistability checks are verified before running. Each graph is initialized at either its already validated patterned equilibrium or its uniform equilibrium. These are the endpoints of the earlier 240-unit frozen chemical continuations, not the raw chemistry at moving time 150.

Each of 16 cells is separately subjected to an instantaneous activator multiplication by 0.5, 0.75, 0.9, 1.1, 1.25, or 1.5. Inhibitor and all other cells are unchanged. The bolus is an external addition/removal, not conserved redistribution; injected amount is recorded. Equal fractional pulses have unequal absolute amounts when starting concentrations or volumes differ. Each trial starts fresh from its equilibrium and runs for 240 chemical time units. A matched unperturbed continuation is shared by all interventions on that equilibrium.

This gives 6 graphs × 2 chemical backgrounds × 16 target cells × 6 pulses = **1,152 pulse trajectories**, plus 12 controls. These are nested interventions on three developmental histories, not independent biological or developmental replicates. Mechanics, geometry, volumes, and polarity do not evolve in this assay.

## Measurements

Changes are measured relative to the matched unperturbed trajectory using logarithmic concentration ratios. Recorded endpoints include target activator peak gain, induced inhibitor peak, time-integrated absolute target activator response, and the maximum volume-weighted chemical response across all other cells. Gains and integrals are divided by the absolute log pulse magnitude, so they are not simply the injected concentration difference. This normalization does not remove all nonlinear dose dependence.

Target recovery is the first sampled time after which both target-species log deviations remain below 10% of the initial activator log displacement. Network recovery uses 10% of the initial volume-weighted two-species log-RMS displacement. At least 24 additional observed time units must follow a claimed recovery. Sampling is every 0.1 units through time 10, every 0.5 through 60, and every 2 through 240. Recovery times are therefore sampling-limited. An unrecovered trajectory is censored, not evidence of permanent memory.

An endpoint is classified as returned only if its log-RMS distance from the control is below 1e-4, its maximum reaction-transport derivative is below 1e-6, and its full chemical Jacobian has negative largest real eigenvalue. A different stable endpoint must satisfy the last two checks while failing the return-distance criterion. Otherwise the outcome is unsettled or numerically unresolved, not silently counted as an identity switch.

## Numerical gates

Every trajectory and control is integrated with DOP853 at relative/absolute tolerances 1e-8/1e-10 and 1e-11/1e-13; maximum log disagreement must be below 1e-5. The most tolerance-sensitive pulse for each equilibrium/graph is independently checked with Radau at 1e-11/1e-13. Positivity and finiteness are required. Results, trajectories, controls, source hashes, and status are saved. A failed numerical gate prevents pooled assessment.

## Continuous exploratory analysis

The descriptive comparison pairs patterned and uniform responses at the same cell position and on the same operator. It changes the entire chemical background, not only the target cell: differences establish a collective-state contribution, not an intrinsic identity of that cell.

For the patterned background, average the normalized activator response integral from −10% and +10% pulses at each cell position. Compare two fixed ridge regressions (penalty 1): geometric context alone (log volume, exposure, log weighted conductance degree, distance from aggregate center) versus the same context plus baseline log activator and log inhibitor. Predict log response with one whole developmental seed held out at a time. Training-only centering/scaling avoids information leakage. No hyperparameter search, p values from pooled cells, or forced clustering is used. Geometry predictors do not encode the complete graph, so improvement cannot prove independence from all context. Three held-out histories provide an exploratory generalization check only.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_response \
  --seed 7 --output outputs/cell-response/seed-7
# Repeat with seeds 8 and 9; choose fresh directories for a repeat study.
python -m embryo.cell_response_summary --output outputs/cell-response
```

This experiment changes no production mechanics or chemical equations. Its response signatures remain dynamical observables, not independently specified cellular functions. Moving-geometry response, state-context exchange, loss of neighbor support, and inheritance remain subsequent questions.

## Completed results

All 1,152 pulse trajectories have completed. Numerical agreement is accepted after the targeted refinement described below. The response classifications are:

| Initial chemical background | Pulse | Returned to original equilibrium | Different stable endpoint | Unsettled at time 240 |
|---|---|---:|---:|---:|
| Patterned | −50% | 84/96 | 12/96 | 0/96 |
| Uniform | −50% | 20/96 | 64/96 | 12/96 |
| Patterned | Each of −25%, −10%, +10%, +25%, +50% | 96/96 per pulse | 0 | 0 |
| Uniform | Each of −25%, −10%, +10%, +25%, +50% | 96/96 per pulse | 0 | 0 |

Totals are 1,064 returns, 76 different stable endpoints, and 12 unsettled trajectories. All observed transitions follow the strongest negative pulse; no transition occurs under the tested positive pulses. This reveals finite-amplitude response asymmetry, not a universal threshold or a measured hysteresis loop. The 12 unsettled cases all belong to the seed-9 feedback-off geometry with initially uniform chemistry. They remain reported rather than being relabeled or excluded.

For ±10% perturbations, median target recovery times are 3.3/3.4 time units on patterned backgrounds versus 8.35/7.7 on uniform backgrounds (negative/positive pulse respectively). Averaging the normalized response integrals from these two pulses at each matched position, the median patterned-to-uniform ratio is 0.479 for seed 7, 0.528 for seed 8, and 0.447 for seed 9. Thus the patterned collective state generally buffers these small local perturbations more strongly, even on the same geometric operator. These numbers are descriptive summaries across paired interventions, not independent-cell confidence estimates.

On patterned backgrounds, adding baseline chemistry to the selected geometry predictors improves pooled out-of-seed R² from −0.204 to 0.224. However, held-out seed-9 R² remains negative (−0.485 with chemistry). This is limited predictive evidence, not a robust context-independent identity classifier. Chemistry, the full network, and the unequal absolute amounts delivered by equal fractional pulses remain intertwined. No clustering or number of cell types has been imposed.

The strongest conclusion is **state-dependent collective responsiveness and local recovery**, with large negative perturbations sometimes moving the network to a different stable chemical endpoint. The result does not establish cell-autonomous memory, inherited identity, or an independent functional phenotype. The uniform-versus-patterned intervention changes every cell's background chemistry; it does not isolate the target cell's state from that of its neighbors.

## Numerical refinement and audit trail

The original runs are preserved in `outputs/cell-response/`. Three background-level gates initially failed: both patterned controls for seed 7 had coarse-to-tight discrepancies above 1e-5, and two −50% uniform-background pulses on the seed-9 feedback-off graph exceeded that tolerance. The remaining stored tight trajectories already passed their comparisons.

`embryo.cell_response_refinement` audits those four affected integrations against both DOP853 and Radau at relative/absolute tolerances 1e-13/1e-15. It retains the original 1e-11/1e-13 trajectories and response measurements rather than replacing simulation results. Maximum discrepancy in these four audits is below 1.46e-9; no threshold, pulse, or horizon is changed. The original failed comparisons remain recorded, and revised classifications are stored separately in `outputs/cell-response-refined/`.

Across the accepted comparisons, maximum DOP853 log discrepancy is 8.96e-6 and maximum independent-reference discrepancy is 1.69e-8, both below 1e-5. All accepted paths are positive and finite. Six focused tests pass, including external-dose accounting, recovery censoring after late rebound, response-integral verification against an analytic decay, and the existing chemical Jacobian/recovery checks.

```bash
python -m embryo.cell_response_refinement \
  --source outputs/cell-response --output outputs/cell-response-refined
python -m embryo.cell_response_summary --output outputs/cell-response-refined
```

The original and refined output directories must be fresh. The refinement's fixed kinetics (β=2, D_a=0.02, D_b=0.4) were checked against all six source configurations before this audit; the utility is specific to this cohort. `assessment.json` retains all descriptive and held-out-seed summaries. `response.png` shows continuous response variation, paired chemical-background effects, and recovery times. Recovery curves exclude censored cases, whose counts are reported above.

An additional endpoint audit of the saved trajectories finds that all 76 alternative stable endpoints remain nonuniform: the 64 transitions from uniform backgrounds have final log-activator SD 0.588–0.697, and the 12 transitions from patterned backgrounds have SD 1.215–1.322. Thus the strong negative pulse can initiate finite-amplitude patterning despite uniform-state linear stability, or change an existing chemical pattern. This assay does not locate the basin boundary or establish that all resulting patterns are distinct attractors.

## Follow-up: chemical-state exchange and response transfer

The new assay in `embryo/cell_response_exchange.py` separates chemical-state relocation from fixed geometric location, then asks whether response properties are closer to the donor or destination. Earlier exchanges on one older no-feedback graph produced both restoration and collective reorganization; those results motivate this assay but are not reused as observations on the newer endpoint graphs.

The protocol uses seeds 7, 8, and 9 and both validated time-150 endpoint graphs (`switch_on`, `keep_off`). These are six geometries from three histories, not six independent histories. Each starts from its validated patterned chemical equilibrium. The minimum- and maximum-activator cells are selected before inspecting exchange outcomes, with ties resolved by cell ID.

Each geometry has four arms:

- Untouched tissue.
- Sham exchange: copy the selected states back into their original positions. This is an implementation control and must agree with the untouched trajectory.
- Exact exchange of both activator and inhibitor concentrations. Unequal volumes can change total amounts; those changes are recorded.
- Conservative exchange using the previous species-specific pair rescaling, preserving each species' total amount while leaving every other cell unchanged. Unequal volumes mean this is not an exact concentration transplant.

Geometry, transport, polarity, and cell IDs stay fixed. Reaction dynamics do not conserve chemical amounts; conservation refers only to the instantaneous intervention. Each arm relaxes for 1200 model-time units. The maximum chemical derivative must be below 1e-6 and the largest real eigenvalue of the complete chemical Jacobian must be negative before interpreting an endpoint as a locally stable chemical state. Unsettled endpoints remain unresolved.

Both selected cells then receive independent 0.9 and 1.1 activator multipliers, each followed for 240 units against an unperturbed continuation of that same relaxed arm. This gives 24 relaxation arms and 96 pulse trials. Pulses are external fractional boluses, with absolute amounts recorded; they are not conserved exchanges. The response sampling and recovery criterion match the existing frozen pulse assay.

For each target, the signed local response is `log(pulsed concentration / matched control concentration) / abs(log(pulse factor))`, retaining activator and inhibitor separately. Its distance to an untouched donor or destination response is the square root of the time-averaged, two-species mean squared difference over 240 units. Comparisons always use the same pulse factor. Both distances and their difference are saved. The descriptive nearest reference is reported only for settled endpoints with donor/destination baseline response separation greater than 0.01 in these normalized units. This exploratory cutoff avoids assigning a nearest reference when the original responses are nearly indistinguishable; it does not define a cell type or establish equivalence to a donor.

State persistence is assessed separately: late pairwise distances over relaxation times 1000–1200 are normalized by initial pair separation, using the earlier 10% return criterion and informative-separation threshold 0.1. Whole-network endpoint distance is also retained. Transferred ordering or response in a reorganized network cannot establish autonomous memory in an individual cell.

Every relaxation, matched control, and pulse trajectory is independently integrated with DOP853 and Radau at relative tolerance 1e-12 and absolute tolerance 1e-14. Maximum inter-solver log discrepancy must be below 1e-5, with finite positive concentrations. Failures remain numerical failures rather than biological outcomes. Eight focused tests passed before preparation. Inputs, sources, saved trajectories, and protocols are hashed; completed arms can be resumed without overwriting the mechanics experiments.

```bash
python -m embryo.cell_response_exchange prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_response_exchange run
```

Output: `outputs/cell-response-exchange/`. One additional process handles these small graph ODEs while the moving-geometry refinement continues. The experiment tests environment-supported chemical states and response behavior, not lineage inheritance, biological function, moving-geometry persistence, or autonomous cell identities.

### Completed exchange-response result

All 24 arms and 96 pulse trials completed. Every relaxation reached a locally stable chemical endpoint: maximum final derivative 4.87e-12 and least negative largest Jacobian eigenvalue -0.07354. Maximum DOP853/Radau log discrepancy was 5.75e-9, below 1e-5. Sham trajectories and responses matched untouched controls exactly. All 96 target pulse responses met the sustained-recovery criterion.

All twelve exchanges (six graphs, two exchange definitions) were classified as collective reorganization, not restoration of the original pair or exact retention of transferred concentrations. The selected cells' activator ordering reversed on all six graphs under exact exchange. Exact and conservative exchanges reached the same endpoints to maximum log discrepancy 4.29e-12.

All 48 exchanged-cell pulse comparisons were closer to the untouched donor response than the destination response. This count is nested within three developmental histories and includes two intervention definitions, two target cells, and two pulse signs; it is not 48 independent biological replicates.

Response transfer was asymmetric. Dividing donor-response distance by the original donor/destination response separation, the cell receiving the low-activator state had discrepancies approximately 0.037–0.119. The cell receiving the high-activator state had discrepancies 0.168–2.947. Thus even a response nearer the donor can be far from both original responses. The strongest supported interpretation is that chemical history redirects collective organization and influences subsequent response, with substantial context dependence. Autonomous identity, unchanged transfer of a full cell phenotype, and mechanics-independent biological memory remain unestablished.

Machine-readable summaries are in `outputs/cell-response-exchange/assessment.json` and `results.json`; each graph retains full relaxation and pulse trajectories, matched controls, solver checks, state outcomes, and signed response distances. The concurrent moving-geometry timestep-refinement batch was left running unchanged.

## Completed experiment: exchange with moving geometry

`embryo/cell_exchange_moving.py` prepares a two-arm seed-7 pilot on the full-coupling time-150 geometry. The first arm begins with a fresh conservative exchange between the same prespecified cells 27 and 20. It asks whether exchange-induced organization can form while geometry evolves. The second begins with the stable reorganized chemistry from the frozen conservative-exchange endpoint, transplanted onto the same starting geometry and polarity. It asks whether that chemical organization survives release into moving mechanics. These are distinct interventions, not two replicates. The transplant deliberately does not carry a donor cell's geometry or polarity.

Both arms run for 60 model-time units at timestep 0.00375 with observation interval 0.15 and no further division. The completed unexchanged patterned control from the moving-response refinement is reused: same source geometry, source clock, polarity, IDs, random streams and preconditioned original chemistry. This avoids another expensive identical control. The fresh exchange conserves amounts at the intervention; the pre-relaxed transplant need not have the same total amount because reaction dynamics acted during frozen relaxation.

The original launch was gated on completion and acceptance of the moving-response timestep refinement. A waiting process ran no mechanics until that evidence passed, then a single worker advanced the two continuations. Inputs and dependency hashes were checked again at launch. This records the original scheduling; both moving exchange arms have since completed.

Matched 60-unit frozen references for the original, freshly exchanged, and pre-relaxed chemistry are checked with DOP853 and Radau. Moving comparisons retain continuous pair and whole-network log distances, activator ordering, chemical contrast, aggregate axis ratio, and transport-operator differences. Late pair distances over elapsed times 36–60 are normalized by the initial fresh-exchange separation, comparing the unexchanged moving control and the conservative exchange of that control. A ratio below 0.1 describes closeness to a reference, not a moving-system attractor. The pre-relaxed reference and frozen continuation are reported separately. Reversed ordering alone does not establish unchanged transferred concentrations.

Every-step numerical quality and sampled boundary checks are inherited from the moving-response runner. Eight focused tests passed before preparation, including prerequisite gating, correct distinction between destination and exchanged states, conservative transplantation, and continuation/restart behavior. The earlier small-pulse refinement does not establish numerical convergence of this larger exchange perturbation: exchange-specific timestep refinement remains necessary before strong quantitative claims.

```bash
python -m embryo.cell_exchange_moving prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_exchange_moving run --wait
```

Both arms completed with passing numerical quality in `outputs/cell-exchange-moving/`. They retain reversed activator ordering and strong final log-activator dispersion (1.49526 fresh; 1.49533 pre-relaxed), but neither is within the 0.1 normalized late-distance criterion of the original destination or transferred concentration reference. Fresh exchange gives destination/donor ratios 1.27261/0.31426; pre-relaxed exchange gives 1.27272/0.31438. This is collective reorganization, not unchanged concentration transfer. The completed next stage tests pulse-response transfer from these moving endpoints. Exchange-formation timestep refinement remains outstanding.

## Moving-endpoint response transfer

The follow-up `embryo/cell_exchange_response_moving.py` tests whether response properties remain donor-like after chemical exchange and subsequent moving development. It uses the same-age time-210 endpoints of the unexchanged, fresh-exchange, and pre-relaxed-exchange branches. Every continuation preserves its native geometry, polarity, chemistry, cell IDs and random streams; there is no further transplantation or equilibration. Geometry and chemistry differ between backgrounds as outcomes of their histories.

Each endpoint has five continuations: an unperturbed moving control plus independent 0.9/1.1 activator pulses in cells 27 and 20. All fifteen runs continue mechanics and chemistry for 60 units at timestep 0.00375, sampled every 0.15 units. The time-150 pulse trajectories cannot substitute for these new time-210 donor/destination references. Pulse amounts are recorded because equal fractional interventions do not imply equal molecular amounts.

The comparison retains signed, two-species local log responses relative to each background's own evolving control, normalized by absolute log pulse. For an exchanged cell, the destination reference is the same ID in the unexchanged background; the donor reference is the other original exchange ID. References use the same pulse sign and the same elapsed observation times. Distances are time-averaged two-species RMS differences over 60 units, so they are not numerically interchangeable with the previous 240-unit frozen response distances. Both distances, their difference, and baseline reference separation are retained. A nearest-reference label is descriptive and omitted when baseline separation is at most 0.01. Recovery and whole-network effects use the established matched-control metrics and 24-unit follow-up rule.

This test concerns response transfer in an interacting moving tissue, not autonomous identity. The three backgrounds are interventions in one developmental history. Initial chemistry may still be drifting, which is why each requires its own moving control. The completed earlier small-pulse timestep check at time 150 does not validate these new endpoint interventions; pulse-specific refinement and more histories remain necessary for quantitative robustness.

Nine focused tests passed before preparation, including correct donor mapping, subtraction of background drift, refusal to use partial or invalid source results, a complete synthetic transferred-response assessment, source-state preservation under pulses, and exact restart continuation. The existing moving runner supplies per-step volume, radius, clipping, and chemical positivity checks and sampled boundary checks. Protocols, sources, starting checkpoints, histories, and source results are hashed. Source endpoints are copied only after completion and validation, and must have matching physical time, timestep, cell IDs, chemistry, and completed cleavage.

```bash
python -m embryo.cell_exchange_response_moving prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_exchange_response_moving run
```

All fifteen continuations completed with passing quality in `outputs/cell-exchange-response-moving/`. All eight exchanged-cell comparisons are nearer the donor response, with closer resemblance for the low-state recipient; the high-state recipient retains appreciable differences. Six targeted control/negative-pulse repetitions at dt=0.001875 subsequently passed response, integral, recovery, raw-chemistry, and classification agreement criteria. See [the complete moving results and refinement](cell_response_moving.md#completed-moving-exchange-response-study). Response transfer and unchanged transferred concentrations are separate outcomes.


### Optimized execution transition (2026-10-02)

The full-horizon optimized-mechanics control/pulse gate passed with exact agreement in sampled trajectories, endpoint fields, and recovery metrics. The moving-endpoint study resumed through `embryo.fast_response_resume` and subsequently the accepted native C++/OpenMP path before completion. Original protocols, source states, interventions, timesteps, and observations were retained. Backend-transition manifests and per-job execution records preserve the validation/checkpoint provenance. This is an archived execution transition, not the current job status.

## Current response definition and moving-history replication

For chemical species $s$ (activator or inhibitor), pulse multiplier $\zeta$, and target cell $i$, the signed response is

$$
R_i^s(t)=\frac{\log[c_i^{s,\mathrm{pulse}}(t)/c_i^{s,\mathrm{control}}(t)]}
{|\log\zeta|},\qquad \zeta\in\{0.9,1.1\}.
$$

Each background has its own unperturbed control. This removes its ongoing chemical drift rather than comparing a moving cell to a presumed equilibrium. The absolute denominator normalizes pulse size; the numerator retains the response sign. A fractional activator pulse is an external amount change, not conserved redistribution. The dose depends on the target's starting concentration and volume.

The distance between a recipient response and an untouched donor/destination reference is

$$
d=\sqrt{\frac{1}{2T}\int_0^T
\sum_{s\in\{a,b\}}(R_i^s(t)-R_{\mathrm{reference}}^s(t))^2\,\mathrm dt}.
$$

Both species and all matched observation times contribute equally; saved samples use trapezoidal integration. The moving horizon is $T=60$, versus $T=240$ in the frozen response assay, so their numerical distances are not interchangeable. Both reference distances and their separation are reported. Separation must exceed 0.01 to resolve a descriptive nearest reference; this does not establish equivalence or define a cell type.

Chemical-state persistence is measured separately with initial-volume-weighted log distances and pair ordering. Actual evolving volumes still set transport capacities and mechanical dilution. Holding diagnostic weights fixed prevents changing weights from masquerading as a changed chemical response. Recovery requires deviation below 10% of initial log displacement throughout the remaining samples, with at least 24 units of subsequent observation; later/unobserved recovery is censored.

The current [moving-history replication](cell_response_moving.md#replication-across-developmental-histories) repeats all three backgrounds on seeds 8 and 9, with fresh/pre-relaxed formation from t=150 to 210 followed by matched pulses through t=270. It comprises 36 mature continuations; seed 7 is the completed historical reference. PyTorch owns resident GPU arrays/matrix operations, and custom CUDA handles mechanics and spatial geometry/polarity. All four full-horizon backend replays and all six new-history starting-context checks passed before scientific execution. The running study's complete response outcomes remain pending.

This protocol is conditional on existing mature patterned basins and explicitly prepared chemistry. It is not fresh zygote-to-identity formation. Histories are the independent developmental units; backgrounds, targets, and pulse signs are nested interventions. Protocol/source/evidence hashes, strict observation alignment, standard checkpoints, and per-step quality screens remain enforced. General-geometry transport closure, new-history time/spatial refinement, cell autonomy, inheritance, and biological function remain separate questions.
