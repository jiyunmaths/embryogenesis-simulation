# Recipient-state matching on shared physical contexts

## Question

The [accepted shared-context assay](phase_carry_shared_context_response_assessment.md) detects response differences in cells 19/29 when transport and volumes are shared. Its chemical starts still differ at the recipient and throughout the network. **Do response differences persist when the recipient's initial activator and inhibitor are identical, while the surrounding chemical states differ?**

Use the same actual t=510 inputs, histories 7/8/9, recipients 22/23, 19/20 and 29/30, and sixteen-cell compartments. For each recipient, test both measured physical contexts: sham and its corresponding redistributed context. A physical context is the exact saved conservative operator together with its matching volume vector and cell-ID order.

Repeat the comparison at two recipient reference states: its measured sham `(a_i,b_i)` and its measured redistributed `(a_i,b_i)`. For each reference, replace both recipient concentrations in each complete chemical start with that same measured pair. Preserve every other cell's two concentrations exactly. This gives two surrounding-network starts, two physical contexts and two matching references per recipient. References are experimental matching choices, not prescribed cell identities.

Matching and challenge occur immediately at t=510. The recipient subsequently evolves freely: there is no clamp, equilibration, reservoir or continuing reset. Its own matched unperturbed continuation supplies the control. Changes in the target state after t=510 are outcomes of the coupled chemical network, not failures of initial matching.

## Amounts, pulse and equations

Matching preserves concentrations outside the recipient, not global chemical amounts. Record the donor amounts, amounts after physical-context reassignment, and amounts after recipient matching, per cell and for each species. Separate volume-reassignment changes from the activator/inhibitor added or removed by target matching. Initial species totals and the chemical distribution of the other fifteen cells remain unmatched; this is a whole-surrounding-network intervention, not a nearest-neighbor, allocation-only or dosage-controlled assay.

Apply one activator-only **10% depletion** to the matched recipient. Inhibitor and all other cells remain unchanged at the pulse. Within each shared physical context, both chemical arms have the same target concentrations and volume, so their removed activator amounts are also identical. The previous fractional/common-amount designs coincide here; do not count a second dose or duplicate integrations as independent evidence. Verify requested versus measured removal within relative error `2e-13`, and retain the actual factor `0.9` and log amplitude `abs(log(0.9))`. Physical-context volume changes amount, and matching reference changes concentration; neither changes the fractional factor.

Continue the unchanged equations to t=570, with `beta=2`, `D_a=0.02` and `D_b=0.55`:

$$
\frac{d a_i}{dt}=\frac{a_i^2}{b_i}-a_i+D_a\sum_j\Delta_{ij}a_j,
$$

$$
\frac{d b_i}{dt}=\beta(a_i^2-b_i)+D_b\sum_j\Delta_{ij}b_j.
$$

Here activator self-enhancement is inhibited by `b_i`; activator generates inhibitor, and both species decay. `Delta=-M^{-1}K`, with `M=diag(V_i)`, exchanges concentrations through reciprocal measured contacts. Fixed volumes imply zero geometric dilution. Exchange conserves each species' volume-weighted amount; reactions and the experimental resets do not. Advance no mechanics, polarity, geometry or divisions. Check conservative operator structure without modifying any saved matrix.

There are 48 factorial arm labels, but **45 distinct chemical–physical backgrounds**: the unchanged sham/sham control is shared by the two recipients within each history. Each target-specific pulse remains separate. This gives **93 logical starts and 93 distinct initial-value problems**, verified from prepared arrays: 45 controls and 48 pulses. All starts are computed anew. Two primary timesteps give 186 production paths; a third adds 93 diagnostics, for 279 GPU paths and 279 separate CPU verifications. DOP853 and Radau each solve all 93 starts, yielding 186 reference trajectories. Exactly unchanged arms provide numerical comparisons with the accepted parent; no parent path is reused as a new scientific outcome.

## Primary response and descriptive outcomes

Sample every 0.15 units throughout the full 60-unit interval. For each arm let `R_C,P,L(t)` be the signed target two-species pulse/own-control log ratio divided by the actual absolute initial log pulse amplitude. `C` denotes the origin of the surrounding chemistry (sham or redistributed), `P` the shared physical context, and `L` the common recipient reference (sham or redistributed).

For each fixed `P,L`, compare the two surrounding starts:

$$
E_{P,L}=\left[\frac{1}{60}\int_0^{60}\frac12\sum_{s\in\{a,b\}}
\left(R_{R,P,L,s}(t)-R_{S,P,L,s}(t)\right)^2dt\right]^{1/2}.
$$

Use trapezoidal observation-grid integration. The primary operational screen remains strictly **`E>0.01`**. Report all 24 comparisons and history counts out of two recipients separately for each matching reference and each physical context. Also report which recipients pass all four comparisons; do not select a favorable reference, collapse reference-dependent outcomes, or infer equivalence below the screen.

A passing difference shows that distinct surrounding chemical states can affect a recipient's response despite identical initial local chemistry and shared physics, under the tested intervention and horizon. It does not establish autonomous identity, immutable memory, direct-neighbor specificity, or biological fate. Loss of a previous difference after matching is not an additive causal fraction: target matching changes the subsequent coupled trajectory.

Compare qualified fractional-pulse effects continuously with the parent under each matching reference/physical context. Describe physical-context contrasts at fixed surrounding start/reference and their signed difference-of-differences. Form the signed waveform before taking its RMS; assign no additional scientific screen to these descriptive contrasts and no mediation percentages.

Report target activator response area, target and volume-weighted network recovery, and full signed waveforms. Recovery retains the inherited definition: all remaining sampled displacement lies within 10% of its initial value, with at least 24 later observed units; otherwise right-censor at t=570. Recovery is toward the arm's own evolving control and is not a scientific positive-result gate.

For unperturbed controls, describe drift from the post-matching start, reappearance of target concentration differences, final activator/inhibitor ratios, late (elapsed 36–60) target and volume-weighted network log separation, and across-cell log-activator contrast. Initial target ratios must be exactly one on every shared context/reference. Late differences and recovery are descriptive; require neither as a positive result. These are immediate post-reset responses, not assumed equilibrium traits.

## Numerical qualification and evidence

Retain float64 PyTorch GPU arrays and matrix operations with unchanged `gpu_backend.gm_step`, checked SSP-RK2 stages, and unchanged NumPy/SciPy verification. Independent backgrounds occupy block-diagonal transport blocks; padding duplicates positive controls, is excluded from reported paths, and creates no exchange between backgrounds. Require every constituent positivity substep count to equal the batch count at each timestep. Check both internal Euler stages and final states; never clip concentrations.

Use `dt=0.00375`, `0.001875` and diagnostic `0.0009375`, with identical prepared states and observations. Use DOP853/Radau at `rtol=1e-13`, `atol=1e-15`, maximum step `0.15`. Require raw and normalized reference disagreement at most `1e-8`; production/reference raw log error at most `0.01`; CPU/PyTorch raw error at most `1e-8`. Resolve order where finest error exceeds 100 times reference disagreement, requiring `1.7–2.3`; report unresolved order otherwise.

Preserve normalized pulse waveform tolerance `0.01`, relative activator-area tolerance `0.02`, matching recovery/censoring and recovery-time difference at most `0.3`. Directly compare signed surrounding-chemical, physical and interaction contrast waveforms and their RMS norms: production limits `0.01`, independent-reference limits `1e-8`. Require primary surrounding-chemical decisions to agree across all timesteps and references. Retain all failures without changing thresholds or dropping targets.

Before production, qualify every actual start through a short CUDA/CPU prefix. After dual-reference qualification, CPU verification may run concurrently with GPU production using four single-thread CPU workers; outputs occupy separate paths. Assess only once both have completed. Reference and scientific kernels remain unchanged.

Pin this protocol, the new runner/tests, inherited scientific dependencies, accepted parent assessment/review/provenance, and every consumed prepared array. Preserve reset and pulse accounting, all reference/CPU/GPU paths, checkpoints and commit hashes, clocks, runtime, timing, status and numerical summaries. Separately reconstruct the saved arrays and reported measures without importing project analysis helpers before accepting a scientific interpretation. Runtime completion alone is not acceptance.

## Execution and scope

Use `python -m embryo.phase_carry_recipient_state_response prepare`, then `run`, with `--output output/studies/phase-carry-recipient-state-response`. Use physical GPU 1 as logical `cuda:0` through `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1`, four CPU workers, and one numerical-library thread per worker. Small chemical outputs stay in the writable repository; read parents in place and copy no 3D fields.

The independent units remain three existing mature histories. Local references, targets, arms, timesteps and solvers add no new developmental histories. Retain historical no-carry developmental inheritance, fixed sixteen-cell compartments, frozen mechanics, unmatched network amounts and prior spatial/conductance limitations. This step updates no ledger or manuscript numerical Results, changes no shared kernel/baseline simulator, and launches no subsequent experiment.
