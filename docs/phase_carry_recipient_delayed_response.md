# Delayed challenge after recipient-state matching

## Question and fixed timing

The [accepted recipient-matching assay](phase_carry_recipient_state_response_assessment.md) equalized the target's activator and inhibitor at t=510, then pulsed it immediately. It detects surrounding-network response effects in 0/2, 2/2 and 2/2 selected recipients in histories 7–9. Local concentration differences reappear; immediate effects in cells 20/30 depend strongly on the initial reset reference despite agreement between their later unpulsed endpoints.

**Do response differences remain after a fixed 60-unit wait following matching?** Challenge at t=570, then observe another 60 units through t=630. This delay is fixed before measuring delayed response outcomes. Retain both original matching references, both measured shared physical contexts and all six selected recipients. Do not select a favorable delay, reference, recipient or physical context. The immediate t=510 assay is accepted comparison evidence, not new replication.

Use the accepted fine unpulsed t=570 whole-network states as the common prepared starts for every new timestep and reference solver. There are 45 backgrounds, including the three shared native sham controls. Copy chemistry only; retain each background's exact saved operator, volumes and cell-ID order. Never take a state from a pulsed parent path. Neither the recipient nor its neighbors are reset at t=570: matching occurred at t=510, and the whole chemical network evolved freely during the wait. Target chemistry and surrounding global amounts can differ at challenge time. Thus delayed responses test behavior on re-established states, not continuing equality of target chemistry or isolated surrounding chemistry at the pulse.

## Source and settling checks

Before launching the new continuations, require every fine control endpoint to agree with both accepted reference endpoints and both other accepted production endpoints within maximum raw concentration log difference `1e-8`. Pin all consumed parent paths and their result files. Every new method starts from the same actual fine endpoint; this qualifies the delayed continuation conditional on that accepted state, rather than a new full-development convergence claim.

For each unpulsed background, record the maximum absolute whole-network concentration log drift relative to its final state over the preceding 24 units (t=546–570), and the maximum absolute instantaneous log-concentration rate `abs(rhs/u)` at t=570. Use the prespecified descriptive settling screen: drift at most `0.001` and rate at most `0.0001` per model time unit. Record both numbers and flags for all backgrounds. This is not a numerical acceptance or positive scientific outcome gate. No unsettled arm is dropped; if a background fails, describe its later response as conditional on remaining baseline drift, not settled behavior. Reassess control drift during the new observation window descriptively. Do not claim mathematically exact equilibrium.

## Challenge and amounts

For each original matching reference and physical context, compare sham versus redistributed surrounding chemical origins at the delayed challenge. Retain each arm's own unpulsed evolving control. Test two activator-only depletion designs:

1. Fractional: remove 10% of that arm's target activator.
2. Common amount: remove `0.1 * min(a_S, a_R) * V_i` from each arm, calculated from the actual delayed pair on the same physical context and original matching reference.

Inhibitor and every non-target concentration stay unchanged at the pulse. Physical contexts stay paired with their original volume vectors. Verify measured versus requested removal within relative error `2e-13`. Preserve actual concentration factors, log amplitudes, target volumes, target amounts and global per-species amounts. Normalize response by each pulse's actual absolute initial log amplitude. Common-amount factors can differ between chemical origins; no dose is called exactly matched across different physical contexts or different matching references.

At least one fractional/common preparation in each pair is identical. Retain both labels and report the exact number of distinct initial-value problems; duplicated dose labels add no biological or numerical independence. There are **45 controls and 96 logical pulses: 141 logical starts**, with 282 primary GPU paths and 141 diagnostic GPU paths, **423 GPU paths and 423 CPU checks** in total. DOP853 and Radau each solve all logical starts, giving 282 reference trajectories. There remain **three existing mature developmental histories**, zero new histories, with two selected recipients nested within each. Reference choices, physical contexts, doses and numerical methods are within-history interventions.

## Equations, responses and decision rules

Advance the unchanged Gierer–Meinhardt reactions with conservative fixed-volume transport:

$$
\frac{d a_i}{dt}=\frac{a_i^2}{b_i}-a_i+0.02\sum_j\Delta_{ij}a_j,
$$

$$
\frac{d b_i}{dt}=2(a_i^2-b_i)+0.55\sum_j\Delta_{ij}b_j.
$$

Here `Delta=-M^{-1}K`, `M=diag(V_i)`. The saved conservative operator and volumes are not modified. No mechanics, polarity, division, changing geometry, reservoir or dilution advances. Reaction and experimental depletion can change chemical amounts; transport conserves volume-weighted amounts.

For each original matching reference, dose and fixed physical context, form the signed target two-species pulse/own-control log response divided by its actual initial log pulse amplitude. The primary effect is the difference between the two chemical origins, with equal species weights and trapezoidal observation-grid integration:

$$
E_{P,L,d}=\left[\frac{1}{60}\int_0^{60}\frac12\sum_{s\in\{a,b\}}
\left(R_{R,P,L,d,s}(t)-R_{S,P,L,d,s}(t)\right)^2dt\right]^{1/2}.
$$

Use the same strict screen **`E>0.01`**, separately for all 48 chemical comparisons. Report history counts out of two recipients for every reference/dose/physical context, plus recipients clearing all eight comparisons. Do not infer equivalence below the screen. Compare delayed fractional results with the accepted immediate fractional result continuously and by screen decisions; changing the delay also changes the whole chemical state, so differences are not additive causal percentages.

Physical-context contrasts and signed difference-of-differences are descriptive with no scientific screen. Form signed waveforms before their RMS; scalar effects do not add. Compare the two original matching references descriptively at fixed chemical origin/physics/dose to assess whether their delayed behavior agrees after their measured state convergence. No new equivalence threshold is introduced.

Retain target response area, signed waveforms and target/volume-weighted network recovery toward each pulse's own evolving control. Recovery means all remaining sampled displacement lies within 10% of initial displacement, with at least 24 later observed units; otherwise censor at the 60-unit observation horizon. Neither recovery nor a positive response contrast is required for numerical acceptance. Describe challenge-time/final target species ratios, late target and network log separation, and control drift. Unequal challenge-time target chemistry remains an explicit limitation of this delayed intervention.

## Numerical qualification and execution

Use float64 PyTorch GPU arrays and matrix operations, the unchanged checked SSP-RK2 chemistry kernel, and the accepted grouping of identical operator/volume/ID triples into independent trajectory columns. Padding duplicates positive controls and is excluded from reported counts. Grouping creates no chemical exchange between columns. Require every constituent positivity substep count to equal the batch count.

Use `dt=0.00375`, `0.001875`, and diagnostic `0.0009375`; sample every 0.15 units and checkpoint every 3 units. Retain tight DOP853/Radau settings (`rtol=1e-13`, `atol=1e-15`, maximum step 0.15), all-stage finite positivity checks and no clipping. Retain inherited raw/normalized reference tolerance `1e-8`, CPU/PyTorch raw tolerance `1e-8`, production raw/response/signed-contrast/RMS tolerance `0.01`, relative activator response-area tolerance `0.02`, and recovery-time tolerance `0.3`. Require matching recovery/censoring. Resolve order when finest error exceeds 100 times reference disagreement, requiring `1.7–2.3`; retain unresolved cases as such. Require primary decisions to agree across all methods/timesteps. Check direct signed physical/interaction contrasts too, without assigning them a positive-result screen.

Run new design tests, a separate literal prepared-input audit and short CPU/CUDA comparisons on every actual start before full production. Dual references must qualify before full GPU and CPU production. Use four single-thread CPU workers overlapping GPU production, logical `cuda:0` on physical GPU 1, and no CPU fallback. Retain source/runtime/parent/input hashes, amount records, references, all production/checkpoint arrays, status and timing. Separately reconstruct saved numerical and scientific measures without project analysis helpers before accepting interpretation. Failed gates remain recorded without changing thresholds or dropping arms.

Execute `python -m embryo.phase_carry_recipient_delayed_response prepare`, then `run`, with output `output/studies/phase-carry-recipient-delayed-response`. Read accepted parent arrays in place; write only small chemical study outputs in the repository, not relocated `outputs`. No 3D fields are copied.

## Scope

This test distinguishes immediate-reset responses from behavior after a fixed waiting interval on frozen physical contexts. It does not freshly equalize target chemistry at the delayed pulse, match surrounding totals, isolate direct neighbors, or establish autonomous/inherited biological identity. Retain historical no-carry developmental inheritance, sixteen-cell resolution and prior conductance/spatial limitations. Update this study's assessment and README/plan reporting after acceptance, without modifying accepted parent evidence, baseline/shared kernels, ledger, manuscript or launching a subsequent experiment.
