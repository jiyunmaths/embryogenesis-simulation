# Initial surrounding-amount control with recipient matching

## Question

The [accepted immediate matching](phase_carry_recipient_state_response_assessment.md) detects surrounding-network effects despite equal initial recipient chemistry. The [accepted delayed challenge](phase_carry_recipient_delayed_response_assessment.md) separates large immediate reset responses from selective later behavior. Both leave surrounding chemical amounts unmatched.

**Can different surrounding chemical allocations produce different states and responses when their initial species amounts and the recipient's initial state are equal?** Test the same three mature histories (7–9), recipients 22/23, 19/20 and 29/30, both measured recipient references and both shared conservative operator/volume contexts. Retain every comparison. These interventions add zero developmental histories.

## Symmetric initial amount matching

Begin with the actual accepted t=510 recipient-matched inputs. For each fixed recipient, matching reference and physical context, the sham (`S`) and redistributed (`R`) chemical vectors already have exactly equal target activator/inhibitor concentrations. Their physical operator, volumes and ID order are also equal. Let `u_s,j` be species `s` in compartment `j` and let `i` denote the recipient.

The surrounding amount excludes the target:

$$
B_{C,s}=\sum_{j\ne i}V_j u_{C,s,j},\qquad C\in\{S,R\},\quad s\in\{a,b\}.
$$

Use a common, symmetric positive budget for each species:

$$
B_s^*=\sqrt{B_{S,s}}\sqrt{B_{R,s}},\qquad
\alpha_{C,s}=\frac{B_s^*}{B_{C,s}}.
$$

Replace only the non-target concentrations once, at t=510:

$$
\widetilde u_{C,s,j}=\alpha_{C,s}u_{C,s,j}\ (j\ne i),\qquad
\widetilde u_{C,s,i}=u_{C,s,i}.
$$

The geometric mean treats the two origins symmetrically. It is a specified experimental budget, not a fitted optimum or a biological law. Each species' volume-weighted allocation fractions `V_j u_s,j / B_s` among non-target cells are preserved. Absolute surrounding concentrations and each non-target cell's activator/inhibitor ratio can change: the species use separate scale factors. Do not claim that complete local chemical states are preserved outside the target.

The target pair and its physical context remain unchanged. Both surrounding amounts and both whole-network amounts are initially equal between chemical origins within each reference/physical-context comparison. Require relative amount discrepancies at most `2e-13` and allocation-fraction preservation to the same relative tolerance. Require finite positive float64 concentrations; do not clip, compensate by changing a selected cell, or modify the operator/volumes. Record pre/post per-cell amounts, surrounding totals, whole-network totals, common budgets, scale factors, target pair and allocation fractions. Retain exact input hashes.

All recipient-specific preparations stay separate: there are **48 background labels**. Even formerly shared sham controls depend on which target is excluded from the budget. Count exact initial-value problems from the arrays; identical labels or dose preparations add no independent evidence. Equalization occurs only initially. The unchanged reactions subsequently change amounts freely; later amount differences are measured outcomes of different initial allocations, not a continuing clamp.

## Two prespecified phases

**Immediate phase:** apply one recipient-only 10% activator depletion at t=510 and observe through t=570. Inhibitor and every other compartment remain unchanged by the pulse. Local states and volumes match within each pair, so fractional and common-amount designs coincide; run only one challenge. Each pulse has its own freely evolving unpulsed control. This phase also supplies the unpulsed waiting trajectories for the delayed phase. There are 48 controls and 48 pulses: **96 logical starts**, 288 GPU paths at three timesteps, 288 CPU checks and 192 DOP853/Radau reference solves, all newly calculated.

**Delayed phase:** wait exactly 60 units and challenge at t=570, then observe through t=630. Use the numerically qualified fine unpulsed endpoints from the immediate phase as common starts for every method. Never use a pulsed source or rematch the target or surrounding amounts at t=570. Require every fine control endpoint to agree with both other production levels and both reference endpoints within maximum raw log discrepancy `1e-8` before preparing the delayed phase. A failed source blocks that phase without dropping the background or changing the gate.

Describe preceding 24-unit raw log drift and instantaneous absolute log-concentration rates for every background. Retain the inherited descriptive settling screen: drift at most `0.001` and rate at most `0.0001` per time unit. This is not a scientific positive-result or numerical acceptance gate. Do not exclude unsettled controls or call their responses settled. Describe drift during the delayed observation too; do not assume exact equilibrium.

At the delayed challenge, target concentrations and network totals can differ. Test both 10% fractional depletion and a common removed activator amount `0.1 * min(a_S,a_R) * V_i`, calculated separately for each original reference/shared-physics pair from the actual delayed states. Record actual factors and initial log amplitudes; require measured/requested removal error at most `2e-13`. There are 48 controls and 96 logical pulses: **144 logical starts**, 432 GPU paths, 432 CPU checks and 288 reference solves. At least one dose preparation per chemical pair coincides; count distinct problems explicitly and retain both labels.

The complete planned test therefore calculates **720 GPU paths, 720 CPU verification paths and 480 reference solves** across the two phases. Distinct problems, copied waiting sources and repeated dose labels are reported separately. The units of independent development remain **three existing histories**, with two selected recipients nested within each, zero new histories and no new moving development.

## Equations and scientific measures

Both phases use the unchanged Gierer–Meinhardt equations and fixed conservative transport:

$$
\frac{d a_i}{dt}=\frac{a_i^2}{b_i}-a_i+0.02\sum_j\Delta_{ij}a_j,
$$

$$
\frac{d b_i}{dt}=2(a_i^2-b_i)+0.55\sum_j\Delta_{ij}b_j,
\qquad \Delta=-M^{-1}K,\quad M=\operatorname{diag}(V_i).
$$

No mechanics, polarity, division, reservoir or geometric dilution advances. Transport conserves volume-weighted species amounts; reaction, initial rescaling and depletion do not.

Retain the target's signed two-species pulse/own-control log response divided by its actual absolute initial log pulse amplitude. The primary effect is the full 60-unit time RMS of the difference between chemical origins, with equal species weights and trapezoidal observation-grid integration. Apply the same strict operational screen **`E>0.01`** separately: 24 immediate comparisons and 48 delayed comparisons. Report history counts out of two recipients for every reference, shared physical context and dose; also report recipients passing every prespecified combination. Decisions must agree across all three timesteps and both reference solvers. Below-screen is not equivalence.

Compare the new immediate fractional effects continuously and by screen decisions with accepted unmatched-amount immediate effects. Compare delayed results with accepted unmatched-amount delayed results at the same dose/reference/physics. Changing initial amount matching changes subsequent coupled trajectories; these comparisons are not additive percentages of allocation versus dosage causation.

Physical contrasts and signed interactions remain descriptive, with no positive-result screen. Form signed waveforms before their RMS. Report target activator response area, target/volume-weighted network recovery, signed waveforms, target species ratios and late target/network log separation. Recovery retains the same sustained 10%-displacement definition and at least 24 later observations; otherwise censor at the observation horizon. No pulse recovery or positive effect is required for numerical acceptance.

For every unpulsed chemical pair, retain volume-weighted whole-network species-amount trajectories and ratios, their initial closure, maximum/late log separation and final values. Initial amount equality is the intervention; subsequent reaction-driven differences are part of its outcome. A retained effect supports sensitivity to **initial allocation at a common initial budget**, conditional on this rescaling, target reference and frozen physical context. It does not establish dose independence at every later time or isolated nearest-neighbor signaling.

## Numerical qualification and execution

Use unchanged float64 PyTorch GPU arrays/matrix operations and checked SSP-RK2 chemistry, grouping identical physical operator/volume/ID triples into independent trajectory columns. Padding duplicates positive controls and is excluded from reported counts. Require batch and constituent positivity substep counts to match. Retain checked internal Euler stages and finite positive outputs, with no clipping or CPU scientific fallback.

Keep `dt=0.00375`, `0.001875`, diagnostic `0.0009375`, observations every 0.15 and checkpoints every 3 units. Retain DOP853/Radau `rtol=1e-13`, `atol=1e-15`, maximum step 0.15. Gates remain: raw/normalized/signed-contrast reference disagreement `1e-8`, CPU/PyTorch raw discrepancy `1e-8`, production raw/normalized response/signed contrast/RMS discrepancy `0.01`, relative activator-area discrepancy `0.02`, matching recovery/censoring and recovery-time discrepancy `0.3`. Resolve order only above 100 times reference noise, requiring `1.7–2.3`; retain unresolved cases without assigning order two.

Before each full phase, run preparation/design/pipeline tests, a separate literal array audit and short CPU/CUDA checks on every actual start. Qualify dual references before full production. Numerically qualify and separately reconstruct the immediate phase before accepting it as a delayed source. Pin scientific sources, this protocol, runtime, accepted parent reports/provenance and all consumed/prepared arrays. Retain failures unchanged, all trajectories/checkpoints, amount accounting, status and timing. Separately reconstruct scientific/numerical measures from saved arrays without project analysis helpers before acceptance.

Run physical GPU 1 as logical `cuda:0`, with up to eight single-thread CPU workers for references/verification alongside GPU production (twelve CPU threads are available). Workers solve independent paths; parallel scheduling changes no equation or step. Keep both small chemical study outputs under writable `output/studies/phase-carry-recipient-amount-response` and `output/studies/phase-carry-recipient-amount-delayed-response`. Read parent data in place and copy no 3D fields.

## Scope

This is one prespecified two-phase frozen-geometry experiment. The immediate phase controls target chemistry, initial surrounding/whole amounts and physics. The delayed phase preserves the initial intervention but allows local and global chemistry to evolve. Retain historical no-carry developmental inheritance, sixteen-cell resolution, three selected mature histories, frozen mechanics and inherited conductance limitations. No autonomous/inherited biological identity, fresh carry-corrected zygote differentiation or population prevalence is established. Accepted parents, baseline/shared kernels, the results ledger and manuscript remain unchanged. Reporting may update after acceptance; no later experiment is launched.
