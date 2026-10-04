# Reporting initiation failures and identifying their moving-geometry cause

This reporting clarification and follow-up design is recorded while the three-history carry study is running, before its final three-history assessment. Some numerical paths and pilot checks have already completed; this is not a claim of a blinded preregistration. It clarifies the existing per-history decision rule without changing the active protocol, source code, thresholds or horizon. The [machine-readable addendum](phase_carry_polarity_reporting.json) records the rule and proposed controls separately from the pinned simulation protocol.

## Reporting rule for the current study

Each history is a replication unit. Coarse/fine timesteps, two mechanical branches and frozen references are nested measurements. Report all three history rows; do not pool cell-level measurements to manufacture replication or interpret a fraction of three histories as a population frequency.

Initiation means SD(log activator) exceeds **0.1 throughout elapsed 216–240**. Both timesteps must agree under the original full-window numerical and physical checks. A missing or failed comparison is **unresolved**, not a noninitiating history.

| Qualified moving outcome within a history | Report |
|---|---|
| chi=0 forms; chi=0.35 does not | Polarity-dependent suppression supported in this history. |
| Neither forms; the matched frozen reference forms | Moving-versus-frozen restriction persists with directional polarity tension disabled; its physical cause is unresolved. This does not establish polarity-specific suppression. |
| Neither forms; frozen-reference matching or accuracy is unverified | No successful moving control; attribution to moving geometry remains unverified. |
| Neither forms; a qualified frozen reference also fails | No formation in this window; no moving-versus-frozen initiation difference established. |
| Both form | No suppression of initiation by this intervention under the stated criterion. Differences in timing/amplitude are separate quantitative claims. |
| Only chi=0.35 forms | Positive-contrast-only initiation; report it even though it opposes the suppression hypothesis. |
| Numerical checks fail or a path is incomplete | Unresolved; exclude from supported mechanistic counts and give the reason. |

If the old outcome repeats, **all three carry comparisons qualify**, and the frozen references are verified against the current starts, use:

> Directional polarity mechanics suppressed initiation in **1 of 3 tested histories** (history 9). In histories 7 and 8, initiation failed even with directional polarity tension disabled, while matched frozen references formed patterns. Thus the moving-versus-frozen restriction extends beyond directional polarity tension; its physical cause in those two histories remains unresolved.

This is the intended meaning of a shared moving-geometry restriction. It is a common observed failure, **not yet evidence that the same physical mechanism causes both failures**. Polarity still evolves at chi=0, and the mature starts retain their previous mechanical history. Do not describe these histories as polarity-free embryos.

If carry changes the outcomes, update the history counts and use the corresponding rows above. Do not force the result into the expected 1/3 account. If any history is unresolved, give supported cases, numerically qualified histories and the three planned histories separately; do not present 1/3 as a fully qualified result.

The old frozen chemistry does not contain moving phase updates, but reference reuse still requires verified provenance and matching starting chemistry, IDs, volumes, transport, kinetics and 240-unit observation window. Check the actual carry run's initial operator. If the existing native-graph reference differs from the exported GPU operator, recompute the small frozen ODE reference on that actual operator with DOP853 and Radau instead of assuming exact matching from the source filename. Require maximum solver log disagreement at most **1e-8**. Neither an instantaneous unstable mode nor a failed endpoint bistability flag substitutes for this formation reference.

## First follow-up: dilution and contact remodeling

Target **histories 7 and 8**, conditional on each retaining a qualified chi=0 moving failure while its matched frozen reference forms. A history that forms with carry is reported as such and does not enter a rescue assay for a failure it no longer exhibits. Reuse the current chi=0 baseline paths; keep the exact mature t=150 initial state, near-uniform chemistry, cell order, random streams, material laws and numerical carry.

Cross two interventions in a complete 2 by 2 design:

| Arm | Mechanical dilution | Contact conductances used for chemistry |
|---|---|---|
| Baseline | On | Evolve with geometry |
| Dilution off | Off | Evolve with geometry |
| Fixed conductances | On | Fixed at the initial matrix |
| Combined | Off | Fixed at the initial matrix |

Mechanics, polarity, current cell volumes and chemical feedback on tension/adhesion continue in **all four arms**. Set chi=0 throughout, retain activity–tension coefficient 0.25 and activity–adhesion coefficient 0.35, and keep beta=2, D_a=0.02 and D_b=0.55. This separates changes after the restart from the earlier conditioning of its geometry.

The primary hypothesis is that evolving contact exchange closes the formation opportunity in the two weakly unstable initial graphs. The competing hypothesis is that volume-driven concentration dilution, alone or together with contact change, prevents formation. These are hypotheses, not inferred mechanisms from the old percentage declines or spectral crossing times.

## Conservative definition of fixed transport

Let G(t) contain symmetric off-diagonal conductances, K(t)=diag(G(t)1)-G(t), and M(t)=diag(V_i(t)). Use the **actual chemical conductances** exported by the carry backend at the restart for G(0).

$$
\Delta_{\rho}(t)=-M(t)^{-1}K_{\rho}(t),\qquad
K_{\rho}(t)=
\begin{cases}
K(t), & \rho=1,\\
K(0), & \rho=0.
\end{cases}
$$

rho=1 retains evolving contacts; rho=0 disables contact remodeling in the chemical exchange law. Both retain current compartment volumes. The identity 1^T M(t) Delta_rho(t)=0 preserves total diffusive amount exchange. Simply keeping the entire initial Delta fixed while cell volumes change would generally violate this identity.

**Fixed conductances do not freeze every transport effect:** changing volumes still change exchange rates per unit volume through M(t). Record this remaining route explicitly. The intervention holds chemical exchange capacity fixed; it does not hold physical geometry fixed.

The four arms correspond to rho=0/1 and sigma=0/1 in:

$$
\frac{da_i}{dt}=\frac{a_i^2}{b_i}-a_i
 +D_a\bigl(\Delta_{\rho}(t)a\bigr)_i
 -\sigma a_i\frac{\dot V_i}{V_i},
$$

$$
\frac{db_i}{dt}=\beta(a_i^2-b_i)
 +D_b\bigl(\Delta_{\rho}(t)b\bigr)_i
 -\sigma b_i\frac{\dot V_i}{V_i}.
$$

sigma=1 is the original dilution rule. In the stepped implementation it multiplies both chemicals by V_old/V_new after mechanics. sigma=0 omits that conversion; reactions and conservative diffusion remain active. Dilution off is a **diagnostic counterfactual**: maintaining concentration as volume changes introduces/removes chemical amount. Record the resulting amount change, rather than reporting it as conservation failure or as a biologically complete replacement model.

For a chemical concentration c just after reaction/transport, the off-arm conversion contributes the known amount change (V_new-V_old)c in each cell. On-arms must conserve amount during volume conversion within **2e-14**; off-arms must match this explicit source term within that tolerance. Chemical reactions have their own amount production/decay in every arm.

## Validation and measurements

- Use resident PyTorch chemistry/matrices and the existing custom CUDA mechanics/geometry/polarity with float64 phase-update carry. Put experimental switches in a separate runner/backend; leave the active study's pinned sources and outputs unchanged.
- Run dt=0.00375 and 0.001875, with a **60-unit numerical pilot followed by the full 240-unit horizon**. A noninitiating pilot is not a scientific failure; the frozen references for histories 7 and 8 form much later. A failed numerical pilot blocks its own long pair under the existing rule.
- The original maximum log-chemical refinement tolerance remains **0.01**, with the original transport, volume, polarity, shape, growth and timing limits. Do not time-align curves to pass. Keep volume error below 5%, equivalent radius at least four grid spacings, no clipping, finite positive states and sampled boundary occupancy below 0.01.
- Validate the altered chemical update against an independent float64 CPU calculation, including conservative flux, explicit dilution-off amount accounting, and exact checkpoint/restart preservation of the fixed matrix and phase residual. Baseline switches must reproduce the unchanged carry path. Save both the geometric contact matrix and the matrix actually used for chemistry.
- Measure sustained formation first; report contrast curves, concentrations, onset, cell volumes/dilution rates, conductances and instantaneous spectra second. Spectra use the **effective chemical operator** in that arm. They remain frozen-state diagnostics, not proofs of full moving-system stability. Any endpoint chemical assay must likewise use the effective operator; in a fixed-conductance arm this differs from the current geometric contact operator.

If both histories remain eligible, this requires **12 new moving paths**: two histories, three new arms and two timesteps. Four completed baseline paths are reused, giving 16 nested paths in the factorial comparison. History 9's existing successful chi=0 carry path is a positive formation benchmark, not a new factorial replicate. No new histories or parameter-map expansion are proposed.

## Mechanistic decision table

Interpret only qualified full-window comparisons within each history:

| Rescue relative to the failed moving baseline | Supported interpretation |
|---|---|
| Dilution off forms; fixed conductances alone does not | Removing explicit dilution is sufficient to restore formation in this tested coupled context. Evolving contact changes may still contribute. |
| Fixed conductances forms; dilution off alone does not | Preventing contact-conductance remodeling is sufficient to restore formation; explicit dilution alone does not account for the failure. |
| Both single interventions form | Each intervention can restore formation; do not assign a unique cause. |
| Neither single intervention forms; combined forms | Joint removal restores formation; a combined restriction is supported. A binary rescue pattern alone is not a quantitative interaction estimate. |
| Combined also fails | Dilution/contact-remodeling ablation is insufficient; the failure remains unexplained by these two interventions. |

Report the combined arm even if a single arm rescues. A combined-arm nonrescue after a single-arm rescue is a nonmonotonic result to investigate, not something to discard. If the two histories have different rescue patterns, report different mechanisms. A shared mechanism is supported only if the same qualified intervention result repeats in both, and then only for these tested contexts.

If the combined arm fails, the next discriminating diagnostic is a **chemical-capacity clamp**: use M(0), K(0) and no dilution while mechanics continues. This removes the remaining volume dependence of chemical exchange. It should reproduce the actual-input frozen chemical trajectory within the original 0.01 log-error limit; failure of that check indicates an implementation/integration mismatch before any new mechanism claim. A qualified rescue relative to the combined arm identifies the changing-capacity route in that comparison, without making the clamp a physical embryo model.

After the volume/contact routes are assessed, compare a **same-history no-feedback geometry at the same mature time t=150** against the current conditioned geometry, with chi=0, the same lineage-indexed near-uniform concentration perturbation and matched continuation parameters. Record inherited polarity, volumes and total amounts; this is a different geometric preparation, not an identical physical start. Do not substitute an earlier t=90 geometry and confound geometry with maturation. Such a test addresses prior geometric conditioning; it does not by itself identify which force caused it or validate fresh carry development from a zygote.

## Status

The reporting rule and follow-up are specified; **no new control runs have been launched**. Finish and qualify the current carry comparisons and frozen-reference matching first. Outcome-level claims remain subject to the [precision audit](results_ledger.md#phase-update-precision-audit), and the old quantitative trajectory measurements are not reused as carry measurements.
