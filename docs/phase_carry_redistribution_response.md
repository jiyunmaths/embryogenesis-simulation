# Delayed responses after amount-preserving chemical redistribution

## Question and controlled history

**Does changing the initial allocation of chemicals across surrounding cells change a recipient's later response to a new challenge, when that original intervention preserves the aggregate's chemical amounts?** This asks about response behavior after reorganization, rather than assigning an identity from concentration alone.

The [completed distribution/dosage factorial](phase_carry_context_dosage_assessment.md) found maintained chemical-state changes in both selected recipients of each of histories 7, 8 and 9 after an amount-preserving redistribution. Initial dosage-only changes caused transients but stayed below the declared late-state screen. The [earlier delayed-response assay](phase_carry_delayed_response_assessment.md) found selective response changes after raw resets; those resets also changed initial chemical amounts. The present assay combines the amount-preserving intervention with the delayed challenge.

Start from the **actual fine unperturbed t=510 endpoints** of the accepted redistribution and sham paths. The original redistribution occurred at t=450: the recipient's chemistry was unchanged, while every other cell received reference concentrations rescaled separately for each species to preserve its initial surrounding amount. It was a one-time intervention, with no reservoir or continuing clamp. **Amounts were matched at t=450, not at the t=510 challenge.** Reactions and co-evolving mechanics can subsequently change amounts, concentrations, volumes, polarity and contact conductances. Measure and report the t=510 target and whole-aggregate amounts in each background; do not describe the challenge backgrounds as amount-matched.

The backgrounds retain their actual phase fields, chemistry, polarity, cell IDs, lineage, inherited float64 update residuals, rounding counters and random streams. Both timestep levels start from the same fine saved endpoint within each background. Redistribution and sham can have different physical states at t=510. They are not equilibrated, frozen, or reset again.

Replication is **three existing mature developmental histories**, with two selected recipients nested in each: 22/23 in history 7, 19/20 in history 8, and 29/30 in history 9. All backgrounds have 16 cells. There are no new zygote histories, new divisions or prescribed identity classes. Model time, concentrations and amounts remain dimensionless.

## Continuing controls and two pulse doses

Continue chemistry, polarity, mechanics, volume-weighted conservative transport and geometric dilution together for **60 model units, t=510–570**. Each distinct background has its own unperturbed moving control. The sham control is shared by both recipients within a history; the two redistributed backgrounds require separate controls because they originated from different recipient-preserving interventions.

At t=510 apply a single negative activator pulse only to the chosen recipient. Its inhibitor and every other cell's concentrations remain unchanged. The pulse preserves the saved physical fields and residual carry; chemistry-dependent constitutive values may immediately change. There is no compensating chemical redistribution or persistent external source.

Use two separate pulse designs:

1. **Fractional depletion:** multiply the recipient's activator concentration by 0.9 in each background. This removes 10% of its own initial activator amount.
2. **Common removed amount:** remove the same feasible activator amount from the sham and redistributed backgrounds, equal to 10% of the smaller target amount at t=510.

For target cell `i` in background `b`, let `V_i,b` be its actual GPU-measured volume and `a_i,b` its activator concentration. Its initial activator amount is `N_i,b = V_i,b * a_i,b`. Define:

$$
\delta N_i = 0.1\min\left(N_{i,\mathrm{sham}},N_{i,\mathrm{redistributed}}\right),
\qquad
f_{i,b}=1-\frac{\delta N_i}{N_{i,b}}.
$$

Here `delta N_i` is the common removal and `f_i,b` is the concentration multiplier in that background. Both resulting concentrations remain positive, and each common-amount pulse removes at most 10% of the target's activator. It can be much smaller in the activator-rich background. Record requested and measured removed amounts, measured starting volumes, starting target amounts, pulse factors and actual initial log amplitudes. Requested and measured removal must agree within the inherited pulse-accounting tolerance, relative error 2e-13.

Fractional pulses control relative depletion; common-amount pulses control the amount removed. They do not generally control the same concentration displacement or log amplitude because the backgrounds can have different volumes and concentrations. Neither design equalizes background amounts at challenge time.

## Logical paths and exact reuse

For each history and timestep there are five sham paths: one unperturbed control and two recipients times two pulse doses. There are six redistributed paths: two recipient-specific backgrounds, each with one control and two pulses. The design therefore contains **66 logical paths**, including 18 controls and 48 pulse paths, at dt=0.00375 and 0.001875.

The earlier delayed assay may supply **18 exact reused sham paths**: its unperturbed controls and fractional pulses at both timesteps. Its common-amount pulses cannot be assumed reusable because their doses were computed from raw-reset endpoints rather than the new redistributed endpoints. The expected workload is **48 new paths and 18 reused paths**, conditional on exact equality of every regenerated preparation and validation of the earlier completed evidence. The frozen protocol records the actual counts. Parameter agreement or numerical closeness is insufficient: array dtypes and values, retiming, carry, rounding, metadata and random streams must match. A changed reused preparation blocks reuse and launch for review; no arm is silently substituted.

Read reused histories, field snapshots and checkpoints in place. Retain their original protocols, preparation, prefix, pilot and full evidence. Before preparing and continuing new paths, estimate storage from actual saved-state sizes and preserve a **2 GiB free-space reserve**. If the reserve cannot be maintained, stop before another continuation and retain committed checkpoints. No protected scientific input is deleted or recompressed as part of this assay.

## Primary response outcome

Separate a **maintained chemical state** from **response behavior**. State describes activator/inhibitor concentrations and their ongoing trajectories without the new challenge. Behavior describes how a pulsed path deviates from its own unperturbed moving background. Large concentration differences do not automatically imply a large response difference.

For species `s` in the target, let `c_s,b,pulse(t)` and `c_s,b,control(t)` be its pulsed and own-control concentrations. The signed normalized response is:

$$
R_{s,b}(t)=
\frac{\ln\!\left[c_{s,b,\mathrm{pulse}}(t)/c_{s,b,\mathrm{control}}(t)\right]}
{|\ln(f_{i,b})|}.
$$

Use the **actual factor for that particular pulse**. A common normalization by `abs(ln(0.9))` would incorrectly shrink the response to a smaller common-amount pulse. Both species are included equally. Normalization controls the initial log displacement, not nonlinear dose dependence or differences in the underlying state.

For each dose separately, the primary effect is the full-window time RMS of the redistribution-minus-sham target response:

$$
E_i=
\sqrt{\frac{1}{60}\int_0^{60}
\frac{1}{2}\sum_{s\in\{a,b\}}
\left[R_{s,\mathrm{redistributed}}(t)-R_{s,\mathrm{sham}}(t)\right]^2\,dt}.
$$

Time here is elapsed since the t=510 pulse. Compute the integral by the existing trapezoidal rule on aligned observations every 0.15 units. The prespecified screen is **strictly greater than 0.01**. Report fractional and common-amount decisions separately, whether each recipient exceeds the screen at both doses, and history-level recipient counts. Preserve below-screen outcomes. The screen is an operational response difference, not a biological identity boundary, equivalence test or statistical significance test. The six-unit pilot qualifies numerics; its effect classification does not replace the full 60-unit scientific outcome.

Also report unnormalized waveform differences, normalized activator-response area, target activator and inhibitor peak gains, effects in other cells, actual pulse doses, and fractional-versus-common-amount normalized waveform differences within each background. State differences between unperturbed backgrounds and their chemical contrast are descriptive companions to the response assay; no new state-success criterion is added. Retain the distinction between persistent baseline differences and recovery from the new pulse.

Recovery is the first sampled time after which all later deviations remain within 10% of their initial displacement, with at least 24 subsequent observed units. Target recovery uses the maximum absolute two-species target log deviation and the actual initial log pulse amplitude. Network recovery uses the initial-volume-weighted two-species log distance and 10% of its own initial value. Recovery is toward each background's continuing control, not toward the original sham state. Failure to recover by the horizon is right-censoring, not proof of permanent memory. Do not require recovery as a scientific positive outcome: agreement of recovered or censored status is a numerical gate.

## Unchanged simulation and numerical acceptance

Use PyTorch GPU arrays and matrix operations with the qualified custom CUDA mechanics and geometry/polarity kernels. Phase fields remain float32 with inherited float64 residual carry, and chemistry remains float64. Parameters remain beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35, chi=0.35, grid 72 cubed, extent 2.24 and interface width 0.085. No previously pinned source, scientific kernel, model parameter or acceptance tolerance is changed.

Before launch, run implementation/regression tests and four full-size native/carry preflight prefixes: the low-amplitude common-amount sham pulses for recipients 20 and 30, at both timesteps. Every other new preparation requires its own context gate; exact reused paths require their immutable accepted gates. A gate checks eight independent chemical steps with maximum log error 1e-11, amount conversion within 2e-14, four exact GPU restart steps, and the inherited strict 0.6-unit native/carry comparison. Pulse and endpoint preparation must preserve positive finite chemistry and exact recipient-only changes.

Within each history, all 22 logical six-unit pilots must pass their paired numerical and physical gates before any new long continuation of that history. Sample every 0.15 units and checkpoint every three units. Coarse/fine observations align exactly, without temporal shifts. Require maximum raw log-chemistry, polarity, relative-axis and relative-transport discrepancies at most 0.01; relative-volume discrepancy at most 0.005; final phase-field discrepancy at most 0.02; agreeing contrast outcomes and loss-time discrepancy at most 0.3. Across-cell contrast remains SD of ln(activator), screened at 0.1 throughout the window and in the final 24 observed units when available.

Each pulse's full-network pulse-minus-control log error, normalized by its actual log amplitude, must be at most 0.01. Relative target activator-area error must be at most 0.02. Target and network recovery presence must agree, with times agreeing within 0.3 when recovered; two censored paths can pass numerical agreement. Response-effect classifications must agree at both timesteps and their RMS discrepancy must be at most 0.01. Full-window refinement remains mandatory after pilot acceptance.

Every-step target-volume error must stay below 5%, equivalent radius at least four grid spacings, clipping zero, and dilution amount-conversion error at most 2e-14. Sampled boundary occupancy must stay below 0.01. Any pilot failure blocks that history's long continuation and stops the coordinator for review. No post hoc threshold relaxation, time realignment or favorable-recipient selection is allowed. Retain the inherited exact time-zero restart wrapper; it can archive and reinitialize only a committed state proven identical to its pinned source, never a changed or evolved state.

## Interpretation and remaining separations

A response effect at both pulse doses would show that an initial redistribution at fixed aggregate species amounts can alter later behavior, as well as maintained concentrations. A below-screen result would retain the state finding while narrowing its extension to this challenge behavior. Neither outcome is prescribed in advance. Report all three histories and both recipients in each, without treating paths, dose arms or timesteps as independent developmental replicates.

Because chemistry, geometry and polarity evolve together between t=450 and t=510, this comparison includes their consequences. It does **not** establish that response differences are independent of challenge-time abundance, volumes or mechanics, or that direct chemical transport is the sole mediator. Additional bulk-dosage pulse arms are unnecessary for this narrow initial-intervention question; they would be needed for a full behavioral distribution-by-dosage factorial. Challenge-time chemical matching and frozen-geometry/replay controls are separate designs for resolving current-state and mechanical mediation, respectively.

This is a finite-horizon test of network-supported, context-sensitive chemical behavior. It does not establish autonomous or inherited identity, irreversible fate commitment, stationary equilibria, network necessity, nearest-neighbor specificity, polarity-specific causation or fresh differentiation from a carry-corrected zygote. These mature histories inherit earlier development and do not replace older ledger claims wholesale.

## Execution and provenance

Use `python -m embryo.phase_carry_redistribution_response prepare`, then `preflight` and `run`, with output `outputs/phase-carry-redistribution-response`. Use the dedicated qualified GTX 1080 Ti with `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1`. The frozen protocol pins this document, the new runner/tests, all inherited scientific sources, qualified runtime/CUDA library hashes, the accepted distribution/dosage review, the earlier delayed review required for reuse, parent endpoint evidence and every prepared start.

Save preparation and dose records, per-path context gates, pilot/full comparisons, physical audits, status and summary separately. Successful preparation, launch or runtime completion is not scientific acceptance. The completed assay requires an independent reconstruction from raw histories and provenance verification before promotion to the results ledger or manuscript.
