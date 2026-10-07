# Delayed moving response after chemical-context reorganization

## Question and starting evidence

Do the different chemical states at the end of a moving surrounding-cell reset also have different responses to a **new, delayed perturbation**? The [completed reset study](phase_carry_network_context_assessment.md) finds maintained state changes in both selected recipients of histories 7, 8 and 9. Its pulses were given immediately at the reset. A descriptive decomposition places 99.30%–99.98% of their response difference within the first six model units. Those pulses do not establish response behavior after reorganization.

This study continues the actual fine **unperturbed t=510 endpoints**, rather than the previously pulsed endpoints or an artificially equilibrated frozen state. Moving backgrounds need not be stationary. Each background retains its chemistry, phase fields, inherited float64 update residuals, rounding counters, polarity, cell IDs, lineage and random streams. Both timestep levels, dt=0.00375 and 0.001875, start from that same fine physical state within each background. Reset and sham backgrounds can now have different physical states; differences include the full consequences of the original coupled intervention.

The unit of replication is **three existing histories**. The two selected recipients per history (22/23, 19/20, 29/30), controls, dose arms and timestep levels are nested interventions. There are no new zygote histories or prescribed identity classes.

## Dose controls and moving arms

For each recipient compare its sham and recipient-preserving reset endpoint. Each background has its own unperturbed moving control. Chemistry, polarity, mechanics, volume-weighted conservative transport and geometric dilution continue together for 60 model units, **t=510–570**.

Two separate activator depletions are applied only to the chosen recipient:

1. **Fractional:** multiply its activator concentration by 0.9, removing 10% of its current activator amount.
2. **Matched amount:** remove the same feasible activator amount from both backgrounds, defined as 10% of the smaller initial target amount.

For target cell i, its measured amount is N_i = V_i a_i, where V_i is its actual GPU-measured volume and a_i is its activator concentration. The common removal is delta_N = 0.1 min(N_i_sham, N_i_reset). The resulting factor in background b is f_b = 1 - delta_N / N_i_b. Therefore each matched pulse is strictly positive and at most a 10% depletion; it can be much smaller in the activator-rich background. Record both requested and actual removed amounts, factors and log amplitudes. Inhibitor and all other cell chemistry remain unchanged at the pulse. There is no redistribution or persistent external reservoir.

The previous absolute dose is not reused: it could empty a strongly depleted target. Fractional pulses control relative displacement; matched pulses control removed amount. They need not control concentration displacement or initial log amplitude because volumes and concentrations differ.

Per history and timestep there are five sham paths (one control plus two recipients times two doses) and six reset paths (two backgrounds times one control plus two doses). The full design contains **66 new paths**, including 18 shared background controls and 48 pulsed paths. Exact path reuse is not assumed even when two dose preparations are very close.

## Response definitions and scientific outcomes

For each species, a response is the signed log concentration ratio between the pulsed path and its own unperturbed background, divided by **abs(ln(f_b)) for that particular pulse**. Activator and inhibitor are included equally. Using a common normalization of abs(ln(0.9)) for the smaller matched pulses would incorrectly make their responses appear weaker.

The primary comparison is the full-window time RMS difference of these two-species target waveforms between reset and sham. The prespecified effect screen is **0.01**, separately for each dose. Report the effect in each recipient and summarize within each history; preserve negative outcomes. Report whether a recipient shifts at both doses, as well as dose-specific results. This screen is an operational response difference, not a biological cell-type boundary or a test of equivalence.

Also record the unnormalized log-waveform difference, actual doses, normalized activator-response area, target/inhibitor peak gains, effects in other cells, and recovery to each background's own continuing control. Compare fractional and matched normalized waveforms within each background descriptively to reveal dose dependence. Normalization does not remove nonlinear dose dependence or all baseline-state differences. A common-amount difference and a fractional difference need not agree.

Recovery is the first sampled time after which all remaining deviations stay within 10% of their initial displacement, with at least 24 later observed units. No recovery by the horizon is right-censoring, not permanent memory. Network recovery uses the initial volume-weighted two-species displacement. Chemical contrast is the across-cell SD of ln(activator); record whether it exceeds 0.1 throughout the window and in the final 24 units. No historical donor-reference classification is added because earlier references have a different model age and physical context.

## Unchanged numerical acceptance

Use the existing PyTorch GPU arrays and matrix operations, custom CUDA mechanics and geometry/polarity kernels, float32 phase fields with inherited float64 residual carry, and float64 chemistry. Model parameters remain beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35, chi=0.35, grid 72 cubed, extent 2.24 and interface width 0.085. No scientific kernel or earlier pinned source is changed.

Before launch, run implementation tests and four full-size native/carry prefixes: the low-amplitude matched-amount sham pulses for recipients 20 and 30, at both timesteps. Every other new context still requires its own gate. Each gate checks eight independent chemical steps (maximum log error 1e-11), amount conversion (2e-14), four exact GPU restart steps, and the existing strict 0.6-unit native/carry comparison.

Within each history run all 22 six-unit pilots before any long continuation of that history. Coarse/fine observations must align exactly, every 0.15 units. Require maximum raw log-chemistry, polarity, relative-axis and relative-transport discrepancies at most 0.01; relative-volume discrepancy at most 0.005; final phase-field discrepancy at most 0.02; agreeing contrast outcomes and loss-time discrepancy at most 0.3. Each pulse's full-network log pulse-minus-control error is normalized by its actual log amplitude and must be at most 0.01. Relative target activator-area error must be at most 0.02. Recovery presence must agree and times must agree within 0.3. Response-effect decisions must agree across timesteps; response-shift RMS discrepancies must be at most 0.01. Both-unrecovered cases can pass numerical agreement while remaining scientifically censored.

Every-step target-volume error must stay below 5%, equivalent radius at least four grid spacings, clipping zero, and dilution amount-conversion error at most 2e-14. Sampled boundary occupancy must stay below 0.01. No time alignment or post hoc threshold relaxation is allowed. Failure of any pilot blocks that history's long continuation. Full-window refinement remains required even after pilot acceptance.

Carry checkpoints are committed every three units. The older runner rejects a duration-zero resume. The new wrapper archives and reinitializes only a committed **exactly source-equivalent time-zero state**, after checking all array dtypes/values, carry, rounding, metadata, chemistry, physical audit and history. Changed states are rejected; nonzero states continue through the unchanged validated runner. Archived records retain hashes. This is operational recovery, not a mechanics/chemistry change.

## Interpretation and next separation

A delayed response shift would extend context sensitivity from immediate reorganization to behavior at actual maintained-state endpoints. No shift would narrow the earlier response finding to the initial challenge period, even if chemical levels remain different. This tests finite-time phenotype resilience, not autonomous cell identity, irreversible gene commitment, inheritance through division, stationary equilibria or differentiation from a fresh carry-corrected zygote.

The original reset changed both spatial distribution and total chemical dosage. Equal pulse doses do **not** control the original reset dosage. An amount-preserving surrounding reset and/or a bulk-dosage-matched rescaling remains the next distinct experiment before attributing the reset effect specifically to spatial organization. Different late geometries also remain part of the coupled context rather than an isolated chemical mediator.

## Execution

Use `python -m embryo.phase_carry_delayed_response prepare`, then `preflight` and `run`, with output `outputs/phase-carry-delayed-response`. The protocol pins inherited sources, inputs, the independently reviewed parent evidence, this document, the new runner/tests, and all prepared starts. The qualified runtime and both CUDA library hashes must match. Status, per-path audits, prefixes, pilot/full paired evidence and a summary are saved separately from scientific acceptance.
