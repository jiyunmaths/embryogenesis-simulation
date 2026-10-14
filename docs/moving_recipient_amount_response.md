# Moving confirmation of the initial chemical-allocation control

## Question and scope

Does the selective delayed response to a chemical challenge survive when the cells, contacts, volumes and polarity evolve after the initially matched chemical-allocation intervention?

The [accepted frozen repeat](compensated_recipient_response_assessment.md) matches the recipient's initial activator/inhibitor pair and each species' surrounding total once. After a fixed wait, cells 19/29 retain detected response differences, while other tested recipients remain below the screen. Frozen physics cannot establish that this behavior survives mechanical reorganization. This test allows that reorganization.

Use the existing histories 7, 8 and 9 and their two previously selected recipients: 22/23, 19/20 and 29/30. Start from the exact accepted t=510 chemical preparations in the immediate initial-amount assay. For this first moving confirmation, use **sham-based recipient matching and the original sham physical context**. The other matching reference and initial physical context remain outside this assay. There are twelve recipient-specific backgrounds, with two chemical allocations per recipient. The two allocations begin with the same complete physical checkpoint, volumes and conservative operator. No new developmental history, division or prescribed identity class is added.

Species budgets match initially, not throughout the experiment. Their surrounding allocation fractions are inherited from the accepted preparation. Changing these fractions also changes cell activity, derived tension/adhesion and subsequent mechanical evolution; those are intended consequences. Geometry may diverge after the shared initial state. This comparison does not separately identify mediation through geometry or isolate nearest neighbors.

## Two stages

1. Reorganization: run every unpulsed background for 60 units, t=510–570, at dt=0.00375 and 0.001875 from the same fine physical start. Chemistry, conservative transport, geometric dilution, polarity and direct mechanical feedback remain active. Require pilot and full timestep checks before using the endpoints.
2. Delayed response: from each actual **new fine** t=570 endpoint, prepare the same start at both timesteps and continue through t=630. Run an own-background control, a −10% recipient activator pulse, and a common removed-amount pulse. The latter removes 10% of the smaller actual recipient activator amount across the paired backgrounds, using their actual differing volumes. Inhibitor and other cells are untouched. No target, budget or geometry is rematched at t=570.

There are 24 reorganization paths and 72 response paths: **96 new moving continuations**, no reused scientific trajectories. There are 48 pulse paths across timesteps, representing 24 fine pulse labels. Report results by three histories, with two nested recipients per history; controls, doses and numerical repetitions are not developmental replicates. Compare the delayed findings descriptively with the corresponding frozen sham-matching/sham-physics cases, rather than the entire frozen factorial.

## Arithmetic and checkpoints

PyTorch keeps chemical arrays and matrix products in float64 on physical GPU 1. Existing custom CUDA mechanics, geometry and polarity calculations remain unchanged, with visible float32 phase fields and inherited float64 phase-update carry. The isolated moving caller retains the original beginning-geometry splitting: chemical stepping, polarity, mechanics, refreshed volumes, dilution. No midpoint or time-order upgrade is introduced.

For the unchanged Gierer–Meinhardt reaction/transport increment, use the qualified compensated SSP-RK2 arithmetic. A concentration is represented numerically by high part `u` and signed correction `e`, approximately `u-e`. At dilution, the ratio `r=V_old/V_new` scales both parts. An error-free product decomposition retains the rounding term from `u*r`; a compensated addition incorporates that term and `-e*r`. This preserves the chemical correction across dilution instead of discarding it. The correction is numerical bookkeeping, not biological memory or a new species.

Initialize chemical corrections to zero only at the explicit t=510 intervention. Preserve them through substeps, dilution, wait-to-response continuation, retiming and checkpoints. At a pulse, apply the existing exact rounded high-concentration preparation and reset only the depleted activator coordinate's chemical correction; retain every other correction and the full physical/phase-carry state. Controls retain all corrections. The pulse's new rounded concentration is the declared intervention, with its removed amount checked explicitly.

Every new checkpoint requires the chemical method/version, correction, phase carry, rounding counters, clocks, protocol/job commitment, history and audit. Missing corrections cannot silently initialize a continuation. The checkpoint schema remains an isolated extension; historical shared sources and accepted evidence remain unchanged.

## Prespecified checks and stopping rules

Unit tests cover multiplication against high-precision decimal arithmetic, repeated volume changes with unequal capacities, CPU/PyTorch agreement, missing/invalid corrections, retiming, source alignment and responses against drifting controls with differing pulse factors. Frozen qualification alone does not authorize moving use.

For **each actual start** at each timestep, first run eight GPU steps with independent NumPy chemistry using the actual evolving GPU operators and volume ratios. Require raw chemical error at most 1e-11 and dilution amount error at most 2e-14. Test four-step interrupted/restored GPU continuation bit for bit, including chemistry correction, phase carry, fields, polarity, geometry and rounding counters. Then compare a 0.6-unit GPU prefix with the native CPU implementation at the unchanged strict backend limits. The native prefix is a short independent physical reference, not full-horizon native equivalence. Up to four independent CPU prefixes use two native threads each alongside one serial GPU process; BLAS threads remain one.

Require pilot checks through six units and full checks through 60 before the next stage. Keep the existing mechanical timestep limits: maximum raw chemical log error 0.01, polarity difference 0.01, relative axis error 0.01, relative volume error 0.005, relative transport error 0.01 and final phase-field error 0.02. Contrast retention decisions must match, with loss-time discrepancy no greater than 0.3. The moving scheme is not assumed second order; no convergence order is assigned from these two levels.

For every pulse against its own continuing control, require normalized response discrepancy 0.01, relative activator response-area discrepancy 0.02 and recovery-time discrepancy 0.3, preserving censoring. Compare signed two-species recipient waveforms between chemical allocations, normalized by each actual pulse's log amplitude. Require agreement of the strict full-window `E>0.01` response screen, signed waveform discrepancy 0.01 and RMS-effect discrepancy 0.01 across timesteps. Below-screen is not equivalence. Recovery uses the existing definitions and continuing controls.

Every step checks finite positive chemistry, cell-volume error below 5%, radius at least four grid intervals, no clipping and dilution closure. Boundary occupancy must stay below 0.01. A failed source, context, quality, pilot or full gate stops the coordinator and retains committed evidence without changing tolerances.

Observe every 0.15 units and checkpoint every three. Preserve pilot field snapshots and final checkpoints, without full-field movies. Reuse the three original physical-source paths by reference; do not copy them per background. Reserve up to 18 GiB for new evidence plus 2 GiB of free storage. Outputs live in `output/studies/moving-recipient-amount-response` inside the writable repository; external historical output storage is read in place.

## Assessment

Reconstruct each initial matched pair and budget, inherited physical source, actual fine response start, pulse/correction intervention, numerical gate, signed waveform and recovery from saved evidence. Report full evolving species amounts and the divergence of geometry and conductances. A coordinator pass remains pending independent saved-array review. Preserve the results ledger and manuscript until evidence is independently assessed. This assay tests context-sensitive behavior in mature moving tissue; fresh development, inherited/autonomous identity, spatial refinement and wider parameter/reference robustness remain separate.
