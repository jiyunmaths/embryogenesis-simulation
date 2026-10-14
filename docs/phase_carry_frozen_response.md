# Frozen challenge-time geometry and delayed chemical responses

## Question

**Are the delayed response differences detected after chemical redistribution maintained when the physical environment stops changing at the challenge time?** The [completed moving assay](phase_carry_redistribution_response_assessment.md) found above-screen responses in cells 19 and 29, with no above-screen recipient in history 7. This control retains all six recipients, both original pulse designs, and each background's own control. Its outcomes are not prescribed by the moving result.

Start from the **actual t=510 states** used in that accepted assay. Freeze each background's measured conservative transport operator and cell volumes, then continue its chemical reactions and exchange to t=570. This removes ongoing mechanical evolution, changing transport and geometric dilution during the challenge window. It tests their combined necessity for the delayed differences. It does not separate these three contributions or undo their effects before t=510.

The independent units remain **three existing mature histories**, 7, 8 and 9, with two selected recipients nested in each: 22/23, 19/20 and 29/30. All contexts contain 16 cells. There are no new zygotes, divisions, imposed identity labels, reservoirs or continuing chemical clamps.

## Actual starts and frozen transport

There are nine backgrounds: one sham and two recipient-specific redistributed backgrounds per history. Use each accepted fine control's first observation at elapsed zero, absolute time t=510. Retain its exact cell order, two-species concentrations, measured volumes and saved transport matrix. Check these against the accepted preparations and source checkpoints; do not replace them with t=570 endpoints, an equilibrium, a common graph, a geometry reconstructed from a coarse grid, or a rounded approximation.

Let `a_i` and `b_i` be activator and inhibitor concentrations in cell `i`. Let `V_i` be its measured t=510 volume and `Delta` its saved conservative exchange operator. The operator is `Delta = -M^{-1} K`, with `M = diag(V_i)` and a reciprocal conductance matrix defining `K`. Keep **both `Delta` and `V_i` constant** within that background. The equations are:

$$
\frac{d a_i}{dt}
= \frac{a_i^2}{b_i}-a_i+D_a\sum_j\Delta_{ij}a_j,
$$

$$
\frac{d b_i}{dt}
= \beta(a_i^2-b_i)+D_b\sum_j\Delta_{ij}b_j.
$$

The term `a_i^2 / b_i` represents activator self-enhancement opposed by inhibitor; `-a_i` is activator loss. Activator produces inhibitor through `beta * a_i^2`, while `-beta * b_i` removes inhibitor. The exchange terms communicate concentrations through the measured contact geometry. Use `beta=2`, `D_a=0.02` and `D_b=0.55`, giving `D_b/D_a=27.5`, exactly as in the moving parent.

Exchange preserves each species' volume-weighted amount because `V^T Delta = 0`. Reactions can change that amount. With fixed volumes, the geometric dilution term is zero; this does **not** make the combined reaction–transport system amount-conserving. Verify positive volumes, finite concentrations, nonnegative off-diagonal transport, zero row sums, volume-weighted conservation and reciprocal conductances without averaging or modifying the saved matrix.

Physical fields and polarity are not advanced. Custom CUDA mechanics and geometry/polarity kernels remain unchanged; they are unnecessary once their measured chemical environment is frozen. This control inherits their previously accepted t=510 states. It does not provide a new spatial-grid, interface-width or conductance-closure validation.

## Pulse designs and controls

Every background has its own unperturbed **chemical continuation on its own fixed operator and volumes**. The sham control is shared by the two recipients within a history. For each selected recipient, pulse sham and its corresponding redistributed background at t=510:

1. **Fractional depletion:** multiply its activator by 0.9, removing 10% of its own activator amount.
2. **Common removed amount:** remove the same feasible amount in the two backgrounds, equal to 10% of their smaller recipient activator amount at t=510.

Reuse the original preparation and dose records exactly. Inhibitor, other cells, the frozen matrix and volumes remain unchanged at the pulse. Verify the requested and measured removal within relative error `2e-13`. Preserve each background's actual pulse factor and initial log amplitude; the common-amount sham pulse can be much weaker than 10%.

The original redistribution preserved species amounts at **t=450**, not at the t=510 challenge. Challenge-time chemistry, volumes, graph and amounts can differ between sham and redistribution. Freezing does not match these differences or equalize the recipient's chemical state. No moving control or moving pulse trajectory can be reused as a frozen result.

There are **33 logical starting cases**: nine controls and 24 pulses. They represent **27 unique initial-value problems**: in six background/recipient pairs, the common amount equals that background's 10% fractional removal, making the two pulse preparations identical. Retain both dose records and paths for transparent design accounting. Coarse and fine integration give **66 primary production paths**, including 18 controls and 48 pulse paths. A third timestep adds **33 numerical diagnostic paths**, for 99 new production paths in total. DOP853 and Radau independently integrate each logical case, giving 66 reference trajectories. None of these paths, doses, solver runs or timesteps adds an independent developmental history.

## Response measures and comparison with moving geometry

Sample every 0.15 model units over the same full 60-unit window. A pulse response is its signed pulse/own-control log concentration ratio, divided by that pulse's actual initial absolute log amplitude. The primary effect is the target's two-species time RMS difference between redistributed and sham normalized responses, weighting the species equally. The operational response screen remains **strictly greater than 0.01**. Below-screen is not an equivalence claim or absence of any effect.

Report all twelve recipient/dose comparisons, with history counts out of two recipients. Preserve the full 60-unit primary measure; early waveforms and six-unit RMS values are descriptive companions. Also report activator response area, target and volume-weighted network recovery, unperturbed chemical levels, and across-cell log-activator contrast. Do not require recovery or maintained contrast as a scientific positive result.

Recovery retains the original definition: first sampled time after which every remaining observation lies within 10% of the initial displacement, with at least 24 subsequent observed units. It is recovery toward that frozen background's own continuing control. An unrecovered trajectory is right-censored at t=570.

Compare the qualified frozen effects continuously with the independently reconstructed moving effects, using identical initial states, doses, observations and normalization. If an effect persists after freezing, continuing physical evolution is not necessary for that effect under the tested intervention and horizon. If the effect falls below screen or changes amplitude, report that difference; a screen crossing alone is not proof that mechanics is the sole mediator. Response amplitudes may increase, signed waveforms may change direction, or differences may become detectable in other recipients. Retain every outcome.

Persistence under freezing still combines inherited chemical and physical context. A later **shared-operator control** would isolate those graph differences; **recipient-state matching** would separately test surrounding-network sensitivity at matched target chemistry. Those are different interventions and are not included here. No interpolated geometry replay is used, so the earlier unresolved source-time-spacing screen is not bypassed.

## Numerical qualification before interpretation

Use the unchanged float64 PyTorch Gierer–Meinhardt SSP-RK2 step for production, with GPU-resident arrays and matrix operations. Freeze the source and runtime versions before launch. Batch independent contexts with a block-diagonal matrix and independent trajectory columns; all inter-context conductances are exactly zero. Padding columns, if needed, duplicate controls and are not reported as scientific paths. Before stepping, require the batched positivity substep count to equal every constituent context's count at each timestep. Compare every unpadded trajectory with unchanged NumPy/SciPy production stepping.

Use `dt=0.00375` and `0.001875` for the primary paired results, and `0.0009375` for numerical refinement. Keep the exact same starting state within each context at all three timesteps. Do not clip concentrations, shift time, change the response screen, selectively drop recipients, or relax a failed tolerance.

For every logical starting case, independently solve the fixed equations with **DOP853 and Radau**, relative tolerance `1e-13`, absolute tolerance `1e-15`, and maximum internal step `0.15`, retaining both trajectories. Require positivity and finiteness, maximum solver-to-solver raw log discrepancy at most `1e-8`, and pulse/own-control normalized disagreement at most `1e-8`. Before protocol freezing, a four-path diagnostic found that unbounded reference steps at `rtol=1e-12`, `atol=1e-14` passed raw agreement but failed this normalized screen for the weakest sham pulses (maximum `3.83e-8`). Tightening integration and capping its step reduced that discrepancy below `1e-11`; the acceptance screen is unchanged. Retain both diagnostic records. These references concern the frozen equations only.

Require all production paths to have maximum raw log discrepancy at most `0.01` against their independent references. At all three timesteps, NumPy/SciPy versus PyTorch maximum raw log discrepancy must be at most `1e-8`. Where the finest production error exceeds 100 times the independent-reference discrepancy, require observed second-order convergence in the interval `1.7–2.3`; otherwise explicitly report that the order is unresolved above reference noise.

For coarse/fine and production/reference pulse comparisons, preserve the existing response gates: full-network normalized response discrepancy at most `0.01`, relative target activator-area discrepancy at most `0.02`, matching recovery/censoring status, and recovery-time difference at most `0.3`. Require agreeing full-window response-effect classifications and RMS-effect discrepancy at most `0.01`. Check the weakest common-amount pulses using their actual amplitudes. Raw concentration agreement alone cannot qualify a small-pulse response.

Save numerical, reference and scientific outcomes separately. A failed solver or refinement gate leaves the corresponding science unqualified; retain its evidence for review. Successful runtime completion is not independent scientific acceptance.

## Execution and evidence

Run `python -m embryo.phase_carry_frozen_response prepare`, then `run`, using `--output output/studies/phase-carry-frozen-response`. This small chemical study stays inside the repository's writable workspace and avoids duplicating 3D phase-field checkpoints. Parent outputs are read in place. Use the qualified dedicated GPU with `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1` and four CPU reference/verification workers, each with one BLAS/PyTorch thread.

The protocol pins this document, the new runner/tests, unchanged chemical and response source dependencies, accepted parent protocol/status/summary, independent review and evidence manifest, exact t=510 source histories/preparations, measured contexts and pulse records. Preserve prepared states, both reference methods, production and CPU verification trajectories, checkpoint/resume clocks, timing, status and summary with hashes. Do not edit parent evidence, prior assessments, the results ledger or manuscript numerical Results to launch this control.

This is a finite-horizon mechanistic test of selected mature chemical phenotypes in three inherited histories. It does not establish autonomous/inherited biological identity, irreversible fate, nearest-neighbor specificity, fresh carry-corrected zygote differentiation, a universal polarity mechanism or population-level prevalence.
