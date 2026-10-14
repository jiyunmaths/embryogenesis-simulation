# Moving confirmation of challenge-time chemistry and amount matching

Prospective protocol, fixed before launch on 2026-10-09.

## Question and scope

Do the allocation-sensitive immediate responses in the accepted challenge-time amount assay survive when chemistry and geometry co-evolve? Use the **primary S recipient-state reference, S surrounding-amount reference and S physical endpoint**. This confirms one predefined slice of the accepted frozen factorial. Reciprocal reference/physical-context extensions are separate follow-ups.

Keep all six previously selected recipients: 22/23 in history 7, 19/20 in history 8, and 29/30 in history 9. These recipients are nested observations; the reporting units are **three existing developmental histories**, with zero new histories. Include cell 23's below-screen comparison. Do not select only positive frozen outcomes.

## Exact initial conditions and controls

Start at t=570 from the actual fine, unpulsed waiting endpoints in `moving-recipient-amount-response`, independently accepted before this study. Source fields are referenced without projection. Shapes, target volumes, polarity, IDs, lineage, random state, material parameters, phase carry and rounding counters are preserved.

Use the exact accepted initial concentrations from `pulse-time-amount-response` for M=S, B=S, P=S. Within each recipient pair, S/R surrounding allocation differs, but recipient activator/inhibitor, surrounding amounts of each species, global amounts and physical fields agree initially. Uniform species-specific surrounding rescaling preserves within-species allocation fractions. Matching is an explicit preparation intervention, not an endogenous mechanism. Quantify these chemical changes using the pinned parent accounting.

For each chemical background, run an unpulsed moving control and a pulse removing 10% of recipient activator, with inhibitor unchanged. Because target concentration and capacity match, the same fractional pulse also removes the same amount in both backgrounds; an extra equal-amount label would duplicate the experiment. Reset the chemical correction to zero only at this explicit rounded-state intervention. Preserve it at all subsequent checkpoints and restarts.

There is no waiting/settling stage: chemistry and mechanics co-evolve immediately from t=570 to t=630. Species amounts are matched only initially; reactions and transport subsequently act freely. Each pulse is compared with **its own continuing moving control**, not a static initial state.

## Design and implementation

12 chemical contexts × pulse/control × two timesteps = **48 new labeled moving paths** and 24 labeled chemical starts. Preserve duplicate labels for transparent recipient comparisons; preparation separately counts distinct initial-value problems. All jobs run with the unchanged compensated SSP-RK2 chemistry, error-free volume dilution, conservative volume-weighted transport, phase-carry CUDA mechanics and polarity kernels. PyTorch holds GPU arrays and performs chemical matrix operations. No shared scientific solver is changed.

Timesteps: 0.00375 and 0.001875. Both start from the same actual fine endpoint fields and exact chemical intervention. The beginning-geometry splitting remains first order as a coupled moving method; two-level agreement does not establish second-order moving convergence. Observe every 0.15 units; checkpoint every 3 units; retain chemistry correction, phase carry and restart metadata. Full duration is 60 units, pilot duration 6 units and native-prefix duration 0.6 units. One serial GPU worker avoids VRAM/throughput contention; up to four two-thread native workers qualify independent starts in parallel. Use the qualified GTX 1080 Ti by UUID.

## Gates before long continuations

1. Independently reconstruct source lineage, accepted input/pulse arrays, amount matching, IDs, actual endpoint capacities and operators. Pin source code, runtime, binaries, protocol, tests and consumed inputs.
2. Qualify every actual start against native CPU mechanics at identical timestep over 0.6 units. Keep the original strict GPU validation limits. Also check eight independently stepped NumPy chemical/dilution steps on evolving GPU geometry and four bitwise checkpoint-restart steps, including chemical correction and physical carry.
3. Independently reconstruct all 48 saved-array prefix comparisons before production. Runtime exact-restart and eight-step checks remain explicitly identified as runtime/source evidence, not independently replayed long trajectories.
4. Complete all 48 six-unit pilots, all 12 timestep-pair reports and all six signed S/R contrast comparisons before any full continuation. Continue the pilot checkpoints to 60 only if every gate passes.

Full moving timestep gates: maximum chemical log difference 0.01; polarity component difference 0.01; relative axis difference 0.01; relative cell-volume difference 0.005; relative operator difference 0.01; final phase-field absolute difference 0.02. Pulse-response difference normalized by actual log pulse ≤0.01; relative activator response AUC error ≤0.02; recovery-time difference ≤0.3, with consistent censoring. Signed S/R response-contrast maximum discrepancy and RMS discrepancy must each be ≤0.01, and both timestep classifications must agree.

Quality: maximum relative volume error <0.05; minimum resolved cell radius ≥4 grid intervals; zero clipping; boundary occupancy <0.01; dilution amount closure ≤2e-14. Initial matching closure ≤2e-13. Reserve 2 GiB storage. Stop and preserve evidence on any failed gate; do not change thresholds after observing outcomes.

## Primary readouts and reporting

For background C=S/R and recipient i, define its normalized response as `r_C(t) = log(u_pulse,C,i(t) / u_control,C,i(t)) / abs(log(0.9))`, retaining both species and their signs. Compare backgrounds through the vector waveform `r_R - r_S`, then calculate its time-and-species RMS over the full 60-unit window. The screen is E>0.01, an operational numerical response criterion rather than a biological identity threshold. Require two-timestep agreement before classifying.

Report each history's two recipients separately, then counts per history; retain the below-screen case. Compare all six moving effects with the exact frozen M=S, B=S, P=S results. Also report pulse/network recovery, evolving species amounts, and geometry/operator changes descriptively. Full completion is provisional until separate saved-array review. Do not promote ledger/manuscript claims at launch.

A positive result supports immediate allocation-sensitive chemical behavior despite live mechanics under the tested primary references. A negative result would identify a limitation of the frozen reset assay. Neither demonstrates stable autonomous identity, inheritance, nearest-neighbor specificity, fresh zygote development or spatial convergence. The historical non-input exported-PDF exception remains recorded separately and is not a scientific input failure.
