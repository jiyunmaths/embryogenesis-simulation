# Frozen t=510 response assay: independent design review

Read-only project review, 8 October 2026. No project files changed. The only integrations performed were the root-authorized four-path independent reference probe saved in /tmp.

## Scientific contrast

For each of the nine actual backgrounds, preserve the accepted t=510 ordered cell IDs, chemistry, saved conservative transport matrix, and measured volumes. Hold that matrix and those volumes fixed through elapsed time 60; run chemical reactions and transport only, without geometric dilution, mechanics, polarity evolution, divisions, or a continuing chemical clamp. Repeat the exact original fractional and common-amount pulse factors and retain an own-background unperturbed frozen control. Initial backgrounds are not equilibrium assumptions: their initial chemical RHS maxima are about 0.000268–0.000339.

Use the same 401 observations, elapsed 0:.15:60, normalization by each pulse's actual absolute log factor, equal two-species target weighting, full-window RMS effect, and strict E > .01 operational screen. Use each background's fixed initial volumes in network metrics. Reaction-driven chemical amounts need not remain constant; conservation is an operator property, not a constant-total-chemical gate.

Persistence of the effect under this clamp demonstrates that post-challenge geometry/transport/volume changes are not required for that effect in this assay. It does not isolate chemical-state causation from inherited geometry, measured volume, graph, or earlier mechanical history. Loss or change of an effect identifies sensitivity to this combined freezing intervention, not a unique mechanical mediator. E_frozen/E_moving is not a fraction mediated: both are norms of differences. Report signed, pulse-normalized response curves and their paired moving-versus-frozen waveform distances as continuous comparisons, alongside each arm's original E screen. The independent developmental sample remains histories 7/8/9, not the number of paths or cells.

## Actual archived starts

Independently checked all 66 original logical paths against the nine fine control starts and source checkpoints. Ordered IDs, measured volumes, saved matrices, and all source chemistry are exact. Pulses match their original factors exactly. Maximum requested/actual dose relative error is 4.218847493575595e-15.

Nine maximal transport exit rates range from 1.3316748802849017 to 1.7051238657809589. Every original matrix and their block diagonal has exactly one unchanged gm_step substep at dt=.00375, .001875, and .0009375. Initial local Jacobian real parts span approximately -3.4194 to -.5475, supporting a cheap reference calculation but not proving global stability. The smallest pulse log amplitude is .005728731320338723, history 9 sham common-amount cell 30.

Artifacts: /tmp/frozen_response_design_audit.py and /tmp/frozen_response_design_audit.json.

## Solver and batching checks

Use per-initial-condition DOP853 and Radau references at rtol=1e-13, atol=1e-15, max_step=.15, with the analytic Jacobian supplied to Radau. Root specified a normalized reference-response gate of 1e-8; retain it unchanged. The first unrestricted-step diagnostic below did not qualify that gate, whereas the second precision-qualified diagnostic did. Preserve both trajectories if assessing pulse-normalized reference errors independently. Check successful completion, all expected observation clocks/shapes, and positive finite concentrations. Require per-path maximum raw log discrepancy <=1e-8; additionally require actual normalized pulse/control response discrepancy <=1e-8 and check effect/recovery decisions rather than relying only on raw error. The worst-case two-path normalized discrepancy implied by a 1e-8 raw gate is 3.4912e-6; the corresponding conservative two-background effect bound is 6.9823e-6, already far below .01. Actual discrepancies are much smaller in the small-dose probe below.

Run the actual unchanged float64 PyTorch gm_step. A 9-operator block diagonal and up to 5 columns per background is algebraically valid: reactions are elementwise and there is no transport between contexts. Fill unused columns with positive duplicate controls, never zeros. Pin the background/path/block/column mapping and reconstruct each initial path exactly. Batch and scalar calculations can differ in floating-point summation order; independent references, a short scalar/batch equality check, and explicit diagonal/off-block checks are appropriate. The measured exit-rate check above ensures batching does not alter the original substep count. CPU PyTorch verifies these chemical equations and algorithm, not CUDA arithmetic or split moving mechanics.

All paths should meet the existing absolute chemistry .01 gate as well as the stricter scientific response checks: normalized network response error .01, relative target activator-response area .02, recovery classification agreement with time discrepancy <=.3, and effect discrepancy <=.01 with matching E>.01 decisions. Compare production directly to independent reference in addition to coarse/fine differences, avoiding shared-bias cancellation. Treat any result whose distance to the effect threshold is comparable to measured numerical uncertainty as unresolved, even if a loose maximum error gate passes.

The third dt=.0009375 is useful and cheap with batching. Expected SSP-RK2 orders should be tested only where errors are above the independent-reference discrepancy floor (existing chemistry_accuracy.py uses finest production error >100 times reference disagreement). Flat near-equilibrium controls may legitimately have unresolvable order; do not require noisy ratios to be second order. Do not relabel a finer diagnostic level as a matched moving run: only the original two timesteps have accepted moving counterparts.

## Independent small-dose probe and feasibility

The root-authorized probe integrated history 8 sham/control and common-amount cell 20, and history 9 sham/control and common-amount cell 30, each with an independently written RHS/Jacobian and both reference solvers. Four one-thread CPU workers completed the eight solves in .535 s excluding source hashing. DOP853 took .0164–.0172 s per path; Radau .1076–.1147 s. Maximum raw dual-solver log discrepancy was 2.2069923668814463e-10; normalized all-cell/species response discrepancies were 3.656910840815449e-8 and 3.8270304601789085e-8. Both methods reported target/network recovery at 2.70 for both pulses.

Independent normalized activator-response areas:
- History 8/cell20: DOP853 1.1457665435445814; Radau 1.145766534773446.
- History 9/cell30: DOP853 1.1478521466920786; Radau 1.1478521382940636.

Thus all 33 independent reference initial paths are plausibly seconds-scale; production stepping and batching overhead dominate. Two original timesteps require 16,000 + 32,000 outer steps per path; a third requires 64,000. Batching removes the 33-fold Python-call overhead without changing the ODE. No expensive phase-field arrays or CUDA mechanics are needed. Runtime for full production stepping should be measured rather than promised from reference timing.

Artifacts: /tmp/frozen_response_reference_probe.py, .json, .npz. JSON pins source histories, initial checkpoints, original protocol, script, and array archive. These are independent golden checks, not the completed full assay.

## Precision qualification before freezing the new protocol

Root then requested the same four-path diagnostic at rtol=1e-13, atol=1e-15, max_step=.15; all original probe evidence is retained. The stricter settings pass the unchanged normalized 1e-8 reference gate. Worst raw dual-solver discrepancy is 3.7259084706427195e-13; normalized all-cell/species response discrepancies are 9.627382224763892e-12 (history8/cell20) and 8.572809129248071e-12 (history9/cell30). Recovery remains 2.70 for both methods and both targets/network metrics.

Updated independent normalized activator areas:
- History8/cell20: DOP853 1.1457665347767725; Radau 1.1457665347771513.
- History9/cell30: DOP853 1.1478521382815625; Radau 1.1478521382771651.

Eight solves finish in .818 seconds wall using four workers. DOP853 takes 6017 RHS evaluations/path and .116–.123 seconds; Radau takes 5658–5868 evaluations, one Jacobian, 598–600 LU factorizations and .284–.307 seconds. Artifacts: /tmp/frozen_response_reference_probe_tight.py, .json, .npz. Recommend this precisely pinned solver scheme for all 33 initial paths; do not widen the normalized 1e-8 gate.
