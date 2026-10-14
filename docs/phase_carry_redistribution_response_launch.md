# Delayed redistribution-response launch — 7 October 2026

The [frozen protocol](phase_carry_redistribution_response.md) launched at **18:07 EDT**, worker PID **229767**, on dedicated physical GPU 1, GTX 1080 Ti (`GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e`). An independent observation at **18:09:47 EDT** confirms the worker and its matching GPU process. **Pilot/full acceptance and scientific outcomes remain pending.** [Live status](../outputs/phase-carry-redistribution-response/status.json).

This follow-up asks whether an initial redistribution at fixed aggregate species amounts changes later challenge-response behavior. It starts from the actual fine **t=510** endpoints of the [accepted state factorial](phase_carry_context_dosage_assessment.md), comparing sham with each recipient's originally amount-preserving redistribution. **Amounts were matched at t=450, not at the t=510 challenge.** Reactions, mechanics, polarity, volumes and conductances subsequently co-evolved, so the challenge backgrounds can differ in abundance and physical state.

Continue coupled chemistry, polarity, mechanics, conservative transport and geometric dilution through **t=570**. Each background has its own unperturbed moving control. Pulse only the selected recipient's activator, using two separate designs: **10% fractional depletion** of its own activator and **common removed amount**, equal to 10% of the smaller sham/redistributed target amount. Both species' pulse-minus-own-control log responses are normalized by that particular pulse's actual initial log amplitude. Dose normalization does not remove nonlinear dose dependence or background-state differences.

The primary effect is the full 60-unit time RMS of the redistribution-minus-sham normalized target response, combining activator and inhibitor equally, screened **strictly above 0.01**. Report both doses separately and retain below-screen outcomes. Recovery is toward each background's continuing control; horizon censoring does not establish permanent memory. The effect screen is operational, without an identity, equivalence or statistical-significance interpretation.

There are **66 logical paths: 48 new and 18 exactly verified reused sham paths**, at dt=0.00375 and 0.001875. They contain 18 controls and 48 pulse paths. The six recipients are nested in **three existing mature developmental histories, with zero new histories**: cells 22/23 in history 7, 19/20 in history 8 and 29/30 in history 9. Reuse retains original evidence and requires exact preparation arrays, dtypes, retiming, residual carry, rounding metadata and random streams; numerical closeness alone is insufficient.

## Verification before launch

- **110 regression cases pass**, including **11 GPU cases**, **37 new assay CPU cases**, and 99 CPU cases overall. The tested runner, tests and fixed protocol document hashes match their implementation record.
- Independent preparation verification accepts **all 66 starts**, **53 scientific source hashes** and **2,412 input hashes**. It reconstructs pulse changes and actual removed amounts, preserves physical/carry metadata, verifies 48 new and 18 exact reused starts, and checks **52 scientific document guards**. Maximum relative pulse-amount error is **4.22e-15**, below 2e-13.
- **Four full 72-cubed native/carry prefixes pass**: the weak common-amount sham pulses for recipients 20 and 30, at both timesteps. Each contains eight independent chemical steps and four exact GPU restart steps. Maximum independent chemical-step log error is **0** and amount-accounting error **2.22e-16**, below 1e-11 and 2e-14 respectively.

| Native/carry prefix quantity | Maximum discrepancy | Unchanged limit |
| --- | ---: | ---: |
| Log chemistry | 4.73757e-07 | 1e-5 |
| Polarity | 8.22998e-08 | 1e-5 |
| Relative axis ratio | 5.47833e-08 | 1e-4 |
| Relative volume | 2.97702e-08 | 1e-5 |
| Relative transport | 7.36999e-07 | 1e-5 |
| Final phase field | 3.09944e-06 | 2e-5 |

These four short prefixes qualify their prepared starts; they do not establish full-horizon backend equivalence. Every remaining new start requires its own context gate, and reused paths retain their immutable accepted evidence. **All 22 logical six-unit pilots within each history must pass numerical and physical gates before its new long continuations.** Full-window refinement remains mandatory. No threshold relaxation, temporal realignment or favorable-recipient selection is authorized.

## Timing and storage

Measured predecessor timings estimate **5.156 hours of moving work**. Allow **6–8 hours overall**, including remaining context checks, checkpointing and paired assessment; interruptions or failed gates can change this estimate. The coordinator uses resident PyTorch arrays/matrices and qualified custom CUDA mechanics/geometry/polarity kernels, with inherited float64 phase-update residual carry and float64 chemistry.

The [separate lossless storage archive](../archive/output-compression/2026-10-07-redistribution-launch/report.md) recovered **4.638616 GiB net**, after actual blob/reference round trips restored a real trajectory and six source checkpoints byte-for-byte to independent files. Those seven restored copies remain in place. Its manifest covers 443 older files and checks 2,482 protected files/status guards; accepted parent evidence and their source/input registries remain available. [Actual restore validation](../archive/output-compression/2026-10-07-redistribution-launch/restore-validation.json).

The frozen preparation estimated **14.908 GiB** for 48 new paths and preserves a **2 GiB free-space reserve**. Exact reuse avoids duplicate trajectories. The worker checks the reserve before another new continuation and retains committed checkpoints if storage blocks progress.

An accepted effect would support later response sensitivity to the original chemical allocation intervention, including its subsequent coupled consequences. It would not establish independence from challenge-time abundance or mechanics, direct chemical transport as the sole mediator, nearest-neighbor specificity, network necessity, autonomous/inherited identity, irreversible fate commitment or fresh carry-corrected zygote development. The three mature histories inherit earlier development. An independent reconstruction of raw histories and provenance is required before promoting any outcome to the ledger or manuscript.

## Evidence

- [Protocol JSON](../outputs/phase-carry-redistribution-response/protocol.json), SHA-256 `5574abd1b40bb711676685611cb84b511bdb342a1c377345bbef5d752d78ccc1`; [fixed protocol document](phase_carry_redistribution_response.md).
- [Implementation verification](../outputs/phase-carry-redistribution-response/implementation-verification.json), [regression log](../outputs/phase-carry-redistribution-response/regression-tests.log), [independent preparation verification](../outputs/phase-carry-redistribution-response/preparation-verification.json), and [full-size preflight](../outputs/phase-carry-redistribution-response/preflight.json).
- [Launch](../outputs/phase-carry-redistribution-response/launch.json), [independent worker/GPU observation](../outputs/phase-carry-redistribution-response/running-verification.json), [live status](../outputs/phase-carry-redistribution-response/status.json), and [worker log](../outputs/phase-carry-redistribution-response/run.log).
- [Preparation verifier](../outputs/phase-carry-redistribution-response/preparation_verification_script.py), [launch script](../outputs/phase-carry-redistribution-response/launch-script.py), [running verifier](../outputs/phase-carry-redistribution-response/running_verification_script.py), and [reporting snapshot](../archive/study-document-snapshots/2026-10-07-before-redistribution-response-launch/manifest.json).
- [Machine-readable launch verification](phase_carry_redistribution_response_launch_verification.json) records exact hashes and numerical maxima for this launch observation.
