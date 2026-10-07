# Distribution/dosage factorial launch — 7 October 2026

The [protocol](phase_carry_context_dosage.md) launched at **07:29 EDT** on the dedicated physical GPU 1, GTX 1080 Ti, worker PID **148455**. The worker and its GPU process were observed independently. Scientific outcomes remain pending.

The [completed surrounding reset](phase_carry_network_context_assessment.md) changed both chemical distribution and bulk amounts. The [completed delayed-response test](phase_carry_delayed_response_assessment.md) controlled pulse doses, but not the original reset dosage. This study separates those two initial interventions using a four-arm factorial: original chemistry, raw reset, amount-preserving reference distribution, and reference-dose original distribution. Activator and inhibitor are rescaled separately over **all nonrecipient cells**, using their measured original volumes. Recipient chemistry and the entire initial physical state stay unchanged.

Start from the same fine exchanged **t=450 state** at dt=0.00375 and 0.001875; continue coupled chemistry, polarity, mechanics, conservative transport and dilution through **t=510**. This first stage tests maintained chemical states without new pulses. The primary screen is the recipient's two-species log-state time RMS in the final 24 units, threshold 0.01. Compare spatial effects at each dosage, dosage effects at each distribution, and their signed interaction before taking its magnitude. Keep below-screen results and report by history.

There are **42 distinct paths: 24 new and 18 exact reused**, nested within **three existing mature histories**, with six selected recipients. These are 12 new recipient/arm preparations at two timesteps. Shared sham and recipient-specific raw-reset controls are reused only after regenerating their preparation temporarily and verifying exact arrays, dtypes, metadata, carry, rounding and random streams, plus completed evidence and field snapshots. Existing trajectories remain untouched and reused field snapshots are read in place. No new zygote history or biological identity class is added.

## Verification before launch

- **55 CPU tests pass:** 19 new checks plus 36 parent regressions; ten GPU cases are skipped in that CPU run. **Four new GPU cases pass separately**, with nontrivial three-cell redistribution and dosage changes at both timesteps.
- The independent preparation review verifies **51 scientific sources, 1,073 inputs and all 42 starts**, including 24 new and 18 reused states. It recomputes the factorial without calling the intervention helper and checks positivity, recipient preservation, per-species distribution ratios, actual amounts, retiming, physical arrays, inherited residual carry, rounding, random streams and configuration.
- Maximum independent concentration-formula log error is **5.83e-16**; relative amount error is **3.84e-16**, below 2e-14.
- **Four full 72-cubed native/carry prefixes pass:** both new arms for recipient 19 in history 8, at both timesteps. Each includes eight independent chemical steps, amount accounting and four exact GPU restart steps. Maximum prefix log-chemistry error is **9.29e-07**, relative-transport error **3.88e-06**, and phase-field error **3.76e-06**, below the unchanged strict limits.

Every remaining new preparation needs its own native/carry gate before its pilot. Each history's 14 six-unit paths must pass physical-quality, raw timestep and state-effect checks before new long continuations. Full-window refinement remains mandatory. The physical pilot gate explicitly blocks an already failed pilot; volume, radius, clipping, boundary and dilution thresholds are unchanged. The existing exact time-zero restart wrapper is reused without modifying its pinned source.

Prelaunch draft protocols and the corresponding source versions are archived under [prelaunch protocol drafts](../outputs/phase-carry-context-dosage/prelaunch-protocol-drafts/). They contain preparation only: no GPU prefix or moving trajectory ran against them. The final protocol adds explicit pilot physical-quality blocking and its implementation checks; all prepared-state bytes and scientific thresholds remain unchanged. Earlier pinned study documents, completed evidence, scientific kernels, ledger and manuscript remain unchanged.

## Timing and storage

Prior comparable path timings estimate **2.82 hours of moving work**, excluding remaining context checks, checkpointing, paired assessment and interruptions. Allow roughly **3–4 hours** initially; update the estimate from actual progress. One coordinator uses the qualified GPU, avoiding contention with the desktop GPU. The runtime retains PyTorch arrays/matrices and custom CUDA mechanics/geometry/polarity with float64 phase-update residual carry.

Preparation reserved an estimated 8.4 GiB for the new study plus 2 GiB headroom. Exact control reuse avoids duplicate trajectories. Before a new continuation, the worker stops if free space is below 2 GiB and retains committed checkpoints. The previous lossless archives remain untouched.

The factorial controls **initial** dosage and distribution. Reactions and moving mechanics can subsequently change both amounts and geometry. A spatial effect at matched amounts supports sensitivity to initial chemical allocation, while a dosage effect supports a bulk contribution. Both effects or a signed interaction suggest a coupled explanation. None alone establishes direct chemical mediation, neighbor-only topology, network necessity, new pulse-response behavior, autonomous/inherited cell identity, or fresh zygote development.

## Evidence

- [Protocol JSON](../outputs/phase-carry-context-dosage/protocol.json), SHA-256 `547d648ebb819e8f892c523328832f57dea7694ccf21b22fcb0848cc3ed02a2a`.
- [Implementation tests](../outputs/phase-carry-context-dosage/implementation-verification.json), [independent preparation verification](../outputs/phase-carry-context-dosage/preparation-verification.json), and [full-size preflight](../outputs/phase-carry-context-dosage/preflight.json).
- [Launch](../outputs/phase-carry-context-dosage/launch.json), [observed worker/GPU](../outputs/phase-carry-context-dosage/running-verification.json), [live status](../outputs/phase-carry-context-dosage/status.json), and [worker log](../outputs/phase-carry-context-dosage/run.log).
- [Preparation verifier](../outputs/phase-carry-context-dosage/prelaunch-verification-script.py), [launch script](../outputs/phase-carry-context-dosage/launch-script.py), and [reporting snapshot](../archive/study-document-snapshots/2026-10-07-before-context-dosage-launch/manifest.json).

Review the state factorial before adding a new endpoint response assay or expanding parameters, histories or cell counts.
