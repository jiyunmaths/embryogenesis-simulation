# Delayed moving-response launch — 6 October 2026

The [delayed-response protocol](phase_carry_delayed_response.md) launched at **17:50 EDT**, using the dedicated physical GPU 1 (GTX 1080 Ti), worker PID **52366**. The worker and its GPU process were independently observed after launch. This starts a new assay; it does not qualify scientific outcomes yet.

The study applies new pulses to the actual fine unperturbed **t=510** sham/reset endpoints of the [completed surrounding-cell reset](phase_carry_network_context_assessment.md), with live chemistry, mechanics, polarity, conservative transport and dilution through **t=570**. Each background has its own continuing control. Compare a -10% fractional activator pulse with a common feasible removal of 10% of the smaller target activator amount. The most unequal pair gives about a 0.57% depletion in sham versus 10% in reset. Normalize every response by its own actual log amplitude.

The design has **66 new paths: 18 controls and 48 pulses**, nested within **three existing histories** and six selected recipient comparisons. Both timesteps start from the identical fine physical preparation within each background. Background physical states can differ after the original reset. There are no new zygote histories, stationary-equilibrium assumptions, or prescribed identity labels.

## Checks completed before launch

- **36 CPU tests pass:** 16 new checks plus 20 parent regressions. The six GPU cases are skipped in that CPU run; the two new GPU cases pass separately on the qualified device.
- **Four full-size 72-cubed native/carry checks pass:** matched-dose sham pulses in recipients 20 and 30 at both timesteps. Every check includes eight independently evaluated chemical steps and four exact GPU restart steps.
- Independent preparation review verifies **49 scientific source hashes, 890 input hashes, all 66 prepared states and nine fine background endpoints**. Array dtypes/values, residual carry, rounding counters, configuration, random streams, retiming, recipient-only pulses, positivity and matched actual removed amounts pass. All 47 inherited scientific sources remain unchanged.
- Maximum full-size prefix errors are log chemistry **4.74e-07**, polarity **8.23e-08**, relative axis **5.7e-08**, relative volume **2.85e-08**, relative transport **7.35e-07**, and final phase field **3.1e-06**, below the unchanged strict native/carry limits. Independent chemical error is **0** and amount-conversion error **2.22e-16**.

The remaining contexts require individual checks during execution. All 22 six-unit pilots within a history must pass its paired timestep gate before its long continuations. Full-window physical and response refinement remains mandatory. Low-amplitude matched pulses use their actual normalization for full-network numerical errors; small raw errors alone cannot pass that gate.

The new wrapper handles the previous duration-zero resume issue only by archiving and restarting an exactly source-equivalent committed initial state. It rejects changed initial states and leaves nonzero continuation validation unchanged. No historical kernel, source, completed evidence, ledger or manuscript is modified.

An initial preflight attempt stopped before any prefix because CUDA's default enumeration selected the other installed GPU. The existing identity gate rejected it. Setting `CUDA_DEVICE_ORDER=PCI_BUS_ID` and `CUDA_VISIBLE_DEVICES=1` restores the recorded physical-device mapping; the two small GPU tests were repeated successfully on the qualified GTX 1080 Ti. Runtime and library hashes were unchanged. [Stopped attempt](../outputs/phase-carry-delayed-response/preflight-device-order-attempt.json).

## Time and evidence

Prior comparable path timings estimate **6.91 hours of moving-path work**, excluding individual context validation, paired assessment, checkpoint overhead and any interruption. This is a planning estimate rather than a deadline. The single qualified GPU runs one coordinator; independent CPU preparation review and GPU prefix checks were performed concurrently.

- [Protocol JSON](../outputs/phase-carry-delayed-response/protocol.json), SHA-256 `d488e3268d154795d8229746ffd4a8adf0c9f689633eff51e68b00a14199442d`.
- [Launch and worker PID](../outputs/phase-carry-delayed-response/launch.json), [runtime verification](../outputs/phase-carry-delayed-response/running-verification.json), [live status](../outputs/phase-carry-delayed-response/status.json), and [worker log](../outputs/phase-carry-delayed-response/run.log).
- [Preparation/prefix verification](../outputs/phase-carry-delayed-response/launch-verification.json) and [full-size preflight](../outputs/phase-carry-delayed-response/preflight.json).
- [Independent verification script](../outputs/phase-carry-delayed-response/prelaunch-verification-script.py) and [launch script](../outputs/phase-carry-delayed-response/launch-script.py).
- [Pre-launch reporting snapshot](../archive/study-document-snapshots/2026-10-06-before-delayed-moving-response-launch/manifest.json).

The next independent mechanism test remains an amount-preserving surrounding reset and/or bulk-dosage-matched rescaling. Matching pulse doses does not control the chemical amounts added or removed by the original reset. A delayed response difference would support behavior sensitivity at the maintained-state endpoints; no difference would narrow the earlier behavior effect to the initial challenge period.
