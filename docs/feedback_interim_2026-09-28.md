# Interim assessment of running feedback studies

Captured 28 September 2026, 15:16 UTC (11:16 EDT). These are archived observations, not final experiment outcomes. Simulations continue unchanged.

## Assessment

The interim evidence strengthens the distinction between pattern initiation and maintenance. Timestep halving closely reproduces the completed seed-7 moving trajectories over the currently available prefixes. A second developmental history generates strong chemical contrast. Component controls point toward polarity-mediated mechanics as a candidate contributor to suppressed initiation, while tension or adhesion alone do not suppress formation in the completed matched controls. None of these observations yet establishes reproducibility of moving-pattern survival across developmental histories.

## 1. Long moving-survival timestep refinement

The fine timestep is 0.00375, compared with 0.0075 for the completed reference. Both start from exactly the same physical time-90 checkpoint. Errors below use matching times and cell IDs over each available fine trajectory, with the fixed initial-volume weighting from the prepared protocol.

| Fine branch | Latest model time / target | Maximum chemical log-RMS discrepancy | Maximum relative axis-ratio discrepancy |
|---|---:|---:|---:|
| switch_on | 111.6 / 150 | 1.82e-05 | 0.000246% |
| keep_off | 121.8 / 150 | 3.15e-05 | 0.000418% |

The chemical tolerance is 0.02 and the shape tolerance is 1%. The worst chemical discrepancy is approximately 635 times smaller than its limit; the worst shape discrepancy is about 2,393 times smaller. This is strong agreement over the available prefixes, not a completed convergence pass. Neither fine branch has reached the prespecified final survival window, time 135–150.

At the common time 111.6, the switched and unchanged arms have log-activator SDs 1.08793 and 1.05189 and initial-cell correlations 0.999671 and 0.999857. Their instantaneous homogeneous-state growth rates are negative (−0.06863 and −0.02960). Existing contrast therefore persists while the instantaneous uniform reference is stabilized, consistent with the completed coarse result. No conclusion about infinite-time stability follows.

## 2. Independent developmental history

Seed 8 has reached time 71.4 of 90, with sixteen cells and no active divisions. Log-activator SD is 0.89841; seed 7 at the same time had 0.85920. The seed-8 statistic first exceeded 0.1 at sampled time 42 and continues to rise. This is early replication of strong chemical heterogeneity from a fresh zygote, not completion of the late-window formation or moving-survival test.

Aggregate axis ratios at time 71.4 are 1.36760 for seed 8 and 1.32949 for seed 7. The independent history differs geometrically, making its future switch test informative. Both are no-feedback histories: this shape difference does not demonstrate chemical-feedback-driven morphogenesis. Cell IDs across embryos are not treated as matched biological states.

The seed-8 paired switch/control will start automatically after eligible completion at time 90. Seed 9 has not started and is queued after seed 8’s paired continuations. There is currently no completed independent replication of the time-90-to-150 survival experiment.

## 3. Matched mechanical-component study

Eight of ten continuations are complete. Formation and persistence polarity-only arms are at time 52.8 of 78. These experiments start from the same mature time-18 geometry. Formation begins near homogeneous chemistry; persistence uses a pattern prepared on that frozen graph. They are not fresh zygote histories or the actual developed time-90 switch experiment.

Completed-arm late-window minima, time 63–78:

| Mechanical arm | Formation: minimum log-activator SD | Persistence: minimum log-activator SD |
|---|---:|---:|
| baseline | 0.12322 | 1.14875 |
| tension | 0.12593 | 1.14786 |
| adhesion | 0.12289 | 1.16040 |
| full | 0.02243 | 1.21569 |

The declared sustained-contrast threshold is 0.1. Baseline, tension-only, and adhesion-only formation pass this screen; full coupling does not. All four completed persistence arms pass. Passing the formation screen does not mean those trajectories have equilibrated: their contrast continues to change over the late window. These data support a state-dependent feedback response, not a claim that every mechanical term always inhibits signaling.

To assess polarity without comparing unequal elapsed times, all five arms are compared at time 52.8:

| Mechanical arm | Formation log-activator SD | Formation homogeneous growth rate | Persistence log-activator SD |
|---|---:|---:|---:|
| baseline | 0.04765 | +0.09430 | 1.11554 |
| tension | 0.04838 | +0.09428 | 1.11463 |
| adhesion | 0.04764 | +0.09424 | 1.12610 |
| polarity | 0.02107 | +0.04785 | 1.16261 |
| full | 0.02135 | +0.04784 | 1.17388 |

At this matched time, polarity-only formation contrast is about 56% lower than baseline and close to full coupling. Its spectral growth is likewise close to full coupling and roughly half the baseline value. Both growth rates remain positive here: the immediate observation is reduced amplification, not yet complete loss of supported unstable modes.

This makes polarity-mediated mechanical action a leading candidate for the combined initiation suppression. It does not establish a unique or dominant cause over the whole developmental trajectory. The polarity runs must reach their final window; a tension-plus-adhesion arm without polarity would be needed to test polarity’s necessity in the full combination. Individual terms can interact nonlinearly. Adhesion’s suppressive short-time response on an already patterned time-90 state should not be generalized to failure of formation on this different time-18 background.

## Numerical health and remaining decisions

All saved audits remain within their declared limits: maximum individual volume error 1.689%, minimum radius 5.097 grid spacings, zero clipping, and maximum sampled boundary occupancy 1.34e-10. Source hashes still match both prepared studies. No failed arm is recorded. These screens do not resolve spatial convergence or the live geometric conductance approximation.

The next decisions await the fine trajectories reaching time 150, the polarity-only component arms reaching time 78, and seed-8/seed-9 developmental and paired-survival completion. Keep the manuscript’s robustness claims pending. The current evidence supports O1 only at the level of collective chemical differentiation and finite-horizon maintenance; it does not establish autonomous or inherited biological cell identity.

## Archived evidence

`outputs/interim-feedback-20260928T151647Z/` contains copied observation histories, statuses, protocols, completed component reports, their SHA-256 manifest, and the computed `assessment.json`. Files were copied once per path while simulations continued; each history is internally consistent, and comparisons use explicitly matched observation times. Status audits can be slightly newer than their corresponding history snapshot.
