# Feedback studies: second assessment, 28 September 2026

Observation snapshot: 20:19 UTC (16:19 EDT). The previous assessment was at 15:16 UTC. Completed and partial results are separated below; running simulations were not altered.

## Main changes

1. All ten matched mechanical-component continuations have completed with passing numerical-quality screens.
2. Seed 8 has completed independent zygote development and both actual time-90-to-150 feedback-switch continuations. Both retain contrast and initial-cell association.
3. Fine-timestep seed-7 feedback-off is complete; feedback-on has reached time 148.8 of 150 in this snapshot.
4. Seed 9 is now developing, at time 15 with sixteen cells. Its small current chemical contrast is an early-stage observation, not a failed patterning outcome.

## Independent moving-survival replication

Late-window minima over time 135–150 for the completed seed-8 continuations:

| Branch | Log-activator SD | Correlation with initial cell activities | Outcome |
|---|---:|---:|---|
| switch_on | 1.17351 | 0.998701 | Contrast and cell association retained |
| keep_off | 1.12224 | 0.999202 | Contrast and cell association retained |

The required minima are 0.1 for contrast and 0.8 for correlation. Both are comfortably exceeded. The final homogeneous-state growth rates are −0.18244 with feedback switched on and −0.09654 with feedback off, with no unstable homogeneous modes on either endpoint graph. Existing nonlinear patterns again persist despite homogeneous spectral stabilization.

Seed 8 began the switch from its own developed geometry: time-90 aggregate axis ratio 1.37398, compared with 1.33501 for seed 7. This is an independent developmental history, not a new chemical perturbation on the original graph. Together, seeds 7 and 8 provide two completed histories with survival in both moving branches. Seed 9 remains in the planned denominator and is unfinished; two successes are not a reliable population-frequency estimate.

This replication concerns maintenance of established chemical differentiation. It does not repeat the full direct-feedback-from-zygote developmental comparison or demonstrate inherited biological cell identity. Frozen-endpoint bistability was established on seed 7 only; seed 8 has not yet received that attractor assay.

## Mechanical components: completed result

All arms start from the same mature time-18 geometry and finish at time 78. Formation starts near uniform chemistry; persistence uses chemistry prepared on the frozen starting graph. Entries are minimum log-activator SD over time 63–78, with a sustained-contrast criterion of >0.1.

| Mechanical arm | Formation | Persistence |
|---|---:|---:|
| baseline | 0.12322 | 1.14875 |
| tension | 0.12593 | 1.14786 |
| adhesion | 0.12289 | 1.16040 |
| polarity | 0.02212 | 1.20381 |
| full | 0.02243 | 1.21569 |

Polarity-only and full feedback fail the formation criterion, while baseline, tension-only, and adhesion-only pass. All five persistence arms retain strong contrast. Polarity-only formation never exceeds log-activator SD 0.02861 over the entire run, and full coupling never exceeds 0.02902. Their first sampled negative homogeneous growth rate occurs at time 65.4. Final growth rates are nearly identical: −0.044791 for polarity-only and −0.044806 for full coupling.

The completed intervention shows that polarity-dependent mechanical action is sufficient to reproduce the qualitative formation suppression on this matched starting geometry. Activity-dependent tension or adhesion alone is insufficient at the tested coefficients. This is stronger than the previous interim similarity, but does not establish that polarity is necessary in the full combination: a tension-plus-adhesion arm with polarity mechanics disabled is missing. Nor does this isolate a unique force-level explanation or establish the result across developmental histories.

Near homogeneous activity, the material response tanh(a−1) is small, and the adhesion modulation multiplies two such responses. Polarity can instead be driven by geometric exposure even at unit activity. This is a plausible explanation for why polarity mechanics alters the formation opportunity more strongly here. It is an inference from the specified equations and intervention outcomes; a quantitative mediation or full coupled stability analysis is still needed.

## Timestep refinement: nearly complete, still a partial assessment

| Fine branch | Available interval | Maximum chemical log-RMS error | Maximum relative axis-ratio error |
|---|---|---:|---:|
| switch_on | 90–148.8 | 5.45e-05 | 0.000381% |
| keep_off | 90–150.0 | 5.34e-05 | 0.000578% |

These are matched-time, matched-cell comparisons with the completed coarse trajectory. Limits are 0.02 for chemical log-RMS error and 1% for shape. The available discrepancies are hundreds to thousands of times smaller than those limits. Fine feedback-off completes with contrast and association retained; feedback-on retains them through time 148.8. The final paired refinement verdict awaits time 150 and the driver’s automatic classification comparison. This is temporal agreement for the original developed state, not full developmental or spatial-grid convergence.

## Numerical quality and scientific implications

Across the recorded studies, including early seed 9, maximum per-cell volume error is below 1.690%, minimum equivalent radius exceeds 5.097 grid spacings, clipping is zero, and sampled boundary occupancy is below 8.53e-10. No failed arm is recorded, and source hashes still match the prepared studies.

The evidence for O1 has strengthened: independently developed aggregates can retain cell-associated chemical differences under co-evolving chemistry and mechanics, and the initiation-versus-maintenance distinction now has a completed mechanical-component intervention behind it. This remains collective chemical differentiation, not proof of autonomous, multivariate, inherited cell identity. No new evidence establishes strong feedback-driven shape organization.

Next assess the completed paired timestep refinement and seed-9 survival when available. A subsequent targeted tension-plus-adhesion control would distinguish polarity sufficiency from necessity. The transport-geometry closure and spatial convergence remain separate unresolved limitations.

## Evidence

Archived under `outputs/interim-feedback-20260928T201929Z/`: copied histories, protocols, statuses, completed component and seed-8 reports, `assessment.json`, and SHA-256 manifest. Snapshot files are captured individually while runs continue; comparisons match physical observation times rather than wall-clock progress. The original 15:16 assessment remains unchanged.
