# Completed feedback studies: assessment

Evidence captured 29 September 2026, 01:23 UTC (28 September, 21:23 EDT). The mechanical-component, polarity-ablation, long survival refinement, and planned developmental-cohort studies are complete with no recorded failures.

## Polarity necessity control

Removing only polarity’s mechanical contribution from full coupling restores formation while retaining activity-dependent tension and adhesion. The restored arm crosses log-activator SD 0.1 at sampled time 61.2, as does baseline, and reaches 0.83349 at time 78. Full coupling reaches only 0.02243. This is not a marginal final-time threshold crossing.

| Mechanical arm | Formation: minimum log-activator SD, time 63–78 | Persistence: minimum log-activator SD, time 63–78 |
|---|---:|---:|
| baseline | 0.12322 | 1.14875 |
| tension | 0.12593 | 1.14786 |
| adhesion | 0.12289 | 1.16040 |
| tension_adhesion | 0.12558 | 1.15970 |
| polarity | 0.02212 | 1.20381 |
| full | 0.02243 | 1.21569 |

The sustained-contrast threshold is 0.1. Polarity-only mechanics reproduces suppression, and removing polar mechanics from full coupling rescues formation. Together these interventions support both sufficiency and necessity for suppression in the tested mechanical combination, on this single starting geometry, chemical perturbation, coefficient set, and horizon. They do not establish necessity across embryogenesis or arbitrary parameters.

The final homogeneous-state growth rate of the restored formation arm is +0.02898, compared with −0.04481 for full coupling and +0.03426 for baseline. This supports a transport-spectrum explanation consistent with the restored chemical contrast. All six prepared-pattern arms retain contrast, reinforcing the distinction between initiation and maintenance. Successful formation trajectories are still evolving; threshold passage is not an equilibrium test.

Polarity dynamics and their activator dependence remain active in every component control. Consequently, these findings do not prove that suppression is independent of chemical regulation. Geometry-driven polarity is a plausible route because its drive remains nonzero at unit activity. Fixing the chemical multiplier in the polarity equation at one is the next direct test of that claim. The component necessity/sufficiency screen itself has not yet received long-run timestep refinement or multiple-geometry replication.

## Three completed developmental histories

All three planned seeds (7, 8, and 9) develop patterned no-feedback starting states and complete both moving continuations through time 150. No history was replaced or excluded. Late-window minima are:

| Seed | Feedback-on contrast | Feedback-off contrast | Feedback-on initial-cell correlation | Feedback-off correlation |
|---|---:|---:|---:|---:|
| 7 | 1.14762 | 1.08683 | 0.998965 | 0.999485 |
| 8 | 1.17351 | 1.12224 | 0.998701 | 0.999202 |
| 9 | 1.25722 | 1.20476 | 0.999587 | 0.999725 |

All six arms exceed both declared criteria throughout time 135–150: chemical contrast >0.1 and initial-cell correlation ≥0.8. All six endpoint graphs have negative maximum homogeneous-state spatial growth rates, yet existing nonlinear chemical patterns remain. The frozen-endpoint attractor coexistence assay remains restricted to seed 7; negative homogeneous growth plus moving contrast does not by itself prove bistability on the other endpoints.

This is a complete three-history pilot of developed-pattern maintenance. It does not estimate a reliable population success probability, repeat direct-feedback development from the zygote, or demonstrate inheritance: the survival challenge occurs after the sixteen-cell cap stops further division.

## Timestep and numerical checks

Both seed-7 fine-timestep continuations pass every declared comparison against the coarse runs. Maximum chemical log-RMS errors are 5.539e-5 (feedback-on) and 5.340e-5 (off), versus a 0.02 limit. Maximum relative shape errors are 0.000383% and 0.000578%, versus 1%. Contrast and association classifications agree. This validates temporal agreement of those continuations, not all developmental trajectories or spatial convergence.

Across the assessed component, survival, and additional development audits, maximum per-cell volume error is 1.690%, minimum equivalent radius is 5.097 grid spacings, clipping is zero, and maximum sampled boundary occupancy is 8.53e-10. Source hashes match the prepared studies. Recomputing late-window minima and correlations from saved histories reproduces the reported classifications.

## Current claim for O1

Within this model and tested parameter regime, polarity-mediated mechanics restricts the formation opportunity for chemical differences, but does not erase established collective chemical organization. Developed patterns survive in all three planned histories, and the original long survival result is insensitive to the tested timestep halving. This is evidence for reproducible finite-horizon, cell-associated chemical differentiation, not autonomous or inherited biological cell identity. General-geometry transport accuracy, spatial convergence, and the independence of polarity-driven suppression from chemical modulation remain open.

## Evidence archive

`outputs/completed-feedback-assessment-20260929T012311Z/` contains copied protocols, statuses, histories, endpoint reports, their SHA-256 manifest, and `assessment.json`. The assessment independently recomputes all component and moving-survival late-window classifications. Earlier interim assessments are preserved unchanged.

## Follow-up: frozen-endpoint replication

After the evidence snapshot above, the unchanged chemical bistability assay was completed on both endpoint graphs of seeds 8 and 9. All 168 additional trajectories passed. Combined with seed 7, all six graphs support locally stable uniform and patterned chemical equilibria across three developmental histories. This supersedes the seed-7-only limitation stated in the original assessment; it does not change the archived snapshot. See [the full assay results and numerical checks](feedback_endpoint_bistability.md).
