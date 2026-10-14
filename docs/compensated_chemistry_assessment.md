# Compensated chemistry qualification

The isolated float64 solver passes the prespecified qualification and a separate review of the saved arrays. Chemical equations, conservative operators, initial conditions, timesteps and acceptance thresholds are unchanged. Thirteen tests cover small accumulated increments, nonlinear second-order convergence, amount conservation, positivity, CPU/PyTorch agreement and checkpoint recovery.

The original delayed control failed because small updates disappeared when added to larger concentrations. Retaining the rounding correction resolves its approximately 1.57e-13 relaxation. Its maximum error against the independent DOP853/Radau references is now 1.44e-15 at each timestep, compared with 1.56e-13 in the original fine calculation. The dual-reference disagreement is 1.55e-15. The original 100-times-noise rule therefore classifies the new control's order as **unresolved**. No convergence order is assigned to a trajectory at the reference-noise floor, and no tolerance is increased.

The full 60-unit qualification covers 36 existing starts: 12 controls and 24 pulses, across three existing histories. It compares 108 GPU paths with 108 independent NumPy paths and 72 previously qualified reference paths read in place. Maximum raw GPU/reference error is 1.20e-6; maximum raw CPU/GPU discrepancy is 8.88e-16. Raw concentration, normalized response, response area, recovery and resolved-order checks pass.

All 144 actual delayed starts also pass short GPU/CPU checks at all three timesteps. Interrupting the GPU checkpoint at 0.3 and restoring concentrations **and corrections** gives exactly the same states, corrections and observation histories as uninterrupted calculation through 0.6. Missing corrections are rejected. The separate review verifies saved positivity, correction bounds, exact inputs, checkpoint commitments and unchanged historical evidence.

Qualification took about 2.2 minutes with eight CPU workers alongside physical GPU 1. It qualifies arithmetic for this frozen chemical assay, not a changed developmental model or the moving-geometry backend. The custom CUDA mechanics and polarity kernels are retained and inactive here. All tests and benchmark trajectories are numerical evidence within three existing histories; none adds a developmental replicate.

The original failed study remains unaccepted and unchanged. The complete delayed-response repeat is prepared independently from the accepted immediate-phase endpoints, with every new start and operator required to equal the original preparation exactly. It needs its own full-window checks and independent review before scientific acceptance.

[Method and prospective qualification](compensated_chemistry.md), [immutable qualification manifest](compensated_chemistry_validation.json), [separate array review](../archive/study-reviews/2026-10-08-compensated-chemistry/independent_review.json), [repeat protocol](compensated_recipient_response.md).
