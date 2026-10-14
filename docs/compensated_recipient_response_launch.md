# Compensated delayed-response repeat: launch record

The independently qualified compensated solver now repeats the complete delayed initial-amount assay. It runs on physical GPU 1 with resident float64 PyTorch arrays and matrix operations, alongside eight single-thread CPU workers. It calculates new DOP853/Radau references before full GPU production; their raw and pulse-normalized agreement must pass first.

There are 144 logical starts, 48 controls and 96 pulse labels, with 120 distinct within-label problems and 118 distinct physical initial-value problems ignoring labels. All three timesteps give 432 GPU paths and 432 CPU verification paths; two independent references give 288 trajectories. These are nested numerical and intervention comparisons in three existing mature histories, not new histories.

The separate preparation audit reconstructs the t=570 starts from accepted unpulsed fine endpoints. All operators, volumes, cell IDs and concentration arrays equal the original failed delayed preparations exactly. It also proves that all 144 starts are covered by the accepted GPU/CPU/restart prefix qualification. That prefix is numerical evidence only; no old scientific trajectory is reused as a result of the new study.

The model and scientific gates are unchanged. The only new retained state is a numerical rounding correction, initialized to zero and required in every checkpoint. The original failed study remains unchanged. The results ledger, manuscript, accepted parents and shared historical kernels are preserved. Acceptance requires a separate full saved-array reconstruction. Moving-geometry confirmation remains the next scientific follow-up after this numerical issue is resolved.

This record describes the original launch. The complete run and separate array review now pass all unchanged gates; the [completed assessment](compensated_recipient_response_assessment.md) reports the accepted results. The coordinator took about 4.8 minutes.

[Protocol](compensated_recipient_response.md), [backend assessment](compensated_chemistry_assessment.md), [preparation audit](../archive/study-launches/2026-10-08-compensated-recipient-response/preparation-audit.json). Status: `output/studies/compensated-recipient-amount-delayed-response/status.json`.
