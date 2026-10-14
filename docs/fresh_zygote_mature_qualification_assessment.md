# Initial mature-clock attempt: scheduling stop

The first attempted mature continuation stopped at **physical t=18, elapsed zero**, before advancing its trajectory. Its **eight-step CPU/GPU actual-start prefix passed**, including exact GPU replay. The stop was a coordinator error: the initial checkpoint is scheduled at elapsed zero, but the original `_steps` helper requires a positive duration. It is not evidence against either candidate timestep or the GPU solver.

The preserved [result](../output/studies/fresh-zygote-mature-qualification/h1101_off_coarse/gpu-result.json), [failed coordinator status](../output/studies/fresh-zygote-mature-qualification/status.json) and [data assessment](fresh_zygote_mature_qualification_assessment.json) document the attempt. No six-unit path or paired comparison completed; the other 23 paths were unexecuted. No long scientific experiment started.

The prefix's maximum relative transport difference was 6.455e-7 and maximum log-concentration difference was 4.712e-10, below the unchanged 1e-5 limits. Initial physical checks passed: volume error 1.137%, resolved radius 8.51 voxels, zero clipping and zero dilution amount error. These checks qualify only that prefix and initial state.

The [separately registered clock fix](fresh_zygote_mature_clock_fix.md) maps elapsed zero to step zero, retains the original source/evidence/failure, and reuses the exact audited start files. All positive clock calculations, solver and kernel files, integration steps and numerical gates remain unchanged. A regression test now exercises the actual coordinator loop through its initial checkpoint and restart. Acceptance still requires completion and independent review.
