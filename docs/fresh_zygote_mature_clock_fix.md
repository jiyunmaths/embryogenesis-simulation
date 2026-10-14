# Mature checkpoint scheduler correction

Registered 12 October 2026 after the [first mature-clock attempt](fresh_zygote_mature_qualification_launch.md) stopped at physical t=18. Its eight-step independent CPU/GPU prefix passed; the continuation itself stopped before a mechanics step because the coordinator called the positive-duration resolver on the elapsed-zero checkpoint. This is a scheduling error, not a failed physical or timestep-accuracy gate. Preserve that attempt's protocol, source, failed status and evidence.

The correction is a separate process-local coordinator adapter: elapsed zero maps to integer step zero after validating dt. Positive durations use the original resolver. Existing solver and kernel files, original coordinator source, chemical/mechanical integration, compensation, physical clocks and thresholds remain unchanged.

Use a fresh output directory, `output/studies/fresh-zygote-mature-clock-fix`. Alias the original 24 independently audited start files read-only rather than duplicating arrays. The original audit compared 432 saved arrays byte for byte; immutable file hashes and exact symlink targets carry that inheritance forward. The new protocol pins the original protocol/audit/failed attempt and the correction's source, regression checks and launch audit.

Reuse only the successful history-1101 feedback-off dt=0.00375 actual-start prefix, whose config, state bytes, solver and eight-step reference limits are identical. Record its original protocol and evidence provenance explicitly. Run the remaining 23 prefixes normally. This reuse does not reuse a six-unit trajectory: none was completed in the failed attempt.

The [original registered design](fresh_zygote_mature_qualification.md) otherwise applies without change: three existing histories, four branches, dt=0.00375/0.001875, physical t=18 to24, 24 new six-unit continuations, twelve comparisons, exact restart at elapsed3 and unchanged spectrum/field/homologous/physical gates. Zero new independent histories; no chemistry reset, feedback switch, automatic finer-clock fallback or 60/240-unit extension.

A regression smoke test exercises the complete coordinator path, including saving the initial checkpoint, advancing to the declared horizon and running the restart branch, using a minimal in-memory simulation. It does not mock away `run_path`. Actual CPU/GPU prefix evidence remains the solver check.

Stop on any failure. Passing execution remains provisional until independent review; this correction does not promote a scientific claim or retroactively accept the original attempt.
