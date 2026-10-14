# Lab notebook

This is the entry point for **dated run status, interim observations, execution decisions and operational history**. The [README](../README.md) presents the model, supported claims and open failures. The [results ledger](results_ledger.md) preserves its separately audited claim snapshot; a notebook entry does not promote a result to scientific acceptance.

## Entries

| Entry date | Record | How to read it |
| --- | --- | --- |
| 2026-10-11, completed at 21:10 America/New_York; independently reviewed 2026-10-12 UTC | [History-1103 completed timestep-refinement assessment](fresh_zygote_time_refinement_assessment.md), [data assessment](fresh_zygote_time_refinement_assessment.json), [independent saved-array review](../archive/study-reviews/2026-10-12-fresh-zygote-time-refinement/independent_review.json) | The dt=0.000234375/0.0001171875 feedback-off pair passes every inherited gate: maximum dense spectral discrepancy 0.4936%, density-field discrepancy 0.1007%. Completed t=18 with sixteen cells and fifteen abscissions in 76.5 minutes. One new numerical path, zero new independent histories. Original failure retained; full cohort remains unqualified. Next: six separately registered finer matched-control paths for history 1103. No new simulation launched. |
| 2026-10-11, launched at 19:53 America/New_York | [History-1103 timestep-refinement launch](fresh_zygote_time_refinement_launch.md), [frozen diagnostic design](fresh_zygote_time_refinement.md) | One new full zygote repeat at dt=0.0001171875 versus the immutable saved dt=0.000234375 reference. Zero new independent histories. Gates/model unchanged; original failed cohort retained. Actual-start prefix passed; independent post-run review is required. |
| 2026-10-11 | [Fresh-zygote pilot launch](fresh_zygote_pilot_launch.md), [frozen design](fresh_zygote_pilot.md) | Three registered independent histories with matched mechanical controls and developmental timesteps. The longer scientific window is separately gated. |
| 2026-10-11, stopped at 17:34 America/New_York | [Fresh-zygote stopped-design assessment](fresh_zygote_pilot_assessment.md), [data assessment](fresh_zygote_pilot_assessment.json), [saved-evidence audit](../archive/study-reviews/2026-10-11-fresh-zygote-stopped/independent_review.json) | Eighteen paths complete; eight of nine pairs pass. History 1103 feedback-off fails the 1% dense spectral gate at t=11.43; six controls remain unexecuted. The cohort is unqualified. Further full-start timestep refinement takes priority over long continuations. |
| 2026-10-11, 11:47 America/New_York | [Fresh-zygote interim assessment](fresh_zygote_pilot_interim.md), [data snapshot](fresh_zygote_pilot_interim.json) | Six paths and three pairs completed at that observation time; all pass their recorded gates. This is a saved interim observation, not live progress or final acceptance. |
| 2026-10-11, 16:13 America/New_York | [Fresh-zygote interim assessment](lab_notebook/2026-10-11-fresh-zygote-interim-161352.md), [data snapshot](lab_notebook/2026-10-11-fresh-zygote-interim-161352.json) | Fifteen paths and seven pairs complete; recorded numerical gates pass. About 4–5 hours of developmental computation remain. The final historical-document check needs explicit reconciliation of authorized manuscript edits. |
| 2026-10-11 | [Matched-development completed assessment](development_control_qualification_assessment.md) | Scoped numerical acceptance on one analytic preparation; zero new independent scientific histories. |
| 2026-10-11, import date | [Status history moved from the README](lab_notebook/2026-10-11-readme-status-history.md) | Original reports, engineering updates and superseded next steps retained as historical text. Import date is not the event date of every paragraph. |
| 2026-10-11 | [README consolidation record](../archive/study-document-snapshots/2026-10-11-before-readme-notebook-consolidation/migration.json) | Verbatim original README and transfer inventory; no scientific result is altered by the reorganization. |

Earlier execution dates, launches, failures and completed reviews remain in the [run-record catalog](lab_notebook/run_records.md). Reports remain at their recorded paths so source/evidence hashes and protocol references remain reproducible. The catalog is a navigation aid, not an acceptance list.

## Live progress

The [history-1103 full-start timestep refinement](fresh_zygote_time_refinement_launch.md) is complete, and its [independent assessment accepts the finer feedback-off pair](fresh_zygote_time_refinement_assessment.md) within the registered short developmental window. The [coordinator status](../output/studies/fresh-zygote-time-refinement/status.json), [trajectory progress](../output/studies/fresh-zygote-time-refinement/h1103_off_development_quarter/gpu-progress.json) and [run log](../archive/study-launches/2026-10-11-fresh-zygote-time-refinement/run.log) remain original execution records; their provisional/pending-review fields are not rewritten. No new matched controls or long scientific continuations have been launched by this review.

For the fresh-zygote study, local live records are:

The developmental coordinator has stopped on the recorded numerical gate. See the [stopped-design assessment](fresh_zygote_pilot_assessment.md) for the retained failure and next diagnostic; the earlier interim entries are dated observations, not the current state.

- [Coordinator status](../output/studies/fresh-zygote-pilot/status.json).
- [Frozen protocol](../output/studies/fresh-zygote-pilot/protocol.json).
- [Process metadata](../archive/study-launches/2026-10-11-fresh-zygote-pilot/process-launch.json) and [execution log](../archive/study-launches/2026-10-11-fresh-zygote-pilot/run.log).
- Per-trajectory `gpu-progress.json`, followed by immutable result and paired-comparison records.

Large study outputs are local artifacts on the data drive; they are not required to browse the public documentation. [Storage layout](output_storage.md).

## Recording future work

1. Add a dated entry or link to a dated launch/interim/assessment record. Include the observation time and timezone; preserve older entries.
2. Record the question, preparation, history IDs, interventions, protocol path/hash, numerical gates and decision. Count histories separately from nested branches, cells and repeats.
3. Separate execution status, numerical acceptance and scientific interpretation. Record failed or uninformative comparisons and unresolved conditions explicitly.
4. Put runtimes, completion forecasts, process/storage notes and interim values here. A running or completed job is not a new README claim.
5. After required checks and review, update the README only if the model, a supported claim or an open failure changes. Curating the hashed ledger/manuscript is a separate evidence update; do not rewrite recorded scientific inputs for editorial consistency.
