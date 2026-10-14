# Lab notebook

This is the entry point for **dated run status, interim observations, execution decisions and operational history**. The [README](../README.md) presents the model, supported claims and open failures. The [results ledger](results_ledger.md) preserves its separately audited claim snapshot; a notebook entry does not promote a result to scientific acceptance.

## Entries

| Entry date | Record | How to read it |
| --- | --- | --- |
| 2026-10-11 | [Fresh-zygote pilot launch](fresh_zygote_pilot_launch.md), [frozen design](fresh_zygote_pilot.md) | Three registered independent histories with matched mechanical controls and developmental timesteps. The longer scientific window is separately gated. |
| 2026-10-11, 11:47 America/New_York | [Fresh-zygote interim assessment](fresh_zygote_pilot_interim.md), [data snapshot](fresh_zygote_pilot_interim.json) | Six paths and three pairs completed at that observation time; all pass their recorded gates. This is a saved interim observation, not live progress or final acceptance. |
| 2026-10-11, 16:13 America/New_York | [Fresh-zygote interim assessment](lab_notebook/2026-10-11-fresh-zygote-interim-161352.md), [data snapshot](lab_notebook/2026-10-11-fresh-zygote-interim-161352.json) | Fifteen paths and seven pairs complete; recorded numerical gates pass. About 4–5 hours of developmental computation remain. The final historical-document check needs explicit reconciliation of authorized manuscript edits. |
| 2026-10-11 | [Matched-development completed assessment](development_control_qualification_assessment.md) | Scoped numerical acceptance on one analytic preparation; zero new independent scientific histories. |
| 2026-10-11, import date | [Status history moved from the README](lab_notebook/2026-10-11-readme-status-history.md) | Original reports, engineering updates and superseded next steps retained as historical text. Import date is not the event date of every paragraph. |
| 2026-10-11 | [README consolidation record](../archive/study-document-snapshots/2026-10-11-before-readme-notebook-consolidation/migration.json) | Verbatim original README and transfer inventory; no scientific result is altered by the reorganization. |

Earlier execution dates, launches, failures and completed reviews remain in the [run-record catalog](lab_notebook/run_records.md). Reports remain at their recorded paths so source/evidence hashes and protocol references remain reproducible. The catalog is a navigation aid, not an acceptance list.

## Live progress

For the fresh-zygote study, local live records are:

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
