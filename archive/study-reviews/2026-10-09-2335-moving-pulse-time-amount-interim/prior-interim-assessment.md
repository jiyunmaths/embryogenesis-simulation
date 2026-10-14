# Moving challenge-time amount control: interim assessment

Snapshot: 2026-10-09 20:07 UTC / 16:07 America/New_York. This is an interim review, not full scientific acceptance. [Prospective protocol](moving_pulse_time_amount_response.md), [launch record](moving_pulse_time_amount_response_launch.md), [literal saved-array reconstruction](../archive/study-reviews/2026-10-09-moving-pulse-time-amount-interim/independent_interim.py), [machine-readable assessment](../archive/study-reviews/2026-10-09-moving-pulse-time-amount-interim/interim.json).

## Progress and numerical checks

The coordinator is running the full continuations. All 48 native/GPU actual-start checks pass, and the separate prefix reconstruction accepts all 48. Qualification took 22.1 minutes. All 48 six-unit pilots, 12 paired-context reports and six signed-contrast comparisons pass before full continuation. This interim review verifies all 83 source and 908 input pins and reconstructs the pilot chemical response differences literally from the saved histories.

Five of 48 paths have reached 60 units; all 48 already reached the six-unit pilot horizon. Recorded simulation work is 19.4% of the full 48-path horizon. The active path at the snapshot is history 7, cell 22, R chemical allocation, coarse timestep, activator pulse. There is one completed full paired chemical background (cell 22, S allocation), but **no complete 60-unit S-versus-R response comparison yet**. The completed full background passes raw chemistry/mechanics and pulse-response refinement: maximum chemical log discrepancy 1.24e-6; normalized pulse-response discrepancy 8.30e-6 against 0.01; recovery-time differences zero on the observation grid.

Across the stable pilots and current/completed path quality records, maximum volume error is 1.12% against 5%; minimum equivalent radius is 5.107 grid intervals against four; maximum boundary occupancy is 1.97e-8 against 0.01; maximum dilution amount mismatch is 4.44e-16 against 2e-14; clipping is zero. These are recorded quality extrema, not an independent replay of every mechanical timestep.

## Early scientific readout

The six-unit pilot retains five above-screen recipients and one below-screen recipient. Both timesteps agree:

| Existing history | Recipient | Six-unit response-contrast RMS E | Above 0.01 |
|---|---:|---:|---|
| 7 | 22 | 0.060653 | Yes |
| 7 | 23 | 0.007139 | No |
| 8 | 19 | 0.121118 | Yes |
| 8 | 20 | 0.296744 | Yes |
| 9 | 29 | 0.125504 | Yes |
| 9 | 30 | 0.301441 | Yes |

The early per-history counts are 1/2, 2/2 and 2/2, matching the accepted frozen assay's full-window classifications under these primary references. Maximum pilot normalized response refinement error is 1.70e-5; maximum signed S/R waveform error is 2.47e-6, both below 0.01. Early context sensitivity therefore survives co-evolving mechanics in these pilots. The pilot integrates over six units whereas the primary result integrates over 60: RMS values cannot be compared directly, and pilot classifications do not accept the full-window outcome. Recovery requires at least 24 subsequent units, so the pilot alone cannot certify recovery.

Completed paths retain strong across-cell chemical contrast. This is maintenance and immediate reset behavior in mature tissue inherited from three existing histories, zero new histories. It does not establish autonomous/inherited cell identity or organization from a fresh zygote.

## Remaining time

The snapshot is about two hours after launch. Finished paths take approximately 7.6 minutes at the coarse timestep and 10.4–10.5 minutes at the fine timestep, including their earlier pilot work. Extrapolating each unfinished path from its measured pilot rate, calibrated against available full continuations and crediting all completed pilot work, gives **5.97 hours remaining**, including a heuristic allowance for initialization, checkpointing and final paired checks. Central predicted finish: **22:06 New York time on October 9** (02:06 UTC October 10). A practical planning window is around 22:00–23:00 New York time if throughput remains similar. This is not a statistical confidence interval or a guarantee. Scientific acceptance requires the separate completed-array review afterward.

No jobs, pinned inputs, solver code, thresholds, ledger claims or manuscript conclusions were changed by this review. All remaining work belongs to this existing study; fresh-zygote development remains the next major priority after division-specific qualification.
