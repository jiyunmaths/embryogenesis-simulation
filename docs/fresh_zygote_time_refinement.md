# Fresh-zygote timestep refinement: history 1103

This is a separately registered numerical diagnosis of the [stopped fresh-zygote pilot](fresh_zygote_pilot_assessment.md). Its history-1103 feedback-off comparison failed the 1% dense equal-count transport-spectrum limit: 1.1836% at t=11.43. Two weak contacts switched across the unchanged hard contact cutoff. Physical and density-field gates passed. Saved-matrix sensitivity calculations suggest a cutoff contribution, but do not qualify a different transport model.

## Question and design

Does another timestep halving resolve the discrepancy when the entire developmental trajectory is repeated from the original zygote?

| Path | Timestep | Steps to t=18 | Role |
| --- | ---: | ---: | --- |
| History 1103, feedback off, original finer path | 0.000234375 | 76,800 | Immutable saved reference; not rerun |
| History 1103, feedback off, new quarter path | 0.0001171875 | 153,600 | New full developmental repeat |

Both use grid 120³, extent 2.24, the original analytic spherical zygote, initial chemistry and polarity, lineage, cell-cycle schedule and all three random states. Only `dt`, `steps` and `save_every` differ. Compare every state array and metadata field after excluding configuration before execution. Preserve visible phase, float64 phase carry and chemical correction. Use PyTorch arrays/matrices and the unchanged custom CUDA mechanics/geometry/polarity kernel and runtime binary.

This is **one additional numerical path within an existing history**, not a new independent biological replicate. The original pilot, failed comparison, six unexecuted controls, thresholds and registered histories remain unchanged. Starting near the failure would omit accumulated cleavage timing differences, so the repeat starts at t=0.

## Unchanged gates and measurements

Inherit the pilot's physical limits and every timestep limit verbatim. Retain the 1% density-field and regular/dense spectral limits, 0.01 chemical distribution distances, 0.005 contrast difference, 1% shape-axis difference, 0.5% total volume/amount difference, 0.01 polarity difference, 0.05 maturity-time difference, and recorded unequal-count limits. Keep all unequal-count dense records and edge-support mismatches visible. Do not exclude t=11.43 or change the contact cutoff (0.02).

Before production, run the unchanged eight-step independent NumPy full-state actual-start prefix and exact correction-aware GPU clone check. Require positive finite chemistry, no clipping, volume error below 5%, effective radius at least four voxels, dilution amount error at most 2e-14, unchanged target inheritance and all original abscission tests. Require sixteen final cells, fifteen abscissions, no active furrows and at least three mature units at t=18. Perform the unchanged eight-step exact GPU restart checks in the live 2→4, 4→8 and 8→16 windows. These prefixes/restarts do not establish full fine CPU equivalence.

Save regular observations every 0.15 units, dense contact/conductance/spectrum records every 0.015 plus event transitions, pressure traces, events/lineage and full checkpoints at 0, 3, 6, 9, 12, 15 and 18. Use the original comparison routines and reconstruction rules with the reused reference. Algebraic no-cutoff reassembly is not an acceptance criterion.

## Stop and interpretation rules

Stop on any prefix, path, restart or timestep failure. Preserve the failed output and do not automatically rerun, reduce the timestep further, change the solver/contact law, substitute seeds or relax thresholds. Require an independent saved-array review after execution. A passing pair is scoped numerical evidence for history 1103 at these two finer timesteps; it does not retroactively accept the original pilot or finish the remaining mechanical controls.

If refinement fails, separately specify and qualify a smoother contact-selection law before pursuing that model change. If it passes, register the remaining matched-control developmental design at accepted timesteps before any long scientific continuation. No 60/240-unit formation, identity or new mechanical-control run is launched by this coordinator.

## Provenance and storage

Freeze the failed pilot protocol/status/comparison, independently reviewed assessment, saved reference report and complete evidence hashes, original software/input pins, new starts and this design before launch. Reuse the reference through a read-only symlink and reject executing reused paths. Do not add mutable manuscript drafts to simulation pins or alter the historical preservation manifest; its authorized editorial reconciliation remains recorded in the original review.

[Coordinator](../embryo/fresh_zygote_time_refinement.py), [failure-path tests](../tests/test_fresh_zygote_time_refinement.py), [independent prelaunch audit](../archive/study-launches/2026-10-11-fresh-zygote-time-refinement/audit_prelaunch.py). Execute serially on the dedicated GTX 1080 Ti. Large artifacts live under `output/studies/fresh-zygote-time-refinement` on the data drive via the existing SHA-verified path alias. The original pilot protocol remains `728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee`.
