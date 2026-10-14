# Fresh-zygote pilot: stopped-design assessment

The coordinator **stopped on a numerical failure at 17:34 EDT on 11 October 2026**, after **18/24 trajectories and 9/12 timestep pairs**. Eight pairs pass; history **1103, feedback off**, fails the dense spectrum gate. Six registered history-1103 control trajectories remain unexecuted. No post-cleavage continuation was started. This is a stopped experiment, not a completed three-history qualification.

[Frozen design](fresh_zygote_pilot.md), [data assessment](fresh_zygote_pilot_assessment.json), [saved-evidence review](../archive/study-reviews/2026-10-11-fresh-zygote-stopped/independent_review.json), [review scope and log](../archive/study-reviews/2026-10-11-fresh-zygote-stopped/review.log). The protocol remains `728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee`. Original cases, failure records, thresholds and unused histories are preserved.

## What passed and what failed

All eighteen individual paths reached t=18 with sixteen cells, fifteen completed abscissions and no live furrows. Their recorded every-step physical, positivity and conservation gates pass. The saved-evidence audit verifies their event/lineage, regular histories, dense graphs, pressure records and evidence hashes; it also checks all actual endpoints and independently reconstructs every density checkpoint in the failed pair. Actual-start prefixes and eight-step live-cleavage restarts remain recorded runtime checks, not new full-horizon CPU/GPU replays.

| History | Branch | Density-field error | Regular spectrum error | Dense spectrum error | Pair |
| --- | --- | ---: | ---: | ---: | --- |
| 1101 | off | 0.3699% | 0.8739% | 0.9196% | Pass |
| 1101 | tension_adhesion | 0.3373% | 0.6386% | 0.6722% | Pass |
| 1101 | polarity | 0.3579% | 0.8733% | 0.9461% | Pass |
| 1101 | full | 0.3271% | 0.6386% | 0.8829% | Pass |
| 1102 | off | 0.4368% | 0.7320% | 0.7512% | Pass |
| 1102 | tension_adhesion | 0.4369% | 0.7312% | 0.7508% | Pass |
| 1102 | polarity | 0.4888% | 0.7318% | 0.9244% | Pass |
| 1102 | full | 0.5027% | 0.7309% | 0.9516% | Pass |
| 1103 | off | 0.2552% | 0.7407% | 1.1836% | **Fail** |

Both density-field and spectral limits are 1%. The failed pair passes every regular-history comparison and its density-field screen; its independently reconstructed field maximum is **0.2552%**. The separate dense equal-count spectrum maximum is **1.1836% at t=11.43**, against the unchanged **1%** limit. This is the only above-limit observation among that pair's comparable dense records. Dense unequal-count observations at t=9.075 and 11.88 remain reported and excluded from equal-dimension spectral comparison.

Across all completed paths, maximum volume error is **1.712%** (5% limit), minimum effective radius **8.497 voxels** (four-voxel limit), maximum boundary occupancy **1.38e-11** (0.01 limit), and dilution amount error at most **4.44e-16** (2e-14 limit). Weighted graph conservation/symmetry closure remains near machine precision. These checks do not override the dense spectral failure.

## Diagnosis: two contacts switch across the cutoff

At t=11.43, both trajectories have the same eight cell IDs, the same parent IDs and eight active furrows. Relative volume discrepancy is only **0.000088%**. The coarser trajectory has **15 edges**, its finer partner **17**. Two contact overlaps lie on opposite sides of the hard cutoff (0.02 times the largest overlap):

| Cell pair | Coarser overlap minus cutoff | Finer overlap minus cutoff | Coarser conductance | Finer conductance |
| --- | ---: | ---: | ---: | ---: |
| 9–14 | -1.18015e-06 | 2.43998e-08 | 0 | 0.00797835 |
| 11–14 | -2.22673e-07 | 7.37755e-07 | 0 | 0.00831347 |

The finer 9–14 edge is just 2.44e-8 above its cutoff, but the hard selection turns it into a finite conductance rather than a gradual onset. Cleavage timing also differs slightly between timesteps; the diagnostic does not assign the contact changes to thresholding alone.

![Recorded failure and explicitly post hoc sensitivity reassembly](../archive/study-reviews/2026-10-11-fresh-zygote-stopped/cutoff-diagnostic.png)

Using the same saved contacts and volumes but an identical contact-edge support in both matrices reduces the peak discrepancy to **0.7216%** with the edge intersection, or **0.7249%** with the union. Reassembling conductances without cutoff gives **0.7303%** at this time and **0.7575%** maximum over equal-count dense observations.

This supports a localized cutoff contribution to the numerical discrepancy rather than a visible volume collapse or broad physical failure. **These are algebraic sensitivity calculations on saved geometries.** No chemical trajectory was rerun with these matrices. They neither change the failed decision nor qualify a no-cutoff model, a smoother contact law or the full moving system.

## Scientific interpretation of the short endpoints

| History | Fine branch | SD of log activator at t=18 | Shape axis ratio | Frozen uniform unstable modes |
| --- | --- | ---: | ---: | ---: |
| 1101 | full | 0.003564 | 1.32700 | 8 |
| 1101 | off | 0.003589 | 1.32446 | 9 |
| 1101 | polarity | 0.003541 | 1.32703 | 8 |
| 1101 | tension_adhesion | 0.003611 | 1.32442 | 9 |
| 1102 | full | 0.001936 | 1.32291 | 8 |
| 1102 | off | 0.001942 | 1.32065 | 9 |
| 1102 | polarity | 0.001923 | 1.32292 | 8 |
| 1102 | tension_adhesion | 0.001956 | 1.32061 | 9 |
| 1103 | off | 0.003092 | 1.28525 | 9 |

These values are descriptive endpoints of the stopped pilot. Chemistry remains close to uniform, with log-activator SD around 0.0019–0.0036. They do not establish persistent differentiation, response memory or cell identity. The short developmental window was designed to qualify cleavage before separately gated long chemical observations.

Shape anisotropy is present in feedback-off controls. Polarity makes a small matched change in the shape and transport spectrum in histories 1101 and 1102, but the stopped cohort provides no accepted demonstration of a large polarity-specific shape effect or identity formation. Nine versus eight positive-growth modes refer to linear chemistry on an actual frozen endpoint; evolving geometry, dilution, polarity dynamics and nonlinear maintenance are excluded from that calculation.

The replication unit is the newly seeded history: three histories were registered, two have their complete matched controls, and the third has only its failed feedback-off pair. Paths, cells, branches and repeats are nested. Eight of nine passing pairs is a numerical inventory, not a population success probability. History 1103 cannot be omitted or replaced to promote the cohort.

## Separate document preservation issue

All **38 source hashes and 382 input hashes** match. The review verifies **1,769 unchanged historical entries** and the original archived bytes for five authorized manuscript edits, using the latest [comment-resolution record](../archive/study-document-snapshots/2026-10-11-authorized-comment-resolution/record.json). Four current files match that record; `sections.md` has further user edits including five new inline comments and does not match the last resolved version. Its observed annotated bytes are preserved in the review folder; further editorial reconciliation remains pending. The review does not edit or resolve these new comments. These documents are not simulation inputs and did **not** cause the observed stop: the coordinator stopped at the numerical pair gate before reaching its final historical-location check. The original preservation manifest is unchanged.

## Next numerical test

Register a separate **history-1103 feedback-off full developmental repeat at dt=0.0001171875**, retaining grid 120, the same analytic zygote/random streams, chemistry, contact cutoff and every acceptance limit. Compare it with the existing dt=0.000234375 path through t=0–18, including the dense t=11.43 neighborhood, physical fields, events and full recorded carry state. A short restart at the offending time would not test the accumulated cleavage-timing error.

The purpose is to determine whether another timestep halving resolves the contact-spectrum failure. If it does not, diagnose the discontinuous contact-selection law under a separately specified model change and repeat its reference, conservation, spectrum and developmental checks. Keep this pilot failed. Complete the remaining matched histories/controls under an accepted prospective design before any 60/240-unit formation or identity follow-up. No new simulation is launched by this review, and the ledger/manuscript scientific claims are unchanged.
