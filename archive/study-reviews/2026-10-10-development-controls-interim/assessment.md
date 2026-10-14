# Interim assessment: matched developmental mechanical controls

Snapshot: **10 October 2026, 20:14 America/New_York** (11 October 00:14 UTC). [Read-only audit](review.py), [machine-readable snapshot](interim-20261011T001445Z.json).

Seven of twelve new production paths are complete: all four tension/adhesion paths and the three spatial polarity paths. The polarity timestep-half path is still active, at t=17.4/18 with sixteen cells, fifteen completed divisions and no active furrows in the 20:14:46 check. The full-combination branch's four paths remain pending. The exact reused feedback-off reference adds no new production path. No full remaining-horizon replays or automatic refinement reports have completed.

PID 470605 is alive on the designated GTX 1080 Ti; GPU utilization was 96% with 8,460 MiB used at the host check. No additional jobs were launched. All 33 pinned sources match the prospective protocol. The storage-resumed polarity-88 run reproduced all 54 earlier observations and the three t=0/3/6 checkpoint arrays exactly, as recorded in [the relocation verification](../../output-relocation/2026-10-10-output/prefix-reproduction.json).

All seven completed production reports pass their runtime gates. A separate read-only audit rechecks report/protocol matching, all event gates, saved observation boundary/connectivity/positive-chemistry checks, aggregate physical audits and small JSON/JSONL evidence hashes. Whole checkpoint hashes and reconstructed field comparisons were not independently repeated in this interim check. Every completed path reaches sixteen cells and fifteen divisions at t=18; the three short restart windows each pass eight exact steps, and the independent eight-step chemical prefix discrepancy is zero. Maximum cell-volume error is 1.617%, minimum equivalent radius is 6.233 voxels, clipping is zero and dilution amount error is at most 4.44e-16.

| Completed comparison | Maximum regular-observation spectrum discrepancy | Dense equal-count spectrum maximum | Relevant spectrum limit |
| --- | ---: | ---: | ---: |
| Tension/adhesion 104→120 | 1.410% | 1.474% | 5% regular spatial; dense spatial diagnostic |
| Polarity 104→120 | 1.452% | 1.530% | 5% regular spatial; dense spatial diagnostic |
| Tension/adhesion 120 timestep halving | 0.765% | 0.765% | 1% regular and dense temporal |

All available label-invariant observation metrics are below their respective prospective limits. This is **partial agreement, not a branch qualification**: spatial field errors, decreasing field-error trends, full-combination comparisons and twelve full restart replays remain unassessed. The dense comparisons retain unequal-count observations as exclusions; they do not compare spectra of different dimensions. Chemical contrast remains small: endpoint log-activator SD is 0.00413–0.00436 across completed paths. These are numerical qualification trajectories from one matched analytic history, not independent histories demonstrating cell identity differentiation.

## Completion estimate

From approximately 20:15 EDT, allow **roughly 5–6 hours**, giving an estimated coordinator completion around **01:15–02:15 EDT on 11 October** if no gate failure interrupts execution. Independent final review is additional.

- Remaining production: approximately **80–90 minutes**, including the short remaining polarity half-step tail and four full-combination paths. Completed paths take 8.8–9.0 minutes at grid 88, 13.4–13.9 at grid 104, 19.9–20.7 at grid 120 and approximately 37–39 minutes for grid-120 time halving.
- Automatic field/refinement assessments: approximately **20–30 minutes**. The preceding feedback-off study used 555.5 seconds for four comparison pairs plus guards; the current study has nine pairs. This scaling is an estimate.
- Twelve full remaining-horizon replays: approximately **3–4 hours**. None has yet supplied a measured runtime. Scaling the measured grid-120 production durations by the remaining cell-count-weighted workload from each saved start gives 49.6–58.8 minutes for three replays per branch, or approximately 3.4–3.9 hours across four branches. Linear and quadratic cell-count weights provide two workload proxies; checkpoint I/O and omitted dense tracing can change the actual time.

The replay phase begins late in development and retains the expensive many-cell portion, so twelve replays cannot be estimated as twelve cheap short checks. Production completion alone should occur around 21:35–21:45 EDT; it is not the study's completion. The resumed coordinator timing excludes the earlier 82.44 minutes of completed tension/adhesion production and the storage interruption. Preserve that distinction when reporting total compute time.

No qualification, accepted ledger claim or scientific history count is promoted by this interim assessment.
