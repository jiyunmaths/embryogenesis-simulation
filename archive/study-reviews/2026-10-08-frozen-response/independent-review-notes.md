# Independent completed frozen-response review

The completed fixed-geometry experiment passes independent reconstruction. The same recipients, history8/cell19 and history9/cell29, remain above the full-window response screen under both original doses. History-level counts remain0/2,1/2,1/2. Post-t510 geometry, transport and volume-dilution evolution are not required for these detected delayed effects at the tested point and horizon, conditional on each inherited chemical state and physical context.

| History/cell | Dose | Moving fine E | Frozen fine E | Frozen minus moving |
|---|---|---:|---:|---:|
|7/22|Fractional|.009352072630|.009387436887|+.000035364256|
|7/22|Common amount|.009577582806|.009612815666|+.000035232860|
|7/23|Fractional|.000051642226|.000050106570|−.000001535657|
|7/23|Common amount|.000100652852|.000099217665|−.000001435188|
|8/19|Fractional|.072663714150|.072162170970|−.000501543180|
|8/19|Common amount|.073146350693|.072643214590|−.000503136103|
|8/20|Fractional|.006403191438|.006364827518|−.000038363920|
|8/20|Common amount|.007979022563|.007939925515|−.000039097048|
|9/29|Fractional|.072964918083|.072470308854|−.000494609229|
|9/29|Common amount|.073469151509|.072974028867|−.000495122643|
|9/30|Fractional|.006645432167|.006612097193|−.000033334974|
|9/30|Common amount|.008251146283|.008217204730|−.000033941553|

All24 fine pulse trajectories recover to their own continuing frozen controls: target recovery2.25–4.35; network recovery2.25–5.40. Coarse, diagnostic and both reference methods agree on recovery times. The minimum across-cell log-activator SD over all fine paths is1.355075, and the minimum positive saved concentration across all reconstructed frozen/moving/reference/CPU paths is.000581542. These descriptive outcomes are not newly introduced acceptance conditions.

Frozen final activator redistribution/sham ratios for cells22,23,19,20,29,30 are1.015326,1.027308,.870825,.073630,.869646,.057085. Thus large maintained concentration deficits in cells20/30 coexist with below-screen response shifts. The fractional activator response-area increase in cells19/29 is34.23%/33.94%. The first six units contain99.769–99.980% of integrated squared response differences; that descriptive early window does not replace the full60-unit primary measure.

## Numerical reconstruction

| Quantity | Maximum/range | Limit |
|---|---:|---:|
|DOP853/Radau raw log discrepancy|4.4855e-12|1e-8|
|Reference normalized response discrepancy|4.2620e-11|1e-8|
|Production/reference raw log error|1.1962e-6|.01|
|CPU/GPU production raw log error|2.4469e-13|1e-8|
|Coarse/fine raw log error|8.9755e-7|.01|
|Fine/diagnostic raw log error|2.2401e-7|.01|
|Resolved convergence orders|2.00086–2.00658|1.7–2.3|
|Production/reference normalized response error|1.1407e-5|.01|
|Relative response-area error across comparisons|7.4750e-6|.02|
|Recovery-time disagreement|0|.3|
|Effect discrepancy across numerical comparisons|2.2781e-7|.01 and same screen decision|

Passed33/33 raw reference cases,99/99 production/reference+CPU verification cases,33/33 resolved order tests,66/66 raw timestep pairs,216/216 response comparison rows and108/108 effect comparison rows. The nine comparison sets comprise dual reference, coarse/fine, fine/diagnostic, and all three production levels against both references. Eight pre-freeze independent golden solver trajectories match the final references exactly. All six duplicate-dose pairs give identical trajectories at every level and method.

The audit uses no project functions and performs no integrations. It verifies99 production paths,99 CPU paths,66 reference trajectories,66 matched moving paths, exact source preparations, frozen matrices/volumes/IDs, pulse accounting,401-point clocks, final checkpoint/history consistency and positive float64 saved states. There are33 logical starting cases and27 unique initial-value problems, with only three independent developmental histories.

All56 scientific source pins and115 input pins pass. In total510 distinct consumed artifacts are unchanged before and after review. There are2415 independently compared recorded values; the largest absolute discrepancy is1.249e-16. Internal Euler stage positivity is a runtime check in the pinned code; those unsaved intermediate states cannot be independently reconstructed from output observations. This distinction is retained in the JSON report. The audit makes no project writes. Manuscript and ledger files are outside this frozen input manifest; their unchanged status requires the root's pre-run reporting snapshot rather than extending the510-file claim.

Production GPU times are21.27s,41.03s and79.90s for coarse/fine/diagnostic, totaling142.20s. The99 CPU verification runs consume361.25 summed worker seconds and the33 dual-reference cases20.10 summed worker seconds; each pool has four workers, so these sums are not total elapsed study time.

## Interpretation boundaries

Freezing removes postchallenge dynamics jointly. Every background still inherits its distinct chemistry, graph and measured volumes at t510; this is not a pure chemical-state intervention, a shared-operator comparison, a test of formation from a zygote, a population estimate, a fate-memory test or an autonomy result. Numerical changes in E are descriptive norm differences, not fractions of a mediated effect. A future shared physical-context control should specify both a common operator and its consistent capacity/volume vector, and state whether concentrations or amounts are preserved during reassignment.

Artifacts: /tmp/frozen_response_independent_review.py, .json and .npz. The JSON includes all reconstructed gates, twelve fine moving/frozen effect pairs, full signed-waveform distances, early-window descriptive measures, late baseline chemical comparisons, recovery/area metrics, exact duplicate groups and input hashes. The NPZ preserves raw-derived signed target waveforms for plotting.
