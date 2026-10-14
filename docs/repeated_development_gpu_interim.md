# Repeated developmental GPU qualification: first interim comparison

The study launched on 10 October 2026 follows naturally scheduled development from one cell to sixteen cells. Sixteen GPU trajectories span four mechanical controls and four discretizations; four complete NumPy references run in parallel. The [prospective protocol](repeated_development_gpu.md) remains unchanged. These are numerical repetitions of one synthetic seed, adding **zero independent scientific histories**.

**The first timestep pair fails qualification.** The coupling-off GPU paths at grid 72, with timesteps 0.00375 and 0.001875, both finish fifteen divisions and reach sixteen cells. Their individual physical, event-conservation and short exact-restart checks pass. Their comparison at the same physical times fails the spatial-density and transport-spectrum limits. Full CPU/GPU and spatial comparisons, other mechanical controls and independent saved-array review remain pending.

| Measurement | Largest discrepancy | Preset limit | First-pair result |
|---|---:|---:|---|
| Occupancy-field relative L2 | 0.0395104 | 0.01 | Fail |
| Activator-density relative L2 | 0.0393918 | 0.01 | Fail |
| Inhibitor-density relative L2 | 0.0393859 | 0.01 | Fail |
| Polarity-density normalized L2 | 0.0030773 | 0.01 | Pass |
| Transport-spectrum relative norm | 0.0503397 | 0.01 | Fail |
| Concentration distribution distance | 0.00065184 | 0.01 | Pass |
| Log-concentration distribution distance | 0.00064773 | 0.01 | Pass |
| Chemical contrast difference | 0.00084023 | 0.005 | Pass |
| Relative tissue-axis difference | 0.00815896 | 0.01 | Pass |
| Relative total-volume difference | 0.00063372 | 0.005 | Pass |
| Relative whole-species amount difference | 0.00013094 | 0.005 | Pass |
| Mean-polarity component difference | 0.00005229 | 0.01 | Pass |
| Maximum cell-count difference | 1 | 1 | Pass |
| Count-mismatched observations | 3 | 30 | Pass |
| Last-abscission time difference | 0.028125 | 0.05 | Pass |

The largest field discrepancy occurs at t=12, during the eight-to-sixteen-cell transition. Endpoint occupancy discrepancy at t=18 is 0.0039555, below the limit, but this does not remove the earlier failure. The spectrum discrepancy peaks at t=11.55, when both paths have eight cells and sixteen contact edges. Equal edge counts do not imply equal conductances or spectra.

Matched parent records show spindle-axis differences below 0.092 degrees and abscission-time differences up to 0.028125. These describe where the paths diverge; they do not yet establish the cause of the failed comparison. No temporal shift or spatial registration was used. Direct mechanical feedback is disabled in this pair, so its failure cannot be assigned to polarity-dependent feedback.

The two paths have maximum individual volume errors of about 1.61%, no mechanical clipping and minimum equivalent cell radii above five grid spacings. All three mixed-count, live-furrow checkpoints pass eight restored/uninterrupted steps. Independent full restart replay remains required for acceptance.

The remaining prespecified runs continue to document the scope of the discrepancy. Fresh independent developmental histories remain gated on complete qualification. If independent review confirms these failures, the next numerical experiment should diagnose repeated-cleavage temporal accuracy before biological interpretation or expansion.

This is a provisional comparison made with the production comparison helper, not independent acceptance. [Raw first-pair comparison](../archive/study-launches/2026-10-10-repeated-development-gpu/first-pair-interim.json), [descriptive diagnosis](../archive/study-launches/2026-10-10-repeated-development-gpu/first-pair-diagnosis.json), [launch and progress record](repeated_development_gpu_launch.md).
