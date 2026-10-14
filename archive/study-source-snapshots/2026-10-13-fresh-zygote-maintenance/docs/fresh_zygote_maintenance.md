# Inherited fresh-history patterns: long moving maintenance

Registered 13 October 2026 UTC before execution, after the [independently accepted formation study](fresh_zygote_formation_assessment.md). Test whether patterns that arose during development persist while mechanics continues to remodel their chemical contact network. This stage tests retention, loss and recovery; it does not impose a new chemical pattern.

## Histories, starts and controls

Continue each of the **24 accepted physical t=78 endpoints**, `gpu-e60.npz`, from the formation study. Use histories **1101, 1102 and1103**, the same four inherited contexts, and the same mature timesteps **0.00375 and0.001875**. There are **three existing independent histories**, four nested contexts and two numerical repeats, with **zero new independent histories**.

| Context | Activity-dependent tension | Activity-dependent adhesion | Polarity tension |
| --- | ---: | ---: | ---: |
| Feedback off | 0 | 0 | 0 |
| Tension/adhesion | 0.25 | 0.35 | 0 |
| Polarity | 0 | 0 | 0.35 |
| Full | 0.25 | 0.35 | 0.35 |

Baseline mechanics and polarity evolution remain active in every context. These branches inherit different developmental preparations; they are not interventions on one common mature shape. Keep conservative transport, dilution, beta=2, D_a=0.02, D_b=0.55, grid120, extent2.24, disabled historical fate labels, sixteen cells and every other coefficient unchanged. No further division occurs at the registered sixteen-cell cap.

Run an **additional180 elapsed units**, physical **t=78–258**. The complete mature trajectory therefore spans **240 units since t=18**, comprising the accepted60-unit formation stage and this180-unit continuation. Save maintenance elapsed time, absolute physical time and total elapsed time since18 separately. No resetting chemistry, feedback switching, replacing histories or repeating the accepted formation window.

Preserve every inherited array byte: visible phase fields and phase carry, chemical states and correction, rounding diagnostic, polarity, targets, cell IDs/parents, division-cycle due times, completed lineage/event audits and random states. Only metadata `config.steps` changes to the declared final horizon t258. Physical time remains78; initial step numbers are20800/41600 and final step numbers68800/137600. Verify starts against the independently reviewed endpoint hashes, not only newly calculated hashes.

## Solver and numerical qualification

Use the unchanged resident compensated DevelopmentGpu solver: **PyTorch GPU arrays and matrices with custom CUDA mechanics, geometry and polarity**. Visible phase fields are float32 with float64 phase residuals; chemical state and its correction are float64. Inherit the complete frozen source/input provenance, model parameters, GPU UUID and runtime binary. New coordinator, design, tests and launch/audit scripts are separately pinned; previous sources, runs, failures and reviews remain immutable.

Each retimed endpoint first receives eight independent NumPy-reference/GPU steps and an exact correction-aware GPU replay. Prefix verification never advances the scientific start. Require unchanged inherited prefix tolerances. During production check every-step volume error <5%, equivalent cell radius ≥4 voxels, zero clipping, positive finite chemistry, dilution amount error ≤2e-14 and unchanged target volume. Require sixteen cells, fifteen completed events, no active furrow and unchanged ancestry/cycle/RNG metadata. Regular observations require boundary occupancy <0.01 and one connected component per cell.

Observe every **0.15** units (**1,201 observations/path**), record conservative networks every **0.015** (**12,001 graphs/path**) and save complete checkpoints at maintenance elapsed **0,36,72,108,144,180**. Require exact eight-step restart replays at elapsed **36 and108**. Compensation is included in all checkpoint/restart states. Sparse full checkpoints constrain storage; dense matrix diagnostics retain network transitions.

Retain all previously registered numerical gates: regular distribution/shape/volume/amount/polarity comparisons, native density errors ≤1%, regular/dense relative spectral errors ≤1%, homologous log concentration error ≤0.01, absolute polarity-component error ≤0.01, relative cell-volume error ≤0.005 and relative conservative-matrix norm error ≤0.01. Require exact cell/ancestry correspondence and no unmatched observations. Record every contact-support mismatch. The maximum difference in frozen-uniform nonzero-mode growth must remain ≤0.001. These growth rates describe chemistry on the instantaneously frozen network, not stability of the full moving patterned system.

## Prespecified maintenance measurements

Contrast is the across-cell SD of **natural-log activator**. Use primary threshold **0.1**, with **0.05 and0.2** sensitivity thresholds. Every accepted start is above all three thresholds. These are screens for chemical differences, not definitions of a biological cell identity.

For each threshold, record:

1. **Continuous retention:** every saved observation over elapsed0–180 remains strictly above threshold. Equality counts as loss. Record the minimum and final contrast.
2. **Loss:** the first saved observation at or below threshold, plus all sampled downward/upward crossings. A short dip is recorded even if contrast later recovers.
3. **Sustained loss:** the first complete24-unit interval during which every saved observation is at or below threshold. Both endpoints count, so this requires161 regular observations. Report the interval start only after the entire interval is observed.
4. **Sustained recovery:** after a sustained-loss interval, the first complete24-unit interval strictly above threshold. A short recovery or a recovery too late to fill24 units is reported among crossings, without a sustained-recovery claim.
5. **Late persistence:** every observation in maintenance elapsed156–180 remains strictly above threshold. This is total elapsed216–240 since t18, physical t234–258. Late persistence after recovery is distinguished from continuous retention.

Both clocks must agree at every threshold on continuous retention, late persistence, sustained loss and sustained recovery. Corresponding first loss, sustained-loss and sustained-recovery start times, when observed, must differ by ≤0.3 units. Report all crossings separately; their full transient counts are descriptive, not an additional acceptance gate. A reproducible loss is a valid scientific outcome. Do not require persistence as a condition for numerical acceptance.

Measure homologous chemistry, polarity and volumes, tissue axis ratio, summed undirected conductance, changing edge support, graph spectrum and supported frozen-uniform growing-mode counts. Conductance changes are relative to each branch's own t78 start; later review may combine the separately frozen t18–78 observations with this continuation while retaining stage boundaries. Do not replace the original formation data or silently relabel maintenance time as development from a new zygote.

## Execution, stopping and interpretation

Execute serial coarse/fine pairs on the dedicated GTX1080Ti, with one CPU numerical thread. Compare a complete pair before advancing to another context. The preceding60-unit study took10.08hours; allow provisionally **28–32hours** for this three-times-longer window, conditional on unchanged speed and all gates passing. Re-estimate after the first completed pair. Require at least100GiB free before launch; the sparse-checkpoint design is expected to use roughly70–80GiB, including complete starts, prefixes, restarts and dense diagnostics. Large outputs stay on the existing data-drive output link.

Stop on any actual-start prefix, physical, exact restart or full-horizon numerical comparison failure. Preserve failure states and unexecuted cases. No automatic timestep fallback, tolerance relaxation, replacement history, rerun or downstream endpoint/identity assay. Completed execution remains **provisional until an independent saved-array review** qualifies this longer horizon.

The study tests long persistence under a continuously changing mechanical context. It does not establish autonomous biochemical memory, commitment, inherited biological identity, universal convergence or large-scale chemically caused shape organization. If qualified, frozen-endpoint coexistence and controlled context/state-response assays become the next separately registered tests. README science claims, manuscript and results ledger remain unchanged by launch; dated runtime status belongs in the lab notebook.
