# Exchange formation and response: timestep refinement

This study halves the mechanical/coupled timestep from **0.00375 to 0.001875**. It follows the completed moving-history replication: all 36 continuations passed quality checks, and all sixteen new-history donor/destination response comparisons favored the donor. With the earlier seed-7 reference, this is 24 comparisons within three developmental histories, not 24 independent replicas.

The next question is whether the exchange-generated organization and its responses depend on timestep. Both stages use identical physical starts at each timestep. The response repeats use the original t=210 endpoints; they do not start from the newly refined formation endpoints. This separates response integration sensitivity from accumulated changes during formation.

The driver was launched on **2026-10-03 at 16:25 EDT**. All four full GPU/native comparisons, ten context checks, and eighteen scientific continuations completed and passed. Seventy-three focused software checks passed before preparation; one opt-in GPU check was skipped in that CPU test invocation. The production gate supplies the required trajectory comparisons.

## Completed assessment

Both histories pass formation and same-state response refinement under every predeclared criterion. All four fresh-exchange target comparisons remain donor-nearer, and all sampled target/network recovery times are unchanged.

| Quantity | Maximum observed discrepancy | Limit |
|---|---:|---:|
| Formation chemical log concentration | 3.071e-4 (about 0.031% concentration) | 0.01 |
| Formation relative transport operator | 1.559e-4 (0.0156%) | 1% |
| Formation relative aggregate axis ratio | 2.591e-6 (0.000259%) | 1% |
| Formation endpoint phase field | 9.835e-5 | 0.02 |
| Normalized response waveform | 2.206e-4 (0.0221% of initial log pulse) | 1% |
| Relative target response integral | 6.156e-4 (0.0616%) | 2% |
| Target/network recovery time | 0 | 0.30 |

The root `outputs/exchange-response-histories-refined/refinement.json` retains every per-history comparison. This supports practical agreement under one timestep halving in the specified mature regime, not a convergence-order proof or full developmental/spatial acceptance. The next [neighbor-context assay](neighbor_context.md) isolates changes in surroundings while preserving the target's initial chemical state.

## Design frozen before execution

| Stage | Repeats | Physical time | Interpretation |
|---|---:|---|---|
| Full CPU/GPU acceptance at the smaller timestep | 4 | 210–270 | Existing accepted native CPU untouched/fresh controls and cell-27 −10% pulses; no new CPU trajectories |
| New-context CPU/GPU prefixes | 10 | 0.6 units each | Six t=150 formation contexts and four t=210 response contexts across seeds 8 and 9 |
| Exchange formation/retention | 6 | 150–210 | Untouched, fresh conservative exchange, and frozen pre-relaxed exchange on each history's original mature geometry |
| Same-state response refinement | 12 | 210–270 | Untouched/fresh backgrounds, each with its control and −10% activator pulses in both selected cells |

The response targets remain IDs 19/20 for seed 8 and 26/30 for seed 9. They were selected by the original initial minimum/maximum activator, before observing responses. All chemistry, geometry, polarity, IDs, lineage, age, and random streams are preserved by retiming. The integer step counter and configured step intervals are rescaled to preserve physical time. Prepared chemical interventions are unchanged. There is no further division, fate switch, A/B classification, or imposed identity count.

The GPU uses resident PyTorch arrays and matrix operations plus the existing custom CUDA mechanics/geometry/polarity kernels. Changing the timestep requires a new gate: the older 0.00375 acceptance is insufficient. The four complete replays reuse immutable accepted CPU trajectories from `outputs/cell-exchange-response-refined/`; a separate path adapter supplies their source/checkpoint layout without modifying them. The unchanged full-horizon validator and gated runner are reused. Its original precision, waveform, recovery, field, and conservation criteria remain fixed.

Only after all four replays pass does the driver validate the ten new starting contexts using matched 0.6-unit native/GPU prefixes and the same backend thresholds. A failed gate or prefix prevents all scientific continuations. One worker uses the dedicated GPU; native prefix checks use four CPU threads.

## Formation acceptance

Every matched sample uses the same cell IDs and physical time, at 0.15-unit intervals. Maximum discrepancies over the complete trajectory must meet:

| Quantity | Limit | Meaning |
|---|---:|---|
| Absolute log-concentration difference | 0.01 | Maximum proportional discrepancy over both species and every cell |
| Absolute polarity-component difference | 0.01 | Maximum vector component discrepancy |
| Relative aggregate axis-ratio difference | 1% | Overall shape elongation |
| Relative per-cell volume difference | 0.5% | Compartment capacity agreement |
| Relative contact-operator difference | 1% | Frobenius norm of the difference divided by the coarse operator norm |
| Endpoint phase-field difference | 0.02 | Maximum absolute occupancy-variable difference across all cell grids |
| Late destination/transferred distance-ratio difference | 0.05 | Change relative to initial exchanged-pair separation |

Contrast retention, informative-pair status, destination/transferred likeness, and reversed pair ordering must also agree with the coarse run. The original contrast threshold (log-activator SD >0.1), state-likeness ratio (<0.1), and final 24-unit window remain unchanged. These are predeclared practical agreement limits, not derived error bounds or a convergence-order estimate.

## Response acceptance

The selected response checks retain the existing criteria:

- Maximum discrepancy in pulse-minus-control log responses below 1% of the initial log pulse, over all sampled cells and both species.
- Target activator response-integral relative discrepancy below 2%.
- Target and network recovery classifications unchanged; recovery-time difference at most 0.30 model units.
- Maximum raw chemical log discrepancy at most 0.01.
- Donor/destination nearest-reference classification unchanged and informative at both timesteps; reference separation must exceed 0.01.

Each pulse is compared with its own moving control. Recovery requires staying within 10% of the initial displacement with at least 24 subsequent units observed. The selected pulses add/remove activator externally; they are not conservative exchanges. Original and refined diagnostic weights use the same measured starting volumes, while evolving volumes continue to govern transport and dilution.

All runs separately require maximum individual volume error <5%, radius ≥4 grid spacings, zero clipping, sampled boundary occupancy <0.01, positive finite chemistry/geometry, and dilution amount error ≤2e-14. Checkpoints are written every three units. Passing quality alone does not establish timestep agreement.

## Operation and outputs

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories_refinement prepare

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories_refinement run

python -m embryo.exchange_response_histories_refinement assess
```

Default output: `outputs/exchange-response-histories-refined/`. Preparation requires a fresh directory. `gpu-validation/` contains the four replay results and canonical CPU reference adapter; `backend-checks/` contains new-context comparisons. Each history has separate formation and response sources, per-job histories, and restart checkpoints. `formation_refinement.json` and response `refinement.json` report each history; the root `refinement.json` combines them.

Source and input hashes are frozen before execution, and the old studies remain unchanged. The root driver lock prevents concurrent drivers against one directory. Restart reuses verified passing gates/prefixes and completed jobs, restoring partial GPU checkpoints without repeating a pulse. A numerical failure stops the batch. A completed batch whose refinement comparisons fail is recorded as completed with `passed=false`, rather than silently accepted or replaced.

This is one timestep halving on two existing histories. Positive/pre-relaxed pulse responses, an end-to-end trajectory beginning with refined development, spatial convergence, identity inheritance, and general-geometry transport accuracy remain outside its scope. Similar fresh/pre-relaxed results do not establish a unique global attractor. [Full original moving protocol](cell_response_moving.md#replication-across-developmental-histories).
