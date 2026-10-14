# Mature timestep qualification after fresh development

Registered 12 October 2026, before running. This is a numerical test using the twelve independently accepted developmental preparations described in [the control-completion assessment](fresh_zygote_control_completion_assessment.md). Histories 1101, 1102 and 1103 each supply four actual t=18 endpoints: feedback off, tension and adhesion, polarity alone, and full feedback. All primary endpoints were developed at dt=0.000234375, on grid120 with extent2.24, and contain sixteen cells with fifteen completed cleavage events.

## Question and scope

Can these actual endpoints be continued accurately with a larger mature timestep? The test compares dt=0.00375 with dt=0.001875 over six elapsed model-time units, from physical t=18 to24. There are 24 new numerical continuations and twelve paired comparisons, within three existing independent developmental histories. These are not 24 biological replicates or new histories.

This short screen qualifies candidate clocks for subsequent, separately registered 60- and 240-unit moving experiments. It cannot establish long-time convergence, differentiation, pattern persistence or a general developmental mechanism. Completion remains provisional until independent review of saved evidence.

## State inheritance and solver

Use the existing compensated resident `DevelopmentGpu` solver: PyTorch arrays and conservative chemical transport, with custom CUDA mechanics and geometry/polarity. No solver or kernel source changes are included. The parent protocol, inherited sources, reference evidence and runtime binary are hashed.

Start from each complete saved `gpu-t18.npz`. Preserve every saved array byte, including visible phase fields, phase carry, rounding counts, chemical compensation, chemical states, polarity, targets, IDs, parents and cell-cycle due times. Preserve physical time, lineage, cleavage audits, all three random states and all other metadata. Only Config dt, integer step horizon, observation stride and metadata step_number change. The horizon is physical t=24; observations remain every0.15 units. No chemical reset, supplied pattern, backend handoff, feedback switch or further division is introduced.

The four mechanical branches retain their developmental coefficients:

| Branch | Chemical tension | Chemical adhesion | Polarity tension |
|---|---:|---:|---:|
| off | 0 | 0 | 0 |
| tension_adhesion | 0.25 | 0.35 | 0 |
| polarity | 0 | 0 | 0.35 |
| full | 0.25 | 0.35 | 0.35 |

These coefficients multiply the model's cell-state response or polarity response; they are not physical material measurements. The off branch retains passive mechanics.

## Registered checks

Each actual endpoint and clock first receives eight independent NumPy-reference versus GPU steps and an exact correction-aware GPU replay, at the unchanged inherited prefix limits. During continuation, apply the inherited every-step physical gates: volume error below5%, resolved radius at least four voxels, zero phase clipping and dilution amount error at most2e-14. Require positive finite chemistry, sixteen cells, fifteen unchanged completed cleavage audits, no active furrow, inherited targets, unchanged lineage/cell-cycle/random states, boundary occupancy below1% and connected tissue at regular observations.

Save dense graph diagnostics every0.015 units (401 records per path), regular observables every0.15 units (41 records), and full checkpoints at elapsed0,3 and6. At elapsed3, require exact eight-step checkpoint replay with all compensation fields restored.

Apply the inherited distribution, shape, amount, polarity and graph-spectrum tolerances without relaxation. Require both regular and dense relative spectral differences at most1%, and each occupancy, chemical-density and polarity-density checkpoint error at most1%. All comparisons use native grid120 and existing density definitions. Report any edge-topology disagreement explicitly.

Because mature cells have fixed ancestry and number, also require exact cell correspondence throughout. Additional homologous gates are maximum absolute log-concentration difference0.01, absolute polarity-component difference0.01, per-cell relative volume difference0.005 and relative transport-matrix norm difference0.01. Cell-count and ancestry mismatch are forbidden, rather than using the developmental allowance for asynchronous cleavage.

## Stop and review policy

Run serially on the dedicated GTX1080Ti, with one CPU numerical thread. Stop on any prefix, physical, restart or paired-comparison failure. Retain failed checkpoints and leave later paths unexecuted. Do not automatically substitute dt=0.0009375, change thresholds, rerun a case or extend to long science. A failure requires diagnosis and a separately registered refinement.

Outputs go to `output/studies/fresh-zygote-mature-qualification`, through the existing data-drive relocation alias. The original failed fresh-zygote pilot remains failed. Qualification acceptance requires an independent saved-array review after completion; the coordinator cannot accept its own result. README claims, ledger claims and manuscript evidence are not revised by this launch.
