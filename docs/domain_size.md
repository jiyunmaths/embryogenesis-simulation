# Fixed-spacing domain-size check

The new conservative model is screened on successively larger domains before interpreting long-time shape. The earlier $t=124.2$ boundary-screen crossing occurred in a **historical normalized-signaling** run; it is not assumed to be the onset time in the new model.

## Larger-domain dashboard preset

```bash
python -m embryo.dashboard --config configs/large_domain.json
```

Use an available port with `--port` if another dashboard is running. This starts a new zygote; it does not resize or overwrite another live session. In an existing dashboard, pause, set **Grid per axis = 56**, **Domain half-width = 2.24**, and **Total steps = 2000**, then Reset. Retain `dt=0.015`, interface width `0.085`, conservative transport, and diffusivities `0.02, 0.4` for the preset experiment. Reset starts the new experiment from a zygote.

| Domain | Grid per axis | Half-width | Voxel spacing | Dense voxel count relative to original |
|---|---:|---:|---:|---:|
| Original | 40 | 1.60 | 0.08 | 1.000 |
| Intermediate | 48 | 1.92 | 0.08 | 1.728 |
| Larger preset | 56 | 2.24 | 0.08 | 2.744 |

The larger preset adds 0.64 model-length units on each side, or eight voxel layers. Interface width, cell size, time step, and diffusion coefficients are unchanged. Increasing only `extent` at fixed grid would coarsen the cells, conflating domain and resolution effects. The additional voxel count increases memory and mechanics work; it does not increase the number of biological cells.

The compact library defaults remain available for quick experiments. The larger-domain settings are explicit in `configs/large_domain.json`; pressing the dashboard's Defaults button restores the compact defaults.

## Controlled protocol

```bash
python -m embryo.domain prepare --output outputs/domain-repeat --until 30
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-repeat --grid 40
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-repeat --grid 48
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-repeat --grid 56
python -m embryo.domain compare --output outputs/domain-repeat
```

The study uses one seed-7 conservative zygote, $D_a=0.02$, $D_b=0.4$, $\Delta t=0.015$, and the original mechanical/cell-cycle parameters. The physical simulation window is $0\le t\le30$, with snapshots every 0.6 time units. Preparation saves a common full checkpoint and its SHA-256 hash before any branch advances. All output directories reject overwrites.

Enlargement pads the phase fields with zero-valued exterior space, placing the original arrays at the center without interpolation. Even grid increments preserve voxel alignment. Existing concentrations, volumes to roundoff, target volumes, cell IDs, fate, polarity, lineage, division state, and random streams are retained. The shared coordinate centers and phase-field entries are preserved exactly in memory. Reflecting mechanics boundaries move outward; the newly introduced space starts empty.

The common source has boundary occupancy approximately $3.01\times10^{-11}$. Padding does not extrapolate omitted diffuse tails, and cannot undo previous boundary effects. For that reason the primary comparison starts from the zygote rather than padding a late, already wall-influenced embryo. `prepare --checkpoint PATH` is available for explicitly labeled continuation studies; its report retains source time, boundary occupancy, and original transport model.

Each branch records actual cell meshes, volume-weighted spectral diagnostics, boundary occupancy, axis ratio, signaling contrast, and an approximate distance from the nearest occupied $\phi\ge0.5$ voxel center to the wall. The latter is a voxel-based clearance estimate, not a reconstructed membrane distance.

Predeclared screens are:

- Maximum sampled boundary occupancy below 0.01 in the largest domain.
- Maximum sampled relative axis-ratio difference below 1% between the two largest domains.

The comparison also reports the final cell-aligned phase-field L2 difference on the common box and activator differences. Shape agreement alone cannot establish convergence of signaling or the full cell geometry. Cell IDs and parents must align before reporting field comparisons. Shared-box field error excludes the exterior added region and is labeled accordingly.

## Completed conservative-model screen

All three seed-7 branches reached $t=30$ with 16 cells. Results are saved in `outputs/domain-conservative`, including `comparison.json`, `comparison.png`, and each branch's mesh playback and full checkpoint.

| Grid | Maximum sampled boundary occupancy | Final axis ratio | Minimum sampled voxel-center clearance |
|---:|---:|---:|---:|
| 40³ | 7.3886e-4 | 1.30005257 | 0.36 |
| 48³ | 1.0705e-7 | 1.30004388 | 0.68 |
| 56³ | 1.9049e-11 | 1.30004406 | 1.00 |

All boxes, including the original, pass the 0.01 wall-occupancy screen over this window. Thus this experiment demonstrates additional clearance and small domain sensitivity, rather than correction of an observed boundary-screen failure in the new model.

Between grids 48 and 56, the maximum sampled relative axis-ratio difference is **1.4086e-7 (0.0000141%)**. Their final common-box phase-field relative L2 difference is **9.4760e-6 (0.0009476%)**, and the maximum absolute activator difference is **5.0671e-8**. Current cell IDs and parents align. For comparison, the original/intermediate pair has a final common-box field difference of 0.0695%. Both predeclared largest-domain/largest-pair screens pass.

The final activator standard deviation is approximately 0.008943 in every run; the early shape is not evidence of a mature signaling pattern or feedback-specific axis. Minimum cell radii remain approximately **3.969 grid spacings**, below the pre-existing four-spacing screen. This is unchanged cell resolution, not a new failure introduced by enlargement.

Use the **56³ / half-width 2.24** preset for the next conservative experiment while retaining boundary monitoring. The measured agreement applies to this seed and time window; the longer continuation still needs its own domain comparison. The enlargement and export workflow passes seven focused tests, with 31 domain/backend tests and 11 dashboard-client tests passing in the combined checks.

## Interpretation limits

This is a single-seed, finite-window domain screen. A low wall occupancy is evidence of clearance, not proof that all boundary effects vanish. The initial $t=30$ result alone did not certify clearance through $t=180$. The completed continuation below repeats those diagnostics and compares the larger domains over the longer window. Voxel spacing and smallest-cell resolution are unchanged, so the earlier spatial-resolution concern is not solved by adding empty space. Spatial/time refinement, independent seeds, and feedback controls remain separate acceptance tests.

## Long-time continuation protocol

Continue the completed domains independently, preserving their distinct $t=30$ states and original common-zygote ancestry:

```bash
python -m embryo.domain extend --source outputs/domain-conservative \
  --output outputs/domain-conservative-extended --until 180
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-conservative-extended --grid 40
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-conservative-extended --grid 48
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.domain run --output outputs/domain-conservative-extended --grid 56
python -m embryo.domain compare --output outputs/domain-conservative-extended
```

`extend` prepares a new study directory; each `run` performs the corresponding continuation. Source data remain unchanged. Preparation records each source checkpoint's SHA-256 hash; the runner verifies it before loading. Each domain retains its original grid, half-width, voxel spacing, diffusivities, mechanics, seed state, and time step. The full $t=0$–180 history and playback are joined without a duplicate $t=30$ frame. Reported runtime covers the new continuation only.

Progress histories are written at each sampled time and `progress_state.npz` is updated every ten samples (six model-time units), plus at completion. The final checkpoint and complete mesh playback are exported at the end. A progress checkpoint preserves the solver state; the complete in-memory playback is only written on completion.

Apply the same predeclared wall-occupancy and largest-pair shape screens over the full extended window. Record the first sampled boundary-threshold crossing, final shared-box field discrepancy, and signaling differences. Do not infer boundary independence from clearance alone, or a feedback-specific axis from elongation alone. This extension does not address the unchanged smallest-cell-resolution shortfall.

### Automatic completion and live status

The current long run writes `outputs/domain-conservative-extended/status.json`. A separate monitor waits for all three final analyses, then generates `comparison.json`, `comparison.png`, and `RESULTS.md`. The status becomes `completed` only after these artifacts have been generated; it becomes `failed` if a watched worker exits without final results or the comparison fails. Until then, partial progress is not a completed domain-convergence result.

```bash
python -m embryo.domain_monitor --output outputs/domain-conservative-extended
```

Optional `--workers GRID:PID ...` enables worker-exit detection for already-running jobs. The monitor does not alter the solver. Source files, checkpoints, and scientific acceptance thresholds remain unchanged.

## Completed continuation through time 180

The conservative domain study has completed in all three boxes. Its independent report is `outputs/domain-conservative-extended/RESULTS.md`.

| Grid | Peak boundary occupancy | Final axis ratio | First sampled boundary-threshold crossing |
|---:|---:|---:|---:|
| 40³ | 0.0218351 | 1.34159229 | 131.4 |
| 48³ | 1.64699e-6 | 1.33848833 | None |
| 56³ | 2.82976e-10 | 1.33849733 | None |

Both declared largest-domain/largest-pair screens pass. Between grids 48 and 56, the maximum relative axis-ratio discrepancy is **6.726e-6 (0.000673%)**, and final shared-box phase-field L2 discrepancy is **6.085e-5 (0.006085%)**. Their final activators differ by at most 1.87e-8. These results support the larger domain through time 180 for this seed and parameter set. The original box fails the boundary screen; this need is now observed directly in the conservative model, rather than inferred from the historical normalized-model result.

The 56³ run's activator standard deviation peaks near **0.02934 at time 55.8**, then falls to **9.30e-8 at time 180**. Its **12 A / 4 B labels** persist. Thus persistent fate labels and geometric elongation coexist with decay of signaling contrast. This is not evidence of a persistent Turing pattern, irreversible biological commitment, or a feedback-specific shape axis. Minimum cell radius remains about 3.969 grid spacings; the domain test does not repair spatial resolution.
