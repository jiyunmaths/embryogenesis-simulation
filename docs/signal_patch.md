# Controlled chemical-patch response and mechanical performance

This experiment asks whether the existing mechanics can turn a maintained chemical pattern into a substantial deformation, and whether that deformation remains after the chemical clamp is removed. It is a **forced-response diagnostic**, not spontaneous signaling or shape symmetry breaking. No target shape, additional mechanical force, growth, or imposed cell displacement is introduced.

The implementation is [embryo/signal_patch.py](../embryo/signal_patch.py). Its numerical and intervention tests are in [tests/test_signal_patch.py](../tests/test_signal_patch.py). The separate [long-time shape experiment](shape_persistence.md#long-time-continuation) keeps autonomous signaling and the original mechanical parameters.

## Equation-preserving kernel optimization

The mechanical kernel now evaluates either the constant-tension Laplacian or the variable-tension flux divergence, instead of computing and discarding the former in polarized cells. Current-step occupancy and volume arrays are reused when calculating centers and volume forces; the contact calculation also reuses its occupancy array. No state-dependent array is cached between steps. Update order, float precision, boundary conditions, mechanical forces, time step, and clipping policy are unchanged.

The frozen [pre-optimization reference](../tests/reference_mechanics.py) tests exact full trajectories through cytokinesis, with and without feedback and polarity, zero directional tension, and stronger fate-dependent tension. A mature grid-40, 16-cell checkpoint also produced an exact 30-step state match. In three alternating single-process timing repeats with BLAS/OpenMP limited to one thread, the median time for 30 steps fell from **2.694 s to 2.016 s**: **1.337 times the throughput**, or **25.2% less runtime**. This is a measured local workload, not an extrapolation to GPU execution or larger populations. The saved timing record is `outputs/kernel-benchmark/analysis.json`.

Reproduce the benchmark on an otherwise idle machine:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/mechanics.py \
  --checkpoint outputs/shape-persistence/full/final_state.npz \
  --output outputs/kernel-repeat
```

The optimized production kernel is also used by the live dashboard after restarting an already running server. This change does not add a new implicit solver or a GPU backend.

## Prescribed chemical input

All five branches start from the same full-feedback checkpoint at $t=90$, preserving phase fields, volumes, cell IDs, fates, polarity, and RNG state. Only the chemical input and the declared mechanical coefficients change. The tissue already has a cleavage history and differentiated state, so this experiment measures the response of a mature aggregate, not development from an unbiased zygote.

A smooth patch is assigned using the initial cell centers $\mathbf x_i^0$, fixed target volumes $V_i^*$, and an explicitly imposed unit vector $\mathbf e$. The default vector is laboratory $+z$. Define

$$
\bar{\mathbf x}^{\,*}=\frac{\sum_i V_i^*\mathbf x_i^0}{\sum_i V_i^*},
\qquad
b_i=\tanh\!\left(\frac{\mathbf e\cdot(\mathbf x_i^0-\bar{\mathbf x}^{\,*})}{w}\right),
\qquad
\bar b=\frac{\sum_i V_i^*b_i}{\sum_i V_i^*}.
$$

With $w=0.35$ and $A=0.8$, the imposed regulators are

$$
a_i^{\mathrm{clamp}}=1+A\frac{b_i-\bar b}{\max_j|b_j-\bar b|},
\qquad h_i^{\mathrm{clamp}}=1.
$$

The target-volume-weighted activator mean is exactly one, and all activities are positive. This is a polar, smooth material patch: values remain attached to the original cell IDs rather than being recalculated as cells move. Cleavage must already be complete. There is no material-axis claim for the spontaneously evolving model; the imposed laboratory direction belongs only to this diagnostic.

From $t=90$ through $t=120$, each regulator update is replaced by this external clamp. Fate, polarity, and mechanics continue evolving with the existing equations. The clamp represents external reservoirs; it is not conservative molecular transport. Graph spectra saved during this phase describe the hypothetical autonomous Gierer–Meinhardt dynamics, which are bypassed during forcing.

At $t=120$, the clamp ends. The original Gierer–Meinhardt reaction and contact-graph exchange resume from the clamped state, without resetting fate, polarity, or geometry, through $t=150$. Thus release removes the external input; it does not disable the internal feedback or erase the resulting chemical pattern. A persistent response after release would still require controls to determine what maintains it.

## Five matched branches

| Case | Chemical clamp | Fate–tension coefficient | Mechanical feedback |
|---|---|---:|---|
| `patch_default` | Spatial patch | 0.25 | Enabled |
| `uniform_default` | $a_i=h_i=1$ | 0.25 | Enabled |
| `patch_strong` | Same spatial patch | 0.75 | Enabled |
| `uniform_strong` | $a_i=h_i=1$ | 0.75 | Enabled |
| `patch_no_feedback` | Same spatial patch | 0.25, inactive | Disabled |

The stronger condition changes only `fate_tension`, in the existing scalar relation

$$
\gamma_i=\gamma_0\left[1+\alpha_\gamma\tanh(f_i)\right].
$$

It does not strengthen the polarity coefficient or pairwise fate adhesion. Both values keep this scalar tension positive. Each patch has a uniform control with the same coupling strength and target-weighted chemical means. Inherited fates are retained even in the uniform controls; chemical uniformity is not an artificial reset to mechanically identical cells. In the feedback-disabled case, fates and polarity remain dynamical readouts but cannot influence mechanics.

This tests the response of the existing force laws; it does not implement planar polarity, explicit active stress, apical constriction of an epithelial shell, or tissue growth.

## Measurements and interpretation

Record every 0.6 time units, retaining the time step 0.015 and spatial grid $40^3$. In addition to the original shape tensor, volume, boundary, connectivity, and clipping diagnostics, measure the ratio along the imposed direction:

$$
R_{\parallel}=\sqrt{\frac{\mathbf e^{\mathsf T}C\mathbf e}
{[\operatorname{tr}(C)-\mathbf e^{\mathsf T}C\mathbf e]/2}}.
$$

This detects elongation or compression along the forcing direction even if the aggregate's longest axis remains elsewhere. At the end of the clamp and release, report patch-minus-control differences in the ordinary axis ratio and $R_{\parallel}$, relative covariance-tensor difference, and angle between measured axes. Differences retain their sign. A change in tensor orientation can matter even without an increased longest-to-shortest ratio.

The comparison chart displays boundary occupancy (maximum diffuse-union occupancy on the outer grid layer) and minimum cell radius in grid spacings alongside the shape and volume curves, with the original thresholds visible. Both the late forced window and the late release window are summarized. Numerical quality is additionally screened across **all sampled times**, not just the final state. The original radius-resolution criterion remains unchanged even if an interesting response appears. These are sampled diagnostics; software correctness, finite-window morphological persistence, mechanical response to external forcing, numerical convergence, and spontaneous symmetry breaking are separate claims.

The pilot is one checkpoint and one imposed direction. Rotation, grid/time refinement, larger-domain checks, varied input width/amplitude, and independent developmental histories remain necessary before generalizing a response.

## Reproduce

Use new output directories; existing results are not overwritten. The checkpoint is produced by the [shape-persistence protocol](shape_persistence.md#reproduce).

```bash
for case in patch_default uniform_default patch_strong uniform_strong patch_no_feedback; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.signal_patch run \
    --checkpoint outputs/shape-persistence/full/final_state.npz \
    --case "$case" --output "outputs/patch-repeat/$case" \
    --clamp-duration 30 --release-duration 30
 done

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.signal_patch compare \
  --runs outputs/patch-repeat/patch_default outputs/patch-repeat/uniform_default \
    outputs/patch-repeat/patch_strong outputs/patch-repeat/uniform_strong \
    outputs/patch-repeat/patch_no_feedback \
  --output outputs/patch-repeat/comparison
```

Independent cases may run in separate CPU processes; keep BLAS/OpenMP thread counts low to avoid oversubscription. Each case writes `protocol.json`, `config.json`, `history.json`, `analysis.json`, `trajectory.json`, `viewer.html`, `lineage.json`, and a final checkpoint. Playback shows the actual simulated surfaces. The comparison produces `response.png`, `clamp_shapes.png`, `final_shapes.png`, and `analysis.json`, with matched camera, spatial limits, and fate-color scale. The imposed chemical axis is labeled separately from measured and first-cleavage axes.

## Completed pilot results

All five branches completed the $t=90$–120 clamp and the autonomous $t=120$–150 release. The initial aggregate axis ratio was 1.322977 in every branch. The imposed chemical patch changes identities and the existing mechanics, but these runs do not develop a conspicuous new geometry.

| Case | Axis ratio at end of clamp | Axis ratio after release | First sampled boundary-screen failure |
|---|---:|---:|---:|
| `patch_default` | 1.328943 | 1.334597 | 128.4 |
| `uniform_default` | 1.326382 | 1.328050 | 130.2 |
| `patch_strong` | 1.329074 | 1.335815 | 118.2 |
| `uniform_strong` | 1.323110 | 1.322361 | 121.2 |
| `patch_no_feedback` | 1.328072 | 1.331981 | None through 150 |

The signed patch-minus-uniform differences are:

| Fate–tension coefficient | End of clamp | After release |
|---|---:|---:|
| 0.25 | +0.002562 | +0.006547 |
| 0.75 | +0.005964 | +0.013453 |

Increasing the coefficient makes the patch/uniform difference larger, but the patched aggregates themselves remain very similar between strengths. The larger contrast is partly due to the stronger uniform control becoming less elongated. Neither stronger coupling nor the longer unforced run establishes large shape symmetry breaking in this regime. An axis ratio or a persistent post-release response alone does not establish that the activator–inhibitor equations select an axis spontaneously.

**No branch passes the whole-run numerical screen.** Every run retains the existing minimum-cell-radius resolution failure. All four feedback-enabled patch/uniform branches additionally cross the boundary threshold; the stronger patch crosses it before the clamp ends. Volume error, clipping, individual-cell connectivity, and graph-connectivity screens pass in this pilot. Boundary occupancy is a local edge diagnostic, not lost tissue volume, and requires a larger-domain comparison before interpreting the late differences. These limitations remain in the JSON reports and are visible in the charts.

The next validation should increase domain extent while retaining spatial resolution, then refine space and time. Extending time again or increasing cell count on the same grid would not resolve the identified limitations. Repeat the patch test on the numerically verified domain before deciding whether the existing force law needs an additional active mechanical mechanism.

Results are saved under `outputs/signal-patch/`. Each case has an offline 3D playback covering both phases. The `comparison` folder contains `response.png`, `clamp_shapes.png`, `final_shapes.png`, and the full signed tensor/axis comparisons in `analysis.json`. The final figures use shared camera, spatial limits, and fate colors; the imposed chemical direction is labeled separately from the measured and first-cleavage axes.
