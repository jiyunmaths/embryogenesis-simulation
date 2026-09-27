# Live simulation dashboard

Start from the repository with its Python dependencies installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard
```

Open <http://127.0.0.1:8765>. Optional arguments are `--port`, `--config path/to/config.json`, and `--benchmark path/to/completed/nonlinear-run`. The installed command `embryo-dashboard` starts the same server. Stop it with Ctrl+C. The Live embryo workspace uses the existing `Simulation.step()` implementation. The Continuum signaling workspace reads completed benchmark artifacts. Neither view introduces alternative dynamics or new runtime dependencies.

## Controls and parameters

1. Edit the configuration in the left panel. Run starts from the single zygote and applies edits made before starting.
2. Pause waits for the current numerical step to finish and records its state. Resume continues the same simulation, including its random generator states. Pausing does not change the trajectory.
3. To change a run after it has started, pause, edit, and Reset. Reset validates the new configuration before replacing the old simulation, then returns to a zygote. It clears recorded history. Reset is also available while running.
4. Defaults changes the form to the model defaults. It does not reset the solver until Reset is clicked, or Run is clicked from the initial ready state.

All `Config` fields are editable. Related numerical constraints still apply: interface width must resolve the grid, the time step must satisfy the mechanical and polarity bounds, and cytokinesis must span enough steps. Invalid settings produce an error without replacing the existing run. Large individual numerical steps may take time; if a control times out, the pause request remains active and the status changes when that step ends.

The controls and progress strip describe the current solver state. The four metric cards, geometry, and chart marker describe the selected recorded frame. The timeline only inspects history; it does not rewind the solver. Follow live selects the newest snapshot as snapshots arrive. Reloading the page reconnects to the same local server session; all tabs share one simulation.

## Reading the visualization

- Drag to orbit and scroll to zoom. Reset camera restores the initial orientation. The cutaway removes sampled surface points above a world-Z plane; it does not reconstruct a closed cross-section.
- Cell identity colors use the continuous regulatory state. The identity-count chart uses the configured thresholds. Neither implies irreversible commitment.
- Activator and inhibitor colors use a fixed scale around the homogeneous activity of one. Small fluctuations are not rescaled to fill the color range. Standard-deviation charts show their measured magnitude over simulated time.
- Polarity arrows show the model's apical–basal polarity vectors. Contact links show existing weighted graph edges, not physical filaments. Surface dots are a rendering approximation and do not affect mechanics.
- Numerical diagnostics flag volume error, boundary occupancy, poorly resolved cells, delayed divisions, and phase-field clipping. These help identify artifacts; absence of a warning does not establish numerical convergence or biological validity.
- Frozen-graph linear stability is computed before a run on reference graphs and during the run on actual contacts. The growing-mode count applies near the homogeneous equilibrium on a frozen graph; it does not by itself predict nonlinear development on a dividing embryo.

The default run lasts `1000 × 0.015 = 15` dimensionless time units. Our [timescale study](timescales.md) found that appreciable signaling growth can require longer. Weak early variation is therefore meaningful and should not be hidden by the visualization. For a longer observation at the same time step, set Total steps to 4000 (time 60); keep the seed and other parameters unchanged for a direct comparison.

## Continuum signaling: nonlinear benchmark playback

Select **Continuum signaling** at the top of the dashboard to inspect the [nonlinear persistence experiment](nonlinear_bridge.md). The server first looks for `outputs/nonlinear-bridge-verified`, then `outputs/nonlinear-bridge`, relative to its working directory. To choose a different completed experiment:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard \
  --benchmark outputs/nonlinear-bridge-verified
```

Restart an already-running dashboard server to load the new backend, then reload the page. With no readable benchmark, this workspace displays the command needed to generate one; Live embryo remains usable. Source files are `analysis.json` and `fields.npz`. Failed scientific acceptance checks are shown as failures, not hidden. Incomplete runs or invalid concentration arrays cannot be played.

- **Recorded run** selects one of the spatial refinements, finer time steps, or equal-diffusivity control. Switching runs preserves the nearest saved physical time.
- **Play / Pause**, **Reset playback**, and the timeline inspect saved concentration snapshots. Playback advances one recorded frame every 0.9 seconds, regardless of the simulated time gap. It does not interpolate concentrations, integrate equations, modify parameters, or reset the embryo solver.
- **Activator / Inhibitor** selects the chemical. Each chemical has one fixed color range across every saved run and time, including the stable control. Early fluctuations are not stretched to appear as a large pattern.
- **3D cutaway** displays the exposed faces of the actual rectangular control volumes using their saved concentrations and graded coordinates. Drag to orbit and scroll to zoom. The slice slider removes compartments above a plane along the selected X, Y, or Z axis.
- **2D cross-section** shows one slab of control-volume averages at its actual unequal spatial widths. Labels give the slab bounds. The unoccupied notch is left empty.
- The history chart shows volume-weighted activator and inhibitor standard deviations for the selected run plus the equal-diffusivity activator control. The marker follows playback; shading identifies the late measurement window.
- The refinement chart switches between spatial and temporal comparisons. It plots the saved **two-species relative differences as percentages**, using conservative averaging for spatial comparison. Undefined relative errors are omitted, not plotted as zero.
- The late-change metric and acceptance checks describe the **whole recorded experiment**, even when the timeline shows an early snapshot. Other metric cards describe the selected run and time. A snapshot absent from scalar history is marked “Not sampled”.

This is a fixed-domain signaling experiment. Its numerical compartments are not biological cells; its colors are concentrations, not identities. It does not show emergent embryo geometry. The explicit separation from Live embryo prevents benchmark playback from being confused with a developing tissue. The header continues to report live embryo status while inspecting the benchmark; switching workspaces does not pause that solver. Returning to Live embryo stops benchmark playback.

Benchmark parameters are recorded, not editable in this view. Use `python -m embryo.nonlinear --help` to generate another experiment and point `--benchmark` at it. Restart after replacing already-loaded files; the server caches their metadata and the latest requested run. Benchmark data is separate from the live embryo's retained-frame download.

## History, downloads, and resource limits

The dashboard records the initial state, every `save_every` steps, and at pause/completion. Display requests never advance the solver or consume random numbers. Retained history is limited to 400 frames or approximately 16 MiB of encoded JSON, retaining at least the newest frame. Earlier frames are pruned, and the timeline/download indicate that history may be incomplete. Actual Python and browser memory use is larger than the encoded JSON size.

Download recorded frames saves the retained surfaces, metrics, contact graphs, parameters, and history-start revision. This is visualization data, not full cell fields or an exact-continuation checkpoint. Closing a browser tab does not stop a running solver; stopping the server discards its in-memory session. Use the batch CLI when you need permanent full experiment outputs and restart checkpoints.

Interactive limits are `grid <= 128`, `max_cells <= 128`, and `grid³ × max_cells <= 16,777,216`. Additional bounds limit excessive signal substeps. These prevent accidental extreme allocations; working arrays and temporaries can still consume substantial memory. Increasing resolution or cell count does not automatically preserve numerical accuracy: use the project's refinement and validation protocols.

## Local API

The server binds to IPv4 loopback and serves only dashboard assets and these JSON routes:

| Route | Operation |
|---|---|
| `GET /api/schema` | Defaults, parameter types, choices, and interactive limits |
| `GET /api/state?after=revision` | Current status and newer retained frames |
| `GET /api/preflight` | Reference graph stability and initial contact graph |
| `GET /api/benchmark` | Recorded nonlinear experiment metadata, histories, comparisons, and checks |
| `GET /api/benchmark/run?id=name` | Exact saved fields and mesh for one allowlisted recorded run |
| `POST /api/run` | Start/resume; optional `{"config": {...}}` |
| `POST /api/pause` | Stop at a step boundary; body `{}` |
| `POST /api/reset` | New zygote with `{"config": {...}}` |

Each reset increments `generation`; frame revisions are only meaningful within that generation. POST requests require `Content-Type: application/json`. Cross-origin and nonlocal Host requests are rejected. The dashboard is designed for one local interactive experiment, not as a shared remote service.

## Development checks

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q
# Optional Node.js 18+ check of controls, resets, history, and network races:
node --test tests/test_dashboard_client.cjs tests/test_benchmark_client.cjs
```

The JavaScript tests use a minimal DOM fixture and the real Python parameter schema. They verify client behavior, including playback, run-switching races, fixed color ranges, and unavailable-data recovery; they do not verify browser layout. Node.js is only needed for these tests, not for running the dashboard.


## Mechanical performance and controlled shape experiments

The dashboard uses the optimized production mechanics kernel; restart an already running server to load the change. The [kernel benchmark](signal_patch.md#equation-preserving-kernel-optimization) verifies the same tested trajectories and measured lower runtime, without changing the time step or force laws.

The [long-time shape continuation](shape_persistence.md#long-time-continuation) and [controlled signal-patch/release protocol](signal_patch.md) are reproducible command-line experiments with per-run offline `viewer.html` playbacks and comparison figures. They are not additional live-dashboard modes. Their chemical clamp is explicitly external forcing and should not be interpreted as spontaneous pattern formation.

## Closed cell surfaces

The live embryo view now defaults to **Show → Cell surfaces**. It draws opaque, lit triangle meshes from each cell's actual $\phi=0.5$ boundary. Orbit, zoom, identity/signaling/lineage colors, polarity arrows, contact-graph overlays, and the world-Z cutaway remain available. **Show → Point samples** restores the previous display. Restart an existing dashboard server and refresh the browser to load the new backend and renderer.

Meshes use full-resolution [Lewiner marching cubes](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes), cropped around each occupied cell for extraction efficiency. The added `scikit-image` dependency is installed by `python -m pip install -e .`. No convex hull, sphere fit, geometric smoothing, or artificial boundary cap changes the shape. Lighting normals average the adjoining triangle faces; the vertices themselves are not moved. Enclosed interfaces have closed triangle meshes; a cell that intersects the computational boundary is reported as open. Moving the cutaway slider intentionally opens the rendered surface. Coarse-grid faceting can remain visible and should not be mistaken for additional physical resolution.

WebGL supplies depth occlusion and interpolated lighting. A shaded-triangle Canvas fallback is available when WebGL cannot initialize. Contact and polarity lines are illustrative overlays. These changes affect visualization only: solver forces, states, metrics, and random streams are unchanged. The actual 16-cell saved-state example contains **8,448 triangles**, with all 16 meshes closed; one measured extraction took approximately **0.06 seconds**. Meshes are generated only when a frame is recorded and reused while orbiting. Larger frame payloads remain subject to the dashboard's existing history-memory limit.

New snapshots and offline viewers contain a `mesh` object per cell (`vertices`, `faces`, outward `normals`, closure status and boundary-edge count), as well as the legacy `points`. Saved viewers inline the shared renderer and need no CDN or server. Older point-only frames continue to display as points: their original phase fields cannot be recovered faithfully from sparse samples.

An existing full-state checkpoint can be inspected as a genuine surface snapshot:

```bash
python -m embryo.surface \
  --checkpoint outputs/signal-patch/patch_strong/final_state.npz \
  --output outputs/cell-surface-repeat
```

Open the resulting `viewer.html`. This exports one saved state without advancing the solver or inventing missing earlier frames. The completed example is `outputs/cell-surface-preview/viewer.html`. Future live runs and normal CLI exports record surface meshes throughout their trajectories.

## Conservative signaling default

Restart the server and select Defaults/Reset to start the conservative live model. `signal_transport` selects `conservative` (new default) or `random_walk` (historical control). In conservative mode, diffusivities have model length²/time units; defaults are 0.02 and 0.4. They are not interchangeable with the historical rates 1 and 20. Preflight uses fixed physical unit-box meshes and the actual initial geometry. Cell volumes, concentrations, conductances, and mass-weighted spectra update during the run. See [equations, conservation, and validation limits](live_transport.md).

## Larger-domain preset

Launch `python -m embryo.dashboard --config configs/large_domain.json` for grid 56 and half-width 2.24, retaining voxel spacing 0.08. Alternatively edit both fields together and Reset to a new zygote. The preset runs 2000 steps to time 30 with conservative signaling; it adds space without refining cells. It does not modify an already running session. [Domain protocol and limits](domain_size.md).
