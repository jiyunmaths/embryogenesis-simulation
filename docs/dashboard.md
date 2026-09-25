# Live simulation dashboard

Start from the repository with its Python dependencies installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard
```

Open <http://127.0.0.1:8765>. Optional arguments are `--port` and `--config path/to/config.json`. The installed command `embryo-dashboard` starts the same server. Stop it with Ctrl+C. The dashboard uses the existing `Simulation.step()` implementation; it introduces no alternative dynamics or new runtime dependencies.

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
| `POST /api/run` | Start/resume; optional `{"config": {...}}` |
| `POST /api/pause` | Stop at a step boundary; body `{}` |
| `POST /api/reset` | New zygote with `{"config": {...}}` |

Each reset increments `generation`; frame revisions are only meaningful within that generation. POST requests require `Content-Type: application/json`. Cross-origin and nonlocal Host requests are rejected. The dashboard is designed for one local interactive experiment, not as a shared remote service.

## Development checks

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q
# Optional Node.js 18+ check of controls, resets, history, and network races:
node --test tests/test_dashboard_client.cjs
```

The JavaScript tests use a minimal DOM fixture and the real Python parameter schema. They verify client behavior, not browser layout. Node.js is only needed for these tests, not for running the dashboard.
