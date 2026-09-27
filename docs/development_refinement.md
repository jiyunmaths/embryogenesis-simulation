# Full developmental trajectory refinement

This study starts each simulation from a single zygote and follows spontaneous shape-selected divisions, conservative signaling, fate dynamics, and mechanical feedback through **t = 90**. It addresses repeated cleavage and evolving geometry absent from the successful short manufactured-state confirmation.

## Fixed physical model and initialization

All cases use seed 7, domain half-width 2.24, interface width 0.085, 16-cell cap, conservative diffusivities 0.02 and 0.4, and the existing default reaction, division, polarity, and mechanical parameters. Signaling and polarity retain the shared 0.02 contact cutoff. Noise amplitudes are unchanged. No spindle axis or chemical patch is imposed.

The radius-0.8 analytic zygote is sampled directly on each voxel grid. Its continuous diffuse occupancy volume is computed by radial quadrature and supplied as the same volume target in every case. The initial chemical state, cycle draw, and random-stream states match. This removes grid-dependent initial target quadrature, without interpolating fields or restarting from a mature coarse embryo. It slightly changes initialization relative to historical runs, so all five cases are new.

| Case | Grid | Time step | Purpose |
|---|---:|---:|---|
| space-56 | 56³ | 0.0075 | Coarse spatial reference |
| space-72 | 72³ | 0.0075 | Intermediate spatial / intermediate temporal case |
| space-88 | 88³ | 0.0075 | Finest spatial case |
| time-0.015 | 72³ | 0.015 | Coarse temporal reference |
| time-0.00375 | 72³ | 0.00375 | Finest temporal case |

The spatial and temporal axes vary independently. Temporal testing on 72³ does not establish temporal accuracy on every finer grid; a further cross-check may be needed if sensitivity remains. The 56³ case is a reference, not a candidate production resolution. Its existing small-cell shortfall is not hidden.

## Observations and comparisons

Metrics and per-cell concentrations, fate, polarity, volumes, centers, IDs, and parent IDs are recorded every 0.6 time units. Every six time units, each branch saves a full checkpoint and a surface frame. These checkpoints cover initial development, signal amplification, and later relaxation. Complete lineage records retain actual division times and axes. A diagnostic subclass records volume and regulator-amount jumps immediately around every abscission; it does not change the equations or division rules.

Matching seeds does not force matching cleavage times, spindle choices, event order, or assignments of random draws to descendants. Numeric cell IDs are therefore not assumed to establish homologous cell identities across resolutions. Comparisons at the same physical times use:

- Volume-weighted Wasserstein distance between activator distributions and between inhibitor distributions; this accommodates differing cell counts without artificial cell matching.
- Unweighted cell-to-cell signal standard deviations.
- A, B, and uncommitted cell fractions.
- Aggregate axis ratio and cell count.
- Exact time of the last abscission, provided the final cell cap is reached.

Spatial signal organization is additionally available through saved centers, contact graphs, signal-dipole diagnostics, surface frames, and checkpoints. **Distribution agreement alone does not establish matching spatial patterns.** This is a developmental observable-refinement screen, not a full-field convergence proof. Saved data allow further spatial analysis if the initial screen succeeds.

## Prospective criteria

For the finest spatial pair and finest temporal pair, the maximum sampled differences must satisfy:

- Relative axis-ratio difference below 1%.
- Each signal distribution distance below 0.01, relative to equilibrium concentration 1.
- Each signal standard-deviation difference below 0.005.
- Each fate-fraction difference at most 1/16.
- Cell-count difference at most one.
- Last-abscission timing difference at most 0.6 model-time units, with completed development.

Shape and signal-distribution discrepancies must not increase from the coarser to finer pair, apart from a 1e-10 numerical floor. These are practical tolerances declared before running; they are not biological calibration or a guaranteed convergence rate.

Each case reports separate quality checks: every-step cell-volume error below 5%, minimum equivalent radius of at least four spacings, and no phase-field clipping; sampled boundary occupancy below 1%; and one connected component per non-dividing cell at phi=0.5. Pinching mothers are excluded from this connectivity test because the progressive ring can produce disconnected lobes before abscission. Relative volume and regulator-amount jumps at each abscission must remain below 1e-6.

Candidate acceptance requires all comparison gates and all quality gates for space-72, space-88, and time-0.00375. Coarse-case failures remain visible. At preparation, the prior single-cleavage screen had an unresolved initial interpolation-calibration failure (0.477% versus 0.25%). A subsequent [separate measurement validation](cleavage_measurement.md) now passes the unchanged thresholds using native-phase reconstruction. Neither the original failed report nor this developmental study’s frozen protocol or dynamics was changed.

## Execution and outputs

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.development_refinement prepare \
  --output outputs/development-refinement-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.development_refinement run \
  --output outputs/development-refinement-repeat
```

The production study is `outputs/development-refinement`. Its immutable protocol records full configurations and source-code hashes before starting. Cases execute sequentially. `status.json` reports the active case, current observation time, and completed cases; `run.log` reports checkpoint progress. Per-case `history.json` and `progress.json` expose observations, quality extrema, and lineage while running. Checkpoints and accumulated data support manual recovery; the runner deliberately refuses to overwrite or automatically restart an existing study.

Each completed case exports a final checkpoint, complete sampled history, lineage and event audits, plus an animated surface viewer spanning development. After all five cases finish, `comparison.json`, `RESULTS.md`, and `comparison.png` are written automatically. No numerical-quality failure is relabeled as success, and no physical parameter is tuned based on an attractive morphology.

These long runs are in progress; no developmental refinement result is claimed yet. Follow-up decisions depend on whether failures come from division timing, signaling, fate, shape, or numerical quality. More seeds, interface-width sensitivity, and causal feedback controls remain separate requirements.
