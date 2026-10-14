# Repeated developmental GPU qualification: prospective protocol

This test follows naturally scheduled cleavage from one analytic zygote toward sixteen cells. It tests the isolated developmental GPU method after the [accepted single-cleavage qualification](development_gpu_qualification_assessment.md). No accepted solver/kernel is changed. These are numerical repetitions of synthetic seed 7, **zero new scientific histories**; fresh independent histories follow only a complete qualification and separate prospective protocol.

## Fixed design

An analytic radius-0.8 phase field is sampled independently on every grid. Continuous radial quadrature supplies one common target volume. Initial activator/inhibitor concentrations are 1.03/0.98, polarity is zero and no identities, spindle directions or spatial chemical patches are prescribed. The existing shape-selected spindle rule chooses a direction within the longest-axis eigenspace; a nearly spherical cell has no unique axis. The original three RNG streams, cycle jitter 0.12 and amount-conserving chemical partition noise 0.001 are matched initially across cases. Matching a seed does not guarantee the same subsequent cleavage order or allocation of draws to descendants.

All paths use domain half-width 2.24, interface width 0.085, conservative transport, diffusivities 0.02/0.55, natural mean division interval 2.0, progressive cytokinesis duration 0.9 and ring strength 8. Follow through t=18, observe every 0.15 and checkpoint every three units. Require fifteen resolved divisions, sixteen final cells, no live furrows and at least three units after the last abscission. No unresolved neck is forcibly divided.

Four mechanical branches retain active chemistry and polarity: coupling off; activity-dependent tension/adhesion only (0.25/0.35); directional polarity tension only (0.35); and their full combination. Each branch has four GPU discretizations:

| Label | Grid | Timestep | Role |
|---|---:|---:|---|
| coarse | 72 cubed | 0.00375 | Temporal reference; complete independent NumPy trajectory |
| fine | 72 cubed | 0.001875 | Temporal candidate and coarsest spatial case |
| space88 | 88 cubed | 0.001875 | Intermediate spatial case |
| space104 | 104 cubed | 0.001875 | Finest spatial case |

Thus sixteen GPU trajectories and four complete coarse NumPy references. Four independent CPU workers may run alongside one serial GPU coordinator on the qualified GTX 1080 Ti UUID. PyTorch owns resident arrays/matrix operations; existing custom CUDA performs mechanics and geometry/polarity. No CPU fallback, altered feedback law or silent correction initialization is allowed.

## Event and physical criteria

Keep the accepted event limits: species-amount and target-volume jumps at most 2e-14; visible occupancy-volume jump at most 1e-6; represented pointwise occupancy jump at most 1e-14; dilution closure at most 2e-14. Every division must meet the existing neck, overlap and daughter-fraction rules. Record each event separately even when multiple mothers complete in one step. Other cells and numerical corrections must survive buffer resizing.

Every step must have relative cell-volume error below 5%, equivalent radius at least four spacings, positive finite concentrations/capacities and no mechanical clipping. At observations, boundary occupancy must be below 0.01 and every nondividing cell must have one component at phi >= 0.5; pinching mothers are excluded. Recorded event/step extrema are runtime evidence, separately from saved-array reconstruction.

Check eight exact restored/uninterrupted steps at the first mixed-count, live-furrow states in the 2-to-4, 4-to-8 and 8-to-16 transitions. Require all three windows, including another mother still dividing, and preserve tensors, corrections, targets, geometry, IDs, parents, lineage, events, clocks and all RNG states. Save the checkpoints for independent full replay. A missing restart window is a reported failure.

## Comparisons without assuming identical cell identities

At the same physical observation times compare volume-weighted activator/inhibitor distributions (Wasserstein distance in concentration and log concentration), chemical standard deviations, whole species amounts, total volume, volume-weighted mean polarity, tissue axis ratio and cell count. Compare eigenvalues of the symmetric volume-weighted transport Laplacian when counts agree; check constant preservation, amount conservation and weighted symmetry separately. Cellwise and matrix comparisons are diagnostics only unless lineage and partition-draw correspondence is established. No time shift, rigid registration, amplitude normalization or cell-ID matching repairs a failed raw comparison.

Backend limits are 1e-5 for chemical distribution/contrast, mean polarity, total volume, whole amounts and relative transport-spectrum norm; relative tissue-axis difference at most 1e-4. Complete NumPy/GPU comparisons require corresponding final parent trees and partition draws, all parent-wise event times within two coarse steps and native checkpoint occupancy/chemical-density/polarity-density discrepancies at most 2e-5. Explicitly allow up to thirty sampled observations to straddle slightly different events (two per division), with count difference at most one. Distribution, shape and whole-amount checks cover every observation, including count mismatches.

When ancestry and partition correspondence is established, also apply the single-cleavage backend limits to corresponding cells: chemical log error, polarity component error, relative cell-volume error and contact-operator Frobenius error at most 1e-5; native checkpoint phase error at most 2e-5. Require matched final descendants. Explicitly report temporary nonhomologous observations/checkpoints; allow at most thirty observations. These conditional checks prevent aggregate agreement from hiding individual-cell or network errors. A preliminary prepared directory was never run and is archived before adding these inherited checks; original acceptance thresholds and solver sources are unchanged.

Temporal limits: each concentration and log-concentration distribution distance at most 0.01; contrast difference at most 0.005; relative axis difference at most 0.01; relative total volume and whole amounts at most 0.005; mean-polarity component difference at most 0.01; relative transport-spectrum norm at most 0.01; raw checkpoint occupancy/chemical density relative L2 at most 0.01 and polarity-density absolute L2 normalized by occupancy L2 at most 0.01. Last-abscission time difference at most 0.05. Allow at most thirty count-mismatched observations and count difference at most one.

Spatial limits retain the same distribution/contrast/axis/volume/amount/polarity and field limits; relative transport-spectrum norm at most 0.05 and last-abscission difference at most 0.6. Both spatial pairs must pass, and finest-pair maximum concentration-distribution, axis and field discrepancies must not increase over the coarser pair, with an absolute 1e-10 numerical floor. These are practical observable/field refinement screens, not a formal convergence order, sharp-interface limit or arbitrary-contact continuum proof. Every candidate path must pass physical gates.

## Field reconstruction and calibration

At t=0,3,6,9,12,15,18 reconstruct individual native phase fields cubically on a common 112-cubed probe, then compute summed occupancy, chemical amount densities (sum of occupancy times cell concentration) and polarity densities. This avoids matching labels and tests spatial organization beyond distributions. Do not clip reconstructed phase values, subtract calibration error or interpolate coarse initial states into fine ones. Report phase ranges.

At the final endpoint also use probes 112/128 with cubic/quintic reconstruction. Every variant must pass the field limits; the spread of discrepancy across variants must be at most 0.001. Independently calibrate all source grids/orders/probes on the actual initial sphere and translated/deformed analytic holdouts before launching trajectories, retaining the established 0.0025 relative occupancy-L2 calibration limit. Native same-grid temporal/backend comparisons need no interpolation.

## Execution and acceptance

`python -m embryo.qualify_repeated_development prepare --output output/studies/repeated-development-gpu-qualification` prepares fresh analytic inputs and fixed source/input pins. `run` requires successful calibration, explicit GPU selection, an unchanged prepared protocol and one coordinator lock. Per-case history, lineage, every-event audits, quality extrema and progress expose partial results; completed/failed cases are not silently rerun. Independent eight-step NumPy chemical/dilution checks accompany each GPU start. CPU references start from their own analytic inputs, not from mature GPU states.

Failures remain recorded, thresholds are not retuned after seeing results, and all evidence from the accepted single-cleavage and mature studies is preserved. Production completion is provisional until a separate literal saved-array review, integrity audit and required restart replay are complete. A positive chemical pattern is neither necessary nor sufficient for numerical acceptance. This qualification does not establish biological identities, inherited differentiation or independent fresh developmental replication.
