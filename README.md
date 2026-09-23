# Embryogenesis in 3D

## Scientific objective

**To what extent can mathematical equations explain the emergence of organized living structure from an initially simple zygote?** This project investigates how small fluctuations in an approximately uniform initial state can develop into reproducible cell identities, spatial signals, and complex, changing three-dimensional shapes.

The central mechanism to investigate is an **explicit activator–inhibitor feedback loop**: an activator promotes its own production and the production of an inhibitor; the inhibitor suppresses activation. We will study how the kinetics, interaction ranges, and coupling to cell mechanics determine whether this feedback amplifies fluctuations into order, suppresses them, or produces unstable patterns. Activator and inhibitor are signaling variables; they are distinct from the two initial cell identities, A and B.

We will evaluate three linked outcomes separately:

| Outcome | Scientific question | Evidence to measure |
|---|---|---|
| Cell identity differentiation | Can initially similar cells acquire distinct, persistent identities through local regulatory dynamics and signaling? | Activity distributions, fate persistence, lineage histories, and responses to perturbations |
| Signaling symmetry breaking | Can nearly uniform signaling develop persistent spatial domains, poles, or axes without a prescribed directional cue? | Pattern onset, spatial correlations, wavelength, domain number, and orientation across independent runs |
| Geometry and shape symmetry breaking | Can signaling and cell behavior produce sustained changes in tissue shape while the embryo divides and rearranges? | Shape anisotropy, axis persistence, cell organization, and later folding or cavity formation where the model supports them |

The intended feedback is reciprocal: **signaling influences identity and cell mechanics; cell movement, division, contacts, and shape alter signaling in return.** The embryo boundary must evolve with these interactions. A spatial identity pattern can exist in an approximately spherical embryo, so signaling and differentiation alone do not establish geometric symmetry breaking.

## How we will assess explanatory power

The aim is to identify the smallest interpretable set of equations that accounts for progressively more developmental organization, and to document where that explanation fails. Two identities and a small cell population are the starting point; more complex structures are later tests of the same approach.

- **Separate emergence from assumptions.** Record which features arise from the dynamics and which are supplied through initial conditions, boundary conditions, cell-cycle rules, or prescribed forces. Begin spontaneous-symmetry-breaking experiments with unbiased fluctuations; label imposed gradients and asymmetries as separate controls.
- **Establish the mechanism.** Analyze steady states and their stability, identify parameter regimes where spatial perturbations grow, and test the necessity of activation, inhibition, signal transport, and mechanical feedback by disabling them individually. An activator–inhibitor loop does not automatically imply a Turing instability or a single developmental axis.
- **Measure robustness and limits.** Compare independent seeds, parameter ranges, perturbations, and numerical resolutions. Report uniform, mixed, fragmented, and failed outcomes as well as organized ones; distinguish reproducible structure from transient noise and grid artifacts.
- **Increase biological specificity only with evidence.** Start with a generic, dimensionless model. Quantitative claims about a particular embryo will require experimental calibration and independent validation. Active signaling, proliferation, and force generation remain explicit biological inputs to the model.

This objective guides subsequent model changes and experiments. Each extension should state the phenomenon it aims to explain, the feedback it introduces, a control that could challenge the proposed explanation, and a measurable success criterion. Visual resemblance alone is insufficient.

## Current starting point

The runnable prototype models one cell dividing into a deformable multicellular aggregate. Gierer–Meinhardt activator and inhibitor activities evolve on its normalized contact graph and drive a downstream bistable fate switch. Apical–basal polarity develops from exposed cortex and neighboring orientations, and changes cortical tension directionally.

**Finite-graph stability is analyzed before 3D evolution.** The solver uses a constant-preserving random-walk normalized Laplacian, checks each discrete eigenvalue against the reaction–diffusion Jacobian, and records changes in spectral modes at every cleavage. A continuous unstable band can contain no supported modes on a small graph. Activator/inhibitor variables are distinct from cell fate; the old independent fate noise and geometry biases are disabled by default. See [graph signaling and polarity](docs/graph_signaling.md) for equations, spectral-gap results, and limitations.

This is an exploratory model in dimensionless units, not a reconstruction of a particular organism. The simulation evolves each cell's shape on a 3D grid. It does not prescribe an embryo outline or assign daughter identities.

## Run

From this directory, using Python 3.10 or newer:

```bash
python -m pip install -e '.[test]'
# Analyze finite reference graphs before a coupled run (no 3D simulation):
OPENBLAS_NUM_THREADS=1 python -m embryo.graph_analysis --output outputs/my-graph-analysis
# Every coupled run also writes a finite-graph preflight before its first step:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo --output outputs/my-run
```

NumPy, SciPy, and Matplotlib are the runtime dependencies. The thread settings avoid excessive BLAS threading overhead for the small contact matrices. Existing output directories are never overwritten.

Open `outputs/my-run/viewer.html` in a browser. It is self-contained and works offline:

- Drag to rotate; scroll to zoom.
- Play or scrub the timeline from the zygote through cleavage and differentiation.
- Move the cutaway plane to inspect internal cells.
- Color by fate, activator, inhibitor, or cell ID; show apical polarity arrows.
- Follow the principal axis ratio and volume diagnostics.

Signal colors are relative to the homogeneous activity of one, not automatically rescaled to exaggerate small variations. The viewer displays samples of the simulated cell surfaces. Its apparent surface texture is a rendering approximation; it does not enter the dynamics.

## Experiments

```bash
# Same seed and parameters, without fate/polarity-dependent mechanical properties:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-feedback --output outputs/no-feedback

# Cleavage and mechanics, with signaling, polarity, and differentiation disabled:
OPENBLAS_NUM_THREADS=1 python -m embryo --mechanics-only --output outputs/mechanics

# Separate signaling/polarity ablations:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-signaling --output outputs/no-signaling
OPENBLAS_NUM_THREADS=1 python -m embryo --no-polarity --output outputs/no-polarity

# More cells and finer geometry (substantially more expensive):
OPENBLAS_NUM_THREADS=1 python -m embryo --grid 64 --max-cells 32 --steps 1200 --output outputs/larger

# Customize coefficients with a JSON file:
OPENBLAS_NUM_THREADS=1 python -m embryo --config examples/default.json --output outputs/custom

# Small paired-seed ablation screen; this is not statistical validation:
OPENBLAS_NUM_THREADS=1 python -m embryo.experiments --seeds 7 11 23 --output outputs/screen

# Deterministic four-cell spatial/time refinement checks:
OPENBLAS_NUM_THREADS=1 python -m embryo.convergence --output outputs/convergence.json
```

`--no-feedback` retains signaling and polarity dynamics but disables their fate/polarity-dependent mechanical effects. `--mechanics-only` also disables signaling, polarity, and differentiation. Mechanical and regulatory RNG streams are separate. Shape-dependent orientation and mechanically gated abscission can change division directions and timing across feedback controls, even at a fixed seed.

## Outputs

| File | Contents |
|---|---|
| `viewer.html` | Offline, interactive 3D playback |
| `summary.png` | Shape snapshots, fate counts, and axis ratio |
| `config.json` | Full experiment parameters |
| `metrics.csv` | Volume errors, shape measures, and fate counts over time |
| `trajectory.json` | Sampled surfaces and activities for playback |
| `graph_preflight.json` | Discrete-mode stability on 4-, 8-, and 16-node reference graphs, computed before 3D stepping |
| `graph_history.json` | Actual weighted contacts, activities, polarity, and frozen-graph stability over time |
| `cleavage_spectra.json` | Before/after graphs, inherited-mode transfer, and signal partition jumps at abscission |
| `lineage.json` | Every cell's birth, parent, and division time |
| `final_state.npz` | Full final fields and RNG states for continuation |
| `diagnostics.json` | Run timing and numerical diagnostics |

The sampled trajectory is for visualization, not full-field reconstruction at earlier times. Current (schema 3) checkpoints can be continued exactly. Earlier model checkpoints are rejected rather than silently resumed with different dynamics:

```python
from embryo import Simulation
sim = Simulation.restore("outputs/my-run/final_state.npz")
for _ in range(100):
    sim.step()
sim.checkpoint("outputs/continued.npz")
```

## What is implemented

- Freely deforming cell interfaces with soft volume constraints and repulsion.
- Interface attraction, optionally dependent on fate similarity.
- Shape-aligned division, progressive equatorial constriction, and mechanically gated abscission.
- Exact final occupancy partitioning, conserved mother volume during cytokinesis, and lineage records.
- Gierer–Meinhardt activator/inhibitor dynamics on a normalized contact graph.
- Exact discrete-mode linear stability, spectral gaps, and cleavage mode-transfer diagnostics.
- A downstream signed fate switch driven by activator activity.
- Apical–basal polarity and directional cortical tension with conservative spatial fluxes.
- Fate-dependent effective surface tension.
- Numerical diagnostics, paired ablations, and reproducible checkpoints.

**The default target is 16 cells on a modest grid.** Resolve and validate this before scaling to the planned 32–64-cell model. Volume constraints are soft between divisions; inspect the measured errors. Surface/interface widths and small daughter cells need resolution studies.

A/B labels are instantaneous activity thresholds. They do **not** establish stable commitment. There is no tensile-stress-based spindle rule, lumen, growth between divisions, extracellular morphogen field, or calibrated gene network yet. Graph activities have normalized exchange rates; they are not a volume-conserving molecular diffusion model. Equatorial contraction is a prescribed ring surrogate, not a resolved actomyosin network. Elongation following cleavage is not proof of a spontaneously selected developmental axis.

## Tests and next steps

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

Tests check shape-oriented spindles, gradual furrowing, volume conservation, daughter connectivity, contact symmetry, single-cell isotropy, exact continuation through active cytokinesis, cell-cap reservations, controls, and a basic time-refinement check. They are software/numerical checks, not biological validation.

The graph preflight is in `outputs/graph-analysis`; the coupled example is `outputs/graph-polarity/viewer.html`, with a matched no-polarity control in `outputs/graph-no-polarity`. Earlier outputs are historical models. In the initial time-15 coupled example, signals remain small and all 16 cells remain uncommitted: a positive frozen-graph growth rate does not guarantee a mature pattern on the available developmental timescale.

Read [the cytokinesis update and validation](docs/cytokinesis.md), [the equations and assumptions](docs/model.md), [the staged implementation plan](docs/plan.md), and [the initial numerical results](docs/results.md).
