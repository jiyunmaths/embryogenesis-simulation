# Persistent shape asymmetry: coupled-model pilot

## Scientific question

Does the activator–inhibitor / fate / polarity feedback produce a persistent embryo-scale shape axis beyond the anisotropy inherited from cleavage and ordinary mechanical relaxation?

This experiment returns to the **coupled deformable-cell model**, rather than prescribing a moving domain. A spherical zygote cleaves progressively; its cells deform under the existing phase-field forces. Signaling follows the normalized contact graph, with downstream fate dynamics and apical–basal polarity. The [conservative continuum benchmarks](continuum_bridge.md) remain separate: their concentration transport has not been substituted into the coupled embryo model.

The implementation is [embryo/shape.py](../embryo/shape.py), with [independent diagnostic and intervention tests](../tests/test_shape.py). No new shape force, imposed signaling gradient, target elongation, or parameter strengthening is introduced in this pilot. It tests the current equations. A negative feedback comparison is a useful result, not a reason to alter criteria after the run.

## Mechanism being tested

The existing model transmits chemical differences to mechanics through two routes:

1. Activator activity biases a continuous two-fate variable. That variable changes each cell's interfacial tension and pairwise attraction when mechanical feedback is enabled.
2. Activator activity changes the strength of the exposed-cortex drive on apical–basal polarity. Polarity changes cortical tension directionally when mechanical feedback and polarity are enabled.

Shape and contacts in turn alter graph exchange and exposure cues. These are proposed feedback routes; the existence of those terms does not show that they are sufficient or necessary for a tissue-scale axis. The equations and assumptions are detailed in [the model documentation](model.md), [graph signaling](graph_signaling.md), and [the README method](../README.md#method-and-mathematical-model).

Legacy fate noise, fate partition noise, exposure bias, and direct neighbor fate inhibition are zero. Small signal partition fluctuations are present during division, with independent mechanical, fate, and signaling random streams. Once the 16-cell cap has been reached, no further cleavage or partition noise occurs. The first cleavage axis is selected without a prescribed spatial direction from the approximately spherical zygote; later spindles follow cell shape.

## Developmental history and interventions

The default pilot uses seed 7, grid $40^3$, extent 1.6, mechanical step 0.015, interface width 0.085, a 16-cell cap, and the existing mechanical and signaling coefficients. The earlier [timescale experiment](timescales.md) already developed this full-feedback zygote through $t=60$. Its full checkpoint supplies exactly the same mature physical state to three branches:

| Run | Initialization | Intervention | Question |
|---|---|---|---|
| Full feedback | Existing mature checkpoint at $t=60$ | Continue all existing couplings | Does its shape persist to $t=90$? |
| Feedback removed | Identical checkpoint | Set `feedback=False`; retain chemical, fate, and polarity states/dynamics | Is continuing mechanical feedback necessary to maintain the existing shape? |
| Polarity tension removed | Identical checkpoint | Set `polarity_tension=0`; retain fate-dependent mechanics and polarity dynamics | Is directional polarity tension necessary for maintenance? |
| Development without feedback | New spherical zygote, same seed and non-intervention parameters | Set `feedback=False` from $t=0$ and run to $t=90$ | Can cleavage and mechanics generate comparable anisotropy without chemical/fate/polarity feedback into mechanics? |

The mature branches preserve all phase fields, targets, cell IDs, regulators, polarity vectors, clocks, and random generator states at intervention. The checkpoint SHA-256 and its original configuration are recorded in `protocol.json` and `analysis.json`.

The developmental control has its **own** evolving geometry and cleavage history. A shared seed does not make those trajectories identical. Signals and fates remain active as passive readouts when `feedback=False`; they cannot change mechanics in that control. Removing feedback from a mature state tests maintenance, while removing it from the zygote tests its developmental contribution. The two experiments cannot be substituted for one another.

## Shape, axis identity, and memory

Shape is measured from the capped diffuse union of the cell occupancies, not only from cell centers:

$$
w(\mathbf x)=\min\!\left(1,\sum_i h(\phi_i(\mathbf x))\right),
\qquad h(\phi)=\phi^2(3-2\phi).
$$

Its centroid and covariance are

$$
\bar{\mathbf x}=\frac{\sum_{\mathbf x}w(\mathbf x)\mathbf x}{\sum_{\mathbf x}w(\mathbf x)},
\qquad C=\frac{\sum_{\mathbf x}w(\mathbf x)
(\mathbf x-\bar{\mathbf x})(\mathbf x-\bar{\mathbf x})^{\mathsf T}}
{\sum_{\mathbf x}w(\mathbf x)}.
$$

With eigenvalues $\lambda_1\leq\lambda_2\leq\lambda_3$, the shape ratio is

$$
R=\sqrt{\lambda_3/\lambda_1}.
$$

A unique long axis is reported only when $(\lambda_3-\lambda_2)/\lambda_3\geq0.05$. An almost spherical or oblate shape can otherwise yield an arbitrary eigenvector; its long-axis direction is recorded as `null`, not treated as measured axis memory.

Axis orientation is unoriented: $\mathbf e$ and $-\mathbf e$ represent the same axis. The angle between two axes is

$$
\theta=\arccos\!\left(\frac{|\mathbf e_1\cdot\mathbf e_2|}
{\|\mathbf e_1\|\|\mathbf e_2\|}\right).
$$

The analysis measures rotation relative to the first sampled late-window axis and alignment with the **recorded first-cleavage axis**. The latter is essential because persistent shape may simply preserve a cleavage-selected direction. A volume-weighted signal dipole and its angle to the shape axis are descriptive additional outputs, not acceptance criteria or proof of causation.

## Predeclared persistence and control criteria

All default runs are compared over the sampled late window $75\leq t\leq90$, with observations every 0.6 time units. A finite-window shape-persistence result requires:

- The late window is covered by at least three observations and begins after the last completed cleavage, with no active cytokinesis during the window.
- Axis ratio stays at least **1.20**.
- The long axis remains identifiable at every sampled time.
- Maximum unoriented rotation from the late-window anchor is at most **15 degrees**.
- The range of axis ratios divided by their mean is at most **5%**.

Activator standard deviation of at least 0.1 throughout the late window is recorded separately. Shape persistence must not be defined by the presence of chemical contrast.

For each causal comparison, the full-feedback axis ratio must exceed the control by at least **0.05 at every matched late sample** to meet the declared feedback-excess screen. Differences are retained with their sign; an equal or more elongated control is not a positive feedback result. The branches are paired by checkpoint hash. The zygote control is explicitly labeled as a different developmental history, with matching non-intervention parameters.

Numerical screens are also separate from morphology: individual volume errors at most 5%, boundary occupancy at most 0.01, phase-field clipping at most 1%, minimum equivalent cell radius at least four grid spacings, each cell connected at $\phi>0.5$, and a connected thresholded contact graph. Union connectivity is recorded at occupancy thresholds 0.1 and 0.5 to expose dependence on the diffuse-interface threshold; these are diagnostics, not proof of a sealed epithelium.

Passing these screens would still not establish convergence, long-time attraction, isotropic axis selection across seeds, or biological validity. The current model's force laws and graph-cutoff sensitivity require their own checks.

## Reproduce

The source mature checkpoint can be reproduced with the timescale capture:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.timescales capture \
  --output outputs/shape-source-baseline --duration 60
```

Then run the three branches with that same checkpoint:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape run \
  --checkpoint outputs/shape-source-baseline/final_state.npz \
  --output outputs/shape-repeat/full --mode full --duration 30

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape run \
  --checkpoint outputs/shape-source-baseline/final_state.npz \
  --output outputs/shape-repeat/no-feedback --mode no_feedback --duration 30

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape run \
  --checkpoint outputs/shape-source-baseline/final_state.npz \
  --output outputs/shape-repeat/no-polarity-tension --mode no_polarity_tension --duration 30

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape run \
  --output outputs/shape-repeat/development-no-feedback --mode no_feedback --duration 90

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape compare \
  --runs outputs/shape-repeat/full outputs/shape-repeat/no-feedback \
    outputs/shape-repeat/no-polarity-tension outputs/shape-repeat/development-no-feedback \
  --output outputs/shape-repeat/comparison
```

Every destination must be new. For checkpoint branches, the saved spatial grid and time step are retained; `--grid`, `--dt`, and `--seed` apply only to runs from a zygote. Resampling a mature state is not part of this protocol.

Each run saves its protocol, configuration, sampled tensor and quality diagnostics, cell surfaces, lineage, and final checkpoint. `viewer.html` provides offline 3D playback; `trajectory.json` contains the same visualization frames. The comparison writes `analysis.json`, `persistence.png`, and `final_shapes.png`. The final-shape figure uses the same camera, axis limits, and fate-color scale for every run, and overlays measured and first-cleavage axes. Morphological criteria can fail without a solver error: the run still saves a complete report. Generated outputs are local and ignored by Git.

## Measured result: persistence without feedback-specific excess

The completed seed-7 experiment is saved in `outputs/shape-persistence`. Its full-feedback source is `outputs/signaling-timescales/baseline/final_state.npz`, whose first development ended at $t=60$ and whose last cleavage was at $t=12.375$. The no-feedback developmental control completed cleavage at $t=12.360$.

| Run | Final axis ratio | Mean late ratio | Maximum late axis rotation | Final angle to first cleavage |
|---|---:|---:|---:|---:|
| Full feedback continuation | 1.32298 | 1.32222 | 0.1272° | 1.3622° |
| Feedback removed at $t=60$ | 1.32784 | 1.32580 | 0.0525° | 1.1259° |
| Polarity tension removed at $t=60$ | 1.32764 | 1.32557 | 0.0267° | 1.1269° |
| Development without feedback | 1.33056 | 1.32852 | 0.0255° | 0.1739° |

All four runs meet the **finite-window shape-persistence** conditions over $t=75$–90. They also retain signaling contrast. However, the mean full-feedback axis-ratio excess is **−0.00358** versus the matched feedback-removed branch, **−0.00335** versus the polarity-tension-removed branch, and **−0.00630** versus development without feedback. All three feedback-excess checks fail. These are paired trajectory differences, not confidence intervals or ensemble effect estimates.

The observation supports persistent elongation in the present equations but **does not demonstrate a feedback-specific shape axis**. Comparable elongation arises when signaling, fate, and polarity cannot affect mechanics, and all final long axes remain close to the first-cleavage direction. That is consistent with inherited cleavage geometry and mechanical relaxation. It does not prove that chemical feedback has no effect in other regimes or under perturbation; it does rule out treating this particular elongated appearance alone as evidence for that explanation.

Every run also fails the predeclared smallest-cell resolution screen: minimum late equivalent radii are approximately **3.969–3.972 grid spacings**, just below the cutoff of four. Other sampled late numerical screens pass: the largest cell-volume error is below **1.284%**, maximum boundary occupancy below **0.00441**, every cell has one connected $\phi>0.5$ component, and the contact graph is connected. The diffuse union has one connected component at occupancy 0.1 but 16 at occupancy 0.5. The view therefore represents interacting diffuse cells, not verification of a sealed epithelial body or lumen.

These are coarse pilot results. No mechanical refinement, new grid-rotation experiment, recovery-after-perturbation test, or seed ensemble has been completed here. Software verification passes **286 tests**; this does not turn a failed numerical or causal screen into a successful scientific conclusion.

Inspect `comparison/persistence.png` for the histories and `comparison/final_shapes.png` for the matched-camera geometry/axis comparison. Each run folder contains a separately labeled `viewer.html` for actual trajectory playback. The comparison JSON preserves the negative excess values and failed checks.

## Limits and decision rule

This is a single-seed controlled pilot. Before a robustness claim, the roadmap requires at least 20 independent seeds per screened parameter set, mechanical grid/time refinement, domain-size checks, and rotational tests. Rotating a covariance tensor in a unit test validates the diagnostic only; it does not test grid bias in the dynamics.

If controls retain comparable elongation, the correct conclusion is that the current experiment has not demonstrated a mechanical-feedback-specific shape axis. The next mechanistic study should distinguish inherited cleavage geometry, packing/relaxation, and chemical control of mechanical stress. Parameter changes should be declared before a new experiment and challenged with the same controls. The conservative continuum-transfer work remains a separate unfinished track.

## Long-time continuation

To distinguish a slowly developing response from an apparent plateau, extend each mature branch from $t=90$ to $t=180$ with its own saved state and **unchanged time step and mechanical parameters**. The full, feedback-removed, and polarity-tension-removed trajectories still descend from the common $t=60$ state. Their $t=90$ checkpoints are distinct because their interventions have already acted for 30 time units. This extension does not add growth or raise the cell cap. The separate developmental control remains a completed $t=90$ experiment unless explicitly extended.

```bash
for branch in full no-feedback no-polarity-tension; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape extend \
    --source "outputs/shape-persistence/$branch" \
    --output "outputs/shape-extended-repeat/$branch" \
    --until 180 --late-duration 30
 done

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.shape compare \
  --runs outputs/shape-extended-repeat/full outputs/shape-extended-repeat/no-feedback \
    outputs/shape-extended-repeat/no-polarity-tension \
  --output outputs/shape-extended-repeat/comparison
```

The `extend` command joins previous and new histories/playbacks without duplicating the boundary frame. It retains the original intervention provenance and additionally records the continuation checkpoint and source analysis hashes. Source results are untouched. Late persistence and feedback-excess criteria are applied over $150\leq t\leq180$; the earlier trajectory remains visible for interpreting slow trends. The saved runtime refers only to the newest continuation.

The [controlled chemical-patch experiment](signal_patch.md) is separate: it imposes a regulator pattern to test mechanical responsiveness, then restores autonomous signaling. It cannot substitute for the unforced long-time comparison.

### Completed continuation results

The three seed-7 branches reached $t=180$ with the optimized, trajectory-equivalent kernel. The late window contains 51 samples over $150\leq t\leq180$.

| Branch | Axis ratio at $t=90$ | Axis ratio at $t=180$ | Late ratio slope per time unit | Final angle to first cleavage |
|---|---:|---:|---:|---:|
| full | 1.322977 | 1.327776 | 6.82e-05 | 2.030° |
| no-feedback | 1.327839 | 1.340992 | 0.000104 | 1.065° |
| no-polarity-tension | 1.327635 | 1.341178 | 0.000111 | 1.037° |

All three pass the finite-window shape-persistence screen. The mean full-minus-control ratio differences are **−0.01296** (feedback removed) and **−0.01305** (polarity tension removed), so neither comparison passes the predeclared +0.05 feedback-excess criterion. This extension does not reveal large or feedback-specific elongation in the tested trajectory.

**The numerical screen does not pass.** Every branch still narrowly fails the four-grid-spacing minimum-radius criterion. In addition, the full-feedback branch first crosses the boundary-occupancy threshold 0.01 at sampled $t=124.2$; its final boundary occupancy is 0.06164. Both ablations remain below that boundary threshold. Here boundary occupancy is the maximum diffuse-union occupancy on the outermost grid layer, not a fraction of total volume lost. Crossing the screen identifies a need for domain-size verification, not by itself proof of a specific boundary artifact. Individual-volume, clipping, cell-connectivity, and contact-connectivity late checks pass. The late full/control difference therefore needs a larger-domain and spatial/time-refinement check before physical interpretation; it is not a validated long-time limit or evidence that feedback can never generate shape.

Results are in `outputs/shape-extended/{full,no-feedback,no-polarity-tension}`. Each playback includes its original $t=60$–90 trajectory plus the continuation. `outputs/shape-extended/comparison` contains tensor/axis histories, matched-camera final shapes, and the signed causal comparisons. Source $t=90$ outputs are unchanged.
