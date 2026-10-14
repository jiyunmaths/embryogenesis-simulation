# Running and reproducing the project

The [README](../README.md) describes the scientific model and claims. The [lab notebook](lab_notebook.md) records execution status; [backend policy](backends.md) defines which numerical regimes are qualified. A demonstration run is not a reproduction of a reviewed scientific experiment.

## Install

From the repository, with Python 3.10 or newer:

```bash
python -m pip install -e '.[test]'
```

Core dependencies are NumPy, SciPy, Matplotlib and scikit-image. Scientific GPU paths additionally require the compatible PyTorch/CUDA environment and compiled kernels described in the backend documentation. They are not automatically enabled by installing the package or opening the dashboard.

## Interactive demonstration

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765). Run/pause/resume advances the simulation; reset applies edited parameters to a new zygote. Cell-surface rendering, polarity/contact views, retained-frame playback and benchmark controls are explained in [the dashboard guide](dashboard.md).

The historical dashboard/core model includes a supplied fate switch and can show A/B labels. Those labels are not used in the current attribute-based identity research. Downloaded visualization frames are not full restart checkpoints.

## Scientific experiments

Use an experiment's frozen protocol, prescribed preparation and gated coordinator. The [notebook report catalog](lab_notebook/run_records.md) links launch records and assessments; each study's method document specifies prerequisites and invocation. Full source checkpoints and local result artifacts are required when the preparation depends on earlier studies.

Do not substitute a dashboard preset, restore into a different simulation class, or change a pinned parameter merely to launch a reviewed experiment. Such changes define a new preparation and need their own numerical checks. Observe the distinction between fresh zygote development and continuation of prepared mature chemistry.

| Compute path | Role |
| --- | --- |
| `AttributeSimulation` / native CPU variants | Attribute-based reference and historical research preparations; supported numerical scope depends on the protocol. |
| `ReferenceDevelopment` / `DevelopmentGpu` | Progressive cleavage with mandatory phase-carry, chemical-correction and developmental checkpoint fields. |
| Resident PyTorch/custom CUDA mature adapters | Gated continuations of compatible mature states, with no unsupported division or silent schema transition. |
| SciPy frozen ODE solvers | Chemical systems on measured fixed operators, with independent integration comparisons. |
| Core `Simulation` | Historical supplied-fate demonstrations and reproduction of their original outputs. |

The qualified GPU architecture uses PyTorch arrays/matrix operations and custom CUDA mechanics/geometry/polarity. GPU choice is machine-specific; use the designated device recorded in the launch metadata. Avoid competing coordinators on the same device. [Backend scope](backends.md), [developmental method](development_gpu_qualification.md), [fresh-zygote protocol](fresh_zygote_pilot.md).

## Outputs and checkpoints

| Artifact | Meaning |
| --- | --- |
| `protocol.json` | Frozen design, thresholds, configuration, source/input hashes and interpretation limits. |
| `status.json`, per-job progress | Live execution state, including stops and failures. Completion is not acceptance. |
| History/network/event records | Observed chemistry, volumes, contacts, spectra, polarity, geometry and cleavage diagnostics. |
| Full `.npz` checkpoint | Physical state, lineage, random state and the correction/schema fields required by that backend. |
| Result/comparison records | Protocol-associated quality gates and paired numerical decisions. |
| Independent assessment | Reviewed evidence with its stated scope; preserves earlier failed comparisons. |
| `viewer.html` / sampled surfaces | Visualization; insufficient for complete state reconstruction or exact continuation. |

The large `outputs` and `output/studies` trees live on the data drive through symlinks. Some protocols retain canonical paths and require the SHA-validated relocation wrapper. Preserve original paths and use the launch record's adapter rather than rewriting hashed paths. [Storage instructions](output_storage.md).

Never overwrite recorded scientific outputs. Use a fresh directory for a new preparation and the study's explicit restart/resumption procedure for an interrupted run; not every coordinator supports resumption. Validate the required checkpoint schema before stepping it.

## Software verification

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q
```

GPU tests are opt-in and need the compatible environment/device; see the corresponding test module and qualified backend protocol before enabling them. Software tests, backend equivalence, spatial/time accuracy and scientific mechanism tests are different decisions. Use checks appropriate to the changed component, and retain all scientific gate failures.

The [verbatim former README](../archive/study-document-snapshots/2026-10-11-before-readme-notebook-consolidation/README.md) preserves earlier CLI examples and output tables as historical documentation. Current experimental parameters and accepted scope come from the recorded protocols and assessments.
