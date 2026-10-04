# Retired backend prototypes

This archive freezes the current source snapshots of three exploratory backend families. They are **retired from new scientific scheduling**:

| Prototype | Why archived | Supported replacement |
|---|---|---|
| Compiled C polarity prototype | Intermediate optimization before the fused native and resident GPU implementations | C++/OpenMP reference or custom CUDA geometry/polarity |
| Host-transfer CUDA mechanics pilot | Isolated precision/timing and partial-step experiments; no accepted resident full-system scientific path | Resident PyTorch/custom-CUDA backend |
| Pure PyTorch snapshot mechanics | Fixed-input component benchmarks and CUDA-graph replays; not evolving coupled trajectories | PyTorch arrays/matrices plus custom CUDA mechanics |

The files under `snapshots/` are byte-for-byte archival copies, not a new installed simulation package. `manifest.json` records each original path, snapshot hash, relevant pinned reports, and whether a recorded version matches this snapshot. Older report versions that differ remain explicit; the archive does not invent missing historical code versions.

Original files remain at their recorded `embryo/` paths when needed for immutable result hashes or historical imports. That retention does not make a prototype scientifically supported. Changing or moving a pinned original would invalidate evidence verification. None is selected by the current mature GPU runner.

The C `FastAttributeSimulation` implementation is **not archived as unused**: `NativeSimulation` inherits its geometry cache and fallback. Likewise `cuda_mechanics.cu` and `cuda_spatial_bench.cu` are active includes of `gpu_kernels.cu`, despite their historical filenames. They are preserved as active dependencies; copies alongside a retired pilot are support snapshots only. Legacy resumption utilities also remain where accepted protocol adapters use them.

Verify this archive with:

```bash
python archive/backends/verify.py
```

This checks both snapshot and retained original bytes. To reproduce a historical prototype, use its recorded original module and exact report version, with the original opt-in tests and device setup. Snapshot timings do not authorize scientific backend use. New moving work must pass the full accepted trajectory and context gates described in [the backend policy](../../docs/backends.md).
