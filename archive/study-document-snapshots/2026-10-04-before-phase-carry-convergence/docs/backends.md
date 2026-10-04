# Backend support and archival policy

The active scientific compute paths are intentionally limited:

| Path | Role and supported scope | Acceptance evidence |
|---|---|---|
| NumPy/SciPy `AttributeSimulation` | Developmental reference, division, and small-grid checks | Original component/development methods; full developmental convergence remains unresolved |
| C++/OpenMP `NativeSimulation` | CPU reference for mature moving mechanics and GPU comparison | Native mechanics/trajectory checks and accepted source transitions |
| Resident PyTorch/custom CUDA `GpuSimulation` | Accepted mature, conservative, direct-coupling continuations with no further division | Four full CPU/GPU control/pulse replays plus each new context's gated prefix; smaller-timestep gates checked separately |
| SciPy DOP853 and Radau | Small frozen chemical systems, independent integration checks, analytic Jacobian | Every neighbor/delayed trajectory checked with both methods; other frozen studies retain their original stated checks |

PyTorch owns GPU arrays and matrix operations. Custom CUDA performs mechanics and spatial geometry/polarity. GPU division, alternative coupling laws, the historical fate/dashboard dynamics, and arbitrary new parameters are not automatically covered by mature-backend acceptance. GPU acceleration is not silently enabled in the dashboard.

## Archived prototypes

[The backend archive](../archive/backends/README.md) contains immutable source snapshots and a manifest for the old C polarity prototype, the host-transfer CUDA mechanics pilot, and the pure-PyTorch fixed-input snapshot mechanics. Their timing/precision results remain historical engineering evidence. They are retired from new scientific scheduling and receive no new correctness claims from the active backend's tests.

Pinned originals remain in `embryo/` to preserve historical imports and protocol hash verification. This is a provenance-preserving archive, not a deletion of recorded evidence. The manifest lists recorded versions that match or differ from each archived snapshot. Original benchmark result directories remain under `outputs/` and are inventoried in the [results ledger](results_ledger.md).

## Dependencies that must stay active

`NativeSimulation` inherits `FastAttributeSimulation`, so `fast_mechanics.py` and `fast_mechanics.c` remain reference dependencies. The GPU library compiles `gpu_kernels.cu`, which includes `cuda_spatial_bench.cu`, which includes `cuda_mechanics.cu`. These shared CUDA files remain active even though their names came from benchmarks. Removing them would break the accepted backend.

Resumption/retiming helpers referenced by accepted protocol adapters also remain. Experimental source code whose bytes are frozen by a scientific protocol must not be rewritten just to rename its backend. A future refactor requires an explicit source transition, old-version retention, and repeat backend acceptance; it cannot alter historical results retroactively.

## Operational rule

For new mature moving experiments, use the accepted GPU runner after verifying the source hashes, full reference gate at that timestep, physical parameter restrictions, and a matched native/GPU prefix for the new context. Use native CPU for reference continuations and NumPy/SciPy for unsupported developmental features. Choose the small CPU ODE path for frozen assays. Do not revive a prototype because an isolated benchmark looks faster.

Backend acceptance, timestep sensitivity, spatial refinement, physical transport accuracy, and scientific mechanism tests answer separate questions. Their ledger entries retain separate decisions, including failures.
