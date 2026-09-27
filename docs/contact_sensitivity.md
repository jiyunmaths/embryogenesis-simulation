# Contact-cutoff sensitivity on archived geometries

The independent 112³ refinement confirmation passed all eleven checks. The next screen tests another numerical choice: discarding weak shell overlaps below a fraction of the largest overlap. It uses the live conservative contact adapter, keeping geometry, cell volumes, interface width, diffusion coefficients, and reaction parameters fixed.

## Protocol

Seven archived states are tested: the boundary-cleared 56³ developmental trajectory at times 30, 60, 90, 120, 150, and 180, plus the refined 112³ manufactured state at time 0.6. The default cutoff is 0.02. Cutoffs 0.01, 0.02, and 0.04 form the acceptance neighborhood; 0, 0.005, and 0.08 are wider stress diagnostics. The protocol and checkpoint hashes are saved before analysis.

Every local comparison must retain the component count and number of unstable spatial modes, change the conductance matrix and sorted spectrum by less than 5% in relative Frobenius/L2 norm, and change each sorted spectral growth rate by less than 0.01 per model-time unit. Constant-preservation and amount-conservation residuals must be below 1e-10 for every cutoff. These are declared numerical sensitivity screens, not biologically calibrated tolerances. Eigenvalues are sorted; this does not track eigenvector identity through crossings.

Reproduce with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.contact_sensitivity \
  --output outputs/contact-cutoff-repeat \
  --checkpoints outputs/domain-conservative/grid-56/final_state.npz \
  outputs/domain-conservative-extended/milestones/t-{60,90,120,150,180}/grid-56/state.npz \
  outputs/refinement-confirmation/space-112/final_state.npz
```

## Results

The overall screen **fails**, because four developmental snapshots exceed the declared growth-rate-change threshold. All other local checks pass, and every source checkpoint remains unchanged.

| Developmental time | Unstable spatial modes | Largest local change in any sorted growth rate | Largest local change in fastest spatial growth rate |
|---:|---:|---:|---:|
| 30 | 5 | 0.02852 | 0.000839 |
| 60 | 1 | 0.08581 | 0.001583 |
| 90 | 0 | 0.01137 | 0.001392 |
| 120 | 0 | 0.05255 | 0.001379 |
| 150 | 0 | 0.00432 | 0.000678 |
| 180 | 0 | 0.00418 | 0.000750 |

All tested graphs remain connected and their unstable-mode counts are unchanged across the entire stress range. Within the local neighborhood, conductance changes are at most 1.10%, and spectral L2 changes are at most 0.90%. The 112³ manufactured state passes every local check and retains 14 unstable modes.

The fastest-rate column is an explanatory follow-up computed from the exported `maximum_spatial_growth`; it does not replace the prospectively specified all-rates check. The larger changes in other spectral rates are not evidence of a changed leading instability. Conversely, unchanged instability counts do not establish quantitative signaling robustness.

Full results are in `outputs/contact-cutoff-sensitivity/{protocol.json,comparison.json,RESULTS.md,status.json}`.

## Scope and next decision

These are homogeneous-equilibrium linearizations on frozen sampled geometries. They are not linearizations around the actual heterogeneous signaling state and omit geometry motion, dilution, and changing eigenvectors. The six developmental snapshots also retain their coarse-cell-resolution limitation. This screen does not validate the physical overlap-to-contact-area closure.

Next, run matched evolving-geometry branches with cutoffs 0.01, 0.02, and 0.04 from the same pre-amplification checkpoint, recording signal contrast, fate labels, axis ratio, and numerical quality. First isolate the signaling cutoff from the polarity-alignment cutoff: the current shared configuration field also affects polarity neighbor averaging. Declare whether a branch changes only signaling or both pathways; a shared-cutoff experiment cannot attribute differences solely to molecular transport. Preserve the frozen-screen failure and compare actual trajectory sensitivity before changing the default cutoff. Cleavage and long-time spatial convergence remain separate requirements.

The [dynamic follow-up](cutoff_dynamics.md) is implemented and running. Its independent polarity-cutoff override preserves historical defaults and isolates the direct signaling-filter intervention.
