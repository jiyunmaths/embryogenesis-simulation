# Reflecting boundaries on skewed meshes

The positive periodic stencil now has a separate reflecting-wall benchmark. **Eight of nine checks pass.** Conservation, positivity, and the finest-grid error tolerance pass, but the declared convergence-order check fails at intermediate shears. This is not yet a second-order boundary treatment across the tested mesh family.

## Boundary construction

The reference domain is the same affine parallelepiped x=xi, y=s xi+eta, z=zeta, now with physical reflecting walls instead of periodic identification. The candidate treatment simply omits each directional link whose endpoint lies outside the logical cube. Retained links keep the existing nonnegative symmetric conductance h times the directional weight.

There is no periodic wrap, prescribed inflow, clipping, or renormalization. With M containing cell volumes and K the retained symmetric graph stiffness, the operator remains −M⁻¹K. Constants and total amount are preserved, and off-diagonal transfer rates stay nonnegative. Those algebraic properties do **not** establish the accuracy of the physical boundary condition; that is the purpose of the independent manufactured solution.

## Nontrivial exact no-flux reference

The logical diffusion tensor is

$$
B=\begin{pmatrix}1&-s&0\\-s&1+s^2&0\\0&0&1\end{pmatrix}.
$$

Define

$$
P(\eta)=\eta^3(1-\eta)^3,\qquad
C(\xi)=\xi-3\xi^2+2\xi^3,
$$

$$
f(\xi,\eta,\zeta)=64\left[P(\eta)+sC(\xi)P'(\eta)\right]\cos(2\pi\zeta).
$$

At both xi walls, C=0 and C′=1, giving f_xi=s f_eta. Therefore the physical conormal flux, proportional to (B grad f)_xi, is zero even though the tangential gradient is generally nonzero. At the eta walls, P=P′=P″=0; at the zeta walls, the sine derivative vanishes. Thus all six physical walls have exact zero normal flux.

The manufactured trajectory is

$$
c=1+0.2e^{-t}f,\qquad
S=0.2e^{-t}\left[-f-\nabla_{\xi}\cdot(B\nabla_{\xi}f)\right].
$$

The zero mean of the cosine ensures zero total forcing and constant total amount. Polynomial antiderivatives and the analytic cosine average give **exact cell averages** for the initial state, forcing, and reference solution. An augmented sparse matrix exponential evolves the discrete forced system, removing a chosen finite time-step error from this spatial study.

Positivity is tested **separately without forcing**, by placing a positive pulse in a cell next to a wall. Manufactured forcing is not used as evidence for source-free positivity.

## Protocol and results

The predeclared study uses grids 8³, 16³, and 32³; shears 0, 0.25, 0.5, and 1; and final time 0.02. It retains the preceding prototype's 0.1% finest-error tolerance, 1.8 minimum observed order, 1e-12 mass/positivity tolerances, and additional checks on the analytic conormal data and discrete invariants. Boundary-cell errors include edges and corners.

| Shear | Finest overall relative L2 error | Finest boundary-cell relative L2 error | Last overall order |
|---:|---:|---:|---:|
| 0 | 0.017796% | 0.017915% | 1.9854 |
| 0.25 | 0.032160% | 0.045685% | 1.5825 |
| 0.5 | 0.038903% | 0.057042% | 1.5198 |
| 1 | 0.028731% | 0.031238% | 1.9913 |

All errors decrease with refinement, and all finest errors are below 0.1%. Source-free wall pulses remain nonnegative and conserve amount. The exact source has zero total integral, and the forced trajectories also conserve amount within tolerance.

The overall result remains **FAIL** because the intermediate-shear orders fall below 1.8. These are observed orders over three grids, not a proof of the eventual asymptotic rate. The result does not show that the limiting no-flux condition is necessarily wrong; it shows that this simple truncation has not met the declared second-order accuracy requirement. Conservation, low final error, or endpoint-shear success must not erase that distinction.

## Reproduction and next decision

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_boundary prepare \
  --output outputs/skew-boundary-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_boundary run \
  --output outputs/skew-boundary-repeat
```

The completed study is `outputs/skew-boundary`. Its frozen protocol records source hashes and criteria; `comparison.json` contains overall, boundary, and interior errors, conservation residuals, and pulse results. `RESULTS.md` retains the failed order check. Regression tests verify missing periodic wraps, nontrivial exact boundary data, constant equilibrium, positive wall/corner pulses, and the orthogonal reference behavior.

Next investigate boundary-local truncation error and construct a positive conservative boundary correction, then repeat the same tests and independently confirm finer-grid behavior. Do not proceed to live embryo integration on the strength of periodic accuracy alone. Nonuniform capacities, irregular faces, phase-field geometry extraction, and changing topology still require their own validation. Live transport and the running developmental-refinement simulation were not changed.

## Corrected-boundary follow-up

The [positive wall-tangential correction](skew_boundary_correction.md) now passes all nine checks, including a second wall family and 64³ confirmation at the previously failing shears. The correction is analytically derived and retains symmetric nonnegative conductances. This separate passing result does not overwrite the original failed boundary screen.
