# Positive wall correction for the skew-mesh prototype

The first reflecting-wall benchmark preserved positivity and conservation but missed its 1.8 convergence-order criterion at intermediate shears. This follow-up derives a boundary correction rather than tuning a parameter to the measured errors. The original failed report remains unchanged.

## Leading boundary defect

Write t=|s| for shear magnitude, and retain the nonnegative directional stencil in the interior. At an xi wall, deleting outgoing diagonal links leaves a leading tangential second-derivative defect. For s≥0, the conormal condition is c_xi=s c_eta. Differentiating along the wall and expanding the boundary-cell operator gives

$$
(L_hc-Lc)_{\xi\text{ wall}}
=-\frac{s(1-s)}{2}\,c_{\eta\eta}+O(h).
$$

The corresponding eta-wall condition is (1+s²)c_eta=s c_xi, giving

$$
(L_hc-Lc)_{\eta\text{ wall}}
=-\frac{s}{2}\left(1-\frac{s}{1+s^2}\right)c_{\xi\xi}+O(h).
$$

The same coefficients use |s| for negative shear. Add tangential conductances only within the boundary-cell layers:

$$
\gamma_{\xi}=\frac{t(1-t)}2,\qquad
\gamma_{\eta}=\frac t2\left(1-\frac{t}{1+s^2}\right),\qquad
G_{\rm added}=h\gamma.
$$

Xi-wall layers gain eta-direction links; eta-wall layers gain xi-direction links. Existing conductances are retained. For |s|≤1, both corrections are nonnegative and symmetric. They add no exterior link and no periodic wrap. Volume capacities remain h³, total amount remains conserved, constants remain stationary, and the diffusion generator retains nonnegative off-diagonal rates. At zero shear the correction vanishes exactly. Zeta walls need no mixed-derivative correction for this shear map.

This boundary-row expansion motivates cancellation of the leading defect. It is not, by itself, a proof of global convergence on arbitrary geometries or corners. The independent refinement checks below test the resulting method.

## Two wall families and finer confirmation

The original xi-wall reference is retained. A second, prospectively specified eta-wall reference uses

$$
f_{\eta}=64\left[P(\xi)+\frac{s}{1+s^2}C(\eta)P'(\xi)\right]\cos(2\pi\zeta),
$$

with P(u)=u³(1−u)³ and C(u)=u−3u²+2u³. This has generally nonzero tangential gradients at the eta walls and exact homogeneous conormal conditions on all walls. It prevents success from depending solely on the original solution's inactive eta-wall derivatives.

Both families use c=1+0.2 exp(−time) f and the exact manufactured source 0.2 exp(−time)(−f−Lf). Exact cell averages are evaluated analytically, and an augmented matrix exponential removes finite-step error from the comparison. Source-free pulses on each skewed wall and at a corner independently test positivity and conservation.

The protocol retains grids 8³, 16³, and 32³; shears 0, 0.25, 0.5, and 1; final time 0.02; the 0.1% error tolerance; the 1.8 minimum order; and the existing mass/positivity tolerances. Before evaluating results, it also specifies **64³ confirmation runs for both wall families at shears 0.25 and 0.5**, where the original method failed. No tolerance or conductance coefficient is fitted after observing these results.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_boundary_correction prepare \
  --output outputs/skew-boundary-correction-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_boundary_correction run \
  --output outputs/skew-boundary-correction-repeat
```

The study is in `outputs/skew-boundary-correction`. Protocol and source hashes are fixed before execution; the final report records both the original-grid and finer-grid checks. The uncorrected boundary report, periodic prototypes, live transport, and developmental runs remain unchanged.

## Scope

The correction is derived for a uniform affine parallelepiped with constant diffusion and |s|≤1. It does not yet establish consistency for changing geometry, unequal cell capacities, curved or irregular walls, or arbitrary corner behavior. A corner pulse tests conservation and positivity, not convergence for general corner singularities. These restrictions remain even if all declared manufactured-solution checks pass.

## Completed outcome

All **nine declared checks pass**. Across all eight shear/wall-family combinations, the 16³-to-32³ orders are 1.985–1.993. The four independent 64³ confirmations have orders 1.997 (rounded), with errors below 0.0052%. Mass, constant preservation, nonnegative rates, and source-free face/corner pulses all pass. No clipping or renormalization was used.

| Shear | Wall family | 64³ relative L2 error | 32³–64³ order |
|---:|---|---:|---:|
| 0.25 | xi | 0.004602% | 1.9968 |
| 0.25 | eta | 0.004561% | 1.9968 |
| 0.5 | xi | 0.005112% | 1.9972 |
| 0.5 | eta | 0.004845% | 1.9974 |

The boundary accuracy issue is resolved for these tested uniform skew meshes by a separate derived correction. The earlier uncorrected failure remains preserved. Twenty-three regression tests pass. Next validate nonuniform capacities and smoothly graded skew geometry before attempting arbitrary live-cell interfaces.
