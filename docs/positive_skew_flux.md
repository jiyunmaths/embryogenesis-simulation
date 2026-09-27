# Positive conservative diffusion on sheared meshes

The centered face-normal correction fixed the nonorthogonal flux error but generated negative concentrations. The new wider-stencil prototype **passes all eleven declared checks**, including nonnegative pulses, while preserving second-order spatial and temporal accuracy on the tested affine mesh family. It remains a standalone prototype; live embryo transport is unchanged.

## Nonnegative directional decomposition

For the same shear map x=xi, y=s xi+eta, z=zeta, the logical diffusion tensor is

$$
B=\begin{pmatrix}1&-s&0\\-s&1+s^2&0\\0&0&1\end{pmatrix}.
$$

Let t=|s| and d=(1,−sign(s),0). For |s|≤1,

$$
B=(1-t)e_xe_x^T+(1+s^2-t)e_ye_y^T+e_ze_z^T+t\,dd^T.
$$

All four weights are nonnegative. At zero shear, the diagonal direction has zero weight and its sign is immaterial. Each directional second difference contributes

$$
(Lc)_i=\sum_d\frac{w_d}{h^2}
\left(c_{i+d}+c_{i-d}-2c_i\right).
$$

This is symmetric conservative exchange with nonnegative off-diagonal rates. Constants are stationary, total amount is conserved, and a nonnegative initial state remains nonnegative under exact discrete evolution. The uniform cell volume is h³, so each virtual directional link has conductance h w_d. These diagonal links are a numerical stencil on this manufactured mesh; they are not assumed to represent direct molecular contact between arbitrary cells.

The implementation rejects |s|>1 rather than silently retaining a negative weight. More general tensors or meshes may require different or more distant neighbor directions; this decomposition is not a universal unstructured-cell method.

## Face flux reconstruction

Each diagonal exchange is divided equally between two grid paths: x-then-y and y-then-x. Accumulating these paths on logical faces gives a single shared flux per face. Its divergence exactly equals the directional operator, so face routing preserves conservation. Affine concentration fields recover the physical face-normal flux on all three face orientations, including negative-shear regression cases.

This is a change to the spatial stencil, not concentration clipping or a smaller time step applied to the failed centered operator. The required diagonal neighbors exist on the periodic reference grid. Their availability and interpretation on a changing cell-contact graph are still unvalidated.

## Positive time integration

The forward-Euler update is a convex combination when

$$
\delta t\le\frac{h^2}{2\sum_d w_d}.
$$

The prototype uses SSP-RK2 assembled directly from these convex-combination Euler stages. Step sizes above the bound are rejected. There is no clipping, concentration floor, or mass renormalization. The nonnegative-rate property provides a structural positivity argument; pulse tests also verify the implementation.

## Completed benchmark

The same periodic geometry, manufactured cell-average solutions, grids 8³/16³/32³/64³, shears 0/0.25/0.5/1, and spatial final time 0.02 are used as for the centered prototype. Original accuracy, flux, conservation, and positivity tolerances are unchanged. Exact Fourier evolution isolates spatial error. An independent three-step-size SSP-RK2 study on 16³ uses final time 0.01 and compares with the exact discrete solution.

| Shear | Finest spatial relative L2 error | Last spatial order | Last temporal order |
|---:|---:|---:|---:|
| 0 | 0.002543% | 1.9984 | 2.0332 |
| 0.25 | 0.002833% | 1.9981 | 2.0321 |
| 0.5 | 0.002993% | 1.9974 | 2.0269 |
| 1 | 0.002744% | 1.9962 | 2.0176 |

All eleven checks pass:

- Nonnegative directional rates and exact affine face fluxes.
- Agreement between directional, reconstructed-face, and Fourier operators.
- Mass drifts below 1e-12 and nonpositive spectra.
- Decreasing spatial errors, last spatial orders above 1.8, and finest errors below 0.1%.
- Nonnegative exact and numerical pulse evolution within the original 1e-12 tolerance.
- Last temporal orders above 1.8.

The exact Fourier pulse minima are between approximately −6e-17 and zero, consistent with roundoff. SSP-RK2 pulse minima are nonnegative without correction. The exact pulse is evaluated at 0.01h², as before; an additional numerical pulse runs for ten steps at 90% of the positivity bound. The previous centered prototype's minimum of −0.00443 is therefore absent for a structural reason, not hidden by a relaxed threshold.

## Reproduction and scope

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.positive_skew_flux prepare \
  --output outputs/positive-skew-flux-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.positive_skew_flux run \
  --output outputs/positive-skew-flux-repeat
```

Results are in `outputs/positive-skew-flux/{protocol.json,comparison.json,RESULTS.md,status.json}`. Source hashes and fixed criteria are recorded before execution. Regression tests cover both shear signs, tensor reconstruction, face routing, nonnegative rates, the time-step bound, and complete reporting.

Next test **reflecting/no-flux boundaries on skewed domains**, followed by irregular geometry and nonuniform capacities. Simply deleting diagonal links at a wall may preserve positivity and conservation while imposing the wrong boundary flux, so boundary consistency needs its own manufactured solution. Validate these before extracting face geometry from phase fields or coupling the prototype to reactions and changing embryo topology. Passing this periodic benchmark does not resolve curvature/gap errors in the existing live contact-area approximation.

## Reflecting-wall follow-up

The [no-flux boundary benchmark](skew_boundary.md) is complete. Simple exterior-link truncation passes conservation, positivity, and the finest-error tolerance but fails the declared convergence-order target at shears 0.25 and 0.5. Boundary correction/refinement is therefore the next task; the periodic passing result does not establish boundary accuracy.
