# Face-normal-aware flux on a sheared 3D mesh

The previous geometric benchmark showed that area divided by Euclidean center distance gives incorrect normal flux on nonorthogonal contacts. This standalone prototype incorporates tangential gradients into the face-normal flux. It **passes seven accuracy/conservation/stability checks but fails positivity**, so it is not ready for live embryo signaling.

## Geometry and flux

Use a periodic parallelepiped obtained from a uniform logical cube by

$$
\mathbf x=F\boldsymbol\xi,\qquad
F=\begin{pmatrix}1&0&0\\s&1&0\\0&0&1\end{pmatrix},\qquad
B=F^{-1}F^{-T}=\begin{pmatrix}1&-s&0\\-s&1+s^2&0\\0&0&1\end{pmatrix}.
$$

The determinant is one, so each cell has volume h³. For physical diffusivity one, the signed transfer into the lower cell through its forward logical-a face is the area-vector contraction with the physical gradient. In logical coordinates it is approximated by

$$
Q_{a,i+1/2}=h^2\left[
 B_{aa}\frac{c_{i+e_a}-c_i}{h}
 +\sum_{b\ne a}B_{ab}\frac{(D_b^0c)_i+(D_b^0c)_{i+e_a}}{2}
\right],
$$

where D⁰ is the centered difference. Cell concentration changes by the signed sum of these shared face fluxes divided by h³. Every face contributes equal and opposite amounts to its two cells; periodic summation therefore conserves total amount.

Tangential derivatives use additional neighbors. This is a multipoint correction, not a replacement of center distance by its normal projection. The flux reduces to ordinary centered diffusion at zero shear. A separate affine-field patch evaluates only interior stencils, avoiding the discontinuous seam that a linear field would have under periodic identification.

## Independent checks

The immutable protocol specifies shears 0, 0.25, 0.5, and 1; grids 8³, 16³, 32³, and 64³; and final time 0.02. Manufactured solutions are sums of three Fourier modes with mixed directions. Initial values and exact references are **cell averages** on the sheared cells, including the sinc averaging factors. The physical diffusion decay of logical mode m is exp(−4π² mᵀBm t).

The discrete operator is evolved exactly in time through its Fourier symbol. Thus the transient comparisons isolate spatial discretization; they do not validate a production time integrator. The independently implemented face-flux divergence is checked against this symbol on a random field.

Predeclared checks require exact affine-field fluxes to 1e-12 per unit area, mass drift below 1e-12, nonpositive discrete eigenvalues, decreasing manufactured-solution errors, last refinement order above 1.8, finest relative L2 error below 0.1%, and no concentration below −1e-12 from a nonnegative pulse.

## Completed results

| Shear | Corrected finest relative L2 error | A/ell finest error | Last refinement order | Positive-pulse minimum |
|---:|---:|---:|---:|---:|
| 0 | 0.002543% | 0.002543% | 1.9984 | −1.08e-17 |
| 0.25 | 0.001215% | 0.722060% | 2.0002 | −0.001076 |
| 0.5 | 0.003647% | 1.434717% | 1.9841 | −0.002226 |
| 1 | 0.014023% | 2.311342% | 1.9721 | −0.004433 |

All affine face tests pass, including tangential gradients that previously generated spurious normal flux. Conservation, agreement between face and spectral operators, spectral stability, and smooth-solution refinement also pass.

The pulse starts at concentration one in one 16³ cell and zero elsewhere, and is evaluated at time 0.01h². Negative concentrations for every nonzero tested shear substantially exceed roundoff. The centered mixed-derivative stencil has negative off-diagonal transfer rates. A symmetric nonpositive spectrum makes this diffusion operator stable in the discrete L2 norm, but does **not** make it positivity preserving.

Because the negative pulse is obtained from the exact discrete exponential, changing the time step cannot make this spatial operator a positivity-preserving generator. Clipping the concentration would conceal the defect and generally break amount conservation. Smooth positive Fourier solutions alone would have missed it.

## Decision and next work

The face-normal correction resolves the demonstrated consistency problem on this mesh family, but cannot be directly substituted into the live positive activator–inhibitor solver. The next prototype should preserve nonnegative transfer rates, for example through a suitable wider stencil on known skewed meshes, and repeat the same linear-field, smooth-refinement, conservation, and pulse tests. Accuracy, positivity, and available neighboring geometry must be assessed together rather than assuming all follow from symmetry.

This benchmark has exact faces, normals, centroids, and volumes. It has periodic boundaries, affine shear, fixed connectivity, and no reactions. Free boundaries, unstructured deformed cells, geometry extraction from phase fields, and topology changes remain unvalidated. The existing embryo model and its running developmental-refinement code were not modified.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_flux prepare \
  --output outputs/skew-flux-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.skew_flux run \
  --output outputs/skew-flux-repeat
```

The completed study is `outputs/skew-flux`. Its `protocol.json`, `comparison.json`, `RESULTS.md`, and `status.json` retain all checks, including the failed positivity decision. Five prototype regression tests cover affine fluxes, face/spectral agreement, conservation, refinement, and detection of the negative pulse; the existing geometric-closure tests also remain in the verification set.

## Positive-stencil follow-up

The [nonnegative directional prototype](positive_skew_flux.md) now passes all eleven checks on the same periodic sheared-mesh family, including exact and SSP-RK2 pulse positivity. This is a separate spatial operator; the centered prototype’s failed positivity result remains unchanged. Boundary and irregular-geometry validation are still required.
