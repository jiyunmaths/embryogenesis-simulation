# Validation of the geometric transport approximation

The live conservative operator preserves constants and molecular amount. This does not by itself establish that its conductances represent the right geometry or flux. This benchmark separates numerical integration from the closure

$$
\widehat A_{ij}=\frac{6\sqrt{2}}{\epsilon}W_{ij},\qquad
G_{ij}=\frac{\widehat A_{ij}}{|\mathbf c_j-\mathbf c_i|}.
$$

**Outcome: the closure passes its flat, orthogonal reference tests but fails the declared general-geometry accuracy checks.** No simulation parameter or live transport code was changed. Ongoing developmental refinement measures numerical sensitivity of this existing approximation; it cannot validate a more accurate physical transport law.

## Reference problems and recorded decisions

The protocol is written before execution in `outputs/geometry-transport-validation/protocol.json`, with a transport source hash. The cases use independent, explicitly constructed geometric references:

1. Complementary planar equilibrium profiles with unit transverse area, evaluated at interface widths 0.02, 0.05, and 0.1. Translational symmetry reduces the 3D integral exactly to one dimension. The inferred area must differ from one by less than 1e-8.
2. Parallel profiles separated by gaps from 0.25 to eight interface widths. For the **sharp direct-contact interpretation**, any positive physical gap has zero shared-face area; leakage is measured against the unit planar-area scale, with a 1% screen. A separate fixed physical gap of 0.1 is tested while reducing interface width.
3. Complementary spherical profiles with radius 0.5, at width/radius ratios 0.05, 0.1, 0.2, and 0.3. Radial quadrature gives the diffuse overlap integral; 3D voxel integration on grids 64³, 96³, and 128³ checks its numerical evaluation. Sharp reference area is exactly 4πR². The declared physical-area tolerance is 1%.
4. Adjacent sheared prisms with known unit face area, volumes, and sharp centroids. Their displacement is d=(1,s,0), while the common face normal is n=(1,0,0). Shears 0, 0.25, 0.5, and 1 are tested with normal, tangential, and mixed linear concentration gradients. With diffusivity one, exact transfer into the left cell is A∇c·n. Flux error must be below 0.01 in these unit tests.

The spherical inner and outer compartments have coincident centroids. Consequently this case tests **area only**; the live adapter correctly rejects their full two-point conductance. No artificial centroid separation is used to disguise that limitation. The prism test supplies exact sharp centroids deliberately to isolate the A/ell approximation from any additional diffuse-centroid error.

## Curvature and width dependence

For radial signed distance s=r−R, the area estimate includes the Jacobian (R+s)². Even with perfectly integrated equilibrium profiles, the diffuse estimator therefore differs from sharp surface area.

| Interface width / radius | Area bias relative to sharp sphere | Finest voxel discrepancy relative to diffuse integral |
|---:|---:|---:|
| 0.05 | 0.1612% | <0.0001% |
| 0.10 | 0.6449% | <0.0001% |
| 0.20 | 2.5797% | 0.0001% |
| 0.30 | 5.8039% | 0.0158% |

The final column includes both voxel integration error and the finite cube's omitted diffuse tail; it is not a pure quadrature-order estimate. All cases pass the 0.1% numerical-reference tolerance, while ratios 0.2 and 0.3 fail the 1% sharp-area tolerance. The small-width curvature bias scales approximately as (epsilon/R)² and decreases when interface width is reduced. Refining voxels at fixed width converges to the biased diffuse integral.

## Gaps and the meaning of communication

| Gap / interface width | Inferred area / unit planar area |
|---:|---:|
| 0 | 1.0000 |
| 0.25 | 0.9876 |
| 0.5 | 0.9515 |
| 1 | 0.8218 |
| 2 | 0.4742 |
| 4 | 0.0780 |
| 8 | 0.000682 |

Diffuse shell overlap remains positive without sharp contact. In this two-cell test, the default relative cutoff cannot suppress it: the only edge is also the largest edge, however weak it becomes. Leakage decreases with interface width at fixed physical gap, but remains a model assumption at any chosen finite width.

This is a failure **if G is interpreted as transport through a shared cell face**. It is not proof that cells cannot communicate across extracellular gaps. Such communication requires a separately specified extracellular or distance-dependent transport law; overlap alone has not been calibrated to that process.

## Nonorthogonal contacts: linear-field consistency failure

For the sheared-prism geometry, a linear concentration field with gradient b gives

$$
Q_{\rm exact}=A\,\mathbf b\cdot\mathbf n,\qquad
Q_{\rm model}=\frac{\widehat A}{|\mathbf d|}\,\mathbf b\cdot\mathbf d.
$$

At shear one, flat-area calibration is effectively exact, yet:

| Gradient | Exact transfer | Model transfer |
|---|---:|---:|
| Normal (1, 0, 0) | 1 | 0.7071 |
| Tangential (0, 1, 0) | 0 | 0.7071 |
| Mixed (1, 1, 0) | 1 | 1.4142 |

The tangential case is decisive: equal-and-opposite graph exchange conserves amount but still gives the wrong flux. This persists independently of voxel quadrature. Replacing Euclidean distance by projected normal distance alone does not solve it—the tangential test would then produce transfer 1 rather than zero. Nonorthogonal bulk diffusion needs a gradient reconstruction or another consistent flux discretization, together with geometric face information.

If concentrations instead represent well-mixed cell compartments linked by phenomenological exchange rates, this continuum linear-field test is not its defining physical requirement. In that interpretation G must be calibrated as a compartment coupling and must not be presented as a generally consistent continuum bulk-diffusion operator.

## Consequence and next implementation decision

Four numerical/reference checks and two width-trend checks pass. All three general-geometry closure checks fail. The source code and active developmental experiment remain unchanged; results are in `outputs/geometry-transport-validation/{comparison.json,RESULTS.md,status.json}`.

Before extending the continuum claim, distinguish communication through sharp contacts, membrane-limited exchange, and extracellular signaling. For the project's bulk-diffusion bridge, the next bounded task is a face-normal-aware flux prototype on known skewed meshes, checked against these linear-field tests and manufactured diffusion solutions. It must also preserve amount, stability, and positivity or explicitly handle their tradeoffs. Then validate face-area/normal extraction from live phase fields, including curvature and gaps, before coupling it to embryo mechanics. Neither cutoff tuning nor mesh refinement alone resolves the demonstrated closure failures.

Reproduce with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.geometry_transport_validation prepare \
  --output outputs/geometry-transport-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.geometry_transport_validation run \
  --output outputs/geometry-transport-repeat
```

## Prototype follow-up

The [centered face-normal-aware prototype](skew_flux.md) now corrects affine fluxes and recovers approximately second-order smooth diffusion on known sheared meshes. Its nonnegative-pulse test fails despite spectral stability and conservation. It is not used in embryo simulations; a positivity-preserving alternative must be tested next.
