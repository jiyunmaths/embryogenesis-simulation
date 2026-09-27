# Conservative transport on an irregular 3D domain

## Question and scope

Does the conservative transport construction still approach the same diffusion equation when the domain is nonconvex and the compartment volumes are unequal?

This extends the [fixed-cube bridge](continuum_bridge.md) to a prescribed three-dimensional L-shaped prism. It tests pure diffusion with an analytic reference and a separate probe of exchange between the prism's arms. The compartments are numerical control volumes, not biological cells. Their faces remain orthogonal, axis-aligned rectangles; this is neither an unstructured mesh nor a representation of curved cell boundaries. The embryo model and default dashboard remain unchanged.

The implementation is in [irregular.py](../embryo/irregular.py), using [masked Cartesian transport](../embryo/transport.py).

## Fixed geometry and unequal compartments

The domain is

$$
\Omega_L=\{(x,y,z)\in[0,L]^3:\ x\le L/2\ \text{or}\ y\le L/2\},
\qquad |\Omega_L|=\frac34L^3.
$$

The missing upper-right quarter extends through every z layer. Homogeneous Neumann conditions apply to all exterior faces, including the reentrant walls. There is no flux into the removed region.

For even resolution $n$, the bounding box has $n$ compartments per axis. Removing the excluded quarter leaves $3n^3/4$ active compartments. Both mesh families resolve the notch exactly at every resolution. The uniform family uses equally spaced edges. The graded family maps $s_k=k/n$ to axis edges

$$
x_k^{(\alpha)}=L\left[s_k+\frac{g_\alpha}{2\pi}\sin(2\pi s_k)\right],
\qquad (g_x,g_y,g_z)=(g,-0.6g,0.4g).
$$

The default is $g=0.35$; the implementation permits $0\le g<0.8$. Endpoint and midpoint edges are set explicitly to $0,L/2,L$. The grading law is fixed during refinement, giving unequal rectangular compartment volumes. The default successive doublings produce nested meshes; arbitrary accepted even resolutions need not be nested. No partially occupied cut cells are introduced.

## Fluxes, capacities, and the masked boundary

For neighboring active compartments, use their true shared face area $A_{ij}$, center separation $\ell_{ij}$, and volume $V_i$:

$$
G_{ij}=\frac{A_{ij}}{\ell_{ij}}=G_{ji},
\qquad
V_i\dot c_i=D\sum_jG_{ij}(c_j-c_i).
$$

With $M=\operatorname{diag}(V_i)$ and $d_i=\sum_jG_{ij}$,

$$
K=\operatorname{diag}(d_i)-G,
\qquad
\Delta_V=-M^{-1}K,
\qquad
\Delta_V\mathbf1=0,
\qquad
\mathbf V^{\mathsf T}\Delta_V=0.
$$

The active conductance submatrix is selected first, and its degree and generator are **rebuilt**. Simply deleting rows and columns from the full-box generator would retain losses toward deleted neighbors, allowing amount to leak through the notch. Rebuilding imposes zero flux at each exposed face and preserves total amount $\sum_iV_ic_i$.

This construction assumes scalar bulk diffusivity on an admissible orthogonal mesh. Its basis is described by [Eymard, Gallouët and Herbin, *Finite Volume Methods*, 2019 author-hosted text](https://raphaeleh.github.io/PUBLI/bookevol.pdf). Membrane permeability, extracellular accessibility, and arbitrary nonorthogonal cell geometry require additional modeling choices.

## Exact smooth diffusion reference

We solve

$$
\partial_t c=D\nabla^2c,\qquad
\nabla c\cdot\mathbf n=0\quad\text{on }\partial\Omega_L.
$$

Products of cosines

$$
\psi_{\mathbf m}(\mathbf x)=\prod_{\alpha\in\{x,y,z\}}
\cos\left(\frac{m_\alpha\pi x_\alpha}{L}\right),
\qquad
\mu_{\mathbf m}=\frac{\pi^2}{L^2}\sum_\alpha m_\alpha^2,
$$

satisfy these boundaries when $m_x$ and $m_y$ are even nonnegative integers and $m_z$ is a nonnegative integer. Even x/y indices make the normal derivative vanish on the reentrant walls at $L/2$. These functions are **only a subset of the L-domain's Neumann modes**. The benchmark does not compute its full spectrum or its complete set of potentially unstable reaction–diffusion modes.

The initial field is

$$
c(\mathbf x,0)=1+0.10\psi_{(2,0,0)}
+0.07\psi_{(0,2,1)}+0.04\psi_{(2,2,2)}.
$$

For each rectangular compartment, the initial and reference values are exact cell averages, not point samples. If $\bar x_{i\alpha}$ and $h_{i\alpha}$ are its center coordinates and widths, then

$$
\overline\psi_{\mathbf m,i}
=\prod_\alpha
\cos\left(\frac{m_\alpha\pi\bar x_{i\alpha}}{L}\right)
\operatorname{sinc}\left(\frac{m_\alpha\pi h_{i\alpha}}{2L}\right),
\qquad
\operatorname{sinc}(u)=\frac{\sin u}{u},\quad\operatorname{sinc}(0)=1.
$$

Each reference mode decays by $\exp(-D\mu_{\mathbf m}t)$. The numerical compartment solution is obtained with the sparse matrix-exponential action `expm_multiply`,

$$
\mathbf c_h(t)=\exp(Dt\Delta_V)\mathbf c_h(0).
$$

This isolates spatial approximation error from explicit time-stepping error, up to the accuracy of the matrix-exponential calculation. On a graded mesh, the exact continuum cell averages need not be eigenvectors of the discrete operator.

The reported error and refinement order are

$$
E_n=\left[\frac{\sum_iV_i(c_{h,i}-\overline c_i)^2}{\sum_iV_i}\right]^{1/2},
\qquad
p_n=\frac{\log(E_n/E_{2n})}{\log2}.
$$

Orders use the refinement parameter $n$; the maximum and minimum actual spacings are also recorded. For refinement ratios other than two, the denominator is $\log(n_{\mathrm{fine}}/n_{\mathrm{coarse}})$. Coarse graded meshes do not have exactly the same spacing ratios as uniform meshes. An order is reported as `null` if either error is at or below the fixed $10^{-13}$ error floor; an undefined final order does not pass the acceptance criterion.

## Default measured results

The default comparison uses $L=1$, $D=0.02$, and final time $t=1$. Values below come from `outputs/irregular-bridge/analysis.json`.

| Bounding-box resolution $n$ | Active compartments | Uniform RMS error | Uniform order from preceding row | Graded RMS error | Graded order from preceding row |
|---:|---:|---:|---:|---:|---:|
| 4 | 48 | 0.00509542 | — | 0.00392438 | — |
| 8 | 384 | 0.00139006 | 1.874 | 0.00183615 | 1.096 |
| 16 | 3,072 | 0.000354716 | 1.970 | 0.000471748 | 1.961 |
| 32 | 24,576 | 0.0000891277 | 1.993 | 0.000118619 | 1.992 |

Both families approach second-order error reduction over their finest refinement intervals. The first graded interval has order **1.10** and is preasymptotic; it is not hidden or classified as second order. The criterion below concerns the final interval only. The graded mesh is not uniformly more accurate than the uniform mesh at equal compartment count.

All eight domains have recorded volume 0.75 and one connected component. At $n=32$, the graded maximum-to-minimum compartment-volume ratio is 4.176. Across the smooth-field comparisons, the largest relative amount drift is $1.45\times10^{-14}$; the largest scaled operator residual is below $4\times10^{-16}$. These observations concern the specified fields and mesh families, not arbitrary data on a reentrant domain.

## A separate pulse checks exchange between arms

The smooth reference modes have zero normal derivative on the entire x/y midpoint planes. They therefore do not, by themselves, expose an erroneous barrier across the interior portions of those planes. A pulse provides an additional connectivity and equilibration check.

On a graded $n=8$ mesh with 384 active compartments, the initial concentration is one in the upper-left arm ($y>L/2$ within the domain) and zero elsewhere. The destination is the lower-right arm ($x>L/2$). Transport must pass through the lower-left junction around the notch.

The pulse is sampled at diffusion times $\tau=Dt/L^2\in\{0,0.05,0.5,5\}$. For the default units these correspond to physical times $0,2.5,25,250$. Initial amount is $1/4$, so the connected closed domain has uniform equilibrium concentration $1/3$ and equilibrium destination amount $1/12$.

| Diffusion time $\tau$ | Amount in destination arm | RMS departure from equilibrium, relative to its initial value |
|---:|---:|---:|
| 0 | 0 | 1 |
| 0.05 | 0.0124102 | 0.660607 |
| 0.5 | 0.0764839 | 0.0497559 |
| 5 | 0.0833333333333 | $5.25\times10^{-13}$ |

The largest pulse amount drift is $4.86\times10^{-14}$, and recorded concentrations remain nonnegative. This is a one-mesh topology/equilibration diagnostic, not an independently converged reference solution for the discontinuous pulse.

## Acceptance checks and reproduction

All ten default checks pass. The fixed acceptance criteria are:

| Check | Criterion |
|---|---|
| Domain volume | Relative error below $10^{-12}$ for every mesh |
| Connectivity | One component in every mesh and in the pulse test |
| Constant and amount invariants | Both scaled operator residuals below $10^{-12}$ |
| Smooth-field amount conservation | Relative drift below $10^{-10}$ |
| Smooth-field positivity | Minimum concentration strictly positive |
| Refinement trend | RMS errors decrease at every refinement in both families |
| Finest refinement order | Final observed order above **1.7** in both families |
| Pulse reaches the other arm | Positive destination amount at $\tau=0.05$ |
| Pulse equilibrates | Final relative equilibrium RMS below $10^{-6}$ |
| Pulse conservation and nonnegativity | Drift below $10^{-10}$ and minimum at least $-10^{-13}$ at each sampled time |

The scaled invariant residuals are $\lVert\Delta_V\mathbf1\rVert_\infty/q_{\max}$ and $\lVert\mathbf V^{\mathsf T}\Delta_V\rVert_\infty/d_{\max}$, where $q_{\max}=\max_i[-(\Delta_V)_{ii}]$ and $d_{\max}=\max_i\sum_jG_{ij}$. These normalizations account for increasing operator magnitude during refinement.

From the repository root, the default reproducer is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.irregular \
  --output outputs/irregular-bridge
```

The recorded output directory already exists in this workspace. The command refuses to overwrite an existing destination; use a fresh path such as `outputs/irregular-bridge-repeat` for a new run. The CLI reports individual checks and exits unsuccessfully if any criterion fails, while retaining the diagnostic outputs.

Outputs are `analysis.json` (parameters, errors, diagnostics and checks), `fields.npz` (mesh edges, masks and sampled fields), `convergence.png`, `fields.png`, and `mesh.png`. Generated outputs remain local and are ignored by Git. The field plots use shared concentration scales and actual unequal compartment edges.

## Interpretation and next tests

This verifies conservative bulk diffusion for selected smooth fields on a fixed nonconvex domain, including unequal capacities and correctly sealed masked faces. The imposed L shape is verification geometry, not emergent shape symmetry breaking. It does not validate nonlinear activator–inhibitor pattern formation, hybrid chemical/cell coupling, moving boundaries, cell identity persistence, or embryo mechanics.

The next useful tests are:

1. Extend irregular-domain verification to nonseparable reference fields and data that probe reentrant-corner behavior; refine the arm pulse rather than using equilibration alone as an accuracy claim.
2. The [weighted spatial spectrum and early GM mode growth](irregular_signaling.md) are now verified using complete discrete spectra. Next, test nonlinear pattern persistence and convergence with the same physical perturbation sampled across meshes. The selected cosine references alone do not determine the full instability spectrum.
3. Specify the biological communication mechanism before coupling an independent chemical mesh to cells. Check source deposition, cell-level sampling, and amount accounting together.
4. For moving or dividing compartments, verify balances for $\mathrm{d}(V_i c_i)/\mathrm{d}t$, transport relative to moving boundaries, dilution, and conservative remapping. Identity and shape outcomes then require their own persistence, perturbation, resolution, and ensemble tests.
