# Discrete-to-continuum bridge: conservative transport

## Question and scope

Can a growing number of spatial compartments approximate the **same physical signaling equation**, rather than change the signaling process as the number increases?

This first bridge is a numerical experiment on a fixed three-dimensional cube. It connects conservative compartment exchange to continuum diffusion and checks the early Gierer–Meinhardt instability. It does not yet convert the deformable embryo into a continuum tissue. Compartment counts below are numerical resolutions, not biological cell counts; the experiment has no division, cell identity dynamics, or shape evolution.

Keeping a fixed domain separates spatial refinement from biological growth and cleavage. It also gives an analytic reference that the actual embryo geometry does not provide. The current embryo simulation and dashboard still use their existing normalized contact-graph model.

## Conservation supplies the missing spatial scale

Let $c_i$ be a concentration averaged over a fixed control volume $V_i$. For two compartments sharing a face of area $A_{ij}$ and with orthogonal center separation $\ell_{ij}$, define the geometric conductance

$$
g_{ij}=\frac{A_{ij}}{\ell_{ij}}=g_{ji}.
$$

With bulk diffusivity $D$, the amount entering $i$ from $j$ per time is $Dg_{ij}(c_j-c_i)$. Therefore

$$
V_i\frac{\mathrm{d}c_i}{\mathrm{d}t}
=V_i R(c_i)+D\sum_jg_{ij}(c_j-c_i).
$$

Reciprocal fluxes cancel when summed over compartments. Without reactions and boundary flux, total amount $\sum_i V_i c_i$ is conserved. Reactions can create or remove amount. This finite-volume construction assumes an orthogonal mesh and scalar bulk diffusivity; it is not a general approximation for arbitrary cell-center geometry. See [Eymard, Gallouët, and Herbin, *Finite Volume Methods*](https://raphaeleh.github.io/PUBLI/bookevol.pdf), especially the sections on admissible meshes, Neumann boundaries, and parabolic equations.

Let $M=\operatorname{diag}(V_i)$, $G=(g_{ij})$, $d_i=\sum_jg_{ij}$, and $K=\operatorname{diag}(d_i)-G$. Then

$$
\Delta_V=-M^{-1}K,
\qquad \Delta_V\mathbf{1}=0,
\qquad \mathbf{V}^{\mathsf T}\Delta_V=0.
$$

The positive symmetric operator $M^{-1/2}KM^{-1/2}$ has the same eigenvalues as $-\Delta_V$. Equivalently, physical modes solve $K\mathbf{u}=\mu M\mathbf{u}$. This replaces a degree-based inner product with a volume-based one.

For the normalized random-walk operator already used by the embryo model,

$$
\Delta_{\mathrm{rw}}=\operatorname{diag}(d_i)^{-1}G-I,
\qquad
\Delta_V=\operatorname{diag}\!\left(\frac{d_i}{V_i}\right)\Delta_{\mathrm{rw}}.
$$

The formula applies on nonisolated compartments; isolated rows are zero in both operators. It makes the required local geometry/capacity factor explicit. Normalized exchange preserves a degree-weighted sum on a fixed graph, which generally differs from molecular amount. That is consistent with its present interpretation as exchange of regulatory activities, but it cannot be reinterpreted as bulk diffusion without addressing the scaling.

In the interior of a uniform cubic mesh with spacing $h$,

$$
V_i=h^3,\quad g_{ij}=h,\quad \frac{d_i}{V_i}=\frac{6}{h^2}.
$$

Thus a fixed normalized exchange rate $\kappa$ corresponds approximately to interior diffusivity $D_{\mathrm{eff}}=\kappa h^2/6$. Halving the compartment spacing while keeping $\kappa$ fixed quarters this effective diffusivity. Boundary compartments have fewer faces, so a single scalar calibration does not exactly match all rows even at the reference resolution.

This benchmark assumes **bulk diffusion**. Membrane-limited communication instead has conductance proportional to permeability times face area, $P_{ij}A_{ij}$; extracellular diffusion requires its own accessible geometry. An inverse-square correction must not be added blindly to the embryo's diffuse contact-overlap weights, which are not measured interface areas.

## Activator–inhibitor equations and mode predictions

On the fixed mesh we retain the existing reduced kinetics:

$$
\frac{\mathrm{d}a_i}{\mathrm{d}t}
=\frac{a_i^2}{b_i}-a_i+D_a(\Delta_V\mathbf{a})_i,
$$

$$
\frac{\mathrm{d}b_i}{\mathrm{d}t}
=\beta(a_i^2-b_i)+D_b(\Delta_V\mathbf{b})_i.
$$

The intended continuum equations on the same fixed domain are

$$
\partial_t a=\frac{a^2}{b}-a+D_a\nabla^2a,
\qquad
\partial_t b=\beta(a^2-b)+D_b\nabla^2b,
$$

with zero normal flux. The equilibrium is $(a,b)=(1,1)$, and

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad r(\mu)=\max\operatorname{Re}\operatorname{eig}
\left[J-\mu\operatorname{diag}(D_a,D_b)\right].
$$

The nonlinear feedback follows the activator–inhibitor mechanism introduced by [Gierer and Meinhardt (1972)](https://www.bio.mpg.de/255219/gierer-and-meinhardt-1972). The particular normalization and numerical benchmark here are our own reduced implementation.

For a cube of side $L$ with Neumann boundaries, the spatial modes are products of cosines with indices $\mathbf{m}=(m_x,m_y,m_z)$:

$$
\mu_{\mathbf{m}}=\frac{\pi^2}{L^2}(m_x^2+m_y^2+m_z^2).
$$

For $n$ compartments along each axis, the actual finite-volume operator has the exact eigenvalues

$$
\mu^h_{\mathbf{m}}=\frac{4}{h^2}\sum_{k\in\{x,y,z\}}
\sin^2\!\left(\frac{m_k\pi}{2n}\right),\qquad 0\le m_k<n.
$$

The implementation checks these eigenpairs against its assembled sparse operator; tests also compare the complete small-cube numerical spectrum. All supported unstable modes can then be counted without a dense large-matrix eigendecomposition. Multiplicities count independent cosine modes, including permutations; the benchmark does not infer a preferred biological axis from an arbitrary basis in a degenerate eigenspace.

## Reproducible protocol

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.continuum \
  --output outputs/continuum-bridge
```

Choose a new output directory for each run. Default resolutions are $n=4,8,16,32$; all use the same unit cube, $\beta=2$, $D_a=0.02$, and $D_b=0.4$. These are illustrative length/time units, not fitted biological values and not a conversion of the embryo's default graph rates. There is no contact pruning on the ideal face-neighbor mesh.

1. **Pure diffusion:** initialize the exact cell averages of $1+0.10\cos(\pi x)\cos(\pi y)+0.06\cos(2\pi z)$. Evolve to time 1 using a sparse matrix exponential and compare with exact continuum cell averages. This suppresses time-stepping error in the spatial convergence test. The even cosine component makes nonconservation of volume-weighted amount by normalized exchange observable.
2. **Fixed-rate comparison:** use the same mesh, initial field, and boundary neighbor graph with normalized exchange. Calibrate once at $n_0=4$ using $\kappa_a=6D_a/h_0^2=1.92$ and hold that rate fixed for every finer mesh. This control asks whether unchanged per-compartment rates represent the same prescribed bulk diffusion.
3. **Spatial instability:** compare all supported finite-volume modes with the continuum band, and inspect specific modes near its edge. For these parameters, the band is $(6.492189,38.507811)$ in inverse-length-squared eigenvalues.
4. **Early growth and decay:** initialize a $10^{-5}$ perturbation in the spatial modes $(1,1,0)$ and $(2,2,0)$, aligned with each discrete chemical eigenvector. Evolve the actual nonlinear GM equations for 0.25 time units and fit the projected activator amplitude. These prescribed perturbations are verification probes, not spontaneous pattern formation. SSP-RK2 uses substeps bounded by the largest local diffusion exit rate plus reaction loss.

Predeclared checks are diffusion mass drift below $10^{-10}$, decreasing spatial error with final observed order above 1.8, measured growth error below $10^{-4}$, and agreement of both the finest-grid unstable-mode count and mode indices with the continuum. They are recorded individually; coarse or altered parameter studies can fail a check without their results being hidden. Passing these checks is not a global convergence proof.

Output files are `analysis.json` (parameters, errors, spectra, mode traces, checks), `convergence.png`, `diffusion.png`, and `diffusion_slices.npz`. Large generated outputs remain ignored by Git. The implementation uses sparse matrices and has no dense operator allocation on the large meshes.

## First results

| Compartments per axis | Total compartments | Conservative diffusion L2 error | Fixed-rate normalized L2 error | Supported unstable modes |
|---:|---:|---:|---:|---:|
| 4 | 64 | 0.002871 | 0.006980 | 10 |
| 8 | 512 | 0.000780 | 0.017603 | 10 |
| 16 | 4,096 | 0.000199 | 0.025589 | 7 |
| 32 | 32,768 | 0.00004994 | 0.027658 | 7 |

The conservative method approaches the same continuum solution with observed orders **1.88, 1.97, and 1.99**. Its largest relative diffusion-only mass drift is $1.47\times10^{-14}$. In contrast, the fixed-rate control increasingly underdiffuses on the refined mesh: its interior effective diffusivity falls from 0.02 to 0.0003125. Its mass drift decreases with refinement, but this does not rescue convergence to the prescribed nonzero diffusivity; the field increasingly retains its initial spatial variation.

The continuum predicts **seven unstable modes**. A near-edge mode, $(2,0,0)$, illustrates why merely increasing cell count is insufficient evidence:

| Resolution | Predicted growth of mode $(2,0,0)$ |
|---|---:|
| $n=4$ | +0.09139 |
| $n=8$ | +0.01505 |
| $n=16$ | −0.00697 |
| $n=32$ | −0.01267 |
| Continuum | −0.01458 |

Coarse discretizations place this mode and its two coordinate permutations inside the unstable band. Finer grids correct their classification and recover seven modes. The two independently measured nonlinear growth/decay probes agree with their **discrete** linear predictions to within $1.45\times10^{-7}$ over all resolutions. At $n=32$, the growing probe has rate 0.209446 versus the continuum value 0.209374; the decaying probe has rate −0.696076 versus −0.700784. Matching a discrete solver to its own linearization and matching that linearization to the continuum are distinct checks. The internal GM time step shrinks with the squared mesh spacing; the growth-fit errors do not establish an independent temporal convergence rate.

The plots use shared concentration limits for the field slices. They depict diffusion of a prescribed field, not emergent tissue organization. All six default checks pass, including agreement of the unstable mode indices as well as their count. This is evidence for the specified fixed-mesh transport bridge, not biological validation.

## What this enables next

The three project outcomes still need distinct transitions:

| Outcome | What this experiment establishes | Remaining bridge |
|---|---|---|
| Signaling symmetry breaking | Conservative exchange and early GM mode growth approach a specified PDE on an ideal fixed mesh | Irregular/moving geometry, nonlinear patterns, noise scaling, and calibrated communication between real cells or through extracellular space |
| Cell identity differentiation | The same GM kinetics can supply a signal in either representation | Cell-level sampling/production, persistent fate memory, and a justified coarse-grained identity description |
| Shape symmetry breaking | No new claim; geometry is held fixed for this experiment | Reciprocal signal–mechanical coupling, free boundaries, constitutive assumptions, and convergence of shape observables |

The next implementation should connect this verified transport to **irregular geometry with known face areas and capacities**, retaining an analytic or independently refined reference. A hybrid model can then keep individual cell mechanics and fate variables while solving chemical fields on an independent spatial mesh; the chemical resolution need not equal the biological cell count.

When compartments move or grow, the balance must be written for amount, $\mathrm{d}(V_i c_i)/\mathrm{d}t$, with fluxes defined relative to the moving boundaries. A material-domain continuum description includes

$$
\partial_t c+\nabla\cdot(c\mathbf{v})=\nabla\cdot(D\nabla c)+R(c).
$$

Volume change introduces dilution through this conservation law. Division must partition amount consistently with the daughters' accessible volumes, and remeshing must not create or destroy signal. Neither effect is represented by silently copying a fixed-domain Laplacian.

Before assigning conductances to the actual embryo, choose and test whether communication is bulk, membrane-limited, extracellular, or an effective regulatory interaction. The current contact-overlap graph alone does not settle this. Contact-cutoff sensitivity, mechanical refinement, independent seeds, and fate-persistence assays remain on the validation plan.
