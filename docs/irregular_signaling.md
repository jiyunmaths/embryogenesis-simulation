# Weighted spectra and early signaling on the irregular domain

## Question and scope

Which spatial modes of the graded L-shaped domain can the Gierer-Meinhardt equations amplify, and does the nonlinear implementation reproduce their predicted early growth?

This follows the [irregular-domain diffusion benchmark](irregular_bridge.md). Geometry remains prescribed and fixed. Compartment counts are numerical resolutions, not biological cell populations. This study computes the **complete discrete spatial spectrum** at each resolution, including multiplicities, and checks small growing and decaying perturbations. It does not establish an exact continuum spectrum, sustained nonlinear patterns, stable identities, or emergent shape.

The implementation is [irregular_signaling.py](../embryo/irregular_signaling.py), using the existing conservative transport and positive SSP-RK2 signaling solver. The coupled embryo simulation and dashboard are separate from this benchmark.

## The correct weighted eigenproblem

Let $M=\operatorname{diag}(V_i)$ contain compartment volumes and $K=\operatorname{diag}(G\mathbf1)-G$ the symmetric conductance stiffness matrix. The transport operator is $\Delta_V=-M^{-1}K$. Spatial modes solve

$$
K\mathbf u_k=\lambda_k M\mathbf u_k,
\qquad \Delta_V\mathbf u_k=-\lambda_k\mathbf u_k.
$$

Diagonalize the symmetric capacity-weighted operator

$$
S=M^{-1/2}KM^{-1/2},
\qquad S\mathbf q_k=\lambda_k\mathbf q_k.
$$

A concentration mode is proportional to $M^{-1/2}\mathbf q_k$, rather than the raw symmetric eigenvector. We normalize it to unit volume-weighted RMS:

$$
\mathbf u_k=\sqrt{\sum_i V_i}\,M^{-1/2}\mathbf q_k,
\qquad \frac{\sum_iV_i u_{k,i}^2}{\sum_iV_i}=1.
$$

Nonconstant modes have zero volume-weighted mean. These capacities differ from the degree weighting in the embryo's normalized contact graph. The finite-volume formulation follows the orthogonal-mesh framework in [Eymard, Gallouët and Herbin, *Finite Volume Methods*](https://raphaeleh.github.io/PUBLI/bookevol.pdf).

## Complete spectrum without a dense 3D eigenproblem

The domain and mesh are products of an L-shaped x/y cross-section and a z interval. Conductance and volume factors give the exact discrete identity

$$
S=S_{xy}\otimes I_z+I_{xy}\otimes S_z.
$$

If $(\alpha_p,\mathbf q_p^{xy})$ and $(\nu_j,\mathbf q_j^z)$ are eigenpairs of the two factors, then

$$
\lambda_{pj}=\alpha_p+\nu_j,
\qquad \mathbf q_{pj}=\mathbf q_p^{xy}\otimes\mathbf q_j^z.
$$

Both factors are diagonalized numerically with a symmetric eigensolver; no uniform-grid cosine formula is substituted for the graded operator. All sums are sorted with their factor indices, retaining multiplicity. Only selected physical 3D modes are reconstructed.

For $n=32$, the full system has 24,576 compartments, but the dense factors have sizes 768 and 32. The CLI therefore limits this implementation to even $4\le n\le32$. This factorization depends on the extruded domain and tensor mesh; it is not a general eigensolver for arbitrary deformable tissue.

Before using the spectrum, the code verifies the factorization against the separately assembled sparse 3D operator. Tests additionally compare every eigenvalue with direct dense diagonalization of the small 3D operator and check the full capacity-weighted mode basis. [SciPy's `eigh` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html) describes the symmetric eigenproblem used for each factor.

## Gierer-Meinhardt linear stability

The concentration equations are

$$
\dot a_i=\frac{a_i^2}{b_i}-a_i+D_a(\Delta_V\mathbf a)_i,
\qquad
\dot b_i=\beta(a_i^2-b_i)+D_b(\Delta_V\mathbf b)_i.
$$

The equilibrium is $(a,b)=(1,1)$. Its local Jacobian and each spatial mode block are

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad B_k=J-\lambda_k\operatorname{diag}(D_a,D_b),
\qquad r_k=\max\operatorname{Re}\operatorname{eig}(B_k).
$$

Local stability requires $\beta>1$. With the default $(\beta,D_a,D_b)=(2,0.02,0.4)$, the local growth rate is $-0.5$. Since the spatial-block trace is negative, stationary instability occurs when its determinant is negative:

$$
\det B_k=\beta+(\beta D_a-D_b)\lambda_k+D_aD_b\lambda_k^2<0.
$$

The candidate eigenvalue band is approximately $(6.49219,38.50781)$. Only actual discrete eigenvalues in this interval are amplified. This band is in the bulk-transport spatial units of the benchmark, not the normalized contact-graph eigenvalue scale.

## Independent early-growth probes

At each resolution, select the fastest growing nonconstant mode and the real, decaying nonconstant mode whose leading rate is closest to zero. Their identities may change with resolution. Selecting modes independently does not track a single eigenvector across refinement, and vectors within degenerate eigenspaces have arbitrary basis and sign.

For a selected spatial mode $\mathbf u_k$ and Euclidean-unit leading chemical eigenvector $(v_a,v_b)$, initialize

$$
a_i(0)=1+\varepsilon v_a u_{k,i},
\qquad b_i(0)=1+\varepsilon v_b u_{k,i},
\qquad \varepsilon=10^{-5}.
$$

The nonlinear solver advances to $t=0.25$, recording at outer intervals of 0.01. Its positivity-limited internal SSP-RK2 steps depend on the actual maximum transport exit rate. Both the complete nonlinear fields and the activator projection are compared with the linear prediction:

$$
A_k(t)=\frac{\sum_iV_i u_{k,i}(a_i(t)-1)}{\sum_iV_i},
\qquad A_k^{\rm lin}(t)=\varepsilon v_a e^{r_kt}.
$$

A least-squares fit of $\log A_k(t)$ estimates growth or decay. The full two-species volume-RMS difference from the linear field is divided by $\varepsilon$, preventing a tiny imposed amplitude from trivially passing an absolute-error check. These are intentionally small perturbations; the experiment does not test saturation or late pattern selection.

Software tests compare the same nonlinear evolution with an independent `solve_ivp` DOP853 integration. They also compare mode-block predictions with all eigenvalues of a directly assembled small 3D two-species Jacobian. Equal diffusivities provide a stable control: all spatial rates are $-0.5-D\lambda_k$ for the default local kinetics.

## Measured default results

All meshes use the previous fixed grading law with $g=0.35$, length $L=1$, and identical kinetics. Results are recorded in `outputs/irregular-signaling/analysis.json`.

| Resolution $n$ | Compartments | First positive eigenvalue | Unstable spatial modes | Maximum predicted growth |
|---:|---:|---:|---:|---:|
| 4 | 48 | 5.115306 | 9 | 0.2096084 |
| 8 | 384 | 5.616789 | 7 | 0.2140618 |
| 16 | 3,072 | 5.799178 | 4 | 0.2149593 |
| 32 | 24,576 | 5.864361 | 4 | 0.2151938 |

The maximum relative change among the first twelve positive eigenvalues is 13.48%, 4.22%, and 1.11% over successive refinements. The last two meshes agree on four unstable modes. This is discrete refinement evidence; the exact full continuum spectrum of this nonconvex domain is not known from these computations. The smooth even-cosine probes in the earlier diffusion benchmark are a subset of spatial modes and omit the first positive mode reported here.

At $n=32$:

| Probe | Spatial eigenvalue | Predicted rate | Measured rate | Absolute rate error |
|---|---:|---:|---:|---:|
| Growing | 15.722499 | 0.215193833 | 0.215193831 | $1.42\times10^{-9}$ |
| Decaying | 39.339584 | -0.012482632 | -0.012482643 | $1.13\times10^{-8}$ |

Across all eight probes, the largest growth-rate error is $5.54\times10^{-7}$ and the largest relative full-field departure from the linear prediction is $5.16\times10^{-6}$. All probe concentrations remain positive. The measured rates verify the nonlinear implementation near the equilibrium against its spatially discrete prediction; they do not establish a continuum nonlinear pattern.

## Acceptance, outputs, and reproduction

All nine default checks pass:

1. Relative factorization residual below $10^{-12}$ against the assembled 3D operator.
2. Locally stable homogeneous kinetics on every mesh.
3. Both growing and decaying probes available on every mesh.
4. Relative physical-mode eigenproblem residual below $10^{-9}$.
5. Absolute measured growth-rate error below $10^{-4}$.
6. Two-species nonlinear-versus-linear field RMS, divided by initial perturbation amplitude, below $10^{-3}$.
7. Positive probe concentrations.
8. Last-refinement relative change below 3% for the first twelve positive eigenvalues.
9. Matching unstable-mode counts on the two finest meshes.

These are benchmark criteria, not rigorous continuum error bounds. A deliberately short refinement sequence can fail the last two criteria while passing growth verification; failed checks remain in the saved report. Parameters without supported unstable modes do not pass this benchmark's requirement for both probe types, even if they correctly describe a stable control.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.irregular_signaling \
  --output outputs/irregular-signaling-repeat
```

A fresh output path is required. The recorded default run is in `outputs/irregular-signaling`. The CLI exits unsuccessfully if an acceptance check fails, while preserving the report. Outputs are:

- `analysis.json`: all discrete eigenvalues, their factor indices and growth rates, selected probe histories, parameters, and checks.
- `modes.npz`: actual mesh edges, masks, volumes, selected physical modes, and final nonlinear concentration fields.
- `spectrum.png`: supported growth rates, low-spectrum refinement, and measured versus predicted early amplitudes.
- `modes.png`: selected physical modes; each panel shows its highest-energy z slab and states the slab bounds. Signs and bases within degenerate eigenspaces are arbitrary.

Generated output directories are ignored by Git. Source, tests, protocol, and these measured results are versioned separately.

The subsequent [nonlinear pattern persistence and convergence experiment](nonlinear_bridge.md) is now complete for one physical initial perturbation shared across meshes. All eleven default criteria pass: the final spatial discrepancy decreases to 1.26%, temporal refinement gives a 0.095% maximum sampled discrepancy for the smaller step pair, and the pattern changes by less than 0.41% over the sampled late window. An equal-diffusivity control returns to uniformity. Initial-condition robustness remains untested. The subsequent [moving-domain amount and dilution test](moving_domain.md) now also passes its eleven checks for prescribed affine motion. Conservative transfer under compartment changes remains next, before coupling to fate variables or mechanical feedback.
