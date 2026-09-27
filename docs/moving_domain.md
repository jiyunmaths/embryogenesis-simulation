# Moving-domain amount balance, dilution, and transport

## Question and scope

Can the chemical transport implementation follow changing geometry without creating or destroying material, and does its solution approach a known moving-domain continuum solution under independent spatial and temporal refinement?

This follows the [fixed-domain nonlinear signaling experiment](nonlinear_bridge.md). The new benchmark evolves **passive chemical amounts on a prescribed moving 3D L-shaped domain**. It tests uniform expansion, anisotropic expansion, and volume-preserving anisotropic deformation. Material compartments move continuously with the domain; their neighbor connectivity stays fixed. There are no reactions, cell identities, cleavage, remeshing, or mechanically generated shape changes in this test.

Implementation: [embryo/moving.py](../embryo/moving.py). Independent verification: [tests/test_moving.py](../tests/test_moving.py). This module does not change the coupled embryo solver or dashboard dynamics.

## Conservation law on a moving domain

For a concentration $c(\mathbf x,t)$ carried by material velocity $\mathbf v$, constant diffusivity $D$, and no reactions,

$$
\partial_t c+\nabla\cdot(c\mathbf v)=D\nabla^2c.
$$

Equivalently,

$$
\frac{\mathrm D c}{\mathrm D t}=D\nabla^2c-c\nabla\cdot\mathbf v.
$$

The final term is dilution during expansion and concentration during compression. Diffusive flux is zero relative to the moving material boundary. Advection and dilution arise from conservation on the evolving material domain; see [Crampin, Gaffney, and Maini (1999), *Reaction and Diffusion on Growing Domains: Scenarios for Robust Pattern Formation*](https://people.maths.ox.ac.uk/maini/PKM%20publications/109.pdf). The three-dimensional affine verification below is our own specialization, rather than a reproduction of that paper's pattern experiments.

The reference domain $\Omega_0$ is the same L-shaped prism as in the [irregular transport benchmark](irregular_bridge.md). Its material coordinates are $\mathbf X$. We prescribe

$$
x_d=s_d(t)X_d,\qquad s_d(t)=e^{g_dt},\qquad
J(t)=\prod_{d=1}^{3}s_d(t),\qquad v_d=g_dx_d.
$$

Thus $\nabla\cdot\mathbf v=\sum_d g_d$ and each compartment volume is $V_i(t)=J(t)V_i(0)$. With no spatial variation, the exact concentration is

$$
c(t)=\frac{c(0)}{J(t)}.
$$

Preserving a uniform concentration during expansion would preserve neither the material amount nor this equation.

## Amount-based finite volumes

The evolved variable is the compartment amount

$$
q_i(t)=V_i(t)c_i(t).
$$

For a face normal to axis $d$, its area scales as $J/s_d$ and its center separation as $s_d$. The finite-volume conductance therefore scales as

$$
G_{ij}(t)=\frac{J(t)}{s_d(t)^2}G_{ij}(0).
$$

Let $K_d$ be the positive stiffness matrix constructed from only the reference faces normal to axis $d$, and $M_0=\operatorname{diag}(V_i(0))$. The amount equation becomes

$$
\dot{\mathbf q}=B(t)\mathbf q,\qquad
B(t)=-D\sum_d s_d(t)^{-2}K_dM_0^{-1}.
$$

Each matrix has zero column sum, so total amount is conserved. Concentrations are recovered by dividing by the **current** capacities. This automatically includes dilution; it does not require a separately discretized dilution source term. All directions are advanced together, with no directional time splitting.

This simplification requires a globally affine, axis-aligned material deformation. It is not a general formula for nonuniform local cell growth, arbitrary mesh motion, or changing contacts. For a mesh that does not follow the material, transport relative to moving faces requires additional advection terms.

## Time integration and positivity

The implementation uses nonautonomous SSP-RK2 (Heun):

$$
\mathbf q^{(1)}=\mathbf q^k+\delta t B(t_k)\mathbf q^k,
$$

$$
\mathbf q^{k+1}=\frac12\mathbf q^k+
\frac12\left[\mathbf q^{(1)}+\delta t B(t_{k+1})\mathbf q^{(1)}\right].
$$

The second stage evaluates the geometry at $t_{k+1}$. A conservative exit-rate bound accounts for both stage times, including the increasing transport rates on compressed axes. The benchmark requires $\delta t$ times that bound to be at most 0.8. It refuses unsafe settings rather than adding hidden substeps or clipping concentrations. Every step is checked for positivity and finiteness.

Amounts are conserved at the algebraic update level. Consequently, an extremely small amount drift verifies conservation but does not alone establish accurate transport. Separate exact-solution and time-refinement tests are essential.

## Exact smooth moving-domain reference

Use the same physical initial concentration at every resolution:

$$
c(\mathbf X,0)=1+\sum_{m=1}^{3}\alpha_m
\prod_{d=1}^{3}\cos(k_{md}\pi X_d/L),
$$

with $(\alpha_m,\mathbf k_m)$ equal to $(0.10,(2,0,0))$, $(0.07,(0,2,1))$, and $(0.04,(2,2,2))$. Even x/y indices enforce no flux on the reentrant walls as well as the outside walls. These selected smooth functions are a verification family, not the complete spectrum of the L-shaped domain.

Define the axis-dependent diffusion clocks

$$
T_d(t)=\int_0^t s_d(u)^{-2}\,\mathrm du=
\begin{cases}
\dfrac{1-e^{-2g_dt}}{2g_d},&g_d\ne0,\\
t,&g_d=0.
\end{cases}
$$

The exact concentration in material coordinates is

$$
c(\mathbf X,t)=\frac{1}{J(t)}\left[
1+\sum_m\alpha_m
\exp\!\left(-\frac{D\pi^2}{L^2}\sum_d k_{md}^2T_d(t)\right)
\prod_d\cos(k_{md}\pi X_d/L)\right].
$$

The finite-volume reference uses exact averages over each moving rectangular compartment, using the reference-cell cosine averages and the dilution factor. It does not compare point samples with compartment averages. Independent tests integrate the physical-space expression with Gaussian quadrature and the diffusion clocks with adaptive quadrature.

Spatial errors are volume-weighted RMS differences divided by the RMS of the exact concentration field. **This normalization uses the full concentration, including its mean**, unlike the contrast-based normalization in the nonlinear persistence experiment. Uniform affine scaling multiplies every capacity by the same Jacobian, so reference-volume and current-volume normalized weights are identical.

## Protocol and acceptance criteria

All quantities are dimensionless. The default experiment uses $D=0.02$, $L=1$, mesh grading 0.35, and final time 2. Each motion advances uniform and patterned initial fields together.

| Motion | Rates $(g_x,g_y,g_z)$ | Physical interpretation |
|---|---|---|
| Isotropic | $(0.10,0.10,0.10)$ | All lengths expand equally |
| Anisotropic | $(0.15,0.05,0.10)$ | Shape proportions change while volume grows |
| Volume preserving | $(0.12,-0.12,0)$ | One axis expands, one contracts, total volume stays fixed |

Spatial refinement uses $n=4,8,16,32$ (48 to 24,576 active compartments), with the same $\delta t=0.002$. Each run saves 21 snapshots. Amount drift, uniform dilution error, and concentration minima are monitored **at every integration step**.

Temporal order is measured independently on $n=8$ with $\delta t=0.04,0.02,0.01$, against adaptive DOP853 integration of the same semidiscrete amount equations with `rtol=1e-12`, `atol=1e-14`. A separate finest-mesh control halves $\delta t$ to 0.001 for every motion, ensuring that the remaining spatial error is much larger than the measured step-halving difference.

The intentionally incorrect control keeps $c=1$ during isotropic expansion and calculates its apparent amount gain. This is an analytic diagnostic of omitted dilution, not a second valid physical model.

Eleven acceptance checks are recorded in the report:

1. Relative total-amount drift below $10^{-10}$ at every step, including fine time controls.
2. Maximum uniform dilution error $|Jc-1|$ below $10^{-11}$.
3. Strictly positive concentrations throughout these positive-initial-data runs.
4. Moving operators and capacities agree with a fresh physical-geometry rebuild to relative residual below $10^{-11}$.
5. Spatial errors decrease at every refinement for all three motions.
6. The last spatial-refinement order exceeds 1.7 for each motion.
7. Final finest spatial RMS error is below 0.1% of the exact full-field RMS.
8. Both temporal orders lie between 1.8 and 2.3 for each motion.
9. Finest temporal-control error against DOP853 is below $10^{-5}$.
10. The finest-mesh step-halving difference is below 1% of that mesh's continuum-reference error.
11. Omitting dilution produces more than 10% artificial amount gain in the default expansion.

The final threshold ensures that this particular diagnostic runs long enough to expose a substantial error; shorter valid experiments can fail it. Failed criteria are retained and the CLI exits unsuccessfully. Thresholds are verification choices, not biological criteria.

## Measured results

The complete default run is saved in `outputs/moving-domain-verified`. **All eleven checks pass.**

| Motion | Relative RMS, $n=4$ | Relative RMS, $n=32$ | Last spatial order | Temporal orders |
|---|---:|---:|---:|---|
| Isotropic | 0.00391205 | 0.000103667 | 1.9984 | 2.0183, 2.0090 |
| Anisotropic | 0.00393373 | 0.000106894 | 1.9976 | 2.0192, 2.0094 |
| Volume preserving | 0.00376082 | 0.000102247 | 1.9984 | 2.0228, 2.0112 |

The largest finest-mesh relative RMS error is **0.01069%**. The first coarse refinement has order approximately 1.21–1.25; the approximately second-order result applies to the better-resolved refinements, not the entire resolution range.

Across all spatial runs and finest-mesh time controls, maximum relative amount drift is $2.16\times10^{-14}$, and maximum uniform dilution error is $3.33\times10^{-16}$. At $t=2$, both expanding motions have volume ratio 1.8221188 and uniform concentration 0.548811636. Volume-preserving deformation keeps the uniform concentration at one even as the shape proportions change. The deliberately incorrect no-dilution control gains **82.21%** apparent chemical amount.

The finest-mesh step-halving differences are $6.78\times10^{-9}$, $6.33\times10^{-9}$, and $7.61\times10^{-9}$ for the three motions. They are less than **0.008% of the corresponding spatial errors**, supporting that the measured spatial refinement is not limited by the chosen time step.

Independent software tests additionally verify pure advection with zero diffusion, the stationary-domain limit, fluxes calculated from rebuilt moving face geometry, positivity rejection without mutation, and second-order agreement with an adaptive solver whose geometry is rebuilt at every RHS evaluation.

## Reproduce and inspect

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.moving \
  --output outputs/moving-domain-repeat

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q
```

Use a fresh output directory; existing experiments are never overwritten. Outputs are:

- `analysis.json`: parameters, limits, histories, convergence errors, time controls, incorrect-control results, and acceptance checks.
- `fields.npz`: reference mesh edges/capacities, masks, times, axis scales, patterned concentrations, and exact final fields. Physical edges at each frame are reference edges multiplied by the corresponding scale. This is verification data, not an embryo restart checkpoint.
- `convergence.png`: uniform dilution, spatial refinement, and independent temporal refinement.
- `balance.png`: conservative amount drift and the incorrect control's artificial gain.
- `moving_fields.png`: a first-material-z-slab view of anisotropic expansion, with fixed physical axes and a shared color scale.

The existing dashboard Continuum signaling view still shows the nonlinear fixed-domain experiment. These moving-domain artifacts have a different schema and are not a new dashboard playback mode.

## Interpretation and next test

This verifies a conservative foundation for transport on prescribed changing geometry. It does not show that activator–inhibitor feedback generates the geometry, or that the established nonlinear pattern persists under growth. The smooth cosine family also does not establish convergence for singular corner fields or arbitrary initial data.

The next bridge test is **conservative transfer when the compartment partition changes**: split/refine, coarsen, and transfer between compatible meshes while preserving total amount and a uniform concentration. Check local transfer accuracy, positivity, and round-trip error before introducing cell-to-field deposition and sampling. That will address changes in numerical topology, which this material-deformation test deliberately holds fixed.
