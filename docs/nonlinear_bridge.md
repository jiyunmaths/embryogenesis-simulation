# Nonlinear signaling persistence and numerical refinement

## Question and scope

Does the early Gierer–Meinhardt instability develop into a persistent spatial pattern, and does that pattern approach the same result when space and time are refined independently?

This experiment follows the [weighted-spectrum study](irregular_signaling.md). It uses a **fixed 3D L-shaped prism**, one prescribed continuous initial perturbation, and conservative concentration transport. The compartments are numerical control volumes, not biological cells. The experiment tests nonlinear signaling; it contains no cleavage, fate switch, mechanical feedback, or moving boundary. It does not establish a unique attractor or robustness across initial conditions.

Implementation: [embryo/nonlinear.py](../embryo/nonlinear.py). Independent verification: [tests/test_nonlinear.py](../tests/test_nonlinear.py). The embryo simulation and dashboard continue to use their existing solvers.

## Physical equations and spatial discretization

For activator $a$ and inhibitor $h$:

$$
\partial_t a = \frac{a^2}{h}-a+D_a\nabla^2a,
\qquad
\partial_t h = \beta(a^2-h)+D_h\nabla^2h.
$$

All external walls, including the reentrant walls, have zero normal flux. The domain is $[0,L]^3$ with the region $x>L/2,\ y>L/2$ removed for every $z$. Defaults are $L=1$, $\beta=2$, $D_a=0.02$, and $D_h=0.4$. The homogeneous equilibrium is $(a,h)=(1,1)$. The activator–inhibitor mechanism follows [Gierer and Meinhardt (1972)](https://www.bio.mpg.de/255219/gierer-and-meinhardt-1972); these dimensionless parameters and this benchmark are our reduced model, without a biological time or length calibration.

The [graded finite-volume mesh](irregular_bridge.md) uses grading $g=0.35$. With compartment capacities $M=\operatorname{diag}(V_i)$ and symmetric face conductances $G_{ij}=A_{ij}/d_{ij}$, define

$$
K=\operatorname{diag}(G\mathbf{1})-G,
\qquad \Delta=-M^{-1}K.
$$

Transport conserves the volume-weighted amount; reactions can produce or remove it. Before the long runs, the code computes the actual complete discrete spectra using the [product-domain method](irregular_signaling.md). The unstable spatial mode counts are 7, 4, and 4 at $n=8,16,32$. A positive linear growth rate motivates the nonlinear test but does not predict its saturated pattern.

## One physical initial perturbation on all meshes

A mesh-specific random draw would change both the physical initial condition and its shortest wavelengths during refinement. Instead, prescribe

$$
q(x,y,z)=\cos(\pi x/L)\cos(\pi z/L)
+0.6\cos(\pi y/L)
+0.3\cos(2\pi x/L)\cos(\pi y/L)-\frac{0.4}{\pi},
$$

$$
a(x,y,z,0)=1+\varepsilon q(x,y,z),
\qquad h(x,y,z,0)=1-0.7\varepsilon q(x,y,z),
\qquad \varepsilon=0.01.
$$

The constant makes the volume mean of $q$ zero on this domain. Each mesh receives the exact rectangular **cell averages** of these functions, with no per-mesh amplitude normalization. For a coordinate interval of width $w$ and center $x_c$,

$$
\overline{\cos(m\pi x/L)}
=\cos(m\pi x_c/L)\operatorname{sinc}\!\left(\frac{mw}{2L}\right),
\qquad \operatorname{sinc}(s)=\frac{\sin(\pi s)}{\pi s}.
$$

The initial function need not satisfy Neumann compatibility on the reentrant walls. The evolution imposes no flux for positive time, so an initial boundary adjustment is part of this experiment. This is controlled amplification of a specified perturbation, not a test of spontaneous selection from an ensemble of disordered states.

## Time integration and independent verification

Long fine-mesh explicit diffusion runs require many small steps. This benchmark uses second-order semi-implicit backward differentiation (SBDF2), with implicit diffusion and explicit reaction extrapolation. For either species $c$, diffusion coefficient $D$, and its reaction component $R$:

$$
\left(\frac32 I-\delta t D\Delta\right)c^{k+1}
=2c^k-\frac12c^{k-1}
+\delta t\left(2R^k-R^{k-1}\right).
$$

The first step uses $(I-\delta t D\Delta)c^1=c^0+\delta t R^0$. Its single-step error is of order $\delta t^2$, consistent with global second-order accuracy. See [Ascher, Ruuth, and Wetton (1995), *Implicit-Explicit Methods for Time-Dependent Partial Differential Equations*](https://epubs.siam.org/doi/10.1137/0732037) for the IMEX framework.

The weighted symmetric operator separates into the L-shaped cross section and the $z$ direction. Transforming into the weighted $z$ eigenbasis reduces each implicit solve to sparse shifted two-dimensional systems. This is an algebraically exact factorization of the discrete three-dimensional solve; it introduces no directional time-splitting approximation.

SBDF2 is **not unconditionally positivity-preserving**. Every computed state is checked for finite, strictly positive concentrations. A failure stops the run and retains its report and previously completed trajectories; no concentration clipping or automatic step-size adjustment is used. Time-step refinement is therefore essential.

The software tests compare the factorized solve against a directly assembled full 3D sparse solve, and the time integrator against independent DOP853 integration on an unequal-capacity graph. Halving the step gives approximately fourfold error reduction in that smooth ODE test. Additional checks cover agreement with the existing SSP-RK2 solver, pure-diffusion amount conservation, the homogeneous equilibrium, independent Gaussian quadrature of the initial fields, conservative restriction, and failure handling.

The discrete amount-balance residual checks the volume sum of the SBDF2 equation after removing its reaction contribution. It is a linear-solve/conservation check, **not** a continuum time-accuracy estimate.

## Comparison protocol and predefined criteria

Spatial refinement uses $n=8,16,32$, corresponding to 384, 3,072, and 24,576 active compartments, all with $\delta t=0.05$ through $t=100$. Independent temporal refinement uses $n=32$ and $\delta t=0.05,0.025,0.0125$. An equal-diffusivity control uses $n=16$, $D_a=D_h=0.02$, and the same initial perturbation; its spectral preflight has no unstable modes.

For spatial comparison, fine concentrations are restricted to coarse compartments by summing their amounts and dividing by the total volume. The nested graded meshes have identical physical boundaries. No field rotation, sign change, translation, or fitted registration is performed.

The relative joint field error is

$$
E(c,\widehat c)=
\frac{\left[\sum_i w_i\{(a_i-\widehat a_i)^2+(h_i-\widehat h_i)^2\}\right]^{1/2}}
{\left[\sum_i w_i\{(\widehat a_i-\langle\widehat a\rangle)^2+
(\widehat h_i-\langle\widehat h\rangle)^2\}\right]^{1/2}},
\qquad w_i=\frac{V_i}{\sum_j V_j}.
$$

The reference is the finer result, restricted where necessary. Errors include differences in means as well as spatial arrangement, and are normalized by reference spatial contrast. Uniform references have undefined relative error and are reported as `null`, not zero. Activator correlations are also volume-weighted.

Full fields are compared at $t=0,20,40,60,75,80,90,100$. Scalar histories are sampled every unit of time. Persistence means activator standard deviation remains at least 0.1 over the sampled late window $75\leq t\leq100$. Approximate stationarity means the joint field changes by at most 1% relative to its $t=75$ pattern over this sampled window. Neither criterion is an asymptotic stability theorem or a recovery-after-perturbation test.

Eleven criteria are recorded before the long runs:

1. Initial joint RMS discrepancy after restriction below $10^{-12}$.
2. Positive concentrations at every integration step in every run.
3. Discrete reaction-adjusted amount-balance residual below $10^{-9}$.
4. Persistent contrast on the finest mesh at all three time steps.
5. Late field change at most 1% for those same three runs.
6. Final spatial error decreases on refinement.
7. Final finest-pair spatial error at most 5%.
8. Final finest-pair activator correlation at least 0.99.
9. Maximum sampled temporal error decreases on time refinement, or both differences are below $10^{-10}$.
10. Finer temporal-pair maximum sampled error at most 0.5%.
11. Equal-diffusivity control has activator standard deviation below $10^{-6}$ throughout the sampled late window.

Failed criteria remain in the output and cause a nonzero CLI exit. These thresholds are practical acceptance choices for this benchmark, not universal biological criteria.

## Measured results

The default verified run is recorded in `outputs/nonlinear-bridge-verified`. **All eleven criteria pass.**

| Resolution | Active compartments | Activator SD at $t=100$ | Maximum late relative field change |
|---|---:|---:|---:|
| $n=8$ | 384 | 0.691659 | 0.329% |
| $n=16$ | 3,072 | 0.684575 | 0.394% |
| $n=32$ | 24,576 | 0.682857 | 0.405% |

| Comparison | Relative field difference | Final activator correlation |
|---|---:|---:|
| $n=8$ versus $16$, final | 5.433% | 0.999131 |
| $n=16$ versus $32$, final | 1.263% | 0.999954 |
| $\delta t=0.05$ versus $0.025$, maximum sampled | 0.3180% | 0.999999999993 |
| $\delta t=0.025$ versus $0.0125$, maximum sampled | 0.09536% | 0.999999999999 |

The ratio of the two maximum sampled temporal differences gives an empirical order of **1.738** for this experiment. It supports decreasing temporal error, but does not demonstrate an asymptotic order of exactly two at these steps. The finest run crosses activator SD 0.1 at approximately $t=13.670$ with $\delta t=0.0125$, compared with $t=13.642$ at $\delta t=0.05$; onset estimates use interpolation between scalar samples.

The equal-diffusivity control ends with activator SD $4.14\times10^{-15}$. Across the long runs, the smallest activator and inhibitor values are approximately 0.03272 and 0.50961, respectively. The largest reaction-adjusted discrete amount-balance residual is $7.90\times10^{-16}$.

This is evidence that the nonlinear spatial pattern in this controlled case persists and is increasingly consistent under numerical refinement. It does not establish that every initial perturbation selects the same pattern, that the pattern recovers after a disturbance, or that the chosen equations explain a particular embryo.

## Reproduce and inspect

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.nonlinear \
  --output outputs/nonlinear-bridge-repeat

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q
```

Use a fresh output path. Generated outputs are ignored by Git; source, tests, and this protocol retain the reproducible method and measured results.

- `analysis.json`: parameters, preflight spectra summaries, scalar histories, errors, criteria, and pass/fail results.
- `fields.npz`: mesh edges, masks, capacities, and saved two-species concentration fields for every run.
- `convergence.png`: contrast history, conservative spatial comparisons, and temporal comparisons.
- `patterns.png`: initial, intermediate, and final activator fields in one fixed near-boundary $z$ slab, with a shared concentration scale. This slice is illustrative; all quantitative comparisons use the complete 3D fields.

## Next scientific step

The [moving-domain amount and dilution benchmark](moving_domain.md) now passes eleven checks for isotropic expansion, anisotropic expansion, and volume-preserving deformation. It verifies inverse-volume dilution, conservation, and approximately second-order transport accuracy for prescribed affine motion. The next bridge test is conservative transfer during changes in the compartment partition, before cell/field coupling or mechanically generated deformation.

Separately, nonlinear robustness still needs multiple physically matched perturbations, parameter variation, and recovery after disturbance. Eventual fate and mechanics experiments must test their own persistence and causal controls. The current fixed-domain result is a foundation for those tests, not their replacement.
