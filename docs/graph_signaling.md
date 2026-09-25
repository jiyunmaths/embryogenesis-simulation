# Discrete graph signaling and apical–basal polarity

This model implements Gierer–Meinhardt activator–inhibitor dynamics directly on the finite normalized contact graph. Spatial instability is decided from its actual eigenvalues, not by assuming that a continuous unstable wavelength is available. The apical–basal polarity module couples signaling and exposed cortex to directional mechanics.

## Normalization and transport assumptions

The existing interface-shell overlap supplies symmetric nonnegative contact weights $W_{ij}$. Self-edges are zero. We remove weights below 2% of the largest current contact (configurable as `graph_contact_cutoff`) and weights below an absolute numerical floor of $10^{-12}$. Threshold sensitivity should be checked; these are diffuse contact proxies, not measured junction areas.

With $d_i=\sum_j W_{ij}$, a connected vertex uses the **random-walk normalized** operator

$$
(\Delta x)_i=\frac{1}{d_i}\sum_j W_{ij}(x_j-x_i).
$$

An isolated vertex has a zero row, so transport is absent and reaction kinetics remain active. Uniform activities are stationary under transport even on an irregular graph. Directly applying the symmetric normalized Laplacian to concentrations would not have this property.

For spectral analysis, use the similar symmetric operator

$$
L_{\mathrm{sym}}=I_{\mathrm{active}}-D^{-1/2}WD^{-1/2},
\qquad
0\le\lambda_k\le2,
\qquad
\operatorname{spec}(\Delta)=\{-\lambda_k\}.
$$

Inverse degrees are zero on isolates. There is one zero eigenvalue per connected component, including isolated vertices. `spectral_gap` is the second-smallest eigenvalue (zero when disconnected, undefined for a single vertex); `smallest_positive_eigenvalue` is reported separately.

Activities are dimensionless signaling states with normalized per-neighbor exchange rates. On a fixed graph this transport preserves degree-weighted activity, **not cell-volume-weighted molecular mass**. Uniform rescaling of all contact weights cancels in the normalization. Weak total contact does not itself slow exchange if its relative graph connections remain. A physically calibrated molecular diffusion model would need contact conductances and cell capacities; it is not implied by this normalization.

## Gierer–Meinhardt equations and equilibrium

Let $a_i$ be activator and $b_i$ inhibitor (`activator` and `inhibitor` in the code):

$$
\begin{aligned}
\dot a_i &= \frac{a_i^2}{b_i}-a_i+D_a(\Delta a)_i,\\
\dot b_i &= \beta(a_i^2-b_i)+D_b(\Delta b)_i.
\end{aligned}
$$

The positive homogeneous equilibrium is $(a_*,b_*)=(1,1)$. Its reaction Jacobian is

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad
\operatorname{tr}J=1-\beta,
\qquad
\det J=\beta.
$$

It is locally asymptotically stable for $\beta>1$. The defaults are $\beta=2$, $D_a=1$, and $D_b=20$. The reaction growth rates at the uniform mode have real part $-0.5$.

For each **discrete** graph eigenvalue, linear perturbations evolve through

$$
M_k=J-\lambda_k\operatorname{diag}(D_a,D_b),
\qquad
r_k=\max\operatorname{Re}\operatorname{eig}(M_k).
$$

A diffusion-driven instability is reported only if the local kinetics are stable and at least one nonzero graph eigenvalue has $r_k>0$. Locally unstable kinetics are not labeled a Turing instability.

The determinant is

$$
\det M_k=\beta+(\beta D_a-D_b)\lambda_k+D_aD_b\lambda_k^2.
$$

For the defaults, its roots delimit the candidate band $(0.129844,0.770156)$. This interval is an envelope in normalized eigenvalue space, not a continuum wavelength prediction. If the graph has no eigenvalues in the interval, this mechanism produces no linear spatial instability.

The numerical signaling update uses positive SSP-RK2 stages with automatic substeps based on the reaction-loss and transport rates. Neither activity is silently clipped to manufacture positivity. A nonfinite or nonpositive inhibitor raises an error. Mechanics and signaling have separate numerical time-step restrictions.

## Spectral gaps and cleavage-induced changes of scale

Before any coupled 3D stepping, the CLI analyzes 4-, 8-, and 16-node cycle and complete reference graphs. The standalone preflight can additionally analyze an actual checkpoint contact graph or a supplied symmetric matrix:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.graph_analysis \
  --checkpoint outputs/progressive/final_state.npz \
  --output outputs/my-graph-preflight

OPENBLAS_NUM_THREADS=1 python -m embryo.graph_analysis \
  --contacts my_contacts.npy --output outputs/my-contact-analysis
```

The checkpoint option reads geometry only, so historical checkpoint geometries can be analyzed without resuming old dynamics. Reference graphs illustrate spectral limitations; they are not substituted for the actual cell contacts in the simulation.

For a cycle of $N$ vertices, the exact spectrum is

$$
\lambda_m=1-\cos(2\pi m/N).
$$

| Graph | Lowest nonzero eigenvalue | Number of unstable modes |
|---|---:|---:|
| Cycle 4 | 1.000000 | 0 |
| Cycle 8 | 0.292893 | 2 |
| Cycle 16 | 0.076120 | 4 |
| Complete 4, 8, or 16 | $N/(N-1)$ | 0 |

The cycle fundamental becomes unstable at eight vertices, then falls below the instability band at sixteen; higher harmonics remain unstable. Modes with multiplicity two count twice but do not imply two independently fixed spatial orientations.

At a refinement $N\to2N$, an inherited cycle harmonic has twice as many vertices per wavelength and a smaller normalized eigenvalue. This is **not automatically a doubling of physical wavelength**. If the same physical domain is discretized more finely, a physical diffusivity would require rates scaled with inverse squared neighbor spacing. The current graph model instead keeps normalized per-neighbor rates fixed and explicitly exposes the resulting scale dependence. Irregular cell graphs have no unique physical wavelength attached to an eigenvalue.

At every actual abscission, the solver:

1. Saves old/new contact matrices, cell IDs, discrete spectra, and growth rates.
2. Constructs a prolongation matrix $P$ that copies the mother's activities to both daughters while retaining all other cells.
3. Projects inherited modes into the new graph eigenbasis.
4. Records the realized partition-noise jump separately from inherited signals.

If $U$ and $U'$ are symmetric Laplacian eigenbases, the mode-transfer matrix is

$$
T=U'^{\mathsf T}D'^{1/2}P D^{-1/2}U.
$$

For this coordinate transformation only, isolated vertices use unit scale. Degenerate eigenvectors can rotate or change sign between computations. We therefore sum $|T|^2$ over each old/new degenerate eigenspace and report normalized energy fractions, rather than arbitrarily matching individual eigenvectors. These fractions describe spectral redistribution, not conserved physical mass or absolute amplitude gain.

`graph_history.json` stores actual graph snapshots. `cleavage_spectra.json` stores all division transitions. The full signaling dynamics automatically evolves on the changing graph; frozen-graph growth rates alone do not solve the time-dependent, mechanically coupled stability problem. Eigenvector mixing, finite time within an unstable regime, and feedback can prevent a mature pattern even when some instantaneous rates are positive.

## Coupling to identity and inheritance

The existing fate variable $f_i$ receives an added bias $g_a(a_i-1)$ in its bistable drift. Activator/inhibitor variables are distinct from the two possible cell identities. Default independent fate noise, neighbor inhibition, exposure bias, and fate partition noise are now zero, so they do not independently manufacture a fate pattern.

The zygote begins at $a=b=1$. Small, bounded, target-volume-balanced activity perturbations are introduced at cleavage (`signal_partition_noise=0.001`), using a separate RNG stream. For lobe fraction $\theta$ and parent activity $x$, a perturbation $\delta$ gives daughters $x+2\delta(1-\theta)$ and $x-2\delta\theta$. Without partition noise, the inherited field is exactly $Px$. This preserves target-volume-weighted activity at cleavage, but does not make subsequent normalized transport molecular-mass-conserving.

## Apical–basal polarity and mechanics

Each cell carries a vector $\mathbf p_i$ whose direction points toward the apical side and whose magnitude is bounded by one. It begins at zero. Using outward diffuse-interface normal $\mathbf n_i$, shell weight $s_i=\phi_i(1-\phi_i)$, and the existing unoccupied-cortex weight $e_i(\mathbf x)$, the geometry cue is

$$
\mathbf q_i=\frac{\int s_i e_i\mathbf n_i\,\mathrm d\mathbf x}{\int s_i\,\mathrm d\mathbf x}.
$$

An isolated spherical cell has zero cue. Contacts can create a net cue toward free cortex without an imposed global axis. Polarity evolves as

$$
\dot{\mathbf p}_i=
\alpha\frac{2a_i}{1+a_i}\mathbf q_i
+\eta(\Delta\mathbf p)_i
-(\mu+\lVert\mathbf p_i\rVert^2)\mathbf p_i.
$$

The discrete update projects magnitudes above one back to one. Both daughters inherit the mother's vector and then adapt to their new contacts. Thus apical polarity can arise from geometry even before an activator pattern develops; it is not by itself evidence of a Turing-selected axis.

Polarity creates a spatial cortical tension field. For regularized radial direction $\widehat{\mathbf r}_i$ from the cell centroid,

$$
\begin{aligned}
\gamma_i(\mathbf x)&=\gamma_i^0[1-\chi\mathbf p_i\cdot\widehat{\mathbf r}_i(\mathbf x)],\\
F_i^{\mathrm{surface}}&=
\epsilon^2\nabla\cdot(\gamma_i\nabla\phi_i)-\gamma_i q'(\phi_i).
\end{aligned}
$$

The apical side has lower effective cortical tension for positive $\chi$; the basal side has higher tension. With $\chi<1$ and $\lVert\mathbf p_i\rVert\le1$, tension stays positive. The isotropic baseline $\gamma_i^0$ retains the existing fate-dependent coefficient.

The implementation computes conservative face fluxes with reflecting boundaries. It includes the spatial gradient of tension; multiplying a Laplacian by $\gamma_i$ alone would omit part of the force. Cell centroids and polarity are frozen during each mechanical substep, then updated. This is a phenomenological active-cortex rule, not a calibrated hydrodynamic force balance, resolved apical protein network, or tensile-stress-based spindle rule. It introduces no preselected embryo outline.

`--no-polarity` removes polarity dynamics and its tension contrast. `--no-feedback` leaves signaling/polarity dynamics active but removes their fate/polarity-dependent mechanical effects. `--no-signaling` holds chemical activities at their inherited values and removes signal partition noise in a new run; geometry can still polarize cells. `--mechanics-only` disables signaling, polarity, and differentiation.

## Verification and observed limits

The graph-only mode tests use small perturbations aligned to a discrete eigenmode and compare fitted exponential rates with the eigenvalues of $M_k$. Default cases match within $4\times10^{-7}$. Tests also cover uniform-state preservation on irregular/disconnected graphs, equal-diffusion stability, positive nonlinear signaling, conservation at inheritance, and basis-independent transfer within degenerate eigenspaces.

Polarity tests verify no grid-selected cue in an isolated sphere, apical cues away from a neighbor, directional mechanical response, and the variable-tension flux term. Checkpoint schema 3 preserves signaling, polarity, three RNG streams, and all graph-event state for exact continuation. Previous model checkpoints cannot resume with identical dynamics; geometry-only spectral analysis remains available.

The initial coupled run uses seed 7, a $40^3$ grid, and time 15. Actual first recorded 4-, 8-, and 16-cell graphs have zero, two, and four unstable spatial modes, respectively. Both the polarity-enabled run and its matched no-polarity control reach 16 cells. At time 15, the polarity-enabled activator standard deviation is only about $0.00137$ and all cells remain below the fate-commitment thresholds. The mean polarity magnitude is about $0.0835$; final axis ratios are 1.2689 with polarity and 1.2703 without. This establishes functioning coupling and finite-mode growth analysis, **not robust fate differentiation or a polarity-driven developmental axis**. Signaling growth times relative to cleavage and mechanical relaxation remain a research variable.

A subsequent [timescale experiment](timescales.md) extends the coupled run to time 60, freezes actual division-event graphs, and compares slower cleavage with controlled graph replay. It observes large later signal differences and quantifies the time needed for amplification. Its exact division-event spectra can differ from the regularly sampled snapshots above.

These are single-seed, coarse-grid demonstrations. Contact threshold sensitivity, evolving-graph transient amplification, longer signaling times, and mechanical/grid refinement remain necessary before biological interpretation.

## Sources

- [Nakao and Mikhailov: Turing patterns in network-organized activator–inhibitor systems](https://www.nature.com/articles/nphys1651) motivates graph-based reaction–diffusion stability.
- [How local spectral gaps regulate the multistability of Turing patterns on graphs](https://journals.plos.org/complexsystems/article?id=10.1371/journal.pcsy.0000044) examines finite spectra and accessible unstable modes.
- [Nissen et al.: Theoretical tool bridging cell polarities with development of robust morphologies](https://elifesciences.org/articles/38407) motivates polarity-dependent mechanical interactions.

The implemented normalized transport and cortical tension rules are explicit modeling choices; this is not a reproduction of these papers' complete models.
