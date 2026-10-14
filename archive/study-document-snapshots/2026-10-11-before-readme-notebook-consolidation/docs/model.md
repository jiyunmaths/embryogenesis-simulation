# Model definition

This document describes the **current attribute-based research model**, its shared phase-field mechanics, and the historical core/fate alternative. The current model uses continuous chemistry, polarity, and cell shape without a supplied identity switch. Historical sections are retained so earlier experiments remain reproducible; their A/B labels are not observations from the new model.

The conservative chemical equations, dilution, and spectral analysis are specified in [live_transport.md](live_transport.md). Historical normalized chemical transport and the current polarity law are distinguished in [graph_signaling.md](graph_signaling.md). Here $a_i,b_i$ are chemical activities; $r_i$ is an instantaneous material response; the downstream $f_i$ appears only in the historical core alternative.

All coordinates, times, activities, energies, and coefficients are dimensionless. The model is generic. Parameters are illustrative rather than measured.

## Current model and experiment map

`AttributeSimulation` in `embryo/attribute_development.py` retains deformable cell fields, Gierer–Meinhardt reactions, conservative volume-weighted exchange, mechanical dilution, apical–basal polarity, and progressive division. It disables downstream fate integration/noise and does not classify A/B. A zero-valued compatibility array remains only for the existing checkpoint/bookkeeping format. Two chemical species do not prescribe two cell identities.

The direct material response is

$$
\begin{aligned}
r_i&=\tanh(a_i-1),\\
\gamma_i^0&=\gamma_0(1+c_\gamma r_i),\\
A_{ij}&=A_0(1+c_A r_i r_j),\qquad A_{ii}=0.
\end{aligned}
$$

The unit reference makes uniform chemistry mechanically unbiased. The tanh bounds material modulation and has no independent evolution or stored memory. Higher activity raises baseline tension; the product gives stronger attraction for same-sign responses. Current contrasts are $c_\gamma=0.25$ and $c_A=0.35$, stored under the legacy names `fate_tension` and `fate_adhesion`. These constitutive assumptions are not inferred chemical identities or calibrated molecular effects. Polarity then makes $\gamma_i(\mathbf x)$ directional, using the existing [polarity equation and tension law](graph_signaling.md#apicalbasal-polarity-and-mechanics).

The no-feedback branch keeps baseline tension/attraction and removes polarity's action on mechanics, while chemistry and polarity dynamics still evolve. Component interventions separately enable tension, adhesion, and polar mechanics. Polarity's chemical activation factor remains active in those controls; removing polar mechanics is not the same as removing chemical modulation of polarity.

| Experiment | Question answered | Starting-state distinction |
|---|---|---|
| Fresh development | Do differences form as a zygote divides? | Analytic zygote, initially unit chemical activities and zero polarity |
| Developed-pattern survival | Does a pattern survive switching mechanics? | Actual no-feedback developmental chemistry at t=90, not a transplanted equilibrium |
| Frozen bistability/recovery | What stable chemical states and basins does one graph support? | Fixed endpoint geometry; independent chemical integration |
| Chemical exchange | Do states return, transfer, or reorganize? | Relocated activator/inhibitor; conservative exchange preserves both species' amounts at the intervention |
| Moving response transfer | Does subsequent behavior follow donor chemistry or destination context? | Mature equilibrated/transplanted chemistry at t=150, moving formation to t=210, then matched pulse/control continuations |

Current research uses a 72³ voxel grid, half-width 2.24, width parameter 0.085, and a sixteen-cell cap. The developmental/component pilots use dt=0.0075; mature response and accepted GPU continuations use dt=0.00375. The targeted seed-7 response refinement uses dt=0.001875 on CPU. No claim of full attribute-development or spatial convergence follows from these selected checks. Common/finite-reservoir chemical assays are separate experiments; a reservoir is not part of the live 3D/GPU state.

Phenotypes are described by log activator, log inhibitor, polarity magnitude, cell axis ratio, and asphericity. Location, exposure, contacts, and lineage are context. Activity-derived material coefficients do not provide independent phenotype evidence. Persistence, response, relocation, context dependence, and eventual inheritance tests are needed before assigning cell identity. See [attribute development](attribute_development.md), [cell response](cell_response.md), and [moving-history replication](cell_response_moving.md#replication-across-developmental-histories).

## Coupled integration and computing paths

For each mature attribute-model step:

1. Measure occupancy volumes, centers, contact shells, and exposed-cortex cues on the current geometry.
2. Construct the conservative chemical operator and separate normalized orientation-averaging operator.
3. Advance reaction/exchange with positivity-restricted SSP-RK2 substeps, keeping geometry fixed.
4. Advance polarity explicitly using updated activator; evaluate activity-dependent material coefficients.
5. Advance phase fields explicitly, including the spatial gradient of directional tension.
6. Apply $a_i\leftarrow a_i V_i^{\mathrm{old}}/V_i^{\mathrm{new}}$ and the same update to $b_i$, conserving chemical amount during mechanical volume change; advance the physical clock.

CPU development additionally completes or initiates progressive cytokinesis and amount-balanced daughter inheritance. **The complete coupling is first-order splitting.** Second-order frozen reaction/exchange integration, exact amount bookkeeping, and positive concentrations do not establish second-order accuracy of the evolving embryo.

| Execution path | Current use |
|---|---|
| Python/NumPy/SciPy `AttributeSimulation` | CPU development, original equations, and reference checkpoints |
| Native C++/OpenMP `NativeSimulation` | Accelerated mature CPU references and refinements; scientific use is scoped to its validation |
| PyTorch/custom-CUDA `GpuSimulation` through the gated adapter | Mature polar direct-feedback conservative continuations at accepted parameters |
| Core `Simulation` / dashboard | Historical supplied-fate model; not the current attribute/GPU research runner |

Every moving production step screens finite positive chemistry, per-cell volume error <5%, equivalent radius ≥4 grid spacings, and zero clipping. Recorded observations screen boundary occupancy <0.01. Mature response observations are every 0.15 units and standard restart checkpoints every three units. Preparation freezes source/input hashes; restart preserves IDs, physical state, clocks, and random streams and does not repeat an intervention.

The [resident GPU validation](#resident-gpu-backend-validation) measures whole trajectories, pulse-minus-control responses, recovery, and endpoint fields against CPU evidence. Its successful backend agreement is separate from timestep/spatial convergence and biological validation. Unsupported regimes or changed physical parameters require additional validation before scientific GPU use.

## Reading the equations

The mechanical rules in `embryo/model.py` are shared with the attribute model; its fate rule belongs only to the core alternative. [Conservative signaling](live_transport.md) and [polarity mechanics](graph_signaling.md#apicalbasal-polarity-and-mechanics) give the complementary rules. Read an equation first as a balance of competing effects, then use its coefficients to ask which effect acts most strongly or most quickly.

| Notation | How to read it |
|---|---|
| $i,j$; $\mathbf{x}$; $t$ | Cell labels; spatial position; time. A field depends on position within the box, whereas a cell-level variable has one value per cell. |
| $\int_\Omega u\,\mathrm d\mathbf{x}$ | Add the quantity $u$ over the 3D computational box; the code approximates this with voxel values times $\Delta x^3$. |
| $\nabla\phi$; $\lVert\nabla\phi\rVert^2$ | Spatial rate of change; its squared magnitude. Large values identify a sharp boundary. |
| $\nabla\cdot(\gamma\nabla\phi)$ | Divergence of a spatial flux; it includes variation of both $\phi$ and $\gamma$. |
| $\dot f$; $\partial_t\phi$ | Time rate of change of a cell variable; time rate of change of a spatial field. |
| $h'(\phi)$; $q'(\phi)$ | Ordinary derivatives of scalar functions. |
| $\delta E/\delta\phi$ | Functional derivative: the first-order change in total energy caused by changing the field locally. |
| $V_i^{\mathrm{target}}$ | Preferred volume, also written $V_i^\star$ in the README; it is distinct from measured $V_i$. |
| $A_{ij}$; $\mathcal A_{ij}$ | Adhesion energy coefficient; geometric contact area. These are different quantities. |

“Dimensionless” means the variables have been scaled by unspecified reference units. Model length, time, and energy still have distinct roles: $\epsilon$ and $\Delta x$ are lengths, diffusion coefficients scale as length squared per time, and a rate scales as inverse time. No conversion to micrometers, seconds, or measured cellular forces has been calibrated.

There are three different kinds of input. **Constitutive choices** specify how cells interact (for example, the adhesion law). **Rate and strength parameters** set their relative importance. **Numerical choices**, such as $\Delta x$ and $\Delta t$, determine how accurately those rules are solved. Changing a constitutive parameter is a different experiment from refining the grid or time step.

## Geometry

Cell $i$ has a diffuse indicator $\phi_i(\mathbf{x}) \in [0,1]$. Its occupancy and volume are

$$
h(\phi) = \phi^2(3-2\phi),
\qquad
V_i = \int_{\Omega} h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

Here $\Omega$ is the computational domain. The initial zygote is a radius-0.8 sphere with a smooth interface. Cell shapes are free to change. The finite computational box uses zero-normal-gradient boundary conditions and must be large enough to avoid affecting the embryo.

### Why occupancy is a cubic

The phase field is an indicator with a diffuse boundary, not a probability distribution or a physical material density. We choose a volume-counting function that satisfies

$$
h(0)=0,\qquad h(1)=1,\qquad h'(0)=h'(1)=0.
$$

Requiring these four conditions of a cubic determines $h(\phi)=3\phi^2-2\phi^3$. It increases monotonically on $[0,1]$, gives half occupancy at $\phi=1/2$, and makes volume insensitive to first-order field changes in uniform interior/exterior regions. Its derivative is

$$
h'(\phi)=6\phi(1-\phi).
$$

Thus volume-restoring forces mainly move the diffuse boundary. The volume integral counts fractional boundary voxels instead of abruptly counting only voxels above a threshold. Other smooth interpolations are possible; this one is a convenient assumption whose consequences still need resolution checks.

### Why the potential is a quartic

Define

$$
q(\phi)=\phi^2(1-\phi)^2,
\qquad q'(\phi)=2\phi(1-\phi)(1-2\phi).
$$

The potential is nonnegative and symmetric under $\phi\mapsto1-\phi$. It has minima at 0 and 1, and a maximum $1/16$ at $1/2$. In isolation, the restoring term $-q'(\phi)$ pushes values below $1/2$ toward zero and values above $1/2$ toward one. It therefore discourages a broad intermediate phase. The gradient cost opposes infinitely abrupt transitions, so the two terms together produce a smooth interface rather than either a step function or a spatially uniform intermediate value.

The squared factors are a simple way to give both pure phases zero energy and zero slope. Their use is a modeling convention, not evidence that cellular cortex energy literally follows a quartic polynomial. Notice also that $q=s^2$ for the interface-shell weight $s=\phi(1-\phi)$: this is why $q_iq_j$ can select contact between two diffuse boundaries.

![Occupancy rises smoothly from zero to one; the potential has minima at zero and one and a maximum at one half.](images/occupancy-potential.svg)

The panels have different vertical scales: occupancy counts volume, whereas the potential contributes to interface energy. They are different functions with different jobs.

### Baseline energy: four competing costs

With material response, polarity, target volumes, cell centroids, and the spatial tension coefficient $\gamma_i(\mathbf{x})$ frozen during the mechanical update, the energy is:

$$
\begin{aligned}
E ={}&
\sum_i \int_{\Omega} \gamma_i(\mathbf{x})
\left[
\frac{\epsilon^2}{2}\lVert \nabla\phi_i \rVert^2
+ q(\phi_i)
\right]\,\mathrm{d}\mathbf{x}
\\
&+ \sum_i \frac{K_V}{2V_i^{\mathrm{target}}}
\left(V_i-V_i^{\mathrm{target}}\right)^2
\\
&+ \frac{R}{2}\sum_{i<j}
\int_{\Omega}\phi_i^2\phi_j^2\,\mathrm{d}\mathbf{x}
\\
&- \sum_{i<j} A_{ij}
\int_{\Omega}q(\phi_i)q(\phi_j)\,\mathrm{d}\mathbf{x}.
\end{aligned}
$$

The equation can be read as $E=E_{\mathrm{interface}}+E_{\mathrm{volume}}+E_{\mathrm{overlap}}+E_{\mathrm{attraction}}$. The first three contributions are nonnegative; attraction is negative. There is no global embryo template in these terms.

| Contribution | Interpretation of the mathematical form |
|---|---|
| Interface | A squared gradient penalizes steep spatial changes in every direction equally. The potential prefers pure phases. Their competition assigns a cost to cell boundary area. An isolated cell with constrained volume consequently tends toward a compact shape. |
| Volume | A quadratic is the simplest smooth penalty with zero value and zero derivative at the preferred volume. It penalizes equal positive and negative volume deviations equally. The normalization by target volume makes the restoring coefficient depend on fractional volume error. |
| Overlap | The product is large when both cells have interior values near one at the same location, and vanishes when either field is zero. It discourages interpenetration without imposing a hard geometric exclusion constraint. |
| Attraction | The two $q$ factors select overlapping interfaces; the minus sign rewards this contact. The term is small deep inside either cell, unlike a bulk-overlap attraction. |

The pair sum $i<j$ avoids counting interactions twice. The factor $1/2$ in the overlap term additionally defines the convention for $R$: differentiating with respect to either cell field gives $R\phi_i\phi_j^2$, not $2R\phi_i\phi_j^2$. Similarly, the volume factor $1/2$ cancels the derivative of the square. These factors set coefficient definitions rather than adding a separate biological mechanism.

For relative volume error $e_i^V=(V_i-V_i^{\mathrm{target}})/V_i^{\mathrm{target}}$, the volume cost becomes

$$
E_{\mathrm{volume},i}=\frac{K_V V_i^{\mathrm{target}}}{2}(e_i^V)^2,
\qquad
\frac{\partial E_{\mathrm{volume},i}}{\partial V_i}=K_V e_i^V.
$$

A larger cell pays proportionally more total energy for the same fractional error, but the derivative with respect to volume has the same magnitude. This is a soft penalty: it permits force balance at a small nonzero error. It is not a solved hydrostatic pressure or a fluid incompressibility equation.

| Symbol | Configuration field / state | Default | Role |
|---|---|---:|---|
| $\epsilon$ | `interface_width` | 0.085 | Diffuse-boundary length scale; also affects interface energy per area. |
| $\gamma_0$ | `surface_tension` | 1 | Baseline coefficient multiplying the interface energy. |
| $K_V$ | `volume_stiffness` | 12 | Restoring strength per relative volume error. |
| $R$ | `repulsion` | 3 | Strength of the bulk-overlap penalty. |
| $A_0$ | `adhesion` | 4 | Baseline strength of diffuse-interface attraction. |
| $V_i^{\mathrm{target}}$ | `sim.target[i]` | Initial measured zygote volume, then partitioned | Preferred cell volume, not a growth rate. |
| $\gamma_i(\mathbf{x}), A_{ij}$ | Computed from chemical response/polarity; legacy core uses fate instead | State-dependent | Local interface cost and pair attraction coefficient. |

Increasing stiffness or repulsion strengthens the corresponding restoring force, but may require a smaller explicit time step. Increasing attraction favors contact, but cannot be interpreted independently of repulsion, tension, and diffuse width. Coefficient magnitudes are not directly comparable: $A_0=4$ versus $R=3$ does not mean adhesion wins, because they multiply different functions and spatial regions.

### How the energy generates shape change

For fixed coefficients and targets, the passive variational force is

$$
\begin{aligned}
-\frac{\delta E}{\delta\phi_i}
={}&\epsilon^2\nabla\cdot(\gamma_i\nabla\phi_i)-\gamma_i q'(\phi_i)\\
&+K_V\frac{V_i^{\mathrm{target}}-V_i}{V_i^{\mathrm{target}}}h'(\phi_i)\\
&-R\phi_i\sum_{j\ne i}\phi_j^2\\
&+q'(\phi_i)\sum_{j\ne i}A_{ij}q(\phi_j).
\end{aligned}
$$

This is the force assembled in the core and attribute mechanical updates, before division forcing, clipping, and volume projection. Their difference is how material coefficients are derived. If a cell is too small, the volume term is positive at its boundary, increasing occupancy; if it is too large, the sign reverses. The repulsion term reduces field values where another cell overlaps. The attraction term's sign depends on the side of the interface through $q'$: it reshapes interfaces to increase their overlap, rather than adding occupancy everywhere.

For variable tension, $\nabla\cdot(\gamma\nabla\phi)=\gamma\nabla^2\phi+\nabla\gamma\cdot\nabla\phi$. Retaining only $\gamma\nabla^2\phi$ would omit the spatial-tension-gradient contribution. The tension field is held fixed in this derivative; the code does not differentiate through its centroid or regulatory dependence.

Mechanics follows overdamped gradient descent with mobility set to one (outside active furrowing):

$$
\frac{\partial\phi_i}{\partial t}
= -\frac{\delta E}{\delta\phi_i}.
$$

A more general notation is $\partial_t\phi_i=-m_\phi\,\delta E/\delta\phi_i$, with positive mobility $m_\phi$. Setting it to one chooses the mechanical time scale. There is no inertial acceleration, velocity field, or fluid momentum equation in this update.

For the ideal continuous passive problem with all coefficients fixed and compatible reflecting boundaries,

$$
\frac{\mathrm dE}{\mathrm dt}
=-\sum_i\int_\Omega\left(\frac{\delta E}{\delta\phi_i}\right)^2\,\mathrm d\mathbf{x}\le0.
$$

This explains the minus sign. It does not prove monotone energy decrease for a finite explicit step, nor for the actual active system with changing coefficients, furrow forces, clipping, and constraints. Energy here is a mechanical construction, not an organism's fitness or metabolic expenditure.

### Interface width is not grid spacing

For an isolated flat interface, constant $\gamma$, and no other forces, stationarity gives

$$
\epsilon^2\phi''=q'(\phi),
\qquad
\frac{\epsilon^2}{2}(\phi')^2=q(\phi).
$$

The second relation follows by multiplying the first by $\phi'$ and using pure-phase limits. Integrating the energy through the interface gives its cost per area:

$$
\sigma_{\mathrm{flat}}
=\gamma\epsilon\int_0^1\sqrt{2q(\phi)}\,\mathrm d\phi
=\frac{\gamma\epsilon\sqrt{2}}{6}.
$$

Consequently, increasing $\epsilon$ at fixed $\gamma$ changes both interface thickness and interface energy. Some other phase-field conventions normalize the gradient and potential differently to hold this cost fixed; those conventions cannot be substituted silently here. $\Delta x$ is instead the numerical voxel spacing. Refining $\Delta x$ at fixed $\epsilon$ resolves the same diffuse model more accurately. Reducing $\epsilon$ changes the model and its contact calibration as well, and requires a separate study.

The production kernel reuses current-state occupancy/volume arrays and evaluates only the spatial operator selected by the tension model; [exact trajectory comparisons and timing](signal_patch.md#equation-preserving-kernel-optimization) document this optimization. The code uses a finite-difference Laplacian and explicit time integration. Values are clipped to $[0,1]$; `clipped_fraction` exposes this numerical safeguard. Check step-size sensitivity, especially after division, rather than assuming clipping makes the method accurate.

The attraction is phenomenological, not a calibrated cadherin or junction model. The volume penalty is soft outside cytokinesis; a dividing mother additionally uses the volume constraint below. Changing chemical response/polarity, cleavage, and stochastic dynamics mean the full simulation is not passive relaxation of a fixed energy.

## Regulatory activity

**Historical core/dashboard only.** The current attribute model does not integrate this equation, use its reporting thresholds, or feed $f$ into material coefficients. Its continuous chemical response is defined above.

Each cell has a signed scalar $f$, a reduced normal form for competition between two possible identities. It is not a concentration, and need not lie in $[-1,1]$. The isolated deterministic switch

$$
\frac{\mathrm{d}f}{\mathrm{d}t}=r_f(f-f^3)
$$

has stable states at $-1$ and $+1$, with an unstable state at $0$.

The switch is the downhill dynamics of a separate local potential $U(f)=f^4/4-f^2/2$: $\dot f=-r_f U'(f)$. Its two minima express **assumed bistability**. They are not derived from the mechanical double-well potential $q(\phi)$, which separates inside from outside rather than A from B. Linearization gives growth rate $r_f$ at zero and decay rate $-2r_f$ at either unbiased stable state. Even a small transient signal bias can therefore choose a fate that persists after the bias disappears; two labels alone do not demonstrate a sustained Turing pattern.

After the configured competence cell count, the discrete update over a time step $\Delta t$ is:

$$
\begin{aligned}
\Delta f_i ={}&
r_f\left[
f_i-f_i^3
-k_n\,\overline{f}_{i,\mathrm{neighbor}}
+b_e(e_i-0.5)+g_a(a_i-1)
\right]\Delta t
\\
&+\sigma_f\sqrt{\Delta t}\,\xi_i,
\qquad \xi_i\sim\mathcal{N}(0,1).
\end{aligned}
$$

| Symbol | Configuration name (default) | Effect |
|---|---|---|
| $r_f$ | `fate_rate` (0.8) | Multiplies the deterministic response; larger values make commitment faster for the same forcing. |
| $g_a$ | `signal_fate_gain` (1) | Bias strength relative to the reference activator level 1. |
| $k_n$ | `neighbor_inhibition` (0) | Optional legacy bias against the contact-weighted mean neighbor fate. |
| $b_e$ | `exposure_bias` (0) | Optional legacy bias from exposure above or below 0.5; distinct from inhibitor $b_i$. |
| $\sigma_f$ | `fate_noise` (0) | Additive noise amplitude. The $\sqrt{\Delta t}$ scaling makes increment variance proportional to elapsed time. |
| Competence count | `competence_cells` (4) | Population gate for the core fate update, supplied as a developmental assumption. |
| Label threshold | `fate_threshold` (0.55) | Reporting threshold; it does not alter the drift equation. |

The activator activity is $a_i$, with homogeneous equilibrium $a_i=1$; $g_a$ is `signal_fate_gain`. The legacy neighbor inhibition, exposure bias, independent fate noise, and fate partition noise are zero by default, so graph signals provide the default bias.

Here $\overline{f}_{i,\mathrm{neighbor}}$ is the contact-weighted neighbor activity, $e_i$ is the exposure proxy, and $\xi_i$ is an independent standard normal draw for each cell and time step.

Neighbor weights are integrals of overlapping $\phi(1-\phi)$ interface shells. The exposure proxy is the shell-weighted mean of

$$
1-\operatorname{clip}\left(
2\sum_{j\ne i}h(\phi_j),\,0,\,1
\right).
$$

It is not an exact free-surface area fraction. The fixed exposure reference $0.5$ is an explicit model assumption, not calibrated biology.

The neighbor term is lateral inhibition: it can favor different activities among contacting cells. The exposure term biases externally exposed cells toward positive activity. Neither guarantees both identities in every run. The frozen activator–inhibitor subsystem has a discrete linear stability analysis, but it does not establish stability of this complete fate–geometry system on a growing graph.

Positive activity is called A and negative activity B. Counts use $f>0.55$ or $f<-0.55$; remaining cells are uncommitted. Stable commitment requires future persistence and perturbation assays.

The core dashboard solver uses the explicit fate update written above. The later [joint fate integration](joint_fate.md) experiment advances deterministic signal and fate variables at shared RK stages to resolve a documented integration sensitivity; its moving pilot has no additional cleavage. That experiment does not silently replace the core solver or validate the optional stochastic terms.

The completed [withdrawal assay](fate_memory.md) tests this distinction directly: short transient activator input seeds later bistable commitment, while a matched relaxing switch loses labels after withdrawal. The subsequent [isolated-switch reversal/noise assay](fate_robustness.md) characterizes reversibility and stochastic sensitivity; robustness within the moving embryo remains a separate test.

## Feedback

**Historical core alternative.** Here material coefficients depend on the supplied $f_i$. Current attribute research replaces $\tanh(f_i)$ with $r_i=\tanh(a_i-1)$, as defined at the start of this document.

When enabled:

$$
\begin{aligned}
\gamma_i^0
&= \gamma_0\left[1+c_{\mathrm{tension}}\tanh(f_i)\right],
\\
A_{ij}
&= A_0\left[1+c_{\mathrm{adhesion}}\tanh(f_i)\tanh(f_j)\right].
\end{aligned}
$$

The feedback coefficients $c_{\mathrm{tension}}$ and $c_{\mathrm{adhesion}}$ correspond to `fate_tension` and `fate_adhesion` in the configuration.

The function $\tanh(f)$ is smooth, preserves the sign of fate, and stays between −1 and 1 even if $f$ goes outside that interval. This bounds tension between $\gamma_0(1-c_{\mathrm{tension}})$ and $\gamma_0(1+c_{\mathrm{tension}})$, and similarly bounds attraction. The same-sign product in $A_{ij}$ increases attraction; opposite signs decrease it. These are saturating constitutive choices, not consequences of the Gierer–Meinhardt equations. Setting a contrast to zero removes that coupling; it does not remove signaling or the fate switch.

These coefficients define the isotropic baseline. When apical–basal polarity is enabled, the spatial tension field and its flux-divergence force are given in [graph_signaling.md](graph_signaling.md#apicalbasal-polarity-and-mechanics). Thus A-like cells have higher baseline effective surface tension, and similarly biased cells have stronger interface attraction. These are hypotheses, not established effects of any named genes. Coefficients are bounded to retain positive tension/attraction. Geometry is measured anew before each regulatory step.

## Cleavage and progressive cytokinesis

The default division rule is Hertwig-style shape alignment. For the mother cell, form its occupancy-weighted covariance:

$$
\begin{aligned}
\mathbf{c}_i &= \frac{1}{V_i}\int_{\Omega}\mathbf{x}\,h(\phi_i)\,\mathrm{d}\mathbf{x},\\
\mathbf{C}_i &= \frac{1}{V_i}\int_{\Omega}
(\mathbf{x}-\mathbf{c}_i)(\mathbf{x}-\mathbf{c}_i)^{\mathsf T}
h(\phi_i)\,\mathrm{d}\mathbf{x}.
\end{aligned}
$$

The spindle follows the largest-eigenvalue eigenvector of $\mathbf{C}_i$; the cleavage plane is **perpendicular** to this direction. Eigenvalues within the configured relative tolerance (`axis_degeneracy`, default 0.03) of the largest form a nearly degenerate subspace. A random vector is projected into that subspace and normalized. Thus a sphere has an isotropic fallback, while an oblate cell selects an axis within its long-axis plane, without privileging grid eigenvectors. `division_orientation="isotropic"` is retained only as an explicit control. No calibrated tensile-stress tensor is currently available, so this rule uses shape, not a claimed stress estimate.

Cell-cycle durations are sampled around a configured mean. Reaching the scheduled time **starts** cytokinesis and reserves a future cell slot. It does not replace the mother, alter its occupancy, or create a new contact. The initial plane offset bisects the mother's occupancy. The axis and offset are held fixed relative to the moving cell centroid throughout the event.

### Mechanical furrow

Let $s$ be signed distance from the cleavage plane, $r_\perp$ the distance from the spindle axis, $t_0$ the onset time, and $T$ the configured constriction duration. The smooth ramp and contracting radius are

$$
\begin{aligned}
u &= \operatorname{clip}\left(\frac{t-t_0}{T},0,1\right),\\
S(u) &= u^2(3-2u),\\
R_{\mathrm{ring}}(t) &= R_0[1-S(u)].
\end{aligned}
$$

The dimensionless progress $u$ runs from zero to one over nominal duration $T$ (`cytokinesis_duration`, 0.9). $S$ is the same cubic smoothstep used for occupancy, but here its argument is progress rather than a phase field: $S(0)=0$, $S(1)=1$, and both endpoint slopes vanish. This avoids an abrupt onset or stop of the prescribed ramp. $R_{\mathrm{ring}}(t)$ is the ring radius, distinct from the repulsion coefficient $R$ in the baseline energy.

Here $R_0$ is the largest radial extent of the mother's $\phi_i\ge 0.5$ interior at onset. A contracting annular potential, confined near the equator, drives furrow ingression:

$$
\begin{aligned}
B(s) &= \exp\left(-\frac{s^2}{2\epsilon^2}\right),\\
H(r_\perp,t) &= \frac12\left[1+\tanh\left(\frac{r_\perp-R_{\mathrm{ring}}(t)}{\epsilon}\right)\right],\\
E_{\mathrm{furrow},i}(t) &= \kappa S(u)\int_{\Omega}
B(s)H(r_\perp,t)h(\phi_i)\,\mathrm{d}\mathbf{x},\\
F_{\mathrm{furrow},i} &= -\kappa S(u)B(s)H(r_\perp,t)h'(\phi_i).
\end{aligned}
$$

Read the product as three spatial/temporal selectors: $S(u)$ gradually turns forcing on; $B(s)$ localizes it near the cleavage plane; and $H$ selects the equatorial material outside the current ring radius. As that radius shrinks, the selected region moves inward. The positive energy cost discourages occupancy there, yielding the negative furrow force. $\kappa$ (`ring_strength`, 8) sets its amplitude. A larger $\kappa$ strengthens the drive but does not guarantee a resolved neck or accurate division at the current time step.

The spatial ring coordinates are held fixed during each mechanical substep and recomputed from the centroid at the next step. This is a **prescribed contracting-ring surrogate**, not a resolved actomyosin network or a calibrated line-tension/fluid model. The potential acts on the diffuse interface through $h'(\phi)=6\phi(1-\phi)$ and displaces cytoplasm from the advancing equatorial furrow. The ring strength is ramped from zero, avoiding a sudden force at onset.

A dividing mother retains its measured onset volume $V_i(t_0)$. After each explicit mechanical update to $\widetilde\phi_i$, a scalar interface-local correction enforces this constraint:

$$
\begin{aligned}
\phi_i^{n+1}
&=\operatorname{clip}\left(
\widetilde\phi_i+\alpha_i h'(\widetilde\phi_i),0,1\right),\\
\int_{\Omega}h(\phi_i^{n+1})\,\mathrm{d}\mathbf{x}
&=V_i(t_0).
\end{aligned}
$$

Here $\alpha_i$ is a per-step correction found by the constraint solve, not the fixed polarity response rate $\alpha$. Multiplication by $h'$ preferentially changes interface values, while the integral condition specifies how much occupancy must be restored. The scalar $\alpha_i$ is solved numerically, with relative volume tolerance $10^{-7}$ before float32 storage. `volume_projection_max` reports the largest correction coefficient per step. This is a numerical incompressibility constraint during cytokinesis, not an inferred hydrostatic pressure. Nondividing cells retain the original soft volume penalty. The measured onset volume can already differ slightly from the target, so the constraint preserves that existing error rather than causing a corrective jump at onset.

### Abscission and daughter inheritance

A completed timer alone cannot trigger abscission. After $T$, the mother remains active until both the maximum phase value in the resolved neck band and the normalized prospective daughter overlap fall below their configured thresholds (`neck_threshold=0.15`, `division_overlap_tolerance=0.002`). The neck band has half-width $\max(\epsilon/2,\Delta x/2)$. Both prospective lobes must carry more than 10% of the mother's volume. Unresolved events continue to constrict; `overdue_divisions` flags events older than $3T$, without forcing a cut.

Only then are daughter occupancies constructed:

$$
\begin{aligned}
w(s)&=\frac12[1+\tanh(s/\epsilon)],\\
h_1&=h(\phi_i)w(s),\qquad h_2=h(\phi_i)[1-w(s)],\\
\phi_1&=h^{-1}(h_1),\qquad \phi_2=h^{-1}(h_2).
\end{aligned}
$$

This final relabeling conserves occupancy pointwise; unlike the previous implementation, it happens after the mother has mechanically formed a thin neck. The overlap criterion uses $\int\phi_1^2\phi_2^2\,\mathrm{d}\mathbf{x}/V_i(t_0)$. The contact graph necessarily changes when two cell IDs appear, but overlapping daughter fields are no longer introduced into a full-width uncontracted mother.

Let $\theta=V_1/(V_1+V_2)$. Daughter target volumes follow the actual lobe fractions, preserving the mother's total target and each lobe's preexisting relative volume error:

$$
V_1^{\mathrm{target}}=\theta V_i^{\mathrm{target}},
\qquad
V_2^{\mathrm{target}}=(1-\theta)V_i^{\mathrm{target}}.
$$

There is no growth between divisions. A partitioning draw $\eta$ gives activities $f_1=f_i+2\eta(1-\theta)$ and $f_2=f_i-2\eta\theta$, preserving target-volume-weighted regulatory activity. This reduces to equal-and-opposite perturbations for equal lobes. Regulatory dynamics need not conserve this activity afterward.

Lineage records distinguish `division_start` from `division` (abscission), and store the spindle axis and neck/overlap diagnostics at completion. Daughter clocks begin at abscission. The cell cap includes reservations for active mothers, including non-power-of-two caps.

Independent mechanical and regulatory random streams remain separate. Coupled controls can nevertheless differ in cleavage orientation and completion timing because shape and mechanical relaxation now determine these events. Full active-event state is checkpointed. Checkpoint schema 3 includes signaling, polarity, and graph-event state. New checkpoints record `signal_transport`; historical schema-3 files lacking it explicitly restore `random_walk`, preserving their original dynamics; earlier checkpoints are rejected explicitly because their subsequent dynamics would differ.

## Shape diagnostics

The embryo's diffuse occupancy is

$$
\rho(\mathbf{x})=\min\left(\sum_i h(\phi_i(\mathbf{x})),\,1\right).
$$

Its volume-weighted covariance has eigenvalues $\lambda_1\le\lambda_2\le\lambda_3$. Define their mean as $\bar{\lambda}=(\lambda_1+\lambda_2+\lambda_3)/3$. The shape measures are

$$
\begin{aligned}
\text{Axis ratio}
&= \sqrt{\frac{\lambda_3}{\lambda_1}},
\\
\text{Asphericity}
&= \frac{3}{2}
\frac{\displaystyle\sum_{k=1}^{3}(\lambda_k-\bar{\lambda})^2}
{\displaystyle\left(\sum_{k=1}^{3}\lambda_k\right)^2}.
\end{aligned}
$$

Covariance eigenvalues measure squared spatial extents along principal directions, so the square root converts their ratio into a length ratio. Using the capped union $\rho$ avoids counting the same location more than once when diffuse cells overlap. Asphericity measures inequality among all three extents and is scale-independent. A sphere has ratio $1$ and asphericity $0$. Neither measure identifies the cause of asymmetry.

Fate separation is the distance between volume-weighted A and B cell centroids divided by embryo radius of gyration. It is undefined when either group is absent and cannot detect concentric inner/outer sorting on its own. Surface area and sphericity are not yet computed.

Boundary occupancy, minimum equivalent cell radius in grid spacings, aggregate/cell volume errors, and clipping are numerical diagnostics. Outputs sample them at `save_every`; unsampled transients may be larger.

See [cytokinesis validation and limitations](cytokinesis.md) for comparisons with the previous cleavage rule.

## References motivating the architecture

- [Maître et al., 2016](https://www.nature.com/articles/nature18958): contractility, positioning, and fate in mouse embryos.
- [Zhu et al., 2017](https://www.nature.com/articles/s41467-017-00977-8): cell polarization and early embryonic symmetry breaking.
- [MorphoSim](https://pmc.ncbi.nlm.nih.gov/articles/PMC9938209/): a multicellular phase-field framework for embryonic morphology.
- [Nissen et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286147/): polarity interactions and emergent 3D tissue morphology.

- [A cytokinetic ring-driven cell rotation achieves Hertwig’s rule (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11194556/): cytokinetic ring mechanics and Hertwig-style orientation in early development.
- [Furrow Constriction in Animal Cell Cytokinesis (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3907238/): furrow constriction, cortical mechanics, and volume conservation.

This implementation is a simplified independent model, not a reproduction of these papers' equations or results.


## Experimental fertilization-inspired polarity cue

The optional `embryo.fertilization_cue.CueSimulation` extends the attribute model without changing the shared production equations or adding fate labels. It tests whether a brief directional bias in the zygote influences subsequent organization. This is a phenomenological cue, not a simulation of sperm penetration, calcium waves, PAR domains, or intracellular transport.

For entry direction $\mathbf n_s$ (a unit vector), the additional polarity drive is

$$
\left.\frac{d\mathbf p}{dt}\right|_{\mathrm{cue}}
=\eta g(t)\mathbf n_s,
\qquad
 g(t)=\begin{cases}
 \dfrac{2}{T}\sin^2\!\left(\dfrac{\pi t}{T}\right), & 0<t<T,\\
 0, & \text{otherwise}.
 \end{cases}
$$

The pulse integrates to one, so $\eta$ is the total imposed polarity drive, not the final polarity magnitude: intrinsic damping and geometry-dependent polarity evolution still act. $T$ is the duration in model-time units. The positive sign points polarity toward the nominal entry site and lowers cortical tension there under the existing tension law. This sign is an experimental modeling assumption. The exact integral of the pulse over each timestep is added after the existing polarity update and before mechanics, a first-order split scheme checked at three timesteps. The extension rejects forcing during division or a polarity magnitude exceeding one.

Durations 0.75 and 1.5 finish before the earliest scheduled division time 1.76. Both daughters inherit the maternal vector through the existing division prolongation; no entry-site direction is reapplied after the pulse. Polarity need not persist: decay and subsequent local cues can erase or redirect it. Chemical concentrations remain well mixed within each cell. At the one-cell stage, there is no intercellular chemical diffusion mode, so this experiment tests a mechanical/polarity route into later signaling asymmetry.

The seed-7 pilot has ten arms: no cue; integrated strengths 0.1 and 0.4 along x; strength 0.4 along -x and the normalized (1,1,1) diagonal; a longer pulse at the same integrated strength; and matched cue/no-cue pairs with polarity-driven tension disabled or chemistry externally clamped to uniform. The clamp also suppresses chemical partition noise and is not a selective molecular inhibition. The zero-polarity-tension pair retains polarity dynamics but uses the same scalar mechanics kernel in both arms, preventing kernel-selection roundoff from becoming an unintended perturbation.

All developmental arms use grid 72³, extent 2.24, timestep 0.00375, a sixteen-cell cap, and horizon 90. A short screen first tests the strong x cue through time 1.5 at timesteps 0.0075, 0.00375, and 0.001875. Coarser endpoints must agree with the finest within 0.002 absolute polarity, 0.002 relative centered-covariance norm, and 0.002 absolute centroid displacement. Standard volume, clipping, boundary, and minimum cell-radius checks apply. The delivered pulse must equal 0.4 and no cleavage may have begun. Passing this pre-cleavage screen does not establish numerical convergence of the full developmental trajectory.

Measurements distinguish translation, centered shape covariance, chemical contrast, and orientation. A shape axis is not reported when the relative gap between its two largest covariance eigenvalues is at most 0.01. A signed chemical dipole relative to the cue direction is reported only when log-activator standard deviation exceeds 0.001 and normalized dipole strength exceeds 0.01. These are diagnostic visibility cutoffs, not identity definitions. A seed-7 directional comparison is an exploratory grid-sensitivity probe, not evidence of population-level axis selection or rotational invariance; those require independent histories and spatial refinement.

The early screen runs in one process. Full developmental runs are queued behind the existing moving-response refinement and moving-exchange study, and do not launch if those prerequisites fail. Five focused implementation tests passed: pulse normalization/shutoff, zero-cue equivalence, matched no-coupling dynamics, exact checkpoint continuation, and undefined axes for a symmetric zygote. Outputs are in `outputs/fertilization-cue/`.

```bash
python -m embryo.fertilization_cue prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.fertilization_cue run
```

The biological motivation is sperm-associated initiation of polarity and asymmetric cortical mechanics in C. elegans ([Cowan and Hyman, 2004](https://pubmed.ncbi.nlm.nih.gov/15343338/); [Munro et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15363415/)). It is not assumed to specify mammalian lineages: experiments in mouse implicate pronuclear configuration and confinement geometry, illustrating why fertilization cues must be distinguished from predetermined identities ([Hiiragi and Solter, 2004](https://pubmed.ncbi.nlm.nih.gov/15254539/); [Kurotaki et al., 2007](https://pubmed.ncbi.nlm.nih.gov/17446354/)).

## Opt-in mechanics optimization and validation

The following sections preserve the sequence of C, C++/OpenMP, and GPU benchmarks and execution transitions. Their timings describe the workloads/load at each stage. For the current scientific path, use the [computing summary](#coupled-integration-and-computing-paths) and [accepted resident GPU backend](#resident-gpu-backend-validation); component benchmarks alone do not authorize scientific migration.

`FastAttributeSimulation` in `embryo/fast_mechanics.py` is an optional backend for the attribute model. It combines polar cortical tension and conservative face-flux calculations in a C99 kernel and caches occupancy, volumes, and centers only within a step while the phase-field array is unchanged. The cache is invalidated when geometry changes and discarded at the end of every step. Chemistry, transport, timestep, material laws, and boundary conditions are unchanged. Active cytokinesis uses the original mechanics/ring/projection code. The production classes and running scientific protocols have not been modified.

The kernel requires a C99 compiler (`cc`) on first use and is compiled with `-O3 -ffp-contract=off`, without fast-math. The compiled shared library is cached under `outputs/mechanics-kernel-cache/` using the C source/platform hash. This prototype currently targets Unix-style shared libraries. The C source is included in package data; no new Python dependency is required. Checkpoints retain the attribute format and must be explicitly restored with the optimized class to select the optimized backend.

A saved 16-cell 72³ state was benchmarked with 48 interleaved original/optimized steps under concurrent production load. Original time averaged 1.5583 seconds/step and optimized time 1.0245 seconds/step: **1.52× faster**, approximately 34% less step time. Cell fields, activator/inhibitor concentrations, polarity, and measured volumes matched exactly over this short trajectory. This is an empirical result on one state, not a guarantee of bitwise equivalence in every configuration. Observation and checkpoint costs are excluded. Results and source hashes are in `outputs/fast-mechanics-benchmark/comparison.json`.

A short concurrency screen measured 0.844, 1.332, 1.660, and 1.850 aggregate benchmark steps/second with one, two, three, and four benchmark workers respectively. Two production simulations remained active in addition. The timings include process startup, checkpoint reading, and warmup, and are not isolated hardware-scaling measurements. Four benchmark workers delivered about 2.19 times one-worker throughput while making individual jobs slower. These results support bounded concurrency and indicate diminishing returns; they do not justify treating twelve hardware threads as twelve independent fast workers.

Six focused tests cover reference flux agreement and conservation, cache invalidation, checkpoint continuation, unchanged division mechanics, completed abscission, and the full validation runner on a short matched-control fixture. An additional scientific gate in `embryo/validate_fast_mechanics.py` repeats a complete time-210 unexchanged control and cell-27 -10% pulse through time 270, comparing against already completed original-kernel trajectories. It checks chemistry, polarity, volumes and axis ratio at every observation, endpoint phase fields, pulse-minus-control response, response integral, and recovery times. It uses one extra worker; original studies continue with their immutable code.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.benchmark_fast_mechanics --output outputs/fast-mechanics-repeat --steps 48 --concurrency
python -m embryo.validate_fast_mechanics prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.validate_fast_mechanics run
```

The full-horizon gate completed and passed: sampled chemistry, polarity, volumes and axis ratio, endpoint phase fields, and matched pulse-response/recovery metrics agreed exactly for both 60-unit continuations. On 2026-10-02 the identity-response study was resumed using `embryo.fast_response_resume`, with the change recorded in a separate `backend-transition.json` and per-job `backend-execution.json` files. Original scientific protocols and completed results remain unchanged; unfinished work resumes from its atomic checkpoint. Nine focused tests passed before this transition, including cross-backend restart. This validates the tested mature control/pulse regime, not all exchanged states or developmental trajectories. The fertilization-cue study remains on the original backend.

### Second optimization: compiled surface-polarity cue

A subsequent profile of the accepted fast-mechanics backend attributed approximately 25% of step time to the surface exposure cue. `FastPolaritySimulation` in `embryo/fast_polarity.py` adds a separate opt-in C kernel for this calculation. It fuses gradient, surface normal, free-surface weighting, and accumulation without constructing full normal-vector fields. It reuses the unchanged occupancy field from the step-local cache. The original fast backend and running scientific source hashes remain unchanged.

Unlike the first optimization, this kernel accumulates float32 per-voxel contributions in float64. This changes reduction order and can introduce small differences from NumPy's reduction, so exact arithmetic equivalence is not claimed. Spatial stencils, boundary one-sided gradients, surface weights, and the symmetric-cue cutoff remain unchanged. C compilation disables fast-math and fused contraction.

Nine focused tests passed, including random and contacting fields, empty/full fields, symmetric cells, checkpoint continuation, clamped chemistry, original mechanics/cleavage checks, and the full validation runner on a short two-cell control/pulse fixture. `embryo.benchmark_fast_polarity` compares against the accepted first optimization; its reported speedup is incremental, not relative to the initial Python/NumPy implementation. `embryo.validate_fast_polarity` separately checks complete saved original control/pulse trajectories with the established tolerances. Results are stored in `outputs/fast-polarity-benchmark/` and `outputs/fast-polarity-validation/`. This second backend is not enabled in either running scientific study before acceptance of its full-horizon gate.

The first surface-cue prototype passed its short accuracy gate but gained only 1.9% in whole-step throughput. A revised nested-loop layout removed per-voxel integer index division and improved the measured incremental speedup to 1.123× (0.9351 to 0.8329 seconds/step over 24 interleaved steps under concurrent load, roughly 11% less step time). Maximum cell-field difference was 5.96e-8; chemical log difference was 1.59e-9. These small differences are not zero. The initial and revised benchmark artifacts are retained separately. Four directly affected tests were rerun after the loop revision and passed. This was an intermediate optimization stage; the current resident GPU path uses its own complete native-CPU reference gate rather than treating this short benchmark as acceptance.


## Threaded C++ backend and optimization shutdown (2026-10-02)

At the user's request, all scientific workers and the dashboard were stopped before further optimization. Completed outputs and atomic checkpoints were retained. The shutdown record is `outputs/optimization-shutdown-20261002T163919Z/`; the older surface-cue validation described above was interrupted, not accepted. No biological study has been restarted by this optimization work.

`embryo.native_mechanics.NativeSimulation` is an opt-in implementation of the attribute model using C++17 and OpenMP. It combines directional surface tension, conservative face fluxes, double-well, volume, repulsion, and adhesion forces in one compiled mechanical update. Independent cell fields run in parallel. The exposure calculation compiles its voxel arithmetic but retains the original NumPy float32 reduction, avoiding the accumulation change in the earlier `fast_polarity` prototype. Occupancy, centers, volumes, and interface weights are reused within a step. Reaction–transport equations, timestep, spatial resolution, and scientific parameters are unchanged. Active cytokinesis and nonpolar mechanics retain the existing backend.

Compilation uses `-O3` with floating-point contraction disabled and without fast-math. A C++ compiler and OpenMP runtime are required on first use; the shared library is cached by source hash. Keep BLAS at one thread to prevent nested thread oversubscription. The `native_threads` execution setting is intentionally separate from model checkpoints and must be set explicitly after restore. The compatibility subclass `embryo.native_cue.NativeCueSimulation` preserves the fertilization experiment's cue, clamp, and checkpoint hooks.

The idle-machine, interleaved 128-step benchmark used the saved mature 16-cell, 72³ state and timestep 0.00375. Compilation, checkpoint IO, and comparison diagnostics were outside the timed steps:

| Backend | Mean seconds per step | Speedup vs reference |
|---|---:|---:|
| Original attribute reference | 0.614 | 1.00× |
| Previously accepted C mechanics | 0.516 | 1.19× |
| New C++/OpenMP, four threads | 0.198 | 3.09× |

All sampled cell-field, chemical, polarity, and volume differences were exactly zero throughout the 128 steps. This is local numerical agreement on this machine and state, **not full developmental or long-horizon response validation**. The earlier 1.52× and 1.12× measurements were made under concurrent production load and must not be multiplied with these isolated measurements. Evidence and code hashes are in `outputs/native-mechanics-benchmark/comparison.json`.

Reproduce the step benchmark in a fresh output directory while other simulations are stopped:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.benchmark_native_mechanics \
  --output outputs/native-mechanics-benchmark-repeat --steps 128
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.benchmark_native_parallel \
  --output outputs/native-parallel-repeat.json --steps 32
```

Use the backend explicitly in a new run or a separately audited continuation:

```python
from embryo.native_mechanics import NativeSimulation
sim = NativeSimulation.restore("path/to/attribute-checkpoint.npz")
sim.native_threads = 4
sim.step()
```

Before replacing the backend of a scientific study, run the separate full control-and-pulse comparison against saved reference trajectories. Its protocol records thread count, source hashes, input hashes, and the existing numerical acceptance criteria. Preparing this validation does not resume any study; running it performs two independent validation replays.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.validate_native_mechanics prepare --threads 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.validate_native_mechanics run
```

The full-horizon native gate has not been run. Existing production protocols and backend transition manifests have not been rewritten to silently select the new implementation.


On the i7-8700 (six physical cores), two reversed-order thread sweeps averaged 0.344, 0.248, 0.208, and 0.192 seconds/step for one, two, four, and six native threads respectively. Use four to six threads for a single mature trajectory; four leaves more headroom for interactive work. This scaling is below ideal because contact/adhesion operations, reductions, and memory traffic remain substantial.

A separate synchronized, 32-step-per-worker screen found:

| Concurrent workers | Native threads per worker | Aggregate steps/second |
|---:|---:|---:|
| 1 | 4 | 4.72 |
| 2 | 2 | 6.59 |
| 3 | 2 | 8.52 |
| 4 | 1 | 8.06 |
| 2 | 1 | 5.26 |

Thus **three independent workers with two native threads each** is the provisional batch allocation for this machine. Keep BLAS at one thread per worker. Use the multiprocessing `spawn` context; do not fork a process that has already initialized OpenMP. The screen excludes restore/warmup/startup time and uses the same mature checkpoint for every worker; actual experiment throughput also includes diagnostics and checkpoint writes. All concurrent copies agreed exactly in final chemistry. Evidence is in `outputs/native-mechanics-benchmark/parallel.json`. More workers are not automatically faster.

The completed regression run passed **66 tests**, covering native mechanics and exposure, reference-model behavior, cleavage through abscission, restart, prescribed chemistry, fertilization cue/clamp hooks, and the validation runner on matched short references. The full native control/pulse protocol is prepared in `outputs/native-mechanics-validation/` but has not been executed. At completion of this optimization session, no project simulation or benchmark process remained running.


### Full native validation accepted and moving-response study resumed

The two 60-unit native validation replays have now completed with passing quality. At every saved observation, chemistry, polarity, volumes, and axis ratios matched the reference exactly; both final phase fields also matched exactly. Normalized response error, relative response-area error, and target/network recovery-time errors were all zero. The perturbed-cell and network recovery times remained 3.15 model units. Evidence: `outputs/native-mechanics-validation/comparison.json` and the two job `result.json` files. This extends acceptance to these mature control/pulse continuations; it does not establish equivalence for every developmental regime.

The next scientific test is the existing moving chemical-exchange response study. Six of fifteen results were retained, and the remaining nine continuations were resumed by `embryo.native_response_resume` with three spawned workers and two native threads each (BLAS one thread). Original scientific protocols and the previous backend manifest remain unchanged; `native-backend-transition.json` separately records validation/code hashes, accepted result hashes, and initial resume-checkpoint hashes. Each migrated run records its own `native-backend-execution.json`. The adapter retains the established observation, quality-control, checkpoint, and final assessment logic. Sixteen focused native/restart tests passed before launch. Completion triggers the existing donor-versus-destination response assessment automatically. The dashboard and fertilization study have not been restarted.


### CUDA precision and mechanics pilot

A separate benchmark used the otherwise idle GTX 1080 Ti (11 GB, sm_61) while all three native response-study workers continued running. It used CUDA 12.0 and GCC 12, targeting sm_61 with fused multiply-add disabled, denormal flushing disabled, and precise division/square root. No production classes or source files were changed. The opt-in implementation is `embryo/cuda_mechanics.cu`, driven by `embryo.cuda_benchmark`. It compiles with the existing nvcc installation; it does not require CuPy.

The pilot evaluates two kernels: spatially varying cortical tension and the mechanical force/update. The mixed-precision mode retains float32 phase fields and reference float32 intermediate products, with float64 tension/flux/force arithmetic. The alternative performs those tension/flux/force calculations in float32, reading the same prepared inputs. Neither mode ports contact construction, geometry reductions, polarity cues, chemistry, or cytokinesis to GPU. These remain on CPU, so this is not a complete GPU simulation backend.

The saved mature 16-cell, 72³ state was used for 100 same-input repetitions timed with CUDA events after warmup. Separate host-call measurements include allocation, transfers, warmup and one measured kernel pair. A further comparison advanced three coupled trajectories (CPU and both GPU modes) for 48 steps at dt=0.00375, covering only 0.18 model-time units.

| Measurement | Reference mixed precision on GPU | Float32 force arithmetic on GPU |
|---|---:|---:|
| GPU-resident mechanical kernel pair | 4.535 ms | 1.232 ms |
| Host call including transfers/allocation | 54.33 ms | 45.60 ms |
| Largest phase-field difference across 48 steps | 0 | 5.96e-8 |
| Largest chemical log difference | 0 | 2.34e-8 |
| Largest relative volume difference | 0 | 2.34e-8 |
| Largest absolute polarity difference | 0 | 5.22e-10 |
| Median complete hybrid step | 0.4695 s | 0.4623 s |

The simultaneous single-thread CPU comparison took 0.5162 s per complete step, giving only about 1.10×–1.12× whole-step throughput in this pilot. Its mechanical kernel alone took 120.8 ms. These host timings were collected under three active production workers and additional low-priority benchmark CPU work; they must not be compared directly with the earlier idle four-thread 0.198 s/step result. GPU resident kernel timings exclude host preprocessing, transfers, allocation, and all other model work. Repeating frozen inputs does not test long-term dynamics; the distinct coupled comparison supplies the limited trajectory evidence.

Six opt-in GPU tests passed, exercising empty, full, and random voxel fields with displaced polarity, nonuniform forces and all stencil boundaries in both precision modes. Exact mixed-precision agreement here does not imply future GPU reductions will preserve CPU summation order. Neither the short trajectory nor the stencil tests establish full developmental or long-duration precision acceptance.

Evidence: `outputs/cuda-precision-pilot/comparison.json`. Reproduce in a fresh directory:

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
python -m embryo.cuda_benchmark --output outputs/cuda-precision-repeat --steps 48 --repeats 100

EMBRYO_CUDA_TESTS=1 CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest tests/test_cuda_mechanics.py -q
```

These measurements motivated retaining reference mixed precision while keeping spatial arrays, contact/adhesion operations, geometry reductions, and surface cues on GPU across steps. That resident architecture has since passed the full gate below. Float32 force arithmetic remains an experimental benchmark mode and is not the accepted scientific force arithmetic.


### Additional GPU benchmarks before migration

`embryo.cuda_spatial_bench` measures the spatial work that remained on CPU in the first pilot. It uses separate CUDA kernels for occupancy/interface arrays and deterministic double-precision reductions of volumes, centers, and exposure cues; cuBLAS computes float32 contacts and float64 adhesion. The underlying per-voxel normal calculation remains float32. It compares all outputs with the reference implementation. These are snapshot microbenchmarks, not an integrated resident GPU timestep.

The saved states were a one-cell cue-screen endpoint at t=1.5 and sixteen-cell states at t=90, 210, and 270. At t=270, resident component timings were:

| Spatial operation | GPU time |
|---|---:|
| Occupancy, interface arrays, occupied field | 0.395 ms |
| Contact matrix | 0.603 ms |
| Adhesion field | 0.933 ms |
| Volumes, centers, exposure cues and reductions | 1.410 ms |
| Sum of separately measured components | 3.341 ms |

The other sixteen-cell states totaled 3.44–4.10 ms. Across the four saved states, maximum relative volume error was below 9e-16, absolute center error below 2.2e-15, absolute cue error below 5.4e-8, and relative contact-matrix Frobenius error below 2.9e-7. Adhesion error stayed below 1.7e-16. The mixed-precision mechanical kernel separately remained exact in phase-field output at all four snapshots. GPU contacts and cue reductions do **not** preserve exact CPU reduction arithmetic, so the earlier exact mechanical-kernel result cannot be extrapolated to a complete GPU backend.

Synthetic, resampled/replicated fields gave spatial-component totals of 1.02 ms for 16 cells on 48³, 7.84 ms for 16 cells on 96³, and 8.21 ms for 32 cells on 72³. These fields test workload scaling; they are not physical resolution or cell-count convergence experiments. Totals exclude mechanics, chemistry, quality checks, checkpointing, transfers and allocation.

Synchronous transfer tests used 20 copies in each direction, after initialization:

| Buffer | Host memory | CPU→GPU | GPU→CPU |
|---|---|---:|---:|
| 24 MiB | Ordinary pageable | 9.15 ms | 5.86 ms |
| 24 MiB | Pinned | 7.49 ms | 4.73 ms |
| 96 MiB | Ordinary pageable | 35.74 ms | 26.61 ms |
| 96 MiB | Pinned | 29.74 ms | 18.88 ms |

Pinned memory reduces the measured round-trip times by roughly 19–22%, but transfers still outweigh resident spatial work. These measurements were made under the active CPU study load, not as a peak PCIe bandwidth test. The 1080 Ti reported a PCIe generation-3, eight-lane link during the check. Keeping spatial arrays on GPU is more consequential than merely changing the host-buffer allocator.

The extended precision assay in `embryo.cuda_precision_extended` replays the first 160 steps (0.6 model-time units) of the saved t=210 control and negative-pulse trajectories in two precision modes. It includes GPU contacts, volume/center/cue reductions, and GPU mechanics. Small chemical systems remain on CPU; post-mechanical volumes use the CPU implementation until the next contact refresh. It is intentionally transfer-heavy and is an accuracy assay, not a performance backend. Limits declared before the run are chemical log error 1e-5, polarity absolute error 1e-5, relative volume error 1e-5, relative axis-ratio error 1e-4, normalized pulse-response error 0.001, zero clipping, and maximum target-volume error below 5%. This short prefix cannot establish full recovery times or long-term pattern preservation.

Nine opt-in GPU tests passed, covering the original mechanical stencils plus spatial reductions and an asymmetric adhesion matrix that detects accidental matrix transposition. Production study source hashes are kept unchanged. Artifacts are `outputs/cuda-spatial-bench/` and `outputs/cuda-precision-extended/`.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
python -m embryo.cuda_spatial_bench --output outputs/cuda-spatial-repeat

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
python -m embryo.cuda_precision_extended --output outputs/cuda-precision-extended-repeat --steps 160
```

All four extended accuracy replays completed and passed. Across both control/pulse pairs and both precision modes, the largest chemical log error was 1.12e-07, polarity error 8.62e-09, relative volume error 2.16e-08, and relative axis-ratio error 3.49e-09. Maximum normalized pulse-response error was 3.61e-07, below the predeclared 0.001 limit. Clipping remained zero and maximum target-volume error was about 1.10%. This supports proceeding to a resident mixed-precision prototype while retaining full-duration control/pulse validation and cleavage checks as requirements before production migration.


### PyTorch GPU benchmark

`embryo.torch_benchmark` compares ordinary PyTorch GPU operations, PyTorch CUDA-graph replay, and the custom CUDA kernels using identical saved one-cell and sixteen-cell 72³ states. It uses the installed PyTorch 2.7.1+cu118 on the GTX 1080 Ti, inference mode, one CPU thread and GPU-resident inputs. TF32 is disabled. There is no Triton/`torch.compile` dependency: CUDA graphs capture and replay ordinary supported GPU operations. The production CPU study continues unchanged.

Each PyTorch timing is the median of three windows of twenty repetitions after warmup, with CUDA events and a synchronized wall-clock cross-check. Custom CUDA operations were rerun in the same session. Host/device copies, first allocation/warmup, and graph capture are excluded. The eager implementation creates temporary tensors during the operations; the custom kernels fuse calculations and preallocate working buffers. This is a comparison of these implementations, not a claim that all possible PyTorch implementations have the same performance.

For the mature sixteen-cell state:

| Operation | PyTorch eager | Custom CUDA |
|---|---:|---:|
| Occupancy/interface arrays | 1.309 ms | 0.393 ms |
| Contact matrix | 0.618 ms | 0.610 ms |
| Adhesion field | 0.939 ms | 0.941 ms |
| Geometry and polarity cues | 8.687 ms | 1.424 ms |
| Mechanical update, reference mixed precision | 17.521 ms | 4.359 ms |
| Mechanical update, float32 force arithmetic | 10.061 ms | 1.241 ms |

The spatial pipeline, excluding mechanics, took 11.564 ms with eager PyTorch and 11.535 ms with CUDA-graph replay: no meaningful improvement on this workload. Graph replay does not fuse tensor operations or remove their intermediate memory traffic. The one-cell pipeline improved from 0.762 to 0.706 ms, where launch overhead has a larger relative effect. Captured and eager outputs matched exactly in both cases.

Mixed-precision mechanical phase-field output matched the CPU kernel exactly at both snapshots. Float32-force maximum phase-field error was 5.96e-8. At sixteen cells, spatial comparisons gave relative volume error 4.44e-16, absolute center error 1.67e-15, absolute cue error 5.36e-8, relative contact-matrix Frobenius error 2.58e-7, and adhesion-field error 1.11e-16. Nine opt-in GPU tests passed, checking empty/full/random mechanical fields, one/two/four-cell spatial reductions, boundary stencils and asymmetric adhesion layout.

These are frozen-input and component benchmarks, not coupled trajectory or full developmental validation. The graph replay repeats fixed inputs and does not advance an embryo. Neither the spatial pipeline nor the sum of component timings is a complete simulation timestep. Timings were measured while three CPU production workers were active; no idle-machine whole-step speedup is asserted.

The results support using PyTorch for GPU array management and matrix operations, with custom CUDA for expensive stencil and geometry/cue operations. Keeping spatial tensors on GPU between operations remains necessary. A direct Python tensor translation is convenient, but it does not match the fused custom kernels here.

This architecture is now implemented and accepted for the mature regime below: PyTorch owns resident arrays and performs contact/adhesion matrix products and small transport-matrix operations; custom CUDA advances mechanics and computes geometry and polarity. Kernel calls share tensor device storage and obey the active CUDA stream, with explicit dtype, layout, device, and lifetime checks. The completed history replication used this gated backend; CPU development and unsupported GPU regimes remain separate. See the [implementation sequence](plan.md#gpu-implementation-design).

Evidence: `outputs/torch-gpu-benchmark/comparison.json`. Reproduce with the compatible installed PyTorch build:

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
python -m embryo.torch_benchmark --output outputs/torch-gpu-repeat --repeats 20

EMBRYO_CUDA_TESTS=1 CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest tests/test_torch_benchmark.py -q
```

### Resident GPU backend validation

`embryo.gpu_backend.GpuSimulation` is a separate, opt-in mature-state backend. PyTorch owns persistent GPU tensors, forms contact and adhesion matrix products, constructs the conservative volume-weighted transport operator, and integrates chemistry with the existing positivity-restricted SSP-RK2 scheme. Stream-aware custom CUDA kernels compute occupancy/interface arrays, excluded-volume fields, volume/center/exposure reductions, polarity evolution, variable-tension mechanical forces, and union-shape moments. No autograd, Triton, or `torch.compile` is used.

Phase fields, occupancy, shell arrays, and contact-product accumulation use float32; shell-squared adhesion buffers, mechanical force arithmetic, geometry accumulation, polarity, transport matrices, and chemistry use float64. The geometry/cue kernel uses deterministic double accumulation, while the native CPU cue has float32 reductions: backend agreement is measured with explicit tolerances rather than assumed bitwise equality. TF32 and fast-math are disabled. Pointer bindings check tensor dtype, shape, contiguity, device, and stream. Tensors retain ownership of device allocations. There are no full-field CPU/GPU round trips within a step; the host receives small diagnostics and explicitly requested standard checkpoints. RNG state and IDs remain intact, and compatibility fate values remain zero.

The accepted environment is PyTorch 2.7.1+cu118 on the dedicated GTX 1080 Ti. The custom library is built with the locally installed `nvcc` and `/usr/bin/g++-12`; the current build targets `sm_61` and disables fused contraction/flush-to-zero. PyTorch is additional to the base CPU dependencies. This is the validated machine/build combination, not a hardware-independent acceptance claim: the adapter checks the device, PyTorch/CUDA version, and compiled-library hash against saved evidence.

The backend currently supports mature polar, direct-feedback attribute embryos with conservative signaling and no active or future cleavage. It rejects unsupported starting regimes. Cue forcing, no-feedback/nonpolar mechanics, and developmental division have not been implemented or validated on GPU.

The `embryo.validate_gpu_backend` gate replays four complete t=210–270 continuations on the 72³ grid at dt=0.00375: moving controls and −10% activator pulses in cell 27 on the unexchanged and fresh-exchange backgrounds. This samples low- and high-activity contexts in the same recipient. Accepted native-CPU trajectories are compared every 0.15 model-time units. Backend and reference sources, inputs, and protocols are hashed; the original scientific studies remain unchanged.

| Agreement measure | Maximum permitted error |
|---|---:|
| Chemical concentration log difference | 1e-5 |
| Absolute polarity difference | 1e-5 |
| Relative cell-volume difference | 1e-5 |
| Relative conservative-operator Frobenius difference | 1e-5 |
| Relative aggregate axis-ratio difference | 1e-4 |
| Final absolute phase-field difference | 2e-5 |
| Pulse/control log waveform difference, normalized by initial log pulse | 0.001 |
| Relative target response-integral difference | 0.002 |
| Recovery time difference, with recovery status unchanged | 0.15 model-time units |
| Relative amount change caused by mechanical dilution alone | 2e-14 |

Every step also requires positive finite chemistry, finite geometry/polarity, volume error below 5%, minimum radius at least four grid spacings, and zero phase-field clipping. Boundary occupancy is screened at every recorded observation and must stay below 0.01. The final assessment independently recomputes trajectory errors, response metrics, recovery, and endpoint field differences.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.validate_gpu_backend prepare

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.validate_gpu_backend run
```

All four full-horizon replays completed and passed. Maximum chemical log difference was 4.92e-7, polarity difference 2.60e-8, relative volume difference 3.46e-8, relative transport-operator difference 3.87e-7, and relative axis-ratio difference 1.63e-8. Every final phase-field comparison had maximum absolute difference 1.19e-7. Volume error remained below 1.101%, minimum radius above 5.1104 grid spacings, clipping zero, boundary occupancy below 5.03e-10, and relative dilution amount error below 4.45e-16.

| Response comparison | Normalized waveform error | Relative response-integral error | Recovery-time differences |
|---|---:|---:|---|
| Unexchanged | 1.56e-6 | 2.59e-5 | Target and network: zero |
| Fresh exchange | 1.12e-6 | 1.05e-6 | Target and network: zero |

Evidence is in `outputs/gpu-backend-validation/comparison.json`, with source/input hashes, complete histories, checkpoints, and per-run audits in the same directory. Twenty-two focused checks passed before preparation; subsequent adapter and PyTorch checks bring the tested set to thirty-nine. Agreement covers the mature regime and parameters above; it does not validate GPU cleavage or the biological interpretation of the transport closure.

`embryo.gpu_response_runner` prepares separate execution directories for mature moving-response protocols, including the current history replication. It refuses to prepare or run until all four validation jobs and response comparisons pass, independently reassesses the evidence without changing accepted files, and admits only validated physical parameters and the GPU/software environment. Seed, age, and horizon can differ; changing grid, timestep, or physical parameters requires another gate. It retains the established history/result schema and restart checkpoints. Adapter acceptance does not establish developmental or biological convergence.

For a prepared child-study directory containing `protocol.json`, `initial_states.npz`, and its source checkpoint:

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.gpu_response_runner prepare \
  --source outputs/NEW-STUDY/CHILD --output outputs/NEW-STUDY-GPU/CHILD

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.gpu_response_runner run \
  --output outputs/NEW-STUDY-GPU/CHILD
```

Replace the example paths with the actual prepared study. This preserves CPU evidence and writes GPU execution to a separate directory. Job clocks must align with observation/checkpoint intervals, initial IDs and compartment masses must match the source, and restarting never repeats the chemical pulse.

A real adapter integration smoke replay also passed: two 0.6-unit control/pulse prefixes agreed with the accepted CPU histories, with normalized response error 3.61e-7. It exercised the actual acceptance gate and standard output/checkpoint path; it is not an additional biological experiment. Evidence: `outputs/gpu-adapter-smoke/comparison.json`. A separate synthetic restart test verifies continuation without reapplying the pulse after loss of the final result file.

The [developmental-history replication](cell_response_moving.md#replication-across-developmental-histories) additionally checks all three t=150 chemical backgrounds on each of seeds 8 and 9 through matched 0.6-unit native/GPU continuations. All six checks passed before formation runs were admitted. Their histories, endpoint fields, discrepancies, and hashes are stored under `outputs/exchange-response-histories/seed-*/prefix/`. These new-context checks supplement full-horizon seed-7 acceptance; they are not full-horizon equivalence or timestep-refinement tests on the new histories.

The [completed exchange timestep study](exchange_response_refinement.md) requires separate acceptance at dt=0.001875. All four full control/negative-pulse GPU replays pass against the existing accepted fine native CPU trajectories; a read-only reference adapter preserves those files. All ten matched 0.6-unit checks pass at the new histories' formation and response starting contexts. The eighteen scientific continuations subsequently pass refinement criteria. Existing dt=0.00375 acceptance is not reused as smaller-timestep validation. Source hashes, precision settings, backend discrepancy thresholds, and quality screens remain fixed. The [neighbor-context assay](neighbor_context.md) uses independent CPU ODE solvers on small frozen graphs; its eventual moving confirmation retains this validated resident GPU architecture.

After full acceptance, a matched 120-step coupled benchmark used the same t=210 state and five warmup steps on both backends. It measured **222.37 ms/step on the four-thread native CPU backend versus 9.97 ms/step on the GTX 1080 Ti: 22.31× throughput**. Both timings include every-step quality audits and observations every 0.15 units. Standard checkpoint writes were measured separately at 1.39 seconds (CPU) and 1.42 seconds (GPU); startup/compilation and a complete assay's I/O are not included in the speedup. Endpoint checks passed. Evidence: `outputs/resident-gpu-coupled-benchmark/comparison.json`.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.benchmark_gpu_backend \
  --steps 120 --threads 4 --output outputs/resident-gpu-coupled-repeat
```
