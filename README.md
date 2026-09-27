# Embryogenesis in 3D

## Scientific objective

**To what extent can mathematical equations explain the emergence of organized living structure from an initially simple zygote?** This project investigates how small fluctuations in an approximately uniform initial state can develop into reproducible cell identities, spatial signals, and complex, changing three-dimensional shapes.

The central mechanism to investigate is an **explicit activator–inhibitor feedback loop**: an activator promotes its own production and the production of an inhibitor; the inhibitor suppresses activation. We will study how the kinetics, interaction ranges, and coupling to cell mechanics determine whether this feedback amplifies fluctuations into order, suppresses them, or produces unstable patterns. Activator and inhibitor are signaling variables; they are distinct from the two initial cell identities, A and B.

We will evaluate three linked outcomes separately:

| Outcome | Scientific question | Evidence to measure |
|---|---|---|
| Cell identity differentiation | Can initially similar cells acquire distinct, persistent identities through local regulatory dynamics and signaling? | Activity distributions, fate persistence, lineage histories, and responses to perturbations |
| Signaling symmetry breaking | Can nearly uniform signaling develop persistent spatial domains, poles, or axes without a prescribed directional cue? | Pattern onset, spatial correlations, wavelength, domain number, and orientation across independent runs |
| Geometry and shape symmetry breaking | Can signaling and cell behavior produce sustained changes in tissue shape while the embryo divides and rearranges? | Shape anisotropy, axis persistence, cell organization, and later folding or cavity formation where the model supports them |

The intended feedback is reciprocal: **signaling influences identity and cell mechanics; cell movement, division, contacts, and shape alter signaling in return.** The embryo boundary must evolve with these interactions.

## How we will assess explanatory power

The aim is to identify the smallest interpretable set of equations that accounts for progressively more developmental organization, and to document where that explanation fails. Two identities and a small cell population are the starting point; more complex structures are later tests of the same approach.

- **Separate emergence from assumptions.** Record which features arise from the dynamics and which are supplied through initial conditions, boundary conditions, cell-cycle rules, or prescribed forces. Begin spontaneous-symmetry-breaking experiments with unbiased fluctuations; label imposed gradients and asymmetries as separate controls.
- **Establish the mechanism.** Analyze steady states and their stability, identify parameter regimes where spatial perturbations grow, and test the necessity of activation, inhibition, signal transport, and mechanical feedback by disabling them individually. An activator–inhibitor loop does not automatically imply a Turing instability or a single developmental axis.
- **Measure robustness and limits.** Compare independent seeds, parameter ranges, perturbations, and numerical resolutions. Report uniform, mixed, fragmented, and failed outcomes as well as organized ones; distinguish reproducible structure from transient noise and grid artifacts.
- **Increase biological specificity only with evidence.** Start with a generic, dimensionless model. Quantitative claims about a particular embryo will require experimental calibration and independent validation. Active signaling, proliferation, and force generation remain explicit biological inputs to the model.

This objective guides subsequent model changes and experiments. Each extension should state the phenomenon it aims to explain, the feedback it introduces, a control that could challenge the proposed explanation, and a measurable success criterion. Visual resemblance alone is insufficient.

## Current starting point

The runnable prototype models one cell dividing into a deformable multicellular aggregate. Gierer–Meinhardt activator and inhibitor concentrations evolve through conservative volume-weighted contact transport and drive a downstream bistable fate switch. Apical–basal polarity develops from exposed cortex and neighboring orientations, and changes cortical tension directionally.

**Finite-graph stability is analyzed before 3D evolution.** The solver uses a conservative volume-weighted Laplacian, checks each discrete eigenvalue against the reaction–diffusion Jacobian, and records changes in spectral modes at every cleavage. A continuous unstable band can contain no supported modes on a small graph. Activator/inhibitor variables are distinct from cell fate; the old independent fate noise and geometry biases are disabled by default. See [live conservative signaling](docs/live_transport.md) for equations, spectral convergence, and limitations, and [polarity mechanics](docs/graph_signaling.md#apicalbasal-polarity-and-mechanics) for the mechanical coupling.

This is an exploratory model in dimensionless units, not a reconstruction of a particular organism. The simulation evolves each cell's shape on a 3D grid. It does not prescribe an embryo outline or assign daughter identities.

## New experiment: identities as emergent cellular phenotypes

The [attribute-based developmental experiment](docs/attribute_development.md) is now running two fresh zygote-to-t=90 branches **without the prescribed bistable fate switch or A/B labels**. It retains activator–inhibitor signaling and compares continuous activity-dependent mechanics with no regulatory mechanical feedback. Chemical activity, polarity strength, and cell shape are recorded as attributes; position, lineage, and exposure are context. No number of identities is prescribed and no clustering is imposed. Persistent differences would initially be candidate phenotypes, requiring perturbation and robustness tests before being called identities.

This experimental branch leaves the historical core/dashboard model unchanged. The method section below describes that historical baseline, including the supplied fate switch, so its assumptions and prior results remain reproducible. The new branch's equations, compatibility details, and protocol are documented separately. Its scientific results are pending.

## Method and mathematical model

Read the equations as a set of coupled rules: **geometry determines contacts; contacts transport signals; signals bias fate and polarity; fate and polarity change mechanics; mechanics changes geometry**. The model does not assume that any one of these links is sufficient to produce an organized embryo.

In the notation below, $i,j$ label cells, $\mathbf{x}$ is a location in the 3D box $\Omega$, a dot means a time derivative, and an integral adds a quantity over the box. A star on $V_i^\star$ denotes a prescribed target, not a measured volume. The occupancy function $h(\phi)$ is unrelated to the inhibitor, which is written $b_i$ here (`inhibitor` and sometimes `h` in the code). Adhesion strength is $A_{ij}$; geometric contact area is $\mathcal A_{ij}$. See [the annotated model](docs/model.md#reading-the-equations) for derivative notation and a term-by-term mechanical derivation.

### State variables and initial conditions

The simulation combines a continuous description of cell shape with a discrete description of signaling between cells. Each cell occupies a deformable region of a common three-dimensional computational domain; the contact network is reconstructed from those regions as they move and divide.

| Variable | Meaning | Initial zygote |
|---|---|---|
| $\phi_i(\mathbf{x},t)$ | Diffuse cell indicator: approximately one inside cell $i$, zero outside | Smooth sphere of radius 0.8 |
| $a_i(t), b_i(t)$ | Nonnegative activator and positive inhibitor activities | Both one |
| $f_i(t)$ | Signed downstream fate variable | Zero, uncommitted |
| $\mathbf{p}_i(t)$ | Apical–basal polarity vector, pointing toward the apical side | Zero |
| $V_i^\star$ | Target cell volume | Measured initial zygote volume |
| Cell ID, parent ID, cycle state | Lineage and division bookkeeping | One founder cell |

All quantities are **dimensionless**. Activities are reduced regulatory variables, not identified genes or measured concentrations. Cell-cycle timing, constitutive laws, and noise amplitudes are supplied assumptions. The model tests their consequences; it does not derive living matter, metabolism, or the cell cycle from chemistry.

The default domain is $[-1.6,1.6]^3$ on a $40^3$ grid, with reflecting mechanical boundaries. The embryo outline is not prescribed. Initially uniform signaling is perturbed by small random partition differences at cleavage; there is no imposed chemical gradient. Random division orientation is used only within geometrically degenerate long-axis subspaces, unless the isotropic control is explicitly selected.

### Deformable cells and mechanical interactions

A phase field lets each cell change shape without treating it as a rigid sphere or prescribing its surface mesh. This representation is motivated by multicellular phase-field work such as [MorphoSim (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9938209/); the energy and solver here are an independent simplified implementation.

Define a smooth occupancy function and an interface potential:

$$
h(\phi)=\phi^2(3-2\phi),
\qquad q(\phi)=\phi^2(1-\phi)^2,
\qquad V_i=\int_\Omega h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

**Occupancy is a smooth volume-counting rule.** It gives $h(0)=0$, $h(1)=1$, and $h(1/2)=1/2$: a grid location in the diffuse boundary contributes a fractional amount to cell volume. The cubic is the lowest-degree polynomial that also has zero slope at both endpoints. Consequently, volume-restoring forces proportional to $h'(\phi)=6\phi(1-\phi)$ act mainly at the boundary, rather than changing the uniform cell interior or empty space. Occupancy is not a probability of cell identity.

**The potential makes inside and outside preferable to an intermediate phase.** The nonnegative quartic $q$ has equal minima at 0 and 1 and a maximum at $1/2$. It penalizes a broad region of intermediate phase; the gradient term below opposes an infinitely sharp boundary. Their competition produces a finite-width interface. These polynomials are simple modeling choices with useful smoothness and symmetry, not unique biological laws.

![Occupancy rises smoothly from zero to one; the potential has minima at zero and one and a maximum at one half.](docs/images/occupancy-potential.svg)

The panels have different vertical scales: occupancy counts volume, whereas the potential contributes to interface energy. They are different functions with different jobs.

With regulatory states and the spatial tension coefficients frozen during each mechanical update, the baseline energy is

$$
\begin{aligned}
E={}&\sum_i\int_\Omega\gamma_i(\mathbf{x})
\left[\frac{\epsilon^2}{2}\lVert\nabla\phi_i\rVert^2+q(\phi_i)\right]\,\mathrm{d}\mathbf{x}\\
&+\sum_i\frac{K_V}{2V_i^\star}(V_i-V_i^\star)^2\\
&+\frac{R}{2}\sum_{i<j}\int_\Omega\phi_i^2\phi_j^2\,\mathrm{d}\mathbf{x}\\
&-\sum_{i<j}A_{ij}\int_\Omega q(\phi_i)q(\phi_j)\,\mathrm{d}\mathbf{x}.
\end{aligned}
$$

“Energy” here means a scalar mechanical cost used to generate shape-restoring forces. Lower values are preferred by the passive shape update. It is not ATP use, a developmental objective, or a claim that an embryo minimizes one fixed energy throughout development.

| Term, in equation order | Why this form? | What competes with it? |
|---|---|---|
| Interface cost | $\lVert\nabla\phi\rVert^2$ penalizes abrupt spatial changes; $q(\phi)$ favors inside/outside values. Together they assign a cost to cell boundary area. | Without a volume constraint, a cell can lower this cost by shrinking. |
| Volume penalty | The square penalizes both swelling and shrinkage and vanishes at the target. Dividing by $V_i^\star$ makes its derivative depend on relative volume error. | Surface and contact forces can sustain a small volume error because this constraint is soft. |
| Overlap repulsion | $\phi_i^2\phi_j^2$ is large where two cell interiors occupy the same region. The positive sign makes overlap costly. | It limits the overlap favored indirectly by cell attraction. |
| Interface attraction | $q(\phi_i)q(\phi_j)$ is appreciable only where both diffuse boundaries meet. The negative sign makes such contact favorable. | It competes with repulsion and interface cost; it does not impose a specific aggregate shape. |

The sum $i<j$ counts each cell pair once. The factors $1/2$ are coefficient conventions that simplify derivatives; they do not represent half a physical interaction. In particular, the attraction coefficient is not a contact area or a molecular binding constant.

| Symbol | Configuration name (default) | Interpretation when increased, with other parameters fixed |
|---|---|---|
| $\epsilon$ | `interface_width` (0.085) | Broadens diffuse interfaces and also changes their energy per area; it is not simply grid resolution. |
| $\gamma_0$ | `surface_tension` (1) | Raises the baseline cost of the cell boundary and the associated restoring force. |
| $K_V$ | `volume_stiffness` (12) | Resists fractional volume changes more strongly; can require smaller explicit time steps. |
| $R$ | `repulsion` (3) | Penalizes interpenetration more strongly. |
| $A_0$ | `adhesion` (4) | Strengthens the preference for overlapping diffuse interfaces. |
| $V_i^\star$ | Stored cell target, not a global configuration constant | Sets the preferred cell volume; daughters share the mother's target and do not grow between divisions. |

These tendencies do not guarantee monotonic changes in whole-embryo shape. With this energy normalization, an isolated flat equilibrium interface with constant $\gamma$ has energy per area $\gamma\epsilon\sqrt{2}/6$. Thus `surface_tension` is a model energy coefficient, not directly a measured cortical tension. The interface attraction is a contact surrogate, not a resolved cadherin or tight-junction model. [Derivation and force signs](docs/model.md#how-the-energy-generates-shape-change).

Mechanical evolution is overdamped, with mobility set to one:

$$
\frac{\partial\phi_i}{\partial t}
=-\frac{\delta E}{\delta\phi_i}+F_{\mathrm{furrow},i}.
$$

The minus sign means the passive field moves in the direction that reduces the instantaneous energy. The functional derivative $\delta E/\delta\phi_i$ asks how the energy changes when the shape field changes locally. “Overdamped” means there is no separate acceleration or momentum equation; mobility converts restoring force into the rate of shape change and is absorbed into the time scale here.

The furrow force is present only during division. Volume constraints are soft for nondividing cells; dividing mothers additionally receive the volume-preserving correction described below. Because signaling changes material properties and cytokinesis supplies active forcing, the complete simulation is not passive relaxation of a single fixed energy.

### Activator–inhibitor signaling on the changing contact graph

For interface shell $s_i=\phi_i(1-\phi_i)$, the raw contact weight is

$$
W_{ij}^{\mathrm{raw}}=\int_\Omega s_i s_j\,\mathrm{d}\mathbf{x},
\qquad W_{ii}=0.
$$

The live simulation defaults to **conservative concentration transport** (`signal_transport="conservative"`). Contacts below 2% of the largest current overlap are removed symmetrically. For complementary flat equilibrium interfaces,

$$
\frac{W_{ij}}{\mathcal A_{ij}}=\int_{-\infty}^{\infty}\phi^2(1-\phi)^2\,\mathrm{d}s
=\frac{\epsilon}{6\sqrt{2}}.
$$

We therefore estimate interface area and conductance by

$$
\widehat{\mathcal A}_{ij}=\frac{6\sqrt{2}}{\epsilon}W_{ij},\qquad
\ell_{ij}=\lVert\mathbf c_j-\mathbf c_i\rVert,\qquad
g_{ij}=\frac{\widehat{\mathcal A}_{ij}}{\ell_{ij}}.
$$

The shell weight $s_i$ is zero in the interior and exterior and largest at the interface, so $W_{ij}$ measures diffuse boundary overlap. The factor $6\sqrt{2}/\epsilon$ converts that overlap volume into an estimated area under the flat-interface assumption. Conductance increases with area and decreases with transport distance: a wider connection permits more exchange, while a longer path permits less. Here $\mathbf c_i$ is the occupancy-weighted cell center. The area estimate is a calibrated closure, not an exact face reconstruction. Gaps, overlapping or curved interfaces, and nonorthogonal center-to-face directions can bias transport. Coincident centers with positive contact are rejected.

With measured cell volumes $M=\operatorname{diag}(V_i)$ and $K=\operatorname{diag}(G\mathbf1)-G$,

$$
\Delta_V=-M^{-1}K,\qquad
(\Delta_V c)_i=\frac{1}{V_i}\sum_jg_{ij}(c_j-c_i).
$$

The difference $c_j-c_i$ sends regulator from higher to lower concentration. Dividing the net incoming amount flux by $V_i$ converts it to a concentration rate; the same incoming amount changes a small cell's concentration more than a large cell's. $M$ stores these volume capacities, and $K$ assembles the equal-and-opposite pairwise fluxes.

Exchange preserves constants and total amount $\sum_i V_i c_i$ on frozen geometry. Isolated cells have zero exchange. The Gierer–Meinhardt equations on moving compartments are

$$
\begin{aligned}
\frac{\mathrm d(V_i a_i)}{\mathrm dt}
 &=V_i\left(\frac{a_i^2}{b_i}-a_i\right)+D_a\sum_jg_{ij}(a_j-a_i),\\
\frac{\mathrm d(V_i b_i)}{\mathrm dt}
 &=V_i\beta(a_i^2-b_i)+D_b\sum_jg_{ij}(b_j-b_i).
\end{aligned}
$$

| Reaction term | Meaning and modeling choice |
|---|---|
| $a_i^2/b_i$ | Activator promotes its own production through the quadratic numerator; inhibitor suppresses that production through the denominator. The reciprocal law is a simplified inhibition rule requiring $b_i>0$. |
| $-a_i$ | Linear activator turnover; its coefficient is set to one by the chosen time scale. |
| $\beta a_i^2$ | Activator induces inhibitor, closing the negative feedback loop. |
| $-\beta b_i$ | Linear inhibitor turnover. The same $\beta$ scales inhibitor production and loss, changing its reaction time scale without shifting the positive uniform equilibrium $(1,1)$. |
| $D_a,D_b$ | Transport strengths (`signal_da`, `signal_dh`). Faster inhibitor transport can oppose broad activation while allowing localized activation; it does not guarantee a pattern on the available graph. |

The powers and coefficients specify a minimal reduced feedback system, not identified biochemical reaction steps. A quadratic production law is one way to supply nonlinear self-amplification; it is not derived from the mere existence of an activator–inhibitor loop. The model must test whether this assumed loop explains the observed outcome.

Thus concentration equations include dilution $-c_i\dot V_i/V_i$. Each step advances reaction/exchange on frozen pre-step geometry with positivity-preserving SSP-RK2 substeps, advances mechanics, then rescales concentrations by $V_i^{old}/V_i^{new}$. This is first-order splitting of the coupled moving problem; the internal RK2 solver does not make the full simulation second order. Abscission partitions measured amounts conservatively. Externally clamped signal experiments instead supply/remove regulators to maintain the imposed concentration.

Defaults are $\beta=2$, $D_a=0.02$, $D_b=0.4$ in model length-squared/time units. These are exploratory diffusivities, not calibrated molecular measurements or a conversion of the old exchange rates. The explicit `random_walk` option preserves the historical model; old checkpoints without the selector restore that option and their recorded coefficients. Polarity neighbor alignment still uses normalized orientation averaging, separately from molecular transport.

### Discrete linear stability before 3D simulation

The positive homogeneous equilibrium is $(a_*,b_*)=(1,1)$, with reaction Jacobian

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad \operatorname{tr}J=1-\beta,
\qquad \det J=\beta.
$$

Each Jacobian entry is the response of one reaction rate to a small change in one activity. The positive upper-left entry represents self-amplification, the negative upper-right entry inhibition, and the positive lower-left entry inhibitor induction. Local stability means a small perturbation in a well-mixed cell decays through their combined action, despite the activator's individual positive feedback. Local kinetics are stable for $\beta>1$. For frozen-geometry spectral analysis, the code diagonalizes

$$
S=M^{-1/2}KM^{-1/2}.
$$

Its eigenvalues are nonnegative, have inverse-length-squared units, and equal the eigenvalues of $-\Delta_V$. Physical right modes are $M^{-1/2}\mathbf q_k$, where $\mathbf q_k$ are orthonormal eigenvectors of $S$. There is no normalized upper bound of two. Each mode has a two-variable linear system

$$
M_k=J-\lambda_k\operatorname{diag}(D_a,D_b),
\qquad r_k=\max\operatorname{Re}\operatorname{eig}(M_k).
$$

A graph mode is a pattern of cell-to-cell activity differences, and $\lambda_k$ measures how strongly transport damps that mode. The zero mode is constant on each connected component. The rate $r_k$ tells whether a small perturbation grows ($r_k>0$) or decays ($r_k<0$); early linear growth is proportional to $\exp(r_k t)$. A mode with small positive growth may need much longer than one cell cycle to become visible.

A diffusion-driven instability requires stable local kinetics and at least one **supported nonzero graph eigenvalue** with $r_k>0$. The determinant is

$$
\det M_k=\beta+(\beta D_a-D_b)\lambda_k+D_aD_b\lambda_k^2.
$$

For the defaults $\beta=2$, $D_a=0.02$, and $D_b=0.4$, the unstable interval is $6.492189<\lambda_k<38.507811$. Its physical units and coefficients stay fixed during refinement; whether discrete eigenvalues enter the interval must converge and is not guaranteed on a coarse mesh. The dashboard and batch preflight use fixed unit-box meshes, followed by actual live-geometry spectra. Every abscission records spectra and volume-weighted mode transfer.

The [calibrated-contact refinement experiment](docs/live_transport.md) uses the same contact adapter on manufactured 3D slabs with unit transverse area. First-mode error decreases from 1.27% at 8 compartments to 0.020% at 64 (observed final order 1.998). The continuum has one unstable mode; 8 compartments incorrectly support two, while 16, 32, and 64 recover one. Measured small-perturbation growth agrees with the discrete prediction within $7.6\times10^{-11}$. This establishes convergence for flat complementary contacts, not arbitrary embryo geometry.

Frozen-geometry spectra describe the signaling subsystem near $(1,1)$. Volume changes introduce dilution; graph changes, mode mixing, and finite growth time require time-dependent analysis. These spectra are not a stability proof for the coupled signaling–fate–mechanics system. Historical normalized-graph experiments remain documented in [graph_signaling.md](docs/graph_signaling.md).

### From signaling to two possible cell identities

Once at least four cells are present, the default downstream fate equation is

$$
\dot f_i=r_f\left[f_i-f_i^3+g_a(a_i-1)\right].
$$

Here $r_f$ (`fate_rate`, 0.8) sets response speed and $g_a$ (`signal_fate_gain`, 1) sets the activator bias. The linear term $+f_i$ amplifies a small signed deviation from zero; the cubic term $-f_i^3$ limits that amplification. Subtracting the reference activity 1 makes homogeneous equilibrium signaling unbiased. A positive bias favors A and a negative bias favors B; it does not assign either identity directly.

Without a signaling bias, this switch has stable states at $f=\pm1$ and an unstable state at zero. Activator above or below its homogeneous value biases the switch toward opposite identities. Independent fate noise, legacy neighbor inhibition, and the direct exposure bias are disabled by default; their optional terms are documented in [the full model definition](docs/model.md).

Cells with $f_i>0.55$ are labeled A, those with $f_i<-0.55$ are labeled B, and the remainder are uncommitted. These labels are thresholds, not proof of irreversible commitment or named biological lineages. Because the fate switch is already bistable by construction, differentiated labels alone would not demonstrate a Turing mechanism; signaling growth must be tested separately.

With mechanical feedback enabled, fate changes baseline tension and interface attraction:

$$
\begin{aligned}
\gamma_i^0&=\gamma_0[1+c_\gamma\tanh(f_i)],\\
A_{ij}&=A_0[1+c_A\tanh(f_i)\tanh(f_j)].
\end{aligned}
$$

The bounded function $\tanh(f)$ prevents unbounded material coefficients when fate leaves the interval $[-1,1]$. The contrast parameters $c_\gamma$ (`fate_tension`, 0.25) and $c_A$ (`fate_adhesion`, 0.35) set the strength of this constitutive coupling. For example, the tension multiplier is bounded between $1-c_\gamma$ and $1+c_\gamma$; these limits are approached only at large fate magnitude. The product in the adhesion law is positive for similar signs and negative for opposite signs. Choosing these forms expresses a hypothesis about identity-dependent mechanics, not a consequence of the signaling equations.

Positive fate has higher baseline tension, and similarly biased cells have stronger attraction. These are explicit constitutive hypotheses, not experimentally calibrated effects of the generic identities A and B.

### Apical–basal polarity and directional mechanics

A cell's exposed cortex supplies a local orientation cue. Define the unoccupied-cortex weight, outward normal, and average cue by

$$
\begin{aligned}
e_i(\mathbf{x})&=1-\operatorname{clip}\left(2\sum_{j\ne i}h(\phi_j),0,1\right),\\
\mathbf{n}_i&=-\frac{\nabla\phi_i}{\max(\lVert\nabla\phi_i\rVert,10^{-10})},\\
\mathbf{q}_i&=\frac{\int_\Omega s_i e_i\mathbf{n}_i\,\mathrm{d}\mathbf{x}}
{\int_\Omega s_i\,\mathrm{d}\mathbf{x}}.
\end{aligned}
$$

An isolated sphere has no net cue. Contacts can create a cue toward free cortex without specifying a global embryo axis. The polarity vector evolves as

$$
\dot{\mathbf{p}}_i=
\alpha\frac{2a_i}{1+a_i}\mathbf{q}_i
+\eta(\Delta_{\mathrm{rw}}\mathbf{p})_i
-(\mu+\lVert\mathbf{p}_i\rVert^2)\mathbf{p}_i.
$$

Here $\alpha$ (`polarity_rate`, 1) controls response to free cortex, $\eta$ (`polarity_alignment`, 0.25) controls averaging toward neighboring orientations, and $\mu$ (`polarity_decay`, 0.5) removes polarity when cues are absent. The activity factor $2a/(1+a)$ equals one at $a=1$ and saturates at two, so signaling modulates rather than indefinitely amplifies the geometric cue. The cubic vector term $-\lVert\mathbf p\rVert^2\mathbf p$ increasingly opposes large polarity. These are response and saturation assumptions, not a resolved polarity-protein network. Each numerical update caps the vector magnitude at one. Daughters inherit the mother's vector and subsequently adapt to their new geometry.

Polarity modifies cortical tension spatially:

$$
\begin{aligned}
\gamma_i(\mathbf{x})&=\gamma_i^0
[1-\chi\mathbf{p}_i\cdot\widehat{\mathbf{r}}_i(\mathbf{x})],\\
\widehat{\mathbf{r}}_i&=\frac{\mathbf{x}-\mathbf{c}_i}
{\sqrt{\lVert\mathbf{x}-\mathbf{c}_i\rVert^2+\epsilon^2}},\\
F_i^{\mathrm{surface}}&=\epsilon^2\nabla\cdot(\gamma_i\nabla\phi_i)-\gamma_i q'(\phi_i).
\end{aligned}
$$

Here $\mathbf{c}_i$ is the cell centroid and $\chi$ (`polarity_tension`, 0.35) is the strength of the directional tension contrast. The dot product selects position relative to the apical direction; the $\epsilon^2$ in the denominator prevents division by zero at the centroid. Positive $\chi$ lowers effective tension on the apical side and raises it on the basal side; $\chi<1$ keeps tension positive. Conservative face fluxes include the spatial gradient of tension. Centroids and polarity are held fixed within each mechanical update.

[Nissen et al. (2018)](https://elifesciences.org/articles/38407) motivates investigating polarity-dependent morphology, but this particular cue and tension law are project-specific assumptions. There is no planar cell polarity or resolved apical protein network. Geometry can polarize cells even with uniform activator, so local polarity is not evidence of a chemically selected global axis.

### Shape-aligned division and progressive cytokinesis

Division uses the largest-eigenvalue direction of the occupancy-weighted cell covariance:

$$
\mathbf{C}_i=\frac{1}{V_i}\int_\Omega
(\mathbf{x}-\mathbf{c}_i)(\mathbf{x}-\mathbf{c}_i)^{\mathsf T}
h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

The spindle follows this long axis; the cleavage plane is perpendicular to it and initially bisects occupancy. Nearly equal maximal eigenvalues define a subspace in which orientation is sampled without privileging grid eigenvectors. This is a Hertwig-style shape rule, not a tensile-stress calculation or a simulation of spindle alignment. Experiments demonstrate a role for ring mechanics in achieving long-axis alignment in particular embryos, but their cell-rotation mechanism is not implemented here. [Middelkoop et al. (2024)](https://pubmed.ncbi.nlm.nih.gov/38870057/)

A cycle timer initiates division while retaining one mother field. With onset $t_0$ and nominal duration $T$, constriction follows

$$
u=\operatorname{clip}((t-t_0)/T,0,1),
\qquad S(u)=u^2(3-2u),
\qquad R_{\mathrm{ring}}(t)=R_0[1-S(u)].
$$

For signed plane distance $z_i$ and distance $r_\perp$ from the spindle axis, the prescribed equatorial force is

$$
F_{\mathrm{furrow},i}=-\kappa S(u)
\exp\left(-\frac{z_i^2}{2\epsilon^2}\right)
\frac{1+\tanh[(r_\perp-R_{\mathrm{ring}})/\epsilon]}{2}
h'(\phi_i).
$$

The smooth ramp starts and ends with zero slope, avoiding an abrupt switch in forcing. The Gaussian confines the force near the cleavage plane; the tanh factor selects tissue outside the shrinking ring; and $h'(\phi_i)$ localizes its action to the cell boundary. The negative sign removes occupancy from that equatorial region, while volume preservation redistributes it into the lobes. $T$ (`cytokinesis_duration`, 0.9) sets the nominal constriction time and $\kappa$ (`ring_strength`, 8) sets force amplitude. Neither alone determines actual abscission time, which also depends on the resolved neck.

This is a contracting-ring surrogate, not a resolved actomyosin network. After each mechanical update, an interface-local scalar correction preserves the mother's measured onset volume. That correction is a numerical constraint, not a hydrostatic pressure solve.

Abscission requires elapsed constriction time, a sufficiently thin resolved neck, low prospective daughter overlap, and two nontrivial lobes. Only then does an occupancy partition replace the mother:

$$
h(\phi_1)=w\,h(\phi_{\mathrm{mother}}),
\qquad h(\phi_2)=(1-w)h(\phi_{\mathrm{mother}}),
\qquad w=\frac{1+\tanh(z_i/\epsilon)}{2}.
$$

Occupancy is conserved pointwise by this construction, up to numerical precision. Daughter target volumes follow their measured lobe fractions; there is no growth between divisions. Signal partition perturbations preserve measured volume-weighted regulator amounts at cleavage; reactions can change these amounts, whereas exchange and mechanical dilution conserve them. Unresolved divisions remain active rather than being forcibly cut. See [cytokinesis details and validation](docs/cytokinesis.md).

### Numerical workflow and default scales

Each time step reconstructs contacts, advances signaling, updates polarity and fate, relaxes cell geometry, and then handles completed or newly scheduled divisions. Signaling uses SSP-RK2 with rate-dependent substeps; polarity, fate, and mechanics use explicit updates. Phase fields are clipped to $[0,1]$, and the fraction clipped is reported. A smaller time step and spatial refinement remain necessary checks, even if the simulation stays numerically finite.

| Parameter group | Default values |
|---|---|
| Domain and time | $40^3$ grid; $\Delta t=0.015$; 1,000 steps; final time 15 |
| Mechanics | $\epsilon=0.085$, $\gamma_0=1$, $K_V=12$, $R=3$, $A_0=4$ |
| Signaling | $\beta=2$, $D_a=0.02$, $D_b=0.4$; partition-noise scale 0.001 |
| Fate | $r_f=0.8$, $g_a=1$; competence at 4 cells |
| Mechanical feedback | $c_\gamma=0.25$, $c_A=0.35$ |
| Polarity | $\alpha=1$, $\eta=0.25$, $\mu=0.5$, $\chi=0.35$ |
| Division | Mean cycle interval 2; nominal constriction duration 0.9; maximum 16 cells |

These values are illustrative. In particular, signaling growth time, cell-cycle time, and mechanical relaxation time must be compared rather than assuming pattern formation finishes before the next cleavage. Complete defaults and validation rules are in [Config](embryo/model.py), with a configurable example in [default.json](examples/default.json).

Shape measurements use the capped aggregate occupancy $\rho=\min(\sum_i h(\phi_i),1)$. If its spatial covariance eigenvalues are $\ell_1\le\ell_2\le\ell_3$, the principal axis ratio is $\sqrt{\ell_3/\ell_1}$, equal to one for a sphere. Signal variation, fate counts, lineage, polarity, graph spectra, and volume errors are recorded alongside shape. Elongation alone cannot identify whether its cause is cleavage, regulatory feedback, or numerical anisotropy.

### What this model can explain, and what remains missing

The implemented system lets us test whether local activation and inhibition, a changing contact network, and polarity-dependent mechanics are sufficient for specific forms of organization. The signal equations do not guarantee two differentiated populations, and local polarity does not guarantee a persistent developmental axis. The current single-seed time-15 example reaches 16 cells with small signal differences and no threshold-classified A/B cells; it is not a validated model of mature differentiation.

**Blastocoel cavitation is not implemented.** Repulsion and cell-volume penalties can leave geometric gaps, but they do not explain accumulation of pressurized extracellular fluid. Na⁺/K⁺-ATPases actively transport ions; aquaporins conduct water passively; a low-leak epithelial barrier permits fluid accumulation. Experimental work also links pump signaling to tight-junction function. [Giannatselis et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21901128/)

A future cavity model needs solute and water balances, barrier permeability, lumen pressure coupled to cell mechanics, and communication between microlumens. Hydraulic opening of contacts and microlumen coarsening can contribute to cavity positioning, as shown in mouse embryos. [Dumortier et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/31371608/) A centered cavity can preserve rotational symmetry; selecting its position is a separate explanatory target. Signaling and polarity could regulate these processes, but cannot substitute for fluid transport equations.

The remaining validation requirements include contact-threshold sensitivity, grid/time refinement, multiple independent seeds, longer signaling times, fate-persistence assays, and comparisons with measured biology. References motivate individual mechanisms; none establishes that this combined implementation reproduces an actual embryo.

## Run

From this directory, using Python 3.10 or newer:

```bash
python -m pip install -e '.[test]'
```

### Live dashboard

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765). The dashboard advances the actual 3D simulation while displaying shaded cell surfaces, signaling variation, cell identities, shape anisotropy, and numerical diagnostics. **Show → Cell surfaces** renders the actual closed cell boundaries; **Point samples** retains the previous display. Restart an existing server and refresh the browser after updating. The [surface-rendering documentation](docs/dashboard.md#closed-cell-surfaces) explains mesh extraction, cutaway behavior, and how to inspect an existing checkpoint.

- **Run** starts the solver; **Pause** stops at a numerical step boundary; **Resume** continues the same state and random streams.
- Edit parameters before running, or pause, edit, then **Reset** to apply them to a new zygote. **Defaults** restores the form values; Reset applies them.
- Orbit, zoom, change cell colors, display polarity/contact links, and cut through the embryo. Scrub retained frames to inspect earlier states; **Follow live** returns to the newest snapshot.
- Download retained visualization frames and their parameters as JSON. This download is not a restart checkpoint. The session is held in memory and ends when the server stops.

The server binds only to this computer and needs no browser dependencies or internet connection. Use `--config examples/default.json` for an initial configuration or `--port 8766` for another port. Reference graph stability is computed before stepping; the dashboard also shows the current contact graph's growing-mode count. See [dashboard controls and limits](docs/dashboard.md).

### Batch simulation and offline playback

```bash
# Analyze finite reference graphs before a coupled run (no 3D simulation):
OPENBLAS_NUM_THREADS=1 python -m embryo.graph_analysis --output outputs/my-graph-analysis
# Every coupled run also writes a finite-graph preflight before its first step:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo --output outputs/my-run
```

NumPy, SciPy, and Matplotlib are the runtime dependencies. The thread settings avoid excessive BLAS threading overhead for the small contact matrices. Existing output directories are never overwritten.

Open `outputs/my-run/viewer.html` in a browser. It is self-contained and works offline:

- Drag to rotate; scroll to zoom.
- Play or scrub the timeline from the zygote through cleavage and differentiation.
- Move the cutaway plane to inspect internal cells.
- Color by fate, activator, inhibitor, or cell ID; show apical polarity arrows.
- Follow the principal axis ratio and volume diagnostics.

Signal colors are relative to the homogeneous activity of one, not automatically rescaled to exaggerate small variations. The viewer displays samples of the simulated cell surfaces. Its apparent surface texture is a rendering approximation; it does not enter the dynamics.

## Experiments

```bash
# Same seed and parameters, without fate/polarity-dependent mechanical properties:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-feedback --output outputs/no-feedback

# Cleavage and mechanics, with signaling, polarity, and differentiation disabled:
OPENBLAS_NUM_THREADS=1 python -m embryo --mechanics-only --output outputs/mechanics

# Separate signaling/polarity ablations:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-signaling --output outputs/no-signaling
OPENBLAS_NUM_THREADS=1 python -m embryo --no-polarity --output outputs/no-polarity

# More cells and finer geometry (substantially more expensive):
OPENBLAS_NUM_THREADS=1 python -m embryo --grid 64 --max-cells 32 --steps 1200 --output outputs/larger

# Customize coefficients with a JSON file:
OPENBLAS_NUM_THREADS=1 python -m embryo --config examples/default.json --output outputs/custom

# Small paired-seed ablation screen; this is not statistical validation:
OPENBLAS_NUM_THREADS=1 python -m embryo.experiments --seeds 7 11 23 --output outputs/screen

# Deterministic four-cell spatial/time refinement checks:
OPENBLAS_NUM_THREADS=1 python -m embryo.convergence --output outputs/convergence.json
```

`--no-feedback` retains signaling and polarity dynamics but disables their fate/polarity-dependent mechanical effects. `--mechanics-only` also disables signaling, polarity, and differentiation. Mechanical and regulatory RNG streams are separate. Shape-dependent orientation and mechanically gated abscission can change division directions and timing across feedback controls, even at a fixed seed.

### Signaling growth versus cleavage timescale

The timescale study compares longer signaling on frozen actual 4-, 8-, and 16-cell graphs, coupled runs with normal and doubled cell-cycle intervals, and replay of an identical contact history with stretched timing. Every-step graph recording allows the replay to be checked against the original coupled signals before interpreting timing effects. See [the protocol, commands, and results](docs/timescales.md).

In the first seed-7 screen, the normal coupled run reaches activator standard deviation 0.1 at time 36.14; doubling the cycle interval delays this to 44.21. Measured from first reaching 16 cells, the delays are nearly equal: 23.76 and 23.34. Thus the original weak time-15 signals primarily reflect insufficient growth time in this case. The frozen 4-cell graph suppresses fluctuations, while the 8- and 16-cell graphs amplify them on different timescales. Multi-seed robustness, mechanical/grid convergence, and fate persistence remain untested by this screen.

### Larger domain at unchanged cell resolution

The new conservative model has a [fixed-spacing domain-size protocol](docs/domain_size.md). The larger dashboard preset uses `grid=56`, `extent=2.24`, preserving the original voxel spacing 0.08 while adding 0.64 length units on each side:

```bash
python -m embryo.dashboard --config configs/large_domain.json
```

This starts a new zygote; the preset runs to $t=30$. The completed conservative-model study compares grids 40, 48, and 56 from the same initial state. Through $t=30$, peak boundary occupancy falls from $7.39\times10^{-4}$ to $1.90\times10^{-11}$; the two larger domains agree within 0.0000141% in axis ratio and 0.000948% in final common-box cell fields. All three boxes pass the boundary screen in this window. Increasing the domain does not improve cell resolution or establish clearance for longer runs. See the protocol for boundary and domain-sensitivity checks.

The $t=180$ domain continuation is complete. The original box crosses the boundary threshold at $t=131.4$; the two larger boxes pass their screens and agree within 0.000673% in axis ratio. The 56³ run retains 12 A / 4 B labels while activator contrast decays almost to zero, so persistent labels do not establish persistent signaling. Full results are in `outputs/domain-conservative-extended/RESULTS.md`. See [continuation and monitoring](docs/domain_size.md#long-time-continuation-protocol).

### Coupled spatial and temporal resolution

The [next resolution screen](docs/coupled_resolution.md) tests grids 56³, 72³, and 88³ at fixed half-width 2.24 and interface width 0.085, with a separate time-step sweep. Every grid samples the same manufactured 16-cell geometry directly. Frozen-geometry spectra agree within 0.000193% between the finer grids. The five coupled cases are complete: 11/12 checks pass, but the finest spatial occupancy-field discrepancy is 1.43%, above the 1% criterion. The completed [projection-error audit](docs/projection_audit.md) finds 1.49% discrepancy even before evolution under the original linear reconstruction; higher-order final comparisons range from 0.24% to 0.32%. This identifies substantial measurement contamination. The original failed result remains unchanged. An [independent 112³ confirmation](docs/refinement_confirmation.md) passed all 11 prospectively specified checks; finest-pair field differences are 0.162–0.207%. This supports short-time consistency for the manufactured state, not developmental convergence.

### Discrete-to-continuum transport bridge

A new fixed-domain benchmark checks conservative compartment exchange against the same continuum diffusion and Gierer–Meinhardt equations at increasing spatial resolution. It uses explicit compartment volumes and face conductances, with a sparse operator that preserves molecular amount under pure diffusion. This is a transport verification experiment; its compartments are numerical elements, not a simulation of tens of thousands of biological cells.

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.continuum --output outputs/continuum-bridge
```

Across 64 to 32,768 compartments, diffusion error falls from 0.002871 to 0.00004994 with approximately second-order convergence. Coarse grids predict ten unstable modes; finer grids recover the continuum prediction of seven. Fixed normalized exchange rates do not approximate the same bulk diffusivity under this refinement. The live simulator now uses conservative transport through a calibrated diffuse-contact adapter; its geometry closure has separate validation requirements. See [equations, protocol, results, and the remaining bridge to changing geometry](docs/continuum_bridge.md).

The next bridge stage now verifies diffusion on a **3D L-shaped domain with unequal compartment volumes**. It uses exact continuum cell averages, a uniform-mesh control, and a separate pulse traveling between the two arms:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.irregular --output outputs/irregular-bridge
```

Across 48 to 24,576 active compartments, graded-mesh error falls from 0.003924 to 0.0001186; the final observed order is 1.992. Pure-diffusion mass drift stays below $1.5\times10^{-14}$ in the refinement runs. The pulse reaches the other arm and approaches the correct volume-weighted equilibrium. All ten default acceptance checks pass. These results establish smooth-solution transport convergence on this fixed nonconvex domain; nonlinear pattern formation, moving geometry, and cell/field coupling remain separate tests. See [the irregular-domain protocol and results](docs/irregular_bridge.md).

The [irregular-domain signaling study](docs/irregular_signaling.md) now computes the complete discrete weighted spectrum and verifies early Gierer-Meinhardt growth:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.irregular_signaling \
  --output outputs/irregular-signaling-repeat
```

Refinement changes the predicted unstable-mode count from **9 to 7 to 4 to 4**. On the finest mesh, the fastest predicted growth rate is 0.215193833 and the nonlinear solver measures 0.215193831. Both growing and decaying probes pass verification, and all nine default checks pass. Agreement of the finest discrete counts is refinement evidence; an exact continuum instability count remains a separate question.

The [nonlinear persistence and convergence study](docs/nonlinear_bridge.md) now follows one physical initial perturbation through $t=100$ on the same graded L-domain:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.nonlinear \
  --output outputs/nonlinear-bridge-repeat
```

All **eleven default criteria pass**. The final two-species spatial difference decreases from 5.43% ($n=8$ versus $16$) to **1.26%** ($n=16$ versus $32$). Independent time refinement gives a **0.095% maximum sampled difference** between the two smaller steps. Late field change stays below 0.41%, while an equal-diffusivity control returns to uniformity. The benchmark uses conservative fine-to-coarse averaging and an independently tested implicit-diffusion solver; it saves concentration snapshots, histories, and comparison figures. This verifies finite-window nonlinear signaling persistence for one prescribed perturbation on fixed geometry. Moving-domain dilution and amount balance are now tested below; initial-condition robustness, fate persistence, and emergent shape remain separate tests.

The dashboard’s **Continuum signaling** workspace now plays the saved concentration fields with 3D cutaways, X/Y/Z cross-sections, activator/inhibitor colors, run selection, and convergence charts. Start with `python -m embryo.dashboard --benchmark outputs/nonlinear-bridge-verified`; an existing server needs restarting. Playback is separate from the live embryo solver. See the [dashboard guide](docs/dashboard.md#continuum-signaling-nonlinear-benchmark-playback).


The [moving-domain transport test](docs/moving_domain.md) now verifies conservation and dilution during prescribed 3D expansion and deformation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.moving \
  --output outputs/moving-domain-repeat
```

All **eleven checks pass**. Isotropic expansion, unequal axis expansion, and volume-preserving stretching conserve total amount to within $2.2\times10^{-14}$ relative drift. Uniform concentration follows inverse volume; spatial and temporal errors approach second-order convergence. The largest finest-mesh transport error is **0.01069%** relative to the exact full concentration field. Omitting dilution instead creates an artificial **82.21% amount gain** in the expanding control. This tests prescribed material motion with fixed connectivity, not emergent shape. Conservative remapping at live-cell abscission is now implemented; general remeshing remains separate.


The [persistent shape pilot](docs/shape_persistence.md) now tests the coupled embryo model through $t=90$, including matched mature-state feedback removal and a separate zygote developed without feedback. All four trajectories retain an identifiable long axis, with final shape ratios **1.323–1.331**. Full feedback does **not** exceed the controls; all final axes remain within **1.37°** of the first-cleavage axis. This is evidence of persistent elongation, not a demonstrated activator–inhibitor-specific shape axis. Every run also narrowly fails the four-grid-spacing smallest-cell-radius screen, so mechanical refinement remains necessary. The experiment exports offline 3D playbacks, axis histories, and a matched-camera comparison; [the protocol](docs/shape_persistence.md#reproduce) reproduces all controls.


The [long-time continuation](docs/shape_persistence.md#long-time-continuation) extends the three mature branches to $t=180$ without increasing the timestep. A separate [controlled chemical-patch experiment](docs/signal_patch.md) tests the response of the existing mechanics to imposed signaling, compares default and stronger fate-dependent tension against matched uniform controls, and then releases the clamp. It is explicitly a forced-response diagnostic, not evidence of spontaneous symmetry breaking.

The completed extension gives final axis ratios **1.328** (full feedback) and **1.341** (both ablations), with no positive feedback-specific excess. The full run crosses the boundary-occupancy screen at $t=124.2$, and all branches retain the smallest-cell resolution failure. Longer time has therefore exposed a domain-size limitation as well as modest shape evolution; [the continuation results](docs/shape_persistence.md#completed-continuation-results) report both.

The completed [patch-and-release pilot](docs/signal_patch.md#completed-pilot-results) also shows modest deformation: patched aggregates remain near axis ratio **1.329** at the end of forcing, at both tested fate–tension strengths. The patch/uniform difference increases with coupling, partly because the stronger uniform control becomes less elongated. All patch runs retain the resolution failure; the four feedback-enabled branches also cross the boundary screen. Larger-domain and refinement checks are therefore the next priority before interpreting later shape differences.

The mechanics kernel now avoids discarded Laplacians and repeated occupancy calculations. It reproduces the previous kernel's tested trajectories exactly; the mature 16-cell benchmark measured **1.34× throughput (25% less runtime)**. This optimization also applies to the live dashboard after restarting its server. The [benchmark and protocol](docs/signal_patch.md#equation-preserving-kernel-optimization) document the scope and reproduction commands.


## Outputs

| File | Contents |
|---|---|
| `viewer.html` | Offline, interactive 3D playback |
| `summary.png` | Shape snapshots, fate counts, and axis ratio |
| `config.json` | Full experiment parameters |
| `metrics.csv` | Volume errors, shape measures, and fate counts over time |
| `trajectory.json` | Sampled surfaces and activities for playback |
| `graph_preflight.json` | Discrete-mode stability on 4-, 8-, and 16-node reference graphs, computed before 3D stepping |
| `graph_history.json` | Actual weighted contacts, activities, polarity, and frozen-graph stability over time |
| `cleavage_spectra.json` | Before/after graphs, inherited-mode transfer, and signal partition jumps at abscission |
| `lineage.json` | Every cell's birth, parent, and division time |
| `final_state.npz` | Full final fields and RNG states for continuation |
| `diagnostics.json` | Run timing and numerical diagnostics |

The sampled trajectory is for visualization, not full-field reconstruction at earlier times. Current (schema 3) checkpoints can be continued exactly. Earlier model checkpoints are rejected rather than silently resumed with different dynamics:

```python
from embryo import Simulation
sim = Simulation.restore("outputs/my-run/final_state.npz")
for _ in range(100):
    sim.step()
sim.checkpoint("outputs/continued.npz")
```

## What is implemented

- Freely deforming cell interfaces with soft volume constraints and repulsion.
- Interface attraction, optionally dependent on fate similarity.
- Shape-aligned division, progressive equatorial constriction, and mechanically gated abscission.
- Exact final occupancy partitioning, conserved mother volume during cytokinesis, and lineage records.
- Gierer–Meinhardt concentrations with conservative contact exchange, mechanical dilution, and amount-preserving cleavage.
- Exact discrete-mode linear stability, spectral gaps, and cleavage mode-transfer diagnostics.
- A downstream signed fate switch driven by activator activity.
- Apical–basal polarity and directional cortical tension with conservative spatial fluxes.
- Fate-dependent effective surface tension.
- Numerical diagnostics, paired ablations, and reproducible checkpoints.

**The default target is 16 cells on a modest grid.** Resolve and validate this before scaling to the planned 32–64-cell model. Volume constraints are soft between divisions; inspect the measured errors. Surface/interface widths and small daughter cells need resolution studies.

A/B labels are instantaneous activity thresholds. They do **not** establish stable commitment. There is no tensile-stress-based spindle rule, lumen, growth between divisions, extracellular morphogen field, or calibrated gene network yet. Conservative exchange uses an approximate diffuse-contact area closure; continuum consistency on arbitrary deformed cell geometries remains unverified. Equatorial contraction is a prescribed ring surrogate, not a resolved actomyosin network. Elongation following cleavage is not proof of a spontaneously selected developmental axis.

## What our current model and results can explain, and what remains missing

**The model currently provides conditional explanations of particular mechanisms of organization, rather than a validated account of embryogenesis as a whole.** Its strongest result concerns signal amplification on a fixed cell-contact network. Cell identity and shape require separate evidence: a signaling pattern, two fate labels, and an elongated aggregate are not interchangeable outcomes.

| Question | What the current evidence supports | What it does not yet establish |
|---|---|---|
| Can small signal differences become persistent patterns? | On one resolved frozen embryo graph, the full activator–inhibitor loop sustains contrast in 20/20 chemical perturbation trials. Removing self-activation, transport, or differential diffusion suppresses persistent contrast. Finite-graph spectra explain which small perturbations can initially grow. | Robust patterning across independently developed embryos, arbitrary moving geometries, or biologically calibrated parameters. The 20 trials share one geometry. |
| Can initially similar cells acquire different identities? | Activator differences can bias the downstream switch toward A or B. Stable signaling ablations still yield both labels in 20/20 trials, whereas removing signal-to-fate forcing leaves the noiseless switch at zero. | That sustained activator–inhibitor patterning is necessary for differentiation. The two stable fate states are built into the switch; their biological meaning and robustness in the coupled moving embryo remain untested. Withdrawal, reversal, and noise responses of the isolated switch are now quantified separately. |
| Can local interactions produce changing, asymmetric shapes? | Deformable interfaces, contact interactions, polarity, and progressive cleavage generate changing 3D geometry without prescribing an aggregate outline. Earlier runs retain elongation and memory of the first-cleavage axis. | That the chemical loop selects a new global axis or causes additional persistent elongation. Earlier shape controls did not demonstrate that effect and had numerical limitations; the completed conservative moving controls also fail the predeclared feedback-specific elongation tests, while passing numerical-quality screens. |

The [joint signal/fate validation](docs/joint_fate.md) strengthens the first two conclusions: against a tightened independent ODE reference, the maximum continuous-fate error is $6.63\times10^{-5}$ at the tested production step, with unchanged cell-wise final labels. This resolves the tested frozen-geometry integration sensitivity, not the accuracy of the entire moving developmental trajectory.

In the [current moving pilot](docs/moving_causal.md), the full-feedback branch has completed $t=18$–78. It ends with axis ratio 1.3285 and 6 A / 10 B cells, but its cell-to-cell activator standard deviation is only 0.0135–0.0244 over the declared late window, below the persistent-contrast criterion of 0.1. All three controls are complete. Full feedback has a lower late mean axis ratio than no mechanical feedback (1.32663 versus 1.32975) and no self-activation (1.32691); it does not achieve the required positive 0.05 excess. This demonstrates why a frozen-graph instability, differentiated labels, and an elongated shape cannot alone establish chemically driven organization under movement; the completed shape comparison does not support additional feedback-specific elongation in this pilot.

Conservative transport explains how exchange and volume changes preserve regulator amount. Manufactured benchmarks establish convergence in their specified geometries. However, the live overlap-to-area approximation fails general-geometry closure screens, and the completed [zygote-to-t=90 refinement study](docs/development_refinement.md) does not pass its combined acceptance criteria: spatial refinement disagrees on signaling, while time refinement disagrees on shape and fate fractions. Numerical conservation and successful component tests therefore do not yet establish a converged developmental prediction.

What remains missing is both numerical and biological:

- **A validated causal link to shape:** refine the completed moving controls' coupled time integration, and repeat across independent developmental geometries and perturbation seeds. Test spatial organization and axis selection alongside elongation.
- **Reliable full-development numerics:** diagnose the observed refinement failures and validate transport on irregular, changing contacts. Successful regular-mesh prototypes are not yet a validated replacement for the live geometric approximation.
- **Evidence for stable cell identity:** extend the completed isolated-switch withdrawal, reversal, and noise assays to the coupled developing embryo, and distinguish fate memory from transient threshold crossings. A/B are generic model states, not identified lineages.
- **Additional developmental mechanisms:** growth between divisions, resolved intracellular/extracellular signaling, tissue fluid mechanics, and a lumen are absent. Cavitation would require an explicit account of fluid accumulation, osmotic/active transport, sealing, and pressure; an empty gap between cells would not demonstrate that mechanism.
- **Biological calibration and prediction:** relate model scales and coefficients to measurements, specify a biological system, and test predictions against independent observations and perturbations.

Order here means amplification and organization of small differences in an already structured dynamical system. Cells, interactions, feedback laws, fate bistability, and division rules are supplied assumptions. The project can test how much organization follows from those assumptions and which couplings are necessary; it does not currently derive the origin of those rules or explain the emergence of life itself.

The [signal-withdrawal experiment](docs/fate_memory.md) now quantifies this built-in fate memory. On the same frozen geometry with twenty paired perturbation seeds, even 0.6 time units of activator input lead to both labels in 20/20 trials by time 120 under the bistable law; no cells were labeled when that short input ended. A matched non-bistable relaxing law returns all cells to the uncommitted state after every tested withdrawal duration. All numerical checks pass (maximum fate timestep discrepancy 0.0009524). This supports transient signal selection followed by intrinsic switch memory, not a requirement for sustained Turing patterns or proof of irreversible biological commitment.

The subsequent [fate-reversal and noise assay](docs/fate_robustness.md) characterizes the isolated switch starting at its ideal equilibria. Opposing bias above the analytic threshold 0.3849 can reverse fate with a sufficiently long pulse; tested switching durations are 24, 12, and 3 for biases 0.4, 0.5, and 0.75. At time 60, additive-noise amplitudes 0.3 and 0.4 produce opposite labels in 42/512 and 167/512 paths. All numerical checks pass. Memory is persistent but reversible in this assumed equation; robustness of the coupled developmental system remains untested.

## Tests and next steps

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

Tests check shape-oriented spindles, gradual furrowing, volume conservation, daughter connectivity, contact symmetry, single-cell isotropy, exact continuation through active cytokinesis, cell-cap reservations, controls, and a basic time-refinement check. They are software/numerical checks, not biological validation.

The graph preflight is in `outputs/graph-analysis`; the coupled example is `outputs/graph-polarity/viewer.html`, with a matched no-polarity control in `outputs/graph-no-polarity`. Earlier outputs are historical models. In the initial time-15 coupled example, signals remain small and all 16 cells remain uncommitted: a positive frozen-graph growth rate does not guarantee a mature pattern on the available developmental timescale.

Read [the cytokinesis update and validation](docs/cytokinesis.md), [the equations and assumptions](docs/model.md), [the staged implementation plan](docs/plan.md), and [the initial numerical results](docs/results.md).

## Manuscript draft

The manuscript describes the earlier normalized live model and has not yet been revised for the conservative coupling. The [Introduction and Methods draft](manuscript/introduction_methods.md) describes the coupled model, signaling-timescale experiments, and conservative transport benchmark with primary references. The [manuscript folder](manuscript/README.md) contains editable source, a self-contained HTML reading copy, and LaTeX/BibTeX files. This is a partial manuscript for development toward an arXiv submission, not a submitted paper.

The subsequent [frozen contact-cutoff screen](docs/contact_sensitivity.md) preserves connectivity and unstable-mode counts across all seven tested geometries, but fails its quantitative spectral-growth sensitivity check at four developmental times. The [paired evolving-geometry experiment](docs/cutoff_dynamics.md) completed t=30 to 90 with all five sensitivity checks passing while polarity filtering was fixed at 0.02. Identity labels match throughout; maximum activator RMS difference is 0.00203 and maximum relative axis-ratio difference is 1.19e-6. Default trajectories retain the historical shared-cutoff behavior; the known cell-resolution shortfall remains explicit.

The next [controlled small-cell cleavage screen](docs/cleavage_resolution.md) compares three grids, three time steps, and axial/oblique division directions. It tests abscission conservation, daughter connectivity, event timing, and shape before full developmental refinement.

The [full developmental refinement study](docs/development_refinement.md) has completed all five fresh zygote-to-t=90 cases: grids 56³/72³/88³ at fixed time step, plus an independent time-step sweep on 72³. The combined acceptance result is **FAIL**: spatial signaling comparisons and temporal shape/fate comparisons exceed their declared limits. Passing component checks do not establish developmental convergence. The coarse 56³ reference also fails the smallest-cell resolution screen.

The cleavage measurement issue is [resolved by calibrated native-phase reconstruction](docs/cleavage_measurement.md): reconstruct phi before applying occupancy, giving at most 0.1503% cubic error across analytic holdouts (limit 0.25%). All 10 revised checks pass; axial and oblique finest-pair evolved differences stay below 0.56%. Original results and ongoing developmental simulations remain unchanged.

The [geometric transport validation](docs/geometric_transport.md) exposes limits beyond numerical conservation: curved-interface area bias reaches 5.80% at width/radius 0.3, a one-width gap retains 82% of planar overlap, and nonorthogonal contacts fail a linear-field flux test. General-geometry closure checks fail despite accurate quadrature and conservation. The current operator is therefore not validated as continuum bulk diffusion on arbitrary embryo geometries; the ongoing refinement runs remain tests of its existing approximation.

The [face-normal-aware skew-mesh prototype](docs/skew_flux.md) now passes linear-flux, conservation, stability, and approximately second-order refinement checks. At shear 1, finest error falls from 2.31% with A/ell to 0.0140%. Its positive-pulse test fails (minimum −0.00443 under exact discrete evolution), so it remains separate from live signaling. A positivity-preserving correction is the next transport task.

A [positive wider-stencil skew-mesh prototype](docs/positive_skew_flux.md) now passes all 11 checks: nonnegative rates and pulse evolution, exact affine face fluxes, mass conservation, and approximately second-order space/time convergence. Worst finest-grid error is 0.0030%. It remains separate from the live model; skewed no-flux boundaries and irregular geometry are next.

The [reflecting-wall skew-mesh test](docs/skew_boundary.md) passes 8/9 checks: positivity, conservation, and finest errors below 0.1% hold, but intermediate-shear refinement orders of 1.58 and 1.52 miss the 1.8 criterion. Boundary accuracy must be improved or independently resolved before irregular-geometry/live integration.

A [derived positive wall-conductance correction](docs/skew_boundary_correction.md) resolves the tested boundary-order failure: all nine checks pass for two independent wall families, with 64³ confirmation orders near 1.997 and errors below 0.0052%. Conservation and positivity remain intact. Nonuniform capacities and graded geometry are next; live transport is unchanged.

The first [causal activator–inhibitor screen](docs/causal_signaling.md) is complete on one frozen resolved embryo graph with 20 paired chemical perturbation seeds. Full feedback sustains signal contrast in 20/20 trials; removing self-activation, transport, or differential diffusion suppresses it. Those controls still produce both fate labels, exposing the separate bistable fate switch as an alternative mechanism. Removing signal-to-fate coupling prevents commitment. The original finer-time checks reproduced labels but failed continuous-fate tolerances. A subsequent [joint signal/fate integrator](docs/joint_fate.md) resolves that sensitivity against a tightened independent ODE reference: maximum fate error is 6.63e-5 at dt=0.0075, below the unchanged 0.05 limit. The corrected trajectories retain the same scientific distinction between persistent patterns and fate labels. Four matched [moving-geometry feedback controls](docs/moving_causal.md) have completed t=18 to t=78. All quality checks pass, but full feedback does not exceed the no-feedback or no-self-activation shape controls, and its persistent-signal check fails. This is not yet a developmental ensemble.


The [geometry replay diagnostic](docs/geometry_replay.md) now tests frozen transport, replayed transport/volumes, omitted dilution, and frozen transport with volume forcing. The six-unit snapshot pilot passes chemical step refinement but fails fidelity to the live trajectory and snapshot-spacing checks; causal attribution is withheld. An exact continuation is capturing conductances and volumes at every mechanical step, with bitwise checkpoint verification, before automatically repeating the replay controls.
