# Model definition

The explicit graph signaling equations, discrete stability analysis, and apical–basal polarity mechanics are specified in [graph_signaling.md](graph_signaling.md). The scalar activity below is the downstream fate variable, distinct from activator and inhibitor.

All coordinates, times, activities, energies, and coefficients are dimensionless. The model is generic. Parameters are illustrative rather than measured.

## Geometry

Cell $i$ has a diffuse indicator $\phi_i(\mathbf{x}) \in [0,1]$. Its occupancy and volume are

$$
h(\phi) = \phi^2(3-2\phi),
\qquad
V_i = \int_{\Omega} h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

Here $\Omega$ is the computational domain. The initial zygote is a radius-0.8 sphere with a smooth interface. Cell shapes are free to change. The finite computational box uses zero-normal-gradient boundary conditions and must be large enough to avoid affecting the embryo.

Define $q(\phi)=\phi^2(1-\phi)^2$. At fixed fate and target volumes, the energy is:

$$
\begin{aligned}
E ={}&
\sum_i \int_{\Omega} \gamma_i
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

Mechanics follows overdamped gradient descent with mobility set to one:

$$
\frac{\partial\phi_i}{\partial t}
= -\frac{\delta E}{\delta\phi_i}.
$$

The code uses a finite-difference Laplacian and explicit time integration. Values are clipped to $[0,1]$; `clipped_fraction` exposes this numerical safeguard. Check step-size sensitivity, especially after division, rather than assuming clipping makes the method accurate.

The attraction is a phenomenological attraction between diffuse interface regions, not a calibrated cadherin or junction model. The volume penalty is soft outside cytokinesis; a dividing mother additionally uses the volume constraint below. Changing fate, cleavage, and stochastic dynamics mean the full simulation is not a passive energy-minimization process.

## Regulatory activity

Each cell has a signed scalar $f$, a reduced normal form for competition between two possible identities. It is not a concentration, and need not lie in $[-1,1]$. The isolated deterministic switch

$$
\frac{\mathrm{d}f}{\mathrm{d}t}=r(f-f^3)
$$

has stable states at $-1$ and $+1$, with an unstable state at $0$.

After the configured competence cell count, the discrete update over a time step $\Delta t$ is:

$$
\begin{aligned}
\Delta f_i ={}&
r\left[
f_i-f_i^3
-k\,\overline{f}_{i,\mathrm{neighbor}}
+b(e_i-0.5)+g_a(a_i-1)
\right]\Delta t
\\
&+\sigma\sqrt{\Delta t}\,\xi_i,
\qquad \xi_i\sim\mathcal{N}(0,1).
\end{aligned}
$$

The activator activity is $a_i$, with homogeneous equilibrium $a_i=1$; $g_a$ is `signal_fate_gain`. The legacy neighbor inhibition, exposure bias, independent fate noise, and fate partition noise are zero by default, so graph signals provide the default bias.

Here $\overline{f}_{i,\mathrm{neighbor}}$ is the contact-weighted neighbor activity, $e_i$ is the exposure proxy, and $\xi_i$ is an independent standard normal draw for each cell and time step.

Neighbor weights are integrals of overlapping $\phi(1-\phi)$ interface shells. The exposure proxy is the shell-weighted mean of

$$
1-\operatorname{clip}\left(
2\sum_{j\ne i}h(\phi_j),\,0,\,1
\right).
$$

It is not an exact free-surface area fraction. The fixed exposure reference $0.5$ is an explicit model assumption, not calibrated biology.

The neighbor term is lateral inhibition: it can favor different activities among contacting cells. The exposure term biases externally exposed cells toward positive activity. Neither guarantees both identities in every run. The homogeneous state does not have a verified pattern-forming dispersion relation on a growing contact graph.

Positive activity is called A and negative activity B. Counts use $f>0.55$ or $f<-0.55$; remaining cells are uncommitted. Stable commitment requires future persistence and perturbation assays.

## Feedback

When enabled:

$$
\begin{aligned}
\gamma_i
&= \gamma_0\left[1+c_{\mathrm{tension}}\tanh(f_i)\right],
\\
A_{ij}
&= A_0\left[1+c_{\mathrm{adhesion}}\tanh(f_i)\tanh(f_j)\right].
\end{aligned}
$$

The feedback coefficients $c_{\mathrm{tension}}$ and $c_{\mathrm{adhesion}}$ correspond to `fate_tension` and `fate_adhesion` in the configuration.

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
R(t) &= R_0[1-S(u)].
\end{aligned}
$$

Here $R_0$ is the largest radial extent of the mother's $\phi_i\ge 0.5$ interior at onset. A contracting annular potential, confined near the equator, drives furrow ingression:

$$
\begin{aligned}
B(s) &= \exp\left(-\frac{s^2}{2\epsilon^2}\right),\\
H(r_\perp,t) &= \frac12\left[1+\tanh\left(\frac{r_\perp-R(t)}{\epsilon}\right)\right],\\
E_{\mathrm{furrow},i}(t) &= \kappa S(u)\int_{\Omega}
B(s)H(r_\perp,t)h(\phi_i)\,\mathrm{d}\mathbf{x},\\
F_{\mathrm{furrow},i} &= -\kappa S(u)B(s)H(r_\perp,t)h'(\phi_i).
\end{aligned}
$$

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

The scalar $\alpha_i$ is solved numerically, with relative volume tolerance $10^{-7}$ before float32 storage. `volume_projection_max` reports the largest correction coefficient per step. This is a numerical incompressibility constraint during cytokinesis, not an inferred hydrostatic pressure. Nondividing cells retain the original soft volume penalty. The measured onset volume can already differ slightly from the target, so the constraint preserves that existing error rather than causing a corrective jump at onset.

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

Independent mechanical and regulatory random streams remain separate. Coupled controls can nevertheless differ in cleavage orientation and completion timing because shape and mechanical relaxation now determine these events. Full active-event state is checkpointed. Checkpoint schema 3 includes signaling, polarity, and graph-event state and resumes exactly; earlier checkpoints are rejected explicitly because their subsequent dynamics would differ.

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

A sphere has ratio $1$ and asphericity $0$. Neither measure identifies the cause of asymmetry.

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

