# Embryogenesis in 3D

## Scientific objective

**How much developmental organization can a small set of mathematical rules explain?** This project studies whether activator–inhibitor chemistry, communication between cells, and polarity-dependent mechanics can turn small fluctuations into persistent differences during development from a single zygote.

Three outcomes are evaluated separately: **cell identity differentiation**, **signaling symmetry breaking**, and **geometry/shape symmetry breaking**. A chemical pattern does not by itself establish a cell identity or explain a new tissue axis. The model is generic and dimensionless; it does not reconstruct a particular species or derive living cells, metabolism, or division rules from first principles.

The central scientific finding so far is **conditional collective chemical organization**: the opportunity to form a pattern differs from the ability to maintain one, and persistent cell-associated chemistry and responses remain sensitive to the surrounding network. Evidence mainly comes from mature tissue inherited from three developmental histories. Organization throughout fresh development remains an open test.

[Current claims](#current-claims) · [Model](#method-and-mathematical-model) · [Open failures](#open-failures-and-missing-evidence) · [Lab notebook](docs/lab_notebook.md)

## Definitions and interpretation

| Term | Meaning in this study |
| --- | --- |
| Cell | A deformable 3D region with well-mixed intracellular activator and inhibitor. Integer cell/parent IDs track individuals and ancestry, not types. |
| Chemical state | Positive activities $(a_i,b_i)$ in cell $i$. Amounts are $V_i a_i$ and $V_i b_i$; neither chemical is an identity label or an identified gene. |
| Cell attributes | Measured chemical activities, polarity magnitude, cell shape and response to specified perturbations. Activity-derived tension/adhesion are constitutive outputs, not independent evidence of differentiation. |
| Cell identity | A distinguishable, relatively persistent profile of attributes and reproducible responses under declared conditions. Chemical contrast, coloring or clustering alone does not establish identity. No number of identities is prescribed. |
| Differentiation | Lasting differences arising between initially similar cells. Chemical differentiation is a narrower claim than biological identity differentiation, function or commitment. |
| Apical–basal polarity | Local directional organization, represented by $\mathbf p_i$. Its direction selects the modeled apical side and its magnitude, capped at one, measures orientation strength. Zero magnitude has no direction. It is not a tissue axis, molecular apical domain or planar polarity. |
| Context and memory | Context includes neighboring chemistry, contacts, position, volume and exposed cortex. A lasting effect of prior conditions may be maintained collectively; persistence in one context does not establish autonomy. |
| Commitment and inheritance | Retention after a defined change in instructive conditions, and transmission through later divisions, respectively. Neither follows from a parent ID, copied variable or mature-state persistence. |
| Symmetry breaking | Emergence of differences from an initially equivalent configuration. Chemical heterogeneity, local polarity and tissue anisotropy are separate measurements; an imposed cue or cleavage axis must be identified as an input. |

## Method and mathematical model

Each cell has a diffuse shape field $\phi_i(\mathbf x,t)$, chemical activities $(a_i,b_i)$, polarity $\mathbf p_i$, and a target volume $V_i^\star$. Shape fields share a 3D computational box; the tissue boundary changes as cells deform, move and divide.

```mermaid
flowchart LR
    G[Cell geometry and volumes] -->|Contacts and exposed cortex| C[Chemical transport and polarity]
    C -->|Material properties| M[Mechanics and progressive cleavage]
    M -->|Deformation and division| G
    G -->|Dilution and chemical partition| C
```

The research model has **no integrated fate switch or A/B classification**. Initial conditions and perturbations are declared per experiment: a spherical zygote with spatially uniform chemistry, small amount-preserving daughter fluctuations, or explicitly prepared mature states. These preparations answer different questions.

### Deformable cells and mechanical interactions

A phase field is approximately one inside a cell and zero outside. Define

$$
h(\phi)=\phi^2(3-2\phi),\qquad
q(\phi)=\phi^2(1-\phi)^2,\qquad
V_i=\int_\Omega h(\phi_i)\,\mathrm d\mathbf x.
$$

**Occupancy $h$ counts volume smoothly**, with fractional contributions at the diffuse boundary. It has $h(0)=0$, $h(1)=1$ and zero slope at both endpoints, so its volume-restoring force mainly acts at interfaces. **The potential $q$ favors inside/outside values** through equal minima at zero and one. Together with the gradient penalty, it produces a finite-width boundary. These polynomials are modeling choices, not identity probabilities or unique biological laws.

With chemical and polarity-dependent coefficients held fixed during a mechanical update, the baseline energy is

$$
\begin{aligned}
E={}&\sum_i\int_\Omega\gamma_i(\mathbf x)
\left[\frac{\epsilon^2}{2}\lVert\nabla\phi_i\rVert^2+q(\phi_i)\right]\,\mathrm d\mathbf x\\
&+\sum_i\frac{K_V}{2V_i^\star}(V_i-V_i^\star)^2\\
&+\frac R2\sum_{i<j}\int_\Omega\phi_i^2\phi_j^2\,\mathrm d\mathbf x
-\sum_{i<j}A_{ij}\int_\Omega q(\phi_i)q(\phi_j)\,\mathrm d\mathbf x.
\end{aligned}
$$

| Term | Purpose and reason for its form |
| --- | --- |
| Interface cost | The gradient term resists abrupt changes; $q$ resists intermediate phase. Their competition creates a smooth surface with a mechanical cost. |
| Volume penalty | A square penalizes both swelling and shrinkage. Division by $V_i^\star$ makes the restoring force depend on fractional volume error. This is a soft constraint. |
| Overlap repulsion | Positive interior overlap costs energy, resisting cell interpenetration. |
| Interface attraction | The negative boundary-overlap term favors contact, competing with repulsion and surface cost. It does not specify an aggregate shape. |

Here $\epsilon$ controls interface width, $\gamma_i$ the surface-cost coefficient, $K_V$ volume stiffness, $R$ repulsion, and $A_{ij}$ attraction. The sum $i<j$ counts pairs once. “Energy” is an instantaneous mechanical cost, not ATP consumption or a fixed objective minimized throughout development.

Mechanical evolution is overdamped:

$$
\frac{\partial\phi_i}{\partial t}=-\frac{\delta E}{\delta\phi_i}+F_{\mathrm{furrow},i}.
$$

The functional derivative turns local energy changes into shape-restoring forces; mobility is absorbed into the time unit. Active furrow forcing and evolving material properties make the complete model more than passive relaxation. [Force derivation and parameter meanings](docs/model.md#baseline-energy-four-competing-costs).

### Activator–inhibitor signaling on the changing contact graph

Cell contacts provide exchange conductances. With interface shell $s_i=\phi_i(1-\phi_i)$,

$$
W_{ij}=\int_\Omega s_i s_j\,\mathrm d\mathbf x,\qquad
\widehat{\mathcal A}_{ij}=\frac{6\sqrt2}{\epsilon}W_{ij},\qquad
g_{ij}=\frac{\widehat{\mathcal A}_{ij}}{\lVert\mathbf c_j-\mathbf c_i\rVert}.
$$

$W_{ij}$ measures diffuse boundary overlap; its conversion to area is calibrated for complementary flat interfaces. Larger area increases exchange and greater center distance reduces it. Small contacts are removed symmetrically using the protocol's cutoff. This is an approximate geometric closure whose general accuracy remains unresolved.

Let $G=(g_{ij})$, $M=\operatorname{diag}(V_i)$ and $K=\operatorname{diag}(G\mathbf1)-G$. The live chemical operator is

$$
\Delta_V=-M^{-1}K,\qquad
(\Delta_V c)_i=\frac1{V_i}\sum_jg_{ij}(c_j-c_i).
$$

This exchanges equal and opposite amounts, preserves uniform concentration, and conserves volume-weighted amount on frozen geometry. It replaces the historical random-walk chemical operator; conservation alone does not prove continuum accuracy.

Gierer–Meinhardt reactions and moving-volume dilution give

$$
\begin{aligned}
\dot a_i&=\frac{a_i^2}{b_i}-a_i+D_a(\Delta_V a)_i-a_i\frac{\dot V_i}{V_i},\\
\dot b_i&=\beta(a_i^2-b_i)+D_b(\Delta_V b)_i-b_i\frac{\dot V_i}{V_i}.
\end{aligned}
$$

$a_i^2/b_i$ supplies self-amplification opposed by inhibition; $-a_i$ supplies turnover. Activation also produces inhibitor through $\beta a_i^2$, closing the negative feedback loop. $\beta$ sets the inhibitor reaction timescale; $D_a,D_b$ set exchange strengths. Faster inhibitor transport can enable spatial patterning, but the available graph modes determine whether it does. Dilution changes concentration when volume changes without creating chemical amount. Reactions still produce and remove chemicals: the full chemical system is not amount-conserving.

### Discrete linear stability before 3D simulation

Uniform chemistry has equilibrium $(a,b)=(1,1)$ and reaction Jacobian

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},\qquad
r_k=\max\operatorname{Re}\operatorname{eig}
\left[J-\lambda_k\operatorname{diag}(D_a,D_b)\right].
$$

$\lambda_k$ are eigenvalues of the symmetric volume-weighted operator $M^{-1/2}KM^{-1/2}$. For $\beta>1$, local uniform kinetics are stable; a diffusion-driven instability additionally requires a supported nonzero mode with $r_k>0$. A continuous instability band may contain no eigenvalue on a small cell graph. Cleavage changes capacities, modes and the time available for growth.

These frozen snapshots analyze chemical initiation near uniformity. They do not establish maintenance of a developed nonlinear pattern or stability of chemistry, polarity and mechanics together. [Transport and stability analysis](docs/live_transport.md).

### Apical–basal polarity and directional mechanics

Activity modulates baseline tension and attraction through a bounded response:

$$
r_i=\tanh(a_i-1),\qquad
\gamma_i^0=\gamma_0(1+c_\gamma r_i),\qquad
A_{ij}=A_0(1+c_A r_i r_j).
$$

The reference value one makes uniform equilibrium chemistry mechanically unbiased; tanh bounds the coefficients. This law has no stored fate, threshold or independent memory.

Let $\mathbf u_i$ be the surface-averaged outward cue from exposed cortex. Polarity evolves as

$$
\dot{\mathbf p}_i=
\alpha\frac{2a_i}{1+a_i}\mathbf u_i
+\eta(\Delta_{\mathrm{orient}}\mathbf p)_i
-(\mu+\lVert\mathbf p_i\rVert^2)\mathbf p_i.
$$

$\alpha$ sets cue response, $\eta$ neighbor orientation averaging, and $\mu$ decay. Activity modulates the cue through a bounded factor; the cubic term limits large polarity. Orientation averaging uses a normalized graph separately from conservative molecular transport. An isolated sphere has no net cue; contacts can establish local polarity even with uniform chemistry.

Directional tension is

$$
\gamma_i(\mathbf x)=\gamma_i^0
\left[1-\chi\mathbf p_i\cdot
\frac{\mathbf x-\mathbf c_i}{\sqrt{\lVert\mathbf x-\mathbf c_i\rVert^2+\epsilon^2}}\right].
$$

Positive $\chi$ lowers the surface-cost coefficient toward the apical side and raises it toward the basal side. It is a constitutive hypothesis, not a resolved polarity-protein mechanism. The mechanical force includes spatial variation in tension. Controls separate activity-dependent tension/adhesion, polarity tension, their combination and baseline mechanics. [Polarity definitions and assumptions](docs/graph_signaling.md#apicalbasal-polarity-and-mechanics).

### Progressive cleavage and numerical method

A prescribed cell-cycle timer initiates division. The spindle follows the longest occupancy-weighted cell axis, with unbiased sampling in nearly degenerate axis subspaces. A smooth contracting-ring surrogate forms an equatorial furrow while retaining the mother field. Abscission waits for elapsed constriction time, a resolved thin neck, small daughter overlap and two connected lobes. Daughter occupancy, target volume and chemical amounts are partitioned conservatively. These are supplied division rules, not a model of spindle machinery or actomyosin biochemistry. [Cleavage model](docs/cytokinesis.md).

Chemical contrast is the across-cell standard deviation of $\log a_i$; its screening threshold is protocol-specific, not an identity definition. Tissue shape is measured from capped aggregate occupancy: the square root of the largest/smallest spatial covariance eigenvalue is one for a sphere and increases with anisotropy. The live studies described here stop at sixteen chemical compartments; a finer voxel grid resolves their shapes, not a many-cell chemical continuum.

The coupled update advances reactions/exchange, polarity, mechanics, volume-induced concentration conversion and division bookkeeping. It is **first-order operator splitting**; second-order frozen chemical substeps do not make the full moving solver second order. Qualified paths retain small phase increments and chemical rounding corrections. PyTorch keeps GPU arrays/matrices resident; custom CUDA computes mechanics, geometry, polarity and furrows. CPU references and independent frozen ODE solvers support numerical checks. [Backend scope](docs/backends.md).

Representative tested coefficients are $\epsilon=0.085$, $\gamma_0=1$, $K_V=12$, $R=3$, $A_0=4$, $\beta=2$, $D_a=0.02$, $D_b=0.55$, $c_\gamma=0.25$, $c_A=0.35$ and $\chi=0.35$. They are dimensionless, exploratory values; controls, timesteps, grids and observation windows belong to individual frozen protocols, not universal accuracy settings.

**Resource availability is assumed permissive.** Oxygen, nutrient supply, ATP and biomass synthesis are not simulated. Cleavage partitions existing target volume rather than producing tissue growth. Cell-cycle timing, reactions and active mechanics do not respond to resource limitation.

## Current claims

The principal mature-tissue evidence uses **three histories (7, 8, 9)**, reused across assays. Cells, branches, preparations, pulses and numerical repeats are nested observations, not additional developmental replicas or population success probabilities. Claims below refer to qualified experiments with their stated initial conditions, coefficients and horizons.

| Finding | Supported claim and scope | Evidence |
| --- | --- | --- |
| Initiation differs from maintenance | Near-uniform chemistry can fail to organize under moving mechanics even when developed chemistry survives. Polarity-specific initiation suppression is demonstrated in **one of three histories**; the other two fail even without directional polarity tension. | [Initiation](docs/phase_carry_polarity_assessment.md), [maintenance](docs/phase_carry_maintenance_assessment.md) |
| Changing communication contributes to failed initiation | Preserving initial chemical conductances restores sustained patterns in all three tested mature histories, including the polarity-associated loss. This supports a conductance-mediated contribution, not network homogenization or complete causal mediation. | [Conductance controls](docs/polarity_conductance_controls_assessment.md) |
| Developed states can persist | Developed chemical contrast survives 240 elapsed units at both tested polarity-tension settings in all three histories. Actual frozen endpoint graphs support local uniform/patterned coexistence. Neither observation proves whole-system stability or cell autonomy. | [Maintenance and endpoint coexistence](docs/phase_carry_maintenance_assessment.md) |
| Chemical exchange transfers response tendencies | Selected exchanged cells respond closer to the chemical donor reference than to their original destination reference in all three histories, while concentrations reorganize. This is response similarity, not exact state transfer or donor equivalence. | [Moving exchange/response](docs/phase_carry_exchange_response_assessment.md) |
| Responses depend on collective context | Surrounding chemical allocation changes selected responses even with recipient chemistry, initial species amounts and physical fields matched. Moving tests retain immediate effects and selective delayed effects; these are finite-window, context-dependent phenotypes. | [Immediate amount-matched response](docs/moving_pulse_time_amount_response_assessment.md), [delayed moving response](docs/moving_recipient_amount_response_assessment.md) |
| Short repeated cleavage is numerically reproducible at a tested resolution | Four matched mechanical branches pass scoped finest spatial/time screens and live-cleavage restart checks at grid 120³ through t=18 on one analytic preparation. This is numerical qualification, not independent biological replication. | [Developmental qualification](docs/development_control_qualification_assessment.md) |

These findings explain conditional organization within the supplied regulatory and mechanical rules. They do **not** establish autonomous or inherited biological cell types, a prescribed number of identities, a chemically selected tissue axis, or differentiation in a living embryo. [Evidence ledger and original-method precision audit](docs/results_ledger.md).

## Open failures and missing evidence

1. **General geometric transport closure fails its accuracy screens.** The overlap-to-area/center-distance approximation passes flat orthogonal references but fails tested curved, separated and nonorthogonal contacts. Conservative bookkeeping and voxel refinement do not remove this modeling error. General tissue-continuum interpretation remains unsupported. [Closure tests](docs/geometric_transport.md).
2. **Developmental numerical acceptance remains scoped.** Historical coarse-grid and timestep comparisons failed; finer acceptance covers one analytic preparation and a short cleavage window. A fresh-history feedback-off control originally failed the 1% dense transport-spectrum timestep gate near a contact-cutoff transition. Another timestep halving passes the unchanged limits for that control, but its three matched mechanical branches still require the same finer comparison. Acceptance of the full fresh-history cohort remains unresolved, and longer observation windows are needed. A finer shape grid does not add chemical compartments. [Retained fresh-zygote failure](docs/fresh_zygote_pilot_assessment.md), [scoped fresh-history refinement pass](docs/fresh_zygote_time_refinement_assessment.md), [retained refinement failures](docs/repeated_development_gpu_assessment.md), [earlier developmental qualification](docs/development_control_qualification_assessment.md).
3. **Older no-carry measurements are not retroactively accepted.** Original-method timing, response amplitudes and spectral trajectories require their own correction-aware rederivation. Later passes apply only to their tested windows and preparations. [Claim/dependency audit](docs/results_ledger.md#phase-update-precision-audit).
4. **The initiation mechanism is not fully isolated.** Conductance interventions identify a contribution; necessity across fresh histories and parameters, chemical modulation of polarity versus geometric cues, and complete coupled stability remain unestablished. [Mechanistic limits](docs/polarity_conductance_controls_assessment.md).
5. **Identity requires further behavioral and developmental evidence.** Persistence and context sensitivity do not establish autonomous memory, commitment, inheritance through later divisions or distinct biological function. Larger independent cohorts and biological calibration remain necessary.
6. **Activator–inhibitor-specific shape organization is unproven.** Cleavage and packing produce anisotropy, including in controls. An additional persistent chemically caused axis, folding or large structural organization has not been demonstrated. [Shape controls](docs/shape_persistence.md).
7. **Cavitation and metabolism are absent.** Geometric gaps do not model a blastocoel. Cavitation would require solute/water transport, barrier permeability and pressure coupled to mechanics; resource limitation would require additional metabolic balances.

## Reading and reproduction

| Purpose | Entry point |
| --- | --- |
| Equations, derivations, parameter meanings and literature motivations | [Detailed model](docs/model.md) |
| Run status, dated observations and report catalog | [Lab notebook](docs/lab_notebook.md) |
| Audited claims, histories, source hashes and retained failures | [Results ledger](docs/results_ledger.md) |
| Priorities and experimental decision rules | [Research plan](docs/plan.md) |
| Installation, demonstrations, checkpoints and scientific reproduction | [Running the project](docs/running.md) |
| Research draft and build/evidence notes | [Manuscript](manuscript/introduction_methods.md), [manuscript guide](manuscript/README.md) |

The README changes when the model, supported claims or open failures change. Launches, progress checks, runtimes, storage operations and interim results belong in the lab notebook. Protocols and numerical acceptance determine scientific claims; completion status and visual resemblance do not.
