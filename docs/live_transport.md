# Conservative signaling in the live embryo

New simulations default to `signal_transport="conservative"`. This couples `transport.conservative_transport` and its positivity-preserving `integrate_gm` solver to measured live-cell geometry. `random_walk` remains an explicit historical control. Old schema-3 checkpoints without the selector restore random-walk transport with their original parameters; no historical trajectory is silently reinterpreted.

## Geometry and equations

For each pair, the code measures the existing shell overlap

$$
W_{ij}=\int \phi_i(1-\phi_i)\phi_j(1-\phi_j)\,\mathrm dV.
$$

For complementary flat interfaces at the equilibrium profile of this model, $|\phi'|=\sqrt{2}\phi(1-\phi)/\epsilon$. Substitution gives

$$
W_{ij}/A_{ij}=\frac{\epsilon}{\sqrt{2}}\int_0^1\phi(1-\phi)\,\mathrm d\phi
=\frac{\epsilon}{6\sqrt{2}}.
$$

After symmetric relative contact filtering (default 0.02, with no absolute floor for conservative transport), define

$$
\widehat A_{ij}=6\sqrt{2}W_{ij}/\epsilon,\quad
\ell_{ij}=|\mathbf c_i-\mathbf c_j|,\quad
g_{ij}=\widehat A_{ij}/\ell_{ij},\quad
V_i=\int h(\phi_i)\,\mathrm dV.
$$

Centers use the same occupancy weighting as volumes. Positive contact between coincident centers is rejected rather than assigned an arbitrary regularized flux. The sparse operator is

$$
K=\operatorname{diag}(G\mathbf 1)-G,\qquad M=\operatorname{diag}(V_i),\qquad
\Delta_V=-M^{-1}K.
$$

It obeys $\Delta_V\mathbf1=0$ and $\mathbf1^TM\Delta_V=0$. Isolated cells have zero exchange. Equal-and-opposite pairwise amount fluxes conserve total amount under diffusion even for unequal volumes. Reactions remain active and can create/remove regulators.

With $R_a=a^2/b-a$ and $R_b=\beta(a^2-b)$, each species $c$ evolves according to

$$
\frac{\mathrm d(V_i c_i)}{\mathrm dt}=V_iR_c+D_c\sum_jg_{ij}(c_j-c_i),\qquad
\dot c_i=R_c+D_c(\Delta_Vc)_i-c_i\frac{\dot V_i}{V_i}.
$$

The step freezes geometry, advances reaction/exchange, advances mechanics, then applies $c_i\leftarrow c_iV_i^{old}/V_i^{new}$. At cleavage, the parent amount is divided between measured daughter volumes, with volume-balanced partition noise. A small volume discrepancy at relabeling cannot create regulator amount. Signals disabled as a control remain frozen; prescribed concentration clamps are external sources/sinks and bypass dilution while active.

The internal SSP-RK2 stages use a maximum loss rate

$$
L=\max(1+D_a r_{\max},\,\beta+D_b r_{\max}),\qquad
r_{\max}=\max_i\frac{\sum_j g_{ij}}{V_i},\qquad \delta t L\le0.2.
$$

This replaces the normalized-graph step restriction. Smaller compartments can require more substeps. Full moving-geometry integration is **first-order operator splitting**, despite second-order integration of frozen reaction/exchange. No higher-order coupled convergence is claimed.

Default diffusivities are $D_a=0.02$, $D_b=0.4$, in model length²/time, with $\beta=2$. They match the exploratory continuum benchmark coefficients; they are not a fitted conversion of the historical rates $1,20$. Biological length/time calibration remains absent. Polarity alignment uses the historical normalized orientation average, separately from concentration transport. `polarity_contact_cutoff=-1` inherits `graph_contact_cutoff` for backward compatibility; setting an explicit value in [0, 1) allows signaling-cutoff interventions with polarity filtering held fixed (see [dynamic cutoff controls](cutoff_dynamics.md)).

## Spectral analysis

For frozen geometry, diagonalize $S=M^{-1/2}KM^{-1/2}$. Its eigenvalues $\lambda_k\ge0$ are those of $-\Delta_V$; right modes are $M^{-1/2}q_k$. Stability uses $J-\lambda_k\operatorname{diag}(D_a,D_b)$, with

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix}.
$$

For these defaults, the physically fixed unstable band is $6.492189<\lambda<38.507811$. The interval formula alone is not evidence of convergence: discrete supported modes must approach their continuum counterparts. Frozen equilibrium spectra also omit dilution during deformation, so they do not certify a moving embryo's stability.

The dashboard and batch runner analyze fixed unit-box reference meshes before stepping and actual cell geometry at recorded frames. Cleavage reports use volume-weighted spectral bases. `embryo.graph_analysis` remains explicitly the historical normalized-cycle experiment. Conservative validation is reproducible with:

```bash
python -m embryo.live_transport_validation --output outputs/live-transport-check
```

## Completed controlled validation

The test constructs smooth phase-field slabs filling a unit-length domain with unit transverse area and reflecting outer walls. Shell overlaps, occupancy volumes, and centers are integrated from the fields, then passed through the **same contact adapter used by the live solver**. Interface width is one sixteenth of slab spacing and is sampled with six points per width. Only longitudinal concentration modes are represented; this is not the complete 3D box spectrum.

| Compartments | First eigenvalue (continuum: 9.8696044) | Relative error | Unstable longitudinal modes |
|---:|---:|---:|---:|
| 8 | 9.7443221 | 1.2694% | 2 |
| 16 | 9.8380551 | 0.3197% | 1 |
| 32 | 9.8616951 | 0.08014% | 1 |
| 64 | 9.8676250 | 0.02006% | 1 |

Final observed order is 1.9985. The continuum has one unstable longitudinal mode, recovered from 16 compartments onward. The corresponding fixed-rate random-walk first eigenvalue instead falls from 0.09903 to 0.001243. At 16 compartments, predicted growth is 0.145412853922 and fitted early nonlinear growth is 0.145412853847 (absolute discrepancy $7.52\times10^{-11}$). Constant preservation and volume-weighted conservation residuals are below $2\times10^{-12}$. All six experiment checks pass. The report is saved in `outputs/live-transport-validation-complete/report.json`.

Regression tests additionally cover unequal-volume flux balance and the full reaction/transport Jacobian, live integration followed by amount-preserving mechanics, cleavage with deliberately mismatched daughter volume, old/new checkpoint selection, and replay of an actual evolving four-cell trajectory including division and dilution.

A coupled 3D smoke run (`outputs/live-conservative-smoke`) completed 600 steps to $t=9$ and reached 16 cells with no active divisions. Its final total-volume error was -1.64%; no cell had crossed the fate-label thresholds. This checks execution through cleavage and export, not mature pattern formation or geometry convergence.

## Limits and next acceptance tests

Conservation is guaranteed by symmetric conductances and positive capacities. **Continuum consistency is not guaranteed by conservation alone.** The overlap calibration assumes complementary equilibrium interfaces; arbitrary cell fields may have gaps, excess overlap, curved contacts, or altered interface profiles. Center-to-center two-point fluxes also need orthogonal geometry for standard finite-volume consistency. Relative cutoff changes can alter connectivity.

Next test voxel refinement at fixed cell geometry and fixed interface width; vary interface width independently; quantify area error on curved/gapped contacts and nonorthogonal flux bias; test cutoff sensitivity; refine coupled time steps; then rerun biological cleavage and shape controls with fixed physical diffusivities and independent seeds. Manufactured compartment refinement is not equivalent to biological cleavage: changing real cell size and shape can legitimately change the physical spectrum. Historical timescale, forced-patch, and shape-persistence results used the older transport model and must be rerun before supporting claims about this new coupling.

## General-geometry validation outcome

The [geometry benchmark](geometric_transport.md) now quantifies curvature bias, diffuse gap leakage, and nonorthogonal linear-field flux errors. It fails all three general-geometry closure screens while passing numerical/reference checks. The calibrated area-over-distance law is not established as a consistent bulk-diffusion discretization on arbitrary deformed contacts. See the benchmark for the distinction between sharp contact, extracellular transport, and phenomenological well-mixed compartment exchange.
