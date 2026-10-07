# Polarity and contact-conductance intervention

This follow-up asks whether preserving chemical contact conductances removes the initiation suppression associated with directional polarity tension. The [completed dilution/contact study](moving_initiation_controls_assessment.md) restored formation in histories 7 and 8 at chi=0 by preserving conductances, with changing volumes and dilution retained. The [matched carry polarity study](phase_carry_polarity_assessment.md) demonstrated original polarity-specific suppression only in history 9: its chi=0 branch formed and its chi=0.35 branch did not. Histories 7 and 8 had neither baseline branch forming. These distinctions govern reporting in the new test.

## Matched design

Reuse the same mature t=150 16-cell geometry, polarity, cell IDs, random state and near-uniform chemical preparation within each existing history. Histories 7, 8 and 9 are the three replication units; no new histories or fresh zygotes are generated.

| Directional tension contrast | Chemical conductances | Dilution and measured volumes |
|---|---|---|
| chi=0 | Evolving | Retained |
| chi=0 | Initial symmetric matrix held fixed | Retained |
| chi=0.35 | Evolving | Retained |
| chi=0.35 | Initial symmetric matrix held fixed | Retained |

Both dt=0.00375 and 0.001875 are included. There are 24 numerical paths nested in three histories. **Sixteen qualified paths are reused read-only; eight new paths are needed:** both fixed-conductance contrasts in history 9, and the positive-contrast fixed-conductance branches in histories 7 and 8. Schedule history 9 first. This is a controlled mechanistic comparison, not an expanded parameter map.

Keep beta=2, D_a=0.02, D_b=0.55, activity tension/adhesion 0.25/0.35, the 72-cubed grid, extent 2.24 and interface parameter 0.085. Polarity dynamics, mechanics and chemical feedback remain active. chi=0 removes directional tension, not polarity or its prior effect on the inherited starting geometry.

## Conservative intervention and compute path

For symmetric conductances G, define K=diag(G 1)-G and M(t)=diag(V_i(t)). The fixed-conductance branch uses

$$
\Delta_{\mathrm{fixed}}(t)=-M(t)^{-1}K(0).
$$

G(0) is the actual exported initial GPU conductance matrix. Holding it fixed preserves exchange capacities while rates per compartment volume still evolve. Geometry and the network used for polarity/mechanics continue to change. Never hold Delta(0) fixed against changing volumes. Mechanical concentration conversion by V_old/V_new remains active in every branch, with zero declared volume amount source. Save geometric and actual chemical matrices separately; spectra and frozen endpoint assays use the chemical matrix.

Use resident PyTorch chemistry and matrices with the unchanged custom CUDA carry mechanics/geometry/polarity kernels. Visible fields/contact accumulation remain float32 and phase-update residuals float64. All previously pinned sources and original evidence remain unchanged; the new runner adds the experiment protocol and orchestration only.

## Checks before and during execution

Verify completed parent studies, final reviews, source/input hashes, exact reused intervention metadata and original endpoint evidence. Initial chemistry, polarity, measured volumes and operators must match the reused starts at both contrasts and timesteps; initial G is identical within a history. New actual contexts require sixteen baseline steps exactly matching the unchanged carry backend, eight independent NumPy chemical/accounting steps for each contact arm, and exact restart checks. Existing zero-contrast history-7/8 gates remain separate qualified reuse.

Every pair runs a 60-unit pilot; only a passing raw timestep comparison permits continuation to 240 elapsed units, ending at physical time 390. Retain observations every 0.15 units and atomic checkpoints every three units. Preserve valid initial-time checkpoints on restart, including carry residuals and the initial conductance matrix.

Formation requires SD(log activator)>0.1 throughout elapsed 216–240. Maximum raw log-concentration discrepancy must be at most 0.01, without time alignment. Retain original polarity, shape, volume, actual/geometric transport, growth, crossing and onset-time limits. Require volume-target error below 5%, minimum equivalent radius at least four grid spacings, zero clipping, boundary occupancy below 0.01 and amount-accounting error at most 2e-14. A failed pilot blocks its own long pair. A failure of formation is an admissible scientific outcome, not a numerical failure.

Run the actual endpoint chemical assays with independent solvers and the original settling/stability criteria. Endpoint preparations and graphs are nested within the histories. Stable patterned endpoints, stable uniform chemistry and local coexistence are separate classifications; none proves full moving-system stability.

## Decision rules recorded before results

First qualify all four intervention pairs within a history, including both timesteps and endpoint checks. Then report the evolving-conductance polarity comparison separately from the fixed-conductance comparison.

| Qualified fixed-conductance outcomes | Interpretation |
|---|---|
| Both chi=0 and chi=0.35 form | Conductance preservation permits formation at both contrasts. In history 9, this removes the observed polarity-associated loss under this intervention. In histories 7/8, it is a context rescue without an original demonstrated polarity-specific loss. |
| chi=0 forms; chi=0.35 does not | Directional tension still suppresses formation with fixed conductances; preserving conductances is insufficient. Changing capacities, dilution and other moving interactions remain candidates. |
| Neither forms | No successful fixed-conductance initiating control; do not attribute failure specifically to polarity. |
| Only chi=0.35 forms | Report the nonmonotonic result explicitly; the original suppression narrative does not describe this control. |
| Any required numerical/physical gate fails | Report that history as unresolved. Preserve its curves and failed gate; do not substitute endpoint-only agreement or loosen a tolerance. |

A successful intervention identifies a tested contribution, not unique or complete causal mediation. Mode spectra are frozen-snapshot diagnostics. The test does not establish network homogenization, biological identities, inheritance, autonomous cell states, fresh carry development, spatial convergence, general geometric conductance closure or a broader backend qualification. Broader ledger carry revalidation remains necessary.

## Execution

Sources: [runner](../embryo/polarity_conductance_controls.py), unchanged [GPU control backend](../embryo/gpu_initiation_controls.py), and [targeted tests](../tests/test_polarity_conductance_controls.py).

Output directory: `outputs/polarity-conductance-controls/`. Preparation freezes sources, inputs, acceptance rules, reused evidence and new-context verification before launch. GPU trajectories run sequentially on the qualified GTX 1080 Ti; two CPU workers run independent endpoint assays alongside them. `status.json` and `run.log` record progress. The [execution verification](polarity_conductance_controls_verification.json) records validation and launch status.
