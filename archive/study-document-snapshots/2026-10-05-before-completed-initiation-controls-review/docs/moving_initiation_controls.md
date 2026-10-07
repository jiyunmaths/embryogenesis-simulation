# Moving dilution and contact-conductance controls

This experiment investigates why histories 7 and 8 fail to initiate chemical patterns with moving geometry even when directional polarity tension is disabled. The [completed carry study](phase_carry_polarity_assessment.md) qualifies both failures and verifies pattern formation on matched frozen initial operators. The [previously specified factorial design](moving_geometry_initiation_controls.md#first-follow-up-dilution-and-contact-remodeling) supplies the controls and decision rules.

Each developmental history is a replication unit. This assay reuses **two existing histories**, four accepted baseline trajectories and the exact mature t=150 starts. It adds twelve numerical paths, not developmental replicas. History 9's existing forming zero-contrast branch is a positive benchmark.

| Arm | Convert concentration after volume changes? | Chemical conductances |
|---|---|---|
| Baseline, reused | Yes | Evolve with geometry |
| Dilution off | No | Evolve with geometry |
| Fixed conductances | Yes | Hold the actual initial symmetric matrix fixed |
| Combined | No | Hold the actual initial symmetric matrix fixed |

Mechanics, polarity, measured volumes and activity-dependent tension/adhesion co-evolve in every arm. All arms use chi=0, beta=2, D_a=0.02, D_b=0.55, activity tension 0.25 and activity adhesion 0.35. The original 72-cubed grid, extent 2.24 and interface parameter 0.085 are retained. Resident PyTorch chemistry/matrices and the existing custom CUDA mechanics/geometry/polarity use the qualified carry kernel.

## Conservative meaning of fixed contacts

G_ij is symmetric exchange conductance, K=diag(G 1)-G, and M=diag(V_i). The fixed-contact arm uses the actual initial GPU matrix G(0), then recalculates the chemical operator using current volumes:

$$
\Delta_{\mathrm{fixed}}(t)=-M(t)^{-1}K(0).
$$

This preserves total diffusive amount exchange because 1^T M(t) Delta_fixed(t)=0. Holding Delta(0) fixed would generally violate conservation as volumes change. This control holds exchange capacities fixed; chemical exchange rates per unit volume still change through M(t). The geometric contact network continues to govern polarity alignment and mechanical evolution. Both geometric and effective chemical matrices are saved, and spectra and frozen endpoint assays use the effective chemical operator.

## Explicit accounting when dilution is removed

The unchanged splitting first advances reaction/transport, then polarity and mechanics, and then converts concentrations using V_old/V_new. The dilution-off controls omit that final concentration conversion. For concentration c after reaction/transport, the volume conversion therefore contributes chemical amount:

$$
\delta Q_i=(V_i^{\mathrm{new}}-V_i^{\mathrm{old}})c_i.
$$

This source can be positive or negative. It is a diagnostic intervention, not a claim that real cells create molecules through volume change. Save the per-step and cumulative sources for both species and each cell. On-arms must conserve amount in the conversion; off-arms must match the declared source. Both total and compartment-level accounting errors must stay below 2e-14. Reactions continue to produce and remove chemical amount in every arm.

## Validation and observation

The new backend is separate from all pinned earlier scientific sources. Its baseline must match the unchanged phase-carry backend exactly for sixteen steps. Each of the four actual history/timestep contexts receives an independent NumPy float64 chemical/accounting check for all four arms and an exact checkpoint/restart check. Checkpoints retain the fixed initial conductances, initial volumes, float64 phase residual, rounding counters and cumulative chemical-amount sources.

Use dt=0.00375 and 0.001875. Each paired context must pass the original **60-unit pilot** before long continuation to **240 elapsed units**, ending at physical time 390. Require raw maximum log-concentration discrepancy at most **0.01**, with the unchanged volume, polarity, shape, effective transport, growth and timing limits. Also check geometric transport agreement against the original transport tolerance. Compare curves at matched times without temporal alignment. A failed pilot blocks its own long pair; nonformation itself is a scientific outcome.

Formation requires SD(log activator) above 0.1 throughout elapsed 216–240. Retain the 5% volume-error limit, minimum radius of four grid spacings, zero clipping and boundary occupancy below 0.01. Endpoint chemical solver/settling acceptance is separate from whether an endpoint supports a pattern.

Every 0.15-unit observation is retained. Full histories are published with atomic state checkpoints every three units, avoiding repeated serialization of a growing history at every observation. Progress status remains available at each observation. This changes file-writing frequency, not equations, timesteps, sampling or the recoverable checkpoint interval.

## Interpretation

Apply the [factorial decision table](moving_geometry_initiation_controls.md#mechanistic-decision-table) only to qualified full-window results. A single-arm rescue supports that intervention's contribution in its tested context. Rescue only in the combined arm supports a joint restriction. Rescue by either single arm is not a unique-cause result. Single-arm rescue followed by combined-arm failure is a nonmonotonic outcome to investigate. Different histories may have different mechanisms.

If the combined arm also fails, changing-volume dependence of chemical exchange remains a candidate. The next diagnostic would clamp M(0), K(0) and disable dilution while mechanics continues, first requiring reproduction of the actual-input frozen chemical trajectory. That diagnostic, and preparation on a different same-time geometry, are subsequent decisions rather than automatically launched arms.

This test does not establish fresh carry development from a zygote, biological cell identities, inheritance, autonomy, spatial convergence or general geometric transport closure. Earlier ledger claims and original failed checks retain their recorded scope.

## Execution

Sources: [GPU control backend](../embryo/gpu_initiation_controls.py), [experiment runner](../embryo/moving_initiation_controls.py), and [targeted tests](../tests/test_moving_initiation_controls.py).

The output directory is `outputs/moving-initiation-controls/`. Its protocol pins sources, inputs, exact initial operators, parameters and acceptance rules. GPU runs are sequential on the qualified GTX 1080 Ti; two CPU processes handle independent frozen endpoint assays alongside moving computation. `status.json` and `run.log` provide progress. The original carry baseline outputs are reused read-only.

Preparation and actual-context validation precede launch. See the [execution verification record](moving_initiation_controls_verification.json) for their completed status and launch details.
