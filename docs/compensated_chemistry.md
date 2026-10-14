# Compensated chemical accumulation: isolated backend qualification

## Why this numerical branch is needed

The [original delayed amount-matching assay](phase_carry_recipient_amount_delayed_response_assessment.md) retains one failed order check. A nearly stationary control has updates smaller than a float64 concentration's representable spacing. At finer timesteps its stored concentrations never change, although independent solvers and higher precision resolve a roughly 1.57e-13 relaxation. The original run remains unaccepted; its data and criteria are preserved.

This branch changes arithmetic, not the chemical equations. Both species follow the same Gierer–Meinhardt reactions and conservative fixed matrix. It uses resident float64 PyTorch arrays and matrix products on physical GPU 1, with independent float64 NumPy verification. Mechanics, polarity, division and dilution do not advance in this frozen assay; their existing custom CUDA kernels are unchanged.

## Update and numerical correction

Let `u` hold the two species' concentration arrays and let `F(u)` be their unchanged reaction–transport rates. For positivity-restricted substep `h`, calculate

$$
v=u+hF(u),\qquad
\delta u=\frac{h}{2}[F(u)+F(v)].
$$

In exact arithmetic, `u + delta_u` is the same second-order SSP-RK2 update as the existing convex average of two forward Euler stages. The new implementation forms the small increment explicitly and retains rounding error when adding it to the larger concentration:

$$
y=\delta u-e,\qquad
u^{new}=u+y,\qquad
e^{new}=(u^{new}-u)-y.
$$

The array `e` stores the signed rounding correction to subtract from the next increment. It is initialized to zero at each exact prepared start and retained through every substep and timestep. It is numerical bookkeeping, not an additional chemical species, identity variable or biological memory. Saved observations use the rounded positive concentration arrays, without replacing them by a projected equilibrium.

The unchanged maximum-loss rule keeps each forward Euler substep at most 0.2 times its loss coefficient. Inspect both Euler stages, the compensated final concentrations and correction for finiteness and positivity; reject failure without clipping. A checkpoint must retain both concentration and correction arrays, their method/version, clock, history and checksum. Missing correction cannot silently reset to zero. Unit tests require interrupted and uninterrupted production to give identical states, corrections and histories.

## Prespecified numerical qualification

Keep the original runtime, equations, timesteps, sources, clocks and numerical gates. Kernel tests cover sub-representable increments against decimal arithmetic, second-order nonlinear convergence against DOP853, conservative exchange on unequal capacities against matrix exponentiation, positivity substeps, CPU/PyTorch agreement and checkpoint retention/rejection.

A full-window numerical benchmark uses 36 existing delayed initial conditions: all controls and both doses for both chemical origins, with original sham recipient matching and sham shared physics, in both selected recipients of each history. These are 12 controls and 24 pulses, including the diagnosed control. Run all three timesteps through 60 units, with observations every 0.15 and checkpoints every 3 units. Compare GPU and NumPy outputs with the already qualified DOP853/Radau reference trajectories read in place from identical original inputs. Reused reference paths are numerical evidence, not new scientific trajectories or new histories.

Retain the original gates: raw production log error 0.01, CPU/GPU error 1e-8, normalized response error 0.01, activator-area discrepancy 0.02, matching recovery/censoring with time discrepancy 0.3, and order 1.7–2.3 only when the finest error exceeds 100 times dual-reference disagreement. No floor is increased. Also require the diagnosed finest control error to lie below that same 100× noise floor and require its stored trajectory to resolve a nonzero relaxation. Near-noise controls remain unresolved rather than being assigned order two.

Before full-window qualification, check all 144 actual delayed starts and all physical contexts at every timestep through 0.6 units. The actual GPU checkpoint is interrupted at 0.3 and resumed with its correction. Compare that result bit for bit with an uninterrupted GPU calculation and against NumPy raw/normalized errors using the same 1e-8/0.01 gates. No CPU scientific fallback is allowed.

Pin module/test/protocol sources, original source/reference arrays, GPU runtime, prepared benchmark inputs, reports, trajectories and checkpoints. A separate saved-array review must reconstruct the qualification decisions before this backend can authorize a scientific repeat. Tests and benchmark paths are numerical validations nested within three existing histories, not additional developmental replicates.

## Subsequent scientific repeat

Only after qualification passes, the [separate repeat protocol](compensated_recipient_response.md) prepares all delayed starts again from the accepted immediate phase. The original failed delayed study is not an accepted source of scientific trajectories. Its starts must agree exactly with the new preparations. No rematching, new amount clamp or source projection occurs. Shared historical kernels, accepted parents, ledger and manuscript remain unchanged.
