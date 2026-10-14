# Initial amount matching: launch record

The [frozen two-phase protocol](phase_carry_recipient_amount_response.md) matches initial recipient chemistry and both surrounding species totals within each shared physical context. Non-target concentrations are rescaled to the symmetric geometric-mean species budget, preserving each species’ spatial allocation fractions. Reactions and transport then evolve freely.

All eleven tests and the independent literal preparation audit pass. All 96 immediate initial-value problems are physically distinct. The actual GPU prefix checks every path at three timesteps against the CPU solver. Nine physical blocks of width16 use48 padding columns; padding is not independent data. Eight single-thread CPU workers verify GPU production.

The immediate phase has48 controls and48 pulses, giving288 GPU paths,288 CPU verification paths and192 independent reference solves. The delayed phase is planned after independent acceptance of the immediate phase:48 controls and96 pulses after a60-unit wait, giving432 GPU/CPU paths and288 references. No rematching at the later challenge. Source qualification can block that phase.

Outputs use the writable `output/studies/phase-carry-recipient-amount-response` and `output/studies/phase-carry-recipient-amount-delayed-response` directories. Accepted parents, shared kernels, ledger and manuscript remain pinned. Both phases are nested within histories7–9, with no new developmental history.

[Launch evidence](../archive/study-launches/2026-10-08-recipient-amount-response/launch-verification.json), [prelaunch snapshot](../archive/study-document-snapshots/2026-10-08-before-recipient-amount-response-launch/manifest.json).
