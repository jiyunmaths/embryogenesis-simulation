# Repeated-cleavage timestep diagnosis: prospective protocol

The [completed repeated-cleavage study](repeated_development_gpu_assessment.md) has four passing coarse CPU/GPU comparisons but failed temporal and spatial screens in every mechanical branch. This experiment diagnoses the shared temporal sensitivity with feedback disabled. It changes no model, CUDA kernel, pressure correction, cleavage rule, stochastic rule or acceptance limit. These are synthetic numerical repetitions of seed 7, **zero new independent scientific histories**.

## Fixed design

Keep the original analytic zygote, grid 72 cubed, domain half-width 2.24, interface width 0.085, chemistry, polarity evolution, conservative transport, natural cycles and progressive cleavage. All mechanical feedback coefficients are zero. Follow to t=18 with standard observations every 0.15 and checkpoints every three units. Require sixteen cells, fifteen resolved events, no active furrows, at least three mature units and every inherited physical/event gate.

| Label | Timestep | Steps | Role |
|---|---:|---:|---|
| coarse | 0.00375 | 4,800 | Exact reproduction plus diagnostics of preserved old path |
| fine | 0.001875 | 9,600 | Exact reproduction plus diagnostics of preserved old path |
| half | 0.0009375 | 19,200 | New temporal candidate |
| quarter | 0.00046875 | 38,400 | New temporal candidate |

Four serial GPU paths use PyTorch resident arrays and the unchanged custom CUDA mechanics/geometry kernels on the designated GTX 1080 Ti UUID. The old two are diagnostic reruns, not independent histories or replacements for failed results. They must reproduce the original history, lineage, events and every saved checkpoint array exactly before finer candidates proceed. The complete existing feedback-off NumPy reference is pinned and retained; no full finer-timestep CPU ladder is claimed.

The earlier four concurrent CPU references took 78–87 minutes each at dt=0.00375, versus 77–82 seconds on GPU. Rough linear extrapolation gives 5–6 hours at dt=0.0009375 and 10–12 hours at 0.00046875 per CPU path, under similar throughput. These are estimates; concurrency and later cell count affect timing. Hardware substitution would not establish timestep/grid accuracy, since the tested CPU/GPU paths agree at the same discretization.

## Diagnostic measurements without changing dynamics

Record a network/furrow snapshot every 0.015 physical units, and additionally whenever the active-parent set or completed-event count changes. Record IDs, parent relationships, volumes, centers, raw phase-contact overlap, geometric conductances, conservative operator spectrum, edge count and raw contact-to-cutoff margins. Record each active mother's smooth ring progress, neck value, overlap, prospective partition fraction and margins to abscission thresholds. These snapshots use the current physical clock, without event alignment or spatial registration.

A guarded read-only wrapper around the original pressure routine records every mother's capacity, represented volume before/after projection and returned correction magnitude. It forwards the original arguments and returns the original result tensors; it performs no input/output writes and restores the function on exit or error. Correction sign and per-iteration Newton residuals are not exposed by the unchanged routine and are not claimed. The existing step audit retains the maximum correction magnitude. The full projection log records every dividing-mother call, rather than subsampling brief pressure effects.

Before production, compare instrumented and uninstrumented stepping at all three preserved live-furrow mixed-count checkpoints (2→4, 4→8, 8→16): eight steps per checkpoint must be bitwise equal, including arrays, phase/chemical corrections, geometry, lineage, targets, events, clocks and all random states. Read-only sampling before and after steps is included. All four production paths also retain the inherited eight-step live-window restart checks and independent initial chemical/dilution prefix checks. Unit tests verify original tensor/input preservation and restoration after failures.

## Unchanged limits and interpretation

Compare all three adjacent timestep pairs at the same 121 original observation times and seven native field checkpoints. Retain each failed pair. Inherit every temporal limit from the repeated-development protocol: distributions 0.01, contrast 0.005, axis ratio 0.01, total volume/whole amounts 0.005, mean polarity 0.01, spectrum 0.01, occupancy/chemical density relative L2 0.01, polarity-density occupancy-normalized L2 0.01 and maturity timing 0.05. Keep count-mismatch allowance and physical/event limits unchanged.

Report error sequences and where differences first grow, not merely final states. The finest pair must pass all inherited temporal limits before claiming it clears this feedback-off time screen. Error decay is descriptive and does not establish a formal convergence order. A failed old pair remains failed even if a finer pair passes. The final result remains provisional pending independent review; this study alone cannot qualify other mechanical controls, spatial accuracy, full restart behavior or new developmental biology.

Use the logs to examine whether furrow deformation/projection differences precede division-clock shifts and conductance differences, or whether contact threshold changes dominate. Temporal association does not prove causality. If finer-step errors stall, specify an isolated update control before changing the method. After a temporal pass, repeat spatial qualification at an appropriately small timestep; fresh independent zygotes remain a later gated test.

All inherited solver/evidence files are pinned. Preserve the previous failed study and completed independent review, reports, accepted evidence, ledger and manuscript. No thresholds are relaxed after inspecting outcomes. The runner refuses reused directories/recorded paths and requires the instrumentation gate before production.

Run `python -m embryo.development_time_diagnosis prepare`, then `gate`, then `run`, with the explicit output directory and designated CUDA device. [Inherited completed assessment](repeated_development_gpu_assessment.md).
