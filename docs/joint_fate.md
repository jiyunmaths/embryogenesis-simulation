# Joint signaling and fate integration

The continuous-fate sensitivity in the frozen causal screen is resolved by advancing the activator, inhibitor, and fate together at every SSP-RK2 stage. The original sequential update held the end-of-step activator fixed during the fate update, introducing a coupling lag. The unstable uncommitted fate state amplifies small forcing errors. Both historical failed reports remain unchanged; acceptance tolerances were not relaxed.

The new experiment module, `embryo/joint_fate.py`, stores deviations $u=a-1$ and $v=h-1$. For the full loop, its reaction terms are algebraically equivalent to Gierer–Meinhardt:

$$
\dot u=\frac{(1+u)(u-v)}{1+v}+D_a\Delta_V u,
\qquad
\dot v=\beta(2u+u^2-v)+D_h\Delta_V v.
$$

The downstream switch is evaluated at the same stages:

$$
\dot z=\rho\left(z-z^3+g u\right).
$$

Working in deviations avoids repeated subtraction of nearly equal production and loss terms around the homogeneous equilibrium. Transport acts on deviations as well, preserving the zero-deviation state without cancellation of a constant background. Positivity checks and rate-dependent substeps remain; the integrator does not clip fate or weaken its bistability.

## Independent reference validation

On the original frozen 72³, 16-cell checkpoint, all twenty paired perturbation seeds were integrated to time 60. Joint SSP-RK2 at steps 0.0075 and 0.00375 was compared against DOP853, with relative/absolute tolerances tightened from 1e-10/1e-13 to 1e-12/1e-15. Both solvers use the same equations; the reference independently checks time integration, not the biological model or geometric closure.

| Intervention | Maximum fate error, dt=0.0075 | Maximum fate error, dt=0.00375 |
|---|---:|---:|
| Full loop | 2.839e-5 | 7.108e-6 |
| No self-activation | 3.086e-5 | 7.728e-6 |
| No transport | 6.628e-5 | 1.659e-5 |
| Equal diffusion | 3.757e-5 | 9.409e-6 |
| No signal-to-fate coupling | 0 | 0 |

All original signal (0.01) and continuous-fate (0.05) limits pass, and final labels agree cell by cell for every tested seed and arm. Tightening the reference changes fate by at most 2.35e-9. Errors decrease approximately fourfold on step halving, consistent with second-order integration on this frozen graph. This does not establish second-order accuracy for the full moving-geometry split update.

The corrected reference retains the earlier scientific distinction: full feedback produces persistent signal contrast in 20/20 trials. No self-activation, no transport, and equal diffusion suppress persistent contrast but still generate both fate labels in 20/20. Removing signal-to-fate coupling leaves fate zero while signaling patterns persist. Thus the downstream bistable switch can retain transient signaling differences without a sustained pattern.

Artifacts are in `outputs/joint-fate-validation`, including the frozen protocol, trajectories, reference comparison, and `RESULTS.md`. Reproduce in a fresh directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.joint_fate \
  --output outputs/joint-fate-repeat
```

The implementation is currently used by the [moving causal experiment](moving_causal.md). The core dashboard/developmental solver is unchanged, preserving the ongoing developmental study's frozen source protocol.
