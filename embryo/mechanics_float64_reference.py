"""Independent NumPy/SciPy mechanics reference; no chemical/polarity evolution.

All reference fields, geometry, coefficients and updates are float64. The CUDA
comparison retains its actual mixed arithmetic. This is a short component
assay, not a reference solution for the full moving chemical system.
"""
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import laplace

from .feedback_long import digest
from .resolution import _steps, write_json


CRITERIA = dict(phase_abs_max=5e-6, phase_relative_l2=2e-6,
                relative_volume_max=2e-6, clipping_max=0.,
                reference_pair_error_ratio_min=1.5)
DTS = (.00375, .001875, .0009375)
HORIZON = .15


def geometry(phi, dx, extent):
    """Occupancy integrals and centroids, with independent double reductions."""
    x = np.asarray(phi, dtype=np.float64)
    h = x*x*(3.-2.*x)
    volume = h.sum(axis=(1, 2, 3))*dx**3
    if np.any(volume <= 0) or not np.isfinite(volume).all():
        raise ValueError('Positive finite reference cell volumes required')
    coordinate = (np.arange(x.shape[1], dtype=np.float64)+.5)*dx-extent
    center = np.stack([np.einsum('cijk,i->c', h, coordinate),
                       np.einsum('cijk,j->c', h, coordinate),
                       np.einsum('cijk,k->c', h, coordinate)], axis=1)*dx**3/volume[:, None]
    return volume, center


def flux_divergence(phi, gamma, dx):
    """Face-average no-flux operator via an independent Laplacian identity.

    div(gamma grad u) = [L(gamma*u) - u*L(gamma) + gamma*L(u)]/2
    for the chosen discrete face averages. Nearest ghosts give zero flux.
    """
    u, g = np.asarray(phi, dtype=np.float64), np.asarray(gamma, dtype=np.float64)
    return (laplace(g*u, mode='nearest')-u*laplace(g, mode='nearest')+
            g*laplace(u, mode='nearest'))/(2.*dx**2)


def coefficients(phi, activator, polarity, target, c):
    x = np.asarray(phi, dtype=np.float64)
    dx = 2*c.extent/c.grid
    volume, center = geometry(x, dx, c.extent)
    response = np.tanh(np.asarray(activator, dtype=np.float64)-1.)
    tension = c.surface_tension*(1.+c.fate_tension*response)
    adhesion = c.adhesion*(1.+c.fate_adhesion*np.outer(response, response))
    np.fill_diagonal(adhesion, 0.)
    shell = x*(1.-x)
    attraction = (adhesion@(shell*shell).reshape(len(x), -1)).reshape(x.shape)
    squares = x*x
    excluded = squares.sum(axis=0)[None]-squares
    vf = c.volume_stiffness*(np.asarray(target, dtype=np.float64)-volume)/target
    gamma = np.empty_like(x)
    coordinate = (np.arange(c.grid, dtype=np.float64)+.5)*dx-c.extent
    for i in range(len(x)):
        if c.polarity_tension == 0:
            gamma[i].fill(tension[i])
        else:
            rx = coordinate[:, None, None]-center[i, 0]
            ry = coordinate[None, :, None]-center[i, 1]
            rz = coordinate[None, None, :]-center[i, 2]
            radius = np.sqrt(rx*rx+ry*ry+rz*rz+c.interface_width**2)
            dot = polarity[i, 0]*rx+polarity[i, 1]*ry+polarity[i, 2]*rz
            gamma[i] = tension[i]*(1.-c.polarity_tension*dot/radius)
    return gamma, vf, excluded, attraction


def force_from_coefficients(phi, gamma, vf, excluded, attraction, dx, width, repulsion):
    """Float64 discrete PDE RHS, independent of CUDA force/update expressions."""
    x = np.asarray(phi, dtype=np.float64)
    rhs = np.empty_like(x)
    for i, u in enumerate(x):
        # q'(u)=2u-6u^2+4u^3; h'(u)=6u(1-u).
        dq = 2.*u-6.*u*u+4.*u*u*u
        rhs[i] = width**2*flux_divergence(u, gamma[i], dx)-gamma[i]*dq
        rhs[i] += vf[i]*6.*u*(1.-u)-repulsion*u*excluded[i]+dq*attraction[i]
    return rhs


def rhs(phi, activator, polarity, target, config):
    return force_from_coefficients(phi, *coefficients(phi, activator, polarity, target, config),
                                   2*config.extent/config.grid, config.interface_width, config.repulsion)


def frozen_gamma_energy(phi, gamma, adhesion, target, dx, width, stiffness, repulsion):
    """Energy for checking the RHS, holding gamma fixed during differentiation.

    Geometry dependence of gamma is an explicit coupling, not differentiated
    by the model. This check therefore freezes gamma, not volumes or overlaps.
    """
    x = np.asarray(phi, dtype=np.float64)
    q = x*x*(1.-x)**2
    energy = np.sum(gamma*q)*dx**3
    for axis in (1, 2, 3):
        left, right = [slice(None)]*4, [slice(None)]*4
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        a, b = tuple(left), tuple(right)
        energy += .25*width**2*dx*np.sum((gamma[a]+gamma[b])*(x[b]-x[a])**2)
    volumes = np.sum(x*x*(3.-2.*x), axis=(1, 2, 3))*dx**3
    energy += np.sum(stiffness/(2*target)*(volumes-target)**2)
    for i in range(len(x)):
        for j in range(i):
            energy += np.sum(.5*repulsion*x[i]**2*x[j]**2-adhesion[i, j]*q[i]*q[j])*dx**3
    return float(energy)


def gpu_mechanics_step(sim):
    """Existing CUDA mechanics only; freeze concentrations and polarity.

    Refresh geometry/overlaps after each update; do not evolve chemistry,
    polarity or dilute concentrations. Not GpuSimulation.step().
    """
    import torch
    from .gpu_precision_control import force
    c = sim.config
    response = torch.tanh(sim.activator-1)
    tension = c.surface_tension*(1+c.fate_tension*response)
    adhesion = c.adhesion*(1+c.fate_adhesion*response[:, None]*response[None, :])
    adhesion.fill_diagonal_(0)
    attraction = (adhesion@sim.shell2.flatten(1)).reshape(sim.phi.shape)
    vf = c.volume_stiffness*(sim.target-sim.geometry[:, 0])/sim.target
    force([sim.phi, sim.excluded, attraction, sim.geometry[:, 1:4].contiguous(),
           sim.polarity, tension, vf, sim.gamma, sim.updated, sim.quality,
           sim.phase_carry, sim.rounding],
          (len(sim.ids), c.grid, int(sim.precision_arm == 'phase_carry'), sim.dx,
           c.extent, c.interface_width, c.polarity_tension, c.repulsion, c.dt))
    sim.phi, sim.updated = sim.updated, sim.phi
    sim.refresh()


def run(root, contexts, criteria=CRITERIA, horizon=HORIZON, dts=DTS):
    """Compare two real geometries at three dt values; cache final field evidence."""
    import torch
    from .attribute_development import AttributeSimulation
    from .gpu_precision_control import PrecisionSimulation
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    rows, field_files, reference_pairs = [], {}, []
    for context in contexts:
        reference_fields = []
        for dt in dts:
            host = AttributeSimulation.restore(context['source'])
            if context.get('chemical_file'):
                with np.load(context['chemical_file']) as z:
                    if not np.array_equal(z['ids'], host.ids): raise ValueError('Reference IDs changed')
                    host.activator, host.inhibitor = z['uniform'].copy()
                host.config.polarity_tension = 0.
            host.config.dt = dt
            ref = host.phi.astype(np.float64)
            sims = {arm:PrecisionSimulation(host, arm) for arm in ('baseline', 'phase_carry')}
            started = time.perf_counter(); max_clipping = 0.; errors = {
                arm:dict(phase_abs_max=0., phase_relative_l2=0., relative_volume_max=0.) for arm in sims}
            first_update = {}
            for step in range(1, _steps(horizon, dt)+1):
                ref += dt*rhs(ref, host.activator, host.polarity, host.target, host.config)
                if not np.isfinite(ref).all(): raise RuntimeError('Nonfinite float64 reference field')
                max_clipping = max(max_clipping, float(np.max(np.maximum(-ref, ref-1.), initial=0.)))
                if max_clipping > criteria['clipping_max']: raise RuntimeError('Float64 reference clipping')
                rv, _ = geometry(ref, host.dx, host.config.extent)
                for arm, sim in sims.items():
                    gpu_mechanics_step(sim)
                    visible = sim.phi.cpu().numpy().astype(np.float64)
                    error = errors[arm]
                    error['phase_abs_max'] = max(error['phase_abs_max'], float(abs(visible-ref).max()))
                    error['phase_relative_l2'] = max(error['phase_relative_l2'], float(
                        np.linalg.norm(visible-ref)/np.linalg.norm(ref)))
                    gv, _ = geometry(visible, host.dx, host.config.extent)
                    error['relative_volume_max'] = max(error['relative_volume_max'], float(abs(gv/rv-1).max()))
                    if int(sim.quality[0]) or int(sim.quality[1]): raise RuntimeError('GPU reference clipping/invalid field')
                    if step == 1:
                        reconstructed = visible+sim.phase_carry.cpu().numpy()
                        first_update[arm] = float(abs(reconstructed-ref).max())
                write_json(root/'status.json', dict(state='running', context=context['key'], dt=dt,
                    step=step, total=_steps(horizon, dt), errors=errors))
            file = root/f"{context['key']}-dt-{dt:g}.npz"
            np.savez_compressed(file, reference=ref, baseline=sims['baseline'].phi.cpu().numpy(),
                phase_carry=sims['phase_carry'].phi.cpu().numpy(),
                residual=sims['phase_carry'].phase_carry.cpu().numpy(), ids=host.ids)
            field_files[str(file.resolve())] = digest(file); reference_fields.append(ref)
            carry_pass = all(errors['phase_carry'][k] <= criteria[k] for k in errors['phase_carry'])
            rows.append(dict(context=context['key'], dt=dt, horizon=horizon, errors=errors,
                first_update_reconstructed_abs_max=first_update, clipping_max=max_clipping,
                phase_carry_pass=bool(carry_pass), wall_seconds=time.perf_counter()-started))
            print(f"float64 mechanics: {context['key']}, dt={dt:g}, carry={errors['phase_carry']}, pass={carry_pass}", flush=True)
            del sims; torch.cuda.empty_cache()
        e = [float(abs(reference_fields[i]-reference_fields[i+1]).max()) for i in (0, 1)]
        ratio = e[0]/e[1] if e[1] > 0 else None
        reference_pairs.append(dict(context=context['key'], adjacent_phase_abs_max=e, error_ratio=ratio,
            observed_pair_order=None if ratio is None else float(np.log2(ratio)),
            passed=ratio is not None and ratio >= criteria['reference_pair_error_ratio_min']))
    result = dict(passed=all(r['phase_carry_pass'] for r in rows) and all(r['passed'] for r in reference_pairs),
        contexts=contexts, criteria=criteria, horizon=horizon, dts=list(dts), results=rows,
        reference_timestep_pairs=reference_pairs, field_sha256=field_files,
        scope='Two saved states from ONE developmental history; six nested mechanics-only comparisons. Freeze activator and polarity, independently recompute evolving float64 fields/volumes/centers/overlaps. CUDA keeps its mixed arithmetic; residual starts at zero in this component assay. No chemical/polarity integration, cleavage, long-horizon reference or spatial convergence.')
    write_json(root/'result.json', result); write_json(root/'status.json', dict(state='completed', passed=result['passed']))
    return result
