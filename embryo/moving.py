"""Passive amount transport on a prescribed affinely moving 3D L-prism.

No reactions, remeshing, cell division, or mechanically generated motion.
The material mesh is stretched independently along its three coordinate axes.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

from .irregular import (MODES, AMPLITUDES, graded_edges, l_prism_mask,
                        cosine_average, _positive, _rms)
from .transport import conservative_transport, masked_cartesian_transport


MOTIONS = {'isotropic': (.1, .1, .1), 'anisotropic': (.15, .05, .1),
           'volume_preserving': (.12, -.12, 0.)}


def _nonnegative(name, value):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or not np.isrealobj(value):
        raise ValueError(f'{name} must be finite and nonnegative')
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')
    return value


def motion(rates, time):
    """Axis scales exp(g*t), Jacobian, and integrated inverse-square scales."""
    if np.iscomplexobj(rates) or np.asarray(rates).dtype == bool:
        raise ValueError('rates must be three finite real numbers')
    rates = np.asarray(rates, dtype=float)
    time = _nonnegative('time', time)
    if rates.shape != (3,) or not np.isfinite(rates).all():
        raise ValueError('rates must be three finite real numbers')
    with np.errstate(over='raise', under='ignore', invalid='raise', divide='raise'):
        try:
            scales = np.exp(rates * time)
            jacobian = float(np.prod(scales))
            clocks = np.array([time if g == 0 else -np.expm1(-2*g*time)/(2*g) for g in rates])
        except FloatingPointError as error:
            raise ValueError('motion exceeds finite numerical range') from error
    if min(scales) <= 0 or not np.isfinite(jacobian) or jacobian <= 0:
        raise ValueError('motion produces degenerate volumes')
    return scales, jacobian, clocks


class MovingTransport:
    """Material-compartment amounts q=V(t)c; dq/dt=B(t)q.

    B(t)=sum_d exp(-2*g_d*t)*B_d, where B_d=-D*K_d*M0^-1.
    Transposed concentration generators have zero column sums. Splitting their
    matrices by edge direction is algebraic; all directions advance together.
    """

    def __init__(self, n, rates, diffusivity=.02, length=1., grading=.35):
        motion(rates, 0.)
        self.rates = np.asarray(rates, dtype=float).copy()
        self.diffusivity = _nonnegative('diffusivity', diffusivity)
        self.edges, self.mask = graded_edges(n, length, grading), l_prism_mask(n)
        self.reference = masked_cartesian_transport(edges=self.edges, mask=self.mask)
        matrix = self.reference.conductance.tocoo()
        centers = self.reference.centers
        self.operators = []
        for axis in range(3):
            selected = centers[matrix.row, axis] != centers[matrix.col, axis]
            conductance = sparse.coo_matrix((matrix.data[selected],
                (matrix.row[selected], matrix.col[selected])), shape=matrix.shape)
            component = conservative_transport(conductance, self.reference.volumes)
            self.operators.append((self.diffusivity * component.delta.T).tocsr())
        self.exit_rates = np.stack([-operator.diagonal() for operator in self.operators])

    def volumes(self, time):
        return motion(self.rates, time)[1] * self.reference.volumes

    def rhs(self, time, amounts):
        scales = motion(self.rates, time)[0]
        result = np.zeros_like(amounts, dtype=float)
        for scale, operator in zip(scales, self.operators):
            result += (operator @ amounts) / scale**2
        return result

    def positivity_number(self, time, dt):
        """Conservative exit-rate bound over both SSP-RK2 stages."""
        a = motion(self.rates, time)[0]**-2
        b = motion(self.rates, time + dt)[0]**-2
        return float(dt * np.max(np.maximum(a, b) @ self.exit_rates))

    def step(self, amounts, time, dt):
        """Nonautonomous SSP-RK2/Heun; refuse unsafe steps, never clip."""
        dt = _positive('dt', dt)
        if np.iscomplexobj(amounts):
            raise ValueError('amounts must be real')
        q = np.asarray(amounts, dtype=float)
        if (q.ndim not in (1, 2) or q.shape[0] != len(self.reference.volumes)
                or not np.isfinite(q).all() or np.any(q < 0)):
            raise ValueError('amounts must be finite, nonnegative, and match the mesh')
        if self.positivity_number(time, dt) > .8:
            raise ValueError('dt exceeds the SSP positivity bound; reduce dt')
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            first = q + dt * self.rhs(time, q)
            result = .5*q + .5*(first + dt*self.rhs(time + dt, first))
        if not np.isfinite(result).all() or np.any(first < 0) or np.any(result < 0):
            raise FloatingPointError('moving transport lost positivity or finiteness')
        return result

    def geometry_residual(self, time):
        """Compare scaled operator with freshly rebuilt physical face geometry."""
        scales, _, _ = motion(self.rates, time)
        physical = masked_cartesian_transport(edges=tuple(s*e for s, e in zip(scales, self.edges)), mask=self.mask)
        rebuilt = self.diffusivity * physical.delta.T
        scaled = sum(operator / s**2 for operator, s in zip(self.operators, scales))
        difference = rebuilt - scaled
        absolute = np.max(abs(difference.data)) if difference.nnz else 0.
        denominator = np.max(abs(rebuilt.data)) if rebuilt.nnz else 1.
        volume_error = np.max(abs(physical.volumes - self.volumes(time))) / np.max(physical.volumes)
        return {'relative_operator_residual': float(absolute / denominator),
                'relative_volume_residual': float(volume_error)}

    def exact(self, time):
        """Exact continuum cell averages, including material dilution."""
        _, jacobian, clocks = motion(self.rates, time)
        result = np.ones(len(self.reference.volumes))
        for amplitude, mode in zip(AMPLITUDES, MODES):
            exponent = -self.diffusivity * np.pi**2 * (np.square(mode) @ clocks) / self.edges[0][-1]**2
            result += amplitude*np.exp(exponent)*cosine_average(self.edges, mode, self.mask)
        return result / jacobian


def _step_count(duration, dt):
    count = round(duration/dt)
    if count < 1 or not np.isclose(count*dt, duration, rtol=1e-11, atol=1e-13):
        raise ValueError('duration must be a positive integer multiple of dt')
    return count


def trajectory(model, duration=2., dt=.002, snapshots=21):
    duration, dt = _positive('duration', duration), _positive('dt', dt)
    steps = _step_count(duration, dt)
    if type(snapshots) is not int or snapshots < 2:
        raise ValueError('snapshots must be an integer >= 2')
    if model.positivity_number(0., duration) * dt/duration > .8:
        raise ValueError('dt exceeds the full-trajectory SSP positivity bound; reduce dt')
    initial = np.column_stack((np.ones(len(model.reference.volumes)), model.exact(0.)))
    amounts = model.reference.volumes[:, None] * initial
    initial_amount = amounts.sum(axis=0)
    weights = model.reference.volumes / model.reference.volumes.sum()
    saved = set(np.rint(np.linspace(0, steps, min(snapshots, steps+1))).astype(int))
    history, fields, times, scales_history = [], [], [], []
    maximum_mass_drift = maximum_uniform_error = 0.
    minimum_concentration = float('inf')
    for step in range(steps+1):
        time = step*dt
        scales, jacobian, _ = motion(model.rates, time)
        concentration = amounts / (jacobian * model.reference.volumes[:, None])
        mass_drift = float(np.max(abs(amounts.sum(axis=0)/initial_amount - 1)))
        uniform_error = float(np.max(abs(jacobian*concentration[:, 0] - 1)))
        maximum_mass_drift = max(maximum_mass_drift, mass_drift)
        maximum_uniform_error = max(maximum_uniform_error, uniform_error)
        minimum_concentration = min(minimum_concentration, float(concentration.min()))
        if step in saved:
            exact = model.exact(time)
            absolute = _rms(concentration[:, 1] - exact, model.reference.volumes)
            relative = absolute / _rms(exact, model.reference.volumes)
            history.append({'time': time, 'volume': float(model.volumes(time).sum()),
                'volume_ratio': jacobian, 'uniform_mean': float(weights @ concentration[:, 0]),
                'expected_uniform': 1/jacobian, 'patterned_amount': float(amounts[:, 1].sum()),
                'relative_amount_drift': mass_drift, 'uniform_relative_error': uniform_error,
                'absolute_rms_error': absolute, 'relative_rms_error': relative})
            fields.append(concentration[:, 1].copy()); times.append(time); scales_history.append(scales)
        if step < steps:
            amounts = model.step(amounts, time, dt)
    report = {'n': model.mask.shape[0], 'compartments': len(weights), 'rates': model.rates.tolist(),
        'duration': duration, 'dt': dt, 'steps': steps, 'diffusivity': model.diffusivity,
        'maximum_relative_amount_drift': maximum_mass_drift,
        'maximum_uniform_relative_error': maximum_uniform_error,
        'minimum_concentration': minimum_concentration,
        'final_relative_rms_error': history[-1]['relative_rms_error'],
        'geometry_check': model.geometry_residual(duration), 'history': history}
    return report, {'time': np.array(times), 'scales': np.array(scales_history),
                    'concentration': np.array(fields), 'exact_final': model.exact(duration)}


def time_refinement(model, duration=2., steps=(.04, .02, .01)):
    """Independent time accuracy against adaptive DOP853 on this fixed mesh."""
    initial = model.reference.volumes * model.exact(0.)
    reference = solve_ivp(model.rhs, (0., duration), initial, method='DOP853', rtol=1e-12, atol=1e-14)
    if not reference.success:
        raise RuntimeError(reference.message)
    exact = reference.y[:, -1] / model.volumes(duration)
    rows = []
    for dt in steps:
        _, fields = trajectory(model, duration, dt, snapshots=2)
        error = _rms(fields['concentration'][-1] - exact, model.reference.volumes) / _rms(exact, model.reference.volumes)
        rows.append({'dt': dt, 'relative_rms_error': error})
    orders = [float(np.log(a['relative_rms_error']/b['relative_rms_error']) / np.log(a['dt']/b['dt']))
              if min(a['relative_rms_error'], b['relative_rms_error']) > 1e-13 else None
              for a, b in zip(rows, rows[1:])]
    return {'n': model.mask.shape[0], 'reference': 'DOP853 rtol=1e-12, atol=1e-14 (amount variables)',
            'runs': rows, 'orders': orders}


def run_benchmark(output, resolutions=(4, 8, 16, 32), duration=2., dt=.002,
                  diffusivity=.02, length=1., grading=.35):
    output = Path(output)
    if output.exists():
        raise FileExistsError('choose a fresh output directory')
    duration, dt = _positive('duration', duration), _positive('dt', dt)
    diffusivity = _positive('diffusivity', diffusivity)
    _step_count(duration, dt)
    if len(resolutions) < 3 or any(type(n) is not int or n < 4 or n > 32 or n % 2 for n in resolutions):
        raise ValueError('supply at least three even resolutions between 4 and 32')
    if any(b != 2*a for a, b in zip(resolutions, resolutions[1:])):
        raise ValueError('use successive factor-two refinements')
    models = {name: [MovingTransport(n, rates, diffusivity, length, grading) for n in resolutions]
              for name, rates in MOTIONS.items()}
    for group in models.values():
        for model in group:
            if model.positivity_number(0., duration) * dt/duration > .8:
                raise ValueError('dt exceeds positivity bound for this benchmark; reduce dt')
    # Time controls are dimensionless fractions of the same duration and on n=8.
    time_steps = (duration/50, duration/100, duration/200)
    time_models = {name: MovingTransport(8, rates, diffusivity, length, grading) for name, rates in MOTIONS.items()}
    for model in time_models.values():
        if model.positivity_number(0., duration) / 50 > .8:
            raise ValueError('temporal-control step exceeds positivity bound; reduce duration or diffusivity')
    output.mkdir(parents=True)
    report = {'schema': 1, 'experiment': 'prescribed-affine-moving-domain',
              'parameters': {'resolutions': list(resolutions), 'duration': duration, 'dt': dt,
                             'diffusivity': diffusivity, 'length': length, 'grading': grading, 'motions': {k: list(v) for k, v in MOTIONS.items()}},
              'scope': 'Passive material transport; prescribed affine deformation with fixed topology. No remapping, reactions, or emergent shape.',
              'acceptance': {'amount_drift_max': 1e-10, 'uniform_dilution_error_max': 1e-11,
                             'geometry_residual_max': 1e-11, 'finest_relative_rms_max': 1e-3,
                             'spatial_order_min': 1.7, 'temporal_order_range': [1.8, 2.3],
                             'temporal_error_max': 1e-5, 'finest_time_change_to_spatial_error_max': .01,
                             'omitted_dilution_drift_min': .1},
              'spatial': {}, 'temporal': {}, 'finest_time_control': {}}
    archive = {}
    def save():
        (output / 'analysis.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for name, group in models.items():
        report['spatial'][name] = []
        for model in group:
            n = model.mask.shape[0]
            print(f'Moving domain {name}, n={n}, t={duration}, dt={dt}', flush=True)
            row, fields = trajectory(model, duration, dt)
            report['spatial'][name].append(row)
            prefix = f'{name}_n{n}'
            archive.update({f'{prefix}_{key}': value for key, value in fields.items()})
            archive[f'n{n}_mask'] = model.mask
            archive[f'n{n}_volumes_initial'] = model.reference.volumes
            for axis, edges in zip('xyz', model.edges):
                archive[f'n{n}_{axis}_edges_initial'] = edges
            save()
        finest_model = group[-1]
        refined, refined_fields = trajectory(finest_model, duration, dt/2)
        difference = _rms(fields['concentration'][-1] - refined_fields['concentration'][-1], finest_model.reference.volumes)
        difference /= _rms(refined_fields['exact_final'], finest_model.reference.volumes)
        refined['relative_step_halving_difference'] = difference
        refined['difference_to_spatial_error_ratio'] = difference / report['spatial'][name][-1]['final_relative_rms_error']
        report['finest_time_control'][name] = refined
        archive.update({f'{name}_n{resolutions[-1]}_half_dt_{key}': value for key, value in refined_fields.items()})
        report['temporal'][name] = time_refinement(time_models[name], duration, time_steps)
        rows = report['spatial'][name]
        for a, b in zip(rows, rows[1:]):
            b['observed_order'] = float(np.log2(a['final_relative_rms_error']/b['final_relative_rms_error']))
        save()
    # Deliberately incorrect uniform-concentration control: expansion without
    # dilution creates apparent chemical amount despite zero reaction/flux.
    times = archive[f'isotropic_n{resolutions[0]}_time']
    volume_ratios = np.exp(sum(MOTIONS['isotropic'])*times)
    report['omitted_dilution_control'] = {'interpretation': 'Intentionally incorrect: c=1 during expansion',
        'time': times.tolist(), 'relative_amount_drift': (volume_ratios-1).tolist()}
    rows = [row for group in report['spatial'].values() for row in group] + list(report['finest_time_control'].values())
    report['checks'] = {
        'amount_conserved_all_steps': all(r['maximum_relative_amount_drift'] < 1e-10 for r in rows),
        'uniform_dilution_matches_inverse_volume': all(r['maximum_uniform_relative_error'] < 1e-11 for r in rows),
        'concentrations_remain_positive': all(r['minimum_concentration'] > 0 for r in rows),
        'moving_operators_match_rebuilt_geometry': all(max(r['geometry_check'].values()) < 1e-11 for r in rows),
        'spatial_errors_decrease': all(all(a['final_relative_rms_error'] > b['final_relative_rms_error'] for a,b in zip(group,group[1:])) for group in report['spatial'].values()),
        'finest_spatial_order_above_1_7': all(group[-1]['observed_order'] > 1.7 for group in report['spatial'].values()),
        'finest_spatial_error_below_0_1_percent': all(group[-1]['final_relative_rms_error'] < 1e-3 for group in report['spatial'].values()),
        'temporal_orders_near_two': all(all(order is not None and 1.8 < order < 2.3 for order in group['orders']) for group in report['temporal'].values()),
        'fine_temporal_error_below_1e_5': all(group['runs'][-1]['relative_rms_error'] < 1e-5 for group in report['temporal'].values()),
        'finest_time_change_small_relative_to_spatial_error': all(r['difference_to_spatial_error_ratio'] < .01 for r in report['finest_time_control'].values()),
        'omitting_dilution_detectably_violates_balance': float(volume_ratios[-1]-1) > .1,
    }
    save()
    np.savez_compressed(output/'fields.npz', **archive)
    plot_report(report, archive, output)
    return report


def plot_report(report, archive, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for name, rows in report['spatial'].items():
        h = rows[-1]['history']; label = name.replace('_', ' ')
        axes[0].plot([r['time'] for r in h], [r['uniform_mean'] for r in h], label=label)
        axes[0].plot([r['time'] for r in h][::4], [r['expected_uniform'] for r in h][::4], 'k.', ms=4)
        axes[1].loglog([r['n'] for r in rows], [r['final_relative_rms_error'] for r in rows], 'o-', label=label)
        t = report['temporal'][name]['runs']
        axes[2].loglog([r['dt'] for r in t], [r['relative_rms_error'] for r in t], 'o-', label=label)
    axes[0].set(title='Uniform concentration follows 1 / volume', xlabel='Time', ylabel='Concentration')
    axes[1].set(title='Moving-mesh spatial convergence', xlabel='Resolution n', ylabel='Relative RMS vs exact continuum solution')
    axes[2].set(title='Independent time refinement (n=8)', xlabel='Time step', ylabel='Relative RMS vs adaptive ODE reference')
    for ax in axes: ax.legend(fontsize=7)
    fig.suptitle('Prescribed moving domain: dilution and transport accuracy')
    fig.savefig(output/'convergence.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    for name, rows in report['spatial'].items():
        h=rows[-1]['history']
        axes[0].plot([r['time'] for r in h], [r['relative_amount_drift'] for r in h], label=name.replace('_',' '))
    bad=report['omitted_dilution_control']
    axes[1].plot(bad['time'],100*np.array(bad['relative_amount_drift']),color='#bd6445')
    axes[0].set(title='Conservative amount evolution',xlabel='Time',ylabel='Relative amount drift')
    axes[0].legend(fontsize=8)
    axes[1].set(title='Incorrect control: dilution omitted',xlabel='Time',ylabel='Artificial amount gain (%)')
    fig.savefig(output/'balance.png',dpi=160);plt.close(fig)
    n=report['parameters']['resolutions'][-1]
    mask=archive[f'n{n}_mask'];prefix=f'anisotropic_n{n}'
    fields=archive[f'{prefix}_concentration'];scales=archive[f'{prefix}_scales'];times=archive[f'{prefix}_time']
    chosen=[0,len(times)//2,len(times)-1]
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    max_x=archive[f'n{n}_x_edges_initial'][-1]*max(scales[:,0]);max_y=archive[f'n{n}_y_edges_initial'][-1]*max(scales[:,1])
    for ax,k in zip(axes,chosen):
        cube=np.full(mask.shape,np.nan);cube[mask]=fields[k]
        x=archive[f'n{n}_x_edges_initial']*scales[k,0];y=archive[f'n{n}_y_edges_initial']*scales[k,1]
        picture=ax.pcolormesh(x,y,np.ma.masked_invalid(cube[:,:,0].T),shading='flat',cmap='viridis',vmin=fields.min(),vmax=fields.max())
        ax.set(aspect='equal',xlim=(0,max_x),ylim=(0,max_y),title=f't={times[k]:g}',xlabel='Physical x',ylabel='Physical y')
    fig.colorbar(picture,ax=axes,label='Passive concentration',shrink=.75)
    fig.suptitle('Prescribed anisotropic expansion · first material z slab · fixed axes and colors')
    fig.savefig(output/'moving_fields.png',dpi=160);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resolutions',type=int,nargs='+',default=[4,8,16,32])
    for name,value in (('duration',2.),('dt',.002),('diffusivity',.02),('length',1.),('grading',.35)):
        parser.add_argument('--'+name,type=float,default=value)
    args=parser.parse_args()
    try:
        report=run_benchmark(**vars(args))
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))
    print(json.dumps(report['checks'],indent=2),flush=True)
    if not all(report['checks'].values()):
        raise SystemExit('Some moving-domain criteria failed; inspect analysis.json.')


if __name__ == '__main__':
    main()
