"""Diagnostic dilution/contact interventions; existing carry kernels unchanged."""
from pathlib import Path

import numpy as np
import torch

from .gpu_backend import call, gm_step
from .gpu_precision_control import PrecisionSimulation, force

# dilution enabled, initial conductances held fixed
ARMS = {'baseline': (True, False), 'dilution-off': (False, False),
        'fixed-conductances': (True, True), 'combined': (False, True)}


def conservative_delta(conductance, volumes):
    """Freeze exchange capacities, never the volume-normalized operator."""
    return (conductance-torch.diag(conductance.sum(dim=1)))/volumes[:, None]


def checkpoint_fields(payload, shape):
    """Reject corrupt scientific state before any GPU allocation."""
    n = shape[0]
    if str(payload['precision_arm']) != 'phase_carry' or str(payload['initiation_arm']) not in ARMS:
        raise ValueError('Invalid initiation checkpoint arm')
    shapes = {'phase_carry': shape, 'initial_conductance': (n, n), 'initial_volume': (n,),
              'amount_source': (2, n), 'last_amount_source': (2, n),
              'last_relative_volume_change': (n,), 'conversion_error': ()}
    for key, expected in shapes.items():
        x = payload[key]
        if x.dtype != np.float64 or x.shape != expected or not np.isfinite(x).all():
            raise ValueError('Invalid initiation checkpoint field: '+key)
    g = payload['initial_conductance']
    if (np.any(g < 0) or np.any(np.diag(g) != 0) or not np.array_equal(g, g.T)
            or np.any(payload['initial_volume'] <= 0) or payload['conversion_error'] < 0):
        raise ValueError('Invalid initial conductances, volumes or accounting error')
    counts = payload['rounding']
    if counts.dtype != np.int32 or counts.shape != (2,) or np.any(counts < 0):
        raise ValueError('Invalid carry rounding counters')
    if ARMS[str(payload['initiation_arm'])][0] and (
            np.any(payload['amount_source'] != 0) or np.any(payload['last_amount_source'] != 0)):
        raise ValueError('Dilution-on checkpoint has an explicit amount source')


class InitiationControlSimulation(PrecisionSimulation):
    def __init__(self, host, arm='baseline', device='cuda:0'):
        if arm not in ARMS:
            raise ValueError('Unknown initiation control')
        self.initiation_arm = arm
        self.dilution, self.fixed_contacts = ARMS[arm]
        super().__init__(host, 'phase_carry', device)
        self.initial_conductance = self.matrices()[3].clone()
        self.initial_volume = self.geometry[:, 0].clone()
        self.amount_source = torch.zeros((2, len(self.ids)), dtype=torch.float64, device=self.device)
        self.last_amount_source = torch.zeros_like(self.amount_source)
        self.last_relative_volume_change = torch.zeros_like(self.initial_volume)
        self._conversion_error = self.initial_volume.new_tensor(0.)

    def chemical_transport(self, geometric_delta, geometric_conductance):
        if not self.fixed_contacts:
            return geometric_delta, geometric_conductance
        return conservative_delta(self.initial_conductance, self.geometry[:, 0]), self.initial_conductance

    @torch.inference_mode()
    def step(self):
        c = self.config
        _, geometric_delta, normalized, conductance = self.matrices()
        delta, _ = self.chemical_transport(geometric_delta, conductance)
        self.activator, self.inhibitor = gm_step(self.activator, self.inhibitor, delta, c.dt,
                                                c.signal_beta, c.signal_da, c.signal_dh)
        aligned = (normalized@self.polarity).contiguous()
        p = torch.empty_like(self.polarity)
        call('gpu_polarity', [self.polarity, self.geometry[:, 4:].contiguous(), self.activator, aligned, p],
             (len(self.ids), c.dt, c.polarity_rate, c.polarity_alignment, c.polarity_decay))
        self.polarity = p
        response = torch.tanh(self.activator-1)
        tensions = c.surface_tension*(1+c.fate_tension*response)
        adhesion = c.adhesion*(1+c.fate_adhesion*response[:, None]*response[None, :])
        adhesion.fill_diagonal_(0)
        attraction = (adhesion@self.shell2.flatten(1)).reshape(self.phi.shape)
        old_volume = self.geometry[:, 0].clone()
        centers = self.geometry[:, 1:4].contiguous()
        vf = c.volume_stiffness*(self.target-old_volume)/self.target
        amounts = torch.stack([old_volume@self.activator, old_volume@self.inhibitor])
        before = torch.stack([self.activator, self.inhibitor])
        force([self.phi, self.excluded, attraction, centers, self.polarity, tensions, vf,
               self.gamma, self.updated, self.quality, self.phase_carry, self.rounding],
              (len(self.ids), c.grid, 1, self.dx, c.extent, c.interface_width,
               c.polarity_tension, c.repulsion, c.dt))
        self.phi, self.updated = self.updated, self.phi
        self.refresh()
        new_volume = self.geometry[:, 0]
        self.last_relative_volume_change = new_volume/old_volume-1
        if self.dilution:
            # Preserve the original arithmetic, including in-place multiplication.
            self.activator *= old_volume/new_volume
            self.inhibitor *= old_volume/new_volume
            self.last_amount_source.zero_()
        else:
            self.last_amount_source = before*(new_volume-old_volume)[None, :]
            self.amount_source += self.last_amount_source
        new_amounts = torch.stack([new_volume@self.activator, new_volume@self.inhibitor])
        expected = amounts+self.last_amount_source.sum(dim=1)
        self._dilution_error = ((new_amounts/expected)-1).abs().max()
        # Also check every compartment, avoiding cancellation in the total.
        after = torch.stack([self.activator, self.inhibitor])*new_volume[None, :]
        expected_cells = before*old_volume[None, :]+self.last_amount_source
        self._conversion_error = ((after-expected_cells)/expected_cells.abs()).abs().max()
        self.step_number += 1
        self.time = self.step_number*c.dt

    def audit(self):
        result = super().audit()
        error = float(self._conversion_error)
        if not np.isfinite(error) or not bool(torch.isfinite(self.amount_source).all()):
            raise FloatingPointError('Nonfinite volume-conversion accounting')
        result['volume_conversion_amount_error'] = max(error, result['dilution_amount_error'])
        return result

    def observe(self, elapsed):
        row = super().observe(elapsed)
        _, geometric_delta, _, geometric_conductance = self.matrices()
        delta, conductance = self.chemical_transport(geometric_delta, geometric_conductance)
        row.update(delta=delta.cpu().numpy().tolist(), geometric_delta=geometric_delta.cpu().numpy().tolist(),
                   chemical_conductance=conductance.cpu().numpy().tolist(),
                   geometric_conductance=geometric_conductance.cpu().numpy().tolist(),
                   cumulative_volume_amount_source=self.amount_source.cpu().numpy().tolist(),
                   last_volume_amount_source=self.last_amount_source.cpu().numpy().tolist(),
                   last_relative_volume_change=self.last_relative_volume_change.cpu().numpy().tolist(),
                   volume_conversion_amount_error=float(self._conversion_error), initiation_arm=self.initiation_arm)
        return row

    def checkpoint(self, path):
        path = Path(path)
        temporary = path.with_name(path.stem+'.initiation.tmp.npz')
        super().checkpoint(temporary)
        with np.load(temporary, allow_pickle=False) as z:
            payload = {k: z[k].copy() for k in z.files}
        payload.update(initiation_arm=np.array(self.initiation_arm))
        for key in ('initial_conductance', 'initial_volume', 'amount_source', 'last_amount_source',
                    'last_relative_volume_change'):
            payload[key] = getattr(self, key).cpu().numpy()
        payload['conversion_error'] = np.array(float(self._conversion_error), dtype=np.float64)
        checkpoint_fields(payload, self.phi.shape)
        np.savez_compressed(temporary, **payload)
        temporary.replace(path)

    @classmethod
    def restore(cls, path, device='cuda:0'):
        from .attribute_development import AttributeSimulation
        with np.load(path, allow_pickle=False) as z:
            payload = {k: z[k].copy() for k in ('precision_arm', 'initiation_arm', 'phase_carry', 'rounding',
                'initial_conductance', 'initial_volume', 'amount_source', 'last_amount_source',
                'last_relative_volume_change', 'conversion_error')}
        host = AttributeSimulation.restore(path)
        checkpoint_fields(payload, host.phi.shape)
        sim = cls(host, str(payload['initiation_arm']), device)
        for key in ('phase_carry', 'rounding', 'initial_conductance', 'initial_volume', 'amount_source',
                    'last_amount_source', 'last_relative_volume_change'):
            getattr(sim, key).copy_(torch.from_numpy(payload[key]).to(sim.device))
        sim._conversion_error = torch.as_tensor(payload['conversion_error'], device=sim.device)
        return sim
