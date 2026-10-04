import json
import numpy as np
import pytest
from test_gpu_backend import source
from embryo.feedback_long import save_checkpoint, restore_checkpoint, digest
from embryo.phase_carry_formation import transfer


def checkpoint(tmp_path):
    sim = source(); sim.time = 156.; sim.step_number = round(156/sim.config.dt)
    job = dict(key='phase_carry_fine', arm='phase_carry')
    history = [dict(elapsed=6.)]; audit = dict(wall_seconds=10.)
    path = tmp_path/'old.npz'; save_checkpoint(sim, path, audit, history, 'old-protocol', job)
    with np.load(path) as z: payload = {k:z[k].copy() for k in z.files}
    payload.update(precision_arm=np.array('phase_carry'), phase_carry=np.full(sim.phi.shape, 1.25e-8),
                   rounding=np.array([50, 2], dtype=np.int32))
    np.savez_compressed(path, **payload)
    return path, job, history, audit


def test_transfer_preserves_every_state_array_and_relabels_only_experiment(tmp_path):
    old, job, history, audit = checkpoint(tmp_path); new = tmp_path/'new.npz'; before = digest(old)
    transfer(old, new, 'old-protocol', 'new-protocol', job, 6.)
    assert digest(old) == before
    with np.load(old) as a, np.load(new) as b:
        for key in a.files:
            if key != 'long_experiment': np.testing.assert_array_equal(a[key], b[key])
        assert json.loads(str(b['long_experiment']))['protocol_hash'] == 'new-protocol'
    host, saved_audit, saved_history = restore_checkpoint(new, 'new-protocol', job)
    assert host.time == 156. and saved_audit == audit and saved_history == history


def test_transfer_rejects_wrong_clock_and_lost_residual_precision(tmp_path):
    old, job, _, _ = checkpoint(tmp_path); new = tmp_path/'new.npz'
    with pytest.raises(ValueError, match='clock'): transfer(old, new, 'old-protocol', 'new-protocol', job, 3.)
    with np.load(old) as z: payload = {k:z[k].copy() for k in z.files}
    payload['phase_carry'] = payload['phase_carry'].astype(np.float32)
    np.savez_compressed(old, **payload)
    with pytest.raises(ValueError, match='carry'): transfer(old, new, 'old-protocol', 'new-protocol', job, 6.)
