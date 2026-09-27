import numpy as np
import pytest

from embryo.geometry_transport_validation import flat_overlap, adapter, spherical_area, flux_patch, prepare, run


def test_flat_calibration_and_gap_leakage_are_distinct():
    epsilon=.05
    assert float(adapter(flat_overlap(epsilon),epsilon).conductance[0,1])==pytest.approx(1.,rel=1e-10)
    gap=adapter(flat_overlap(epsilon,epsilon),epsilon)
    assert 0<float(gap.conductance[0,1])<1
    assert gap.conductance.nnz==2  # Relative cutoff cannot reject the sole weak edge.


def test_spherical_area_bias_has_expected_small_width_scaling():
    radius=.5;exact=4*np.pi*radius**2
    errors=[spherical_area(radius,e)/exact-1 for e in (.025,.05)]
    assert all(x>0 for x in errors)
    assert errors[1]/errors[0]==pytest.approx(4.,rel=1e-4)


def test_nonorthogonal_flux_detects_tangential_contamination_despite_conservation():
    orthogonal=flux_patch(0,[0,1,0])
    skew=flux_patch(1,[0,1,0])
    assert orthogonal['estimated_flux']==0
    assert skew['exact_flux']==0
    assert skew['estimated_flux']==pytest.approx(1/np.sqrt(2),rel=1e-10)
    assert skew['projected_distance_flux']==1
    assert skew['mass_residual']<1e-12


def test_complete_validation_reports_physical_failures_without_overwrite(tmp_path):
    output=tmp_path/'validation';prepare(output)
    original=(output/'protocol.json').read_bytes()
    result=run(output)
    assert (output/'protocol.json').read_bytes()==original
    assert all(result['numerical_checks'].values())
    assert not result['passed']
    assert not any(result['physical_closure_checks'].values())
    with pytest.raises(ValueError):run(output)
    with pytest.raises(FileExistsError):prepare(output)
