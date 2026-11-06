"""Restricted analytic Gaussian selected-root diagnostic contracts."""

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.stats import chi

from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_region import (
    selected_gaussian_root_test,
)


def test_two_rows_have_unselected_chi_law():
    result = selected_gaussian_root_test(np.array([[0.0, 1.0], [2.0, 3.0]]), sigma=2.0)
    assert result.intervals == ((0.0, np.inf),)
    assert result.p_value == pytest.approx(chi.sf(result.radius, 2))


def test_three_rows_have_exact_merge_boundary_and_replay():
    x = np.array([[0.0], [1.0], [4.0]])
    result = selected_gaussian_root_test(x, sigma=1.0)
    # First merge is (0,1); its fixed distance 1 must beat distance (1,2).
    # Path coordinates: (0.5-r/sqrt(6), 1.5-r/sqrt(6), 1+2r/sqrt(6)).
    boundary = np.sqrt(6.0) / 2.0
    assert result.intervals[0][0] == pytest.approx(boundary)
    assert result.intervals[0][1] == np.inf
    assert result.p_value == pytest.approx(chi.sf(result.radius, 1) / chi.sf(boundary, 1))
    for radius in [0.1, boundary * 0.9, boundary * 1.1, result.radius, 20.0]:
        candidate = result.offset + radius * result.slope
        z = linkage(pdist(candidate, metric="sqeuclidean"), method="average")
        selected = tuple(z[0, :2].astype(int)) == (0, 1)
        inside = any(lo < radius < hi for lo, hi in result.intervals)
        assert selected == inside
    np.testing.assert_allclose(result.offset + result.radius * result.slope, x)
    np.testing.assert_allclose(result.offset.T @ result.contrast, 0.0, atol=1e-14)


def test_extreme_tail_selected_mass_does_not_underflow():
    result = selected_gaussian_root_test(np.array([[0.0], [100.0], [400.0]]), sigma=1.0)
    assert np.isfinite(result.log_selection_probability)
    assert 0.0 <= result.p_value <= 1.0


@pytest.mark.parametrize(
    "x,sigma",
    [
        ([[0.0], [1.0]], 0.0),
        ([[0.0], [1.0]], np.nan),
        ([[0.0], [0.0]], 1.0),
        ([[0.0], [1.0], [2.0]], 1.0),
    ],
)
def test_invalid_scale_degeneracy_and_observed_ties_are_explicit(x, sigma):
    with pytest.raises(ValueError):
        selected_gaussian_root_test(np.array(x), sigma=sigma)


def test_four_rows_replay_complete_merge_sequence_over_path():
    x = np.array([[0.0, 0.0], [0.3, 1.0], [3.0, 0.0], [4.0, 2.0]])
    result = selected_gaussian_root_test(x, sigma=1.0)
    for radius in np.linspace(0.01, 12.0, 121):
        candidate = result.offset + radius * result.slope
        rebuilt = linkage(pdist(candidate, metric="sqeuclidean"), method="average")
        selected = np.array_equal(rebuilt[:, :2], result.linkage_matrix[:, :2])
        inside = any(lo < radius < hi for lo, hi in result.intervals)
        assert selected == inside


def test_quadratic_solver_preserves_disconnected_and_narrow_cells():
    from fractions import Fraction

    from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_region import (
        _quadratic_region,
    )

    assert _quadratic_region(tuple(map(Fraction, [-2, 3, -1]))) == [(0.0, 1.0), (2.0, np.inf)]
    delta = Fraction(1, 10**10)
    # (r-1)*(r-1-delta) <= 0.
    assert _quadratic_region((1 + delta, -2 - delta, Fraction(1))) == [(1.0, 1.0 + float(delta))]


def test_narrow_far_tail_interval_has_positive_log_mass():
    from scipy.integrate import quad

    from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_region import (
        _log_chi_mass,
    )

    lo, hi = 40.0, 40.0 + 1e-10
    log_mass = _log_chi_mass(lo, hi, 2)
    scale = chi.logpdf(lo, 2)
    reference = scale + np.log(quad(lambda r: np.exp(chi.logpdf(r, 2) - scale), lo, hi)[0])
    assert log_mass == pytest.approx(reference, abs=1e-9)
