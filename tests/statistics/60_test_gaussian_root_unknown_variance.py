"""Scalar selected-root inference with unknown common Gaussian variance."""

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.stats import f, kstest

from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_unknown_variance import (
    _log_half_t_mass,
    selected_scalar_gaussian_root_test,
)


def test_three_rows_have_explicit_truncated_f_reference():
    result = selected_scalar_gaussian_root_test([-1.0, 0.0, 3.0])
    np.testing.assert_allclose(result.intervals, [(np.sqrt(3), np.inf)])
    np.testing.assert_allclose(result.statistic, 49 / 3)
    np.testing.assert_allclose(result.p_value, f.sf(49 / 3, 1, 1) / f.sf(3, 1, 1))


def test_half_t_integral_reduces_to_f_without_selection_and_handles_narrow_interval():
    np.testing.assert_allclose(_log_half_t_mass(0, np.inf, 2), 0, atol=1e-15)
    np.testing.assert_allclose(np.exp(_log_half_t_mass(2, np.inf, 2)), f.sf(8, 1, 2))
    lo, hi = 1.0, 1.0 + 1e-12
    # For nu=1 the density is the half-Cauchy density, 2/(pi*(1+h*h)).
    expected = (hi - lo) / np.pi
    np.testing.assert_allclose(np.exp(_log_half_t_mass(lo, hi, 1)), expected, rtol=1e-10)


def test_two_leaves_per_child_can_detect_strong_separation_without_variance_plugin():
    result = selected_scalar_gaussian_root_test([-100.1, -99.91, 99.9, 100.12])
    assert result.residual_degrees_of_freedom == 2
    assert result.p_value < 0.001
    # Conditioning raises the tail above the unselected fixed-label F tail.
    assert result.p_value > f.sf(result.statistic, 1, 2)


def test_scale_location_and_path_replay():
    values = np.array([-2.1, -1.3, 0.6, 2.7, 3.2, 3.6])
    result = selected_scalar_gaussian_root_test(values)
    transformed = selected_scalar_gaussian_root_test(7 * values + 13)
    np.testing.assert_allclose(transformed.p_value, result.p_value, rtol=1e-11)
    np.testing.assert_allclose(result.reconstruct(result.ratio), values, atol=1e-13)
    for h in np.geomspace(0.01, 100, 40):
        replay = linkage(
            pdist(result.reconstruct(h)[:, None], metric="sqeuclidean"), method="average"
        )
        agrees = np.array_equal(replay[:, :2], result.linkage_matrix[:, :2])
        assert agrees == any(lo < h < hi for lo, hi in result.intervals)


def test_selected_root_null_control_for_four_rows():
    rng = np.random.default_rng(190)
    pvalues = [selected_scalar_gaussian_root_test(rng.normal(size=4)).p_value for _ in range(200)]
    assert kstest(pvalues, "uniform").pvalue > 0.01
    assert 0.04 <= np.mean(np.array(pvalues) < 0.1) <= 0.17


@pytest.mark.parametrize("values", [[1, 2], [1, np.nan, 3], [[1, 2], [3, 4], [5, 6]]])
def test_invalid_scalar_data_rejected(values):
    with pytest.raises(ValueError, match="finite scalar"):
        selected_scalar_gaussian_root_test(values)


def test_zero_residual_and_tied_merge_rejected():
    with pytest.raises(ValueError, match="residual"):
        selected_scalar_gaussian_root_test([-1, -1, 2])
    with pytest.raises(ValueError, match="tie"):
        selected_scalar_gaussian_root_test([-1, 0, 1])
