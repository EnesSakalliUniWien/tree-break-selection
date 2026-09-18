"""Contracts for the nuisance-derived Gaussian selected radial prototype."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import log_ndtr
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    RadialSelectionInterval,
    build_gaussian_selected_projection_path,
    compute_selected_gaussian_tail,
)


def _gaussian_path_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    feature_matrix = np.array(
        [
            [2.0, -1.0, 0.5],
            [0.5, 1.5, -0.5],
            [-1.0, 0.25, 2.0],
            [1.25, -0.75, 1.0],
            [-0.5, 2.0, 0.25],
        ],
        dtype=np.float64,
    )
    contrast = np.array([1.0, 1.0, -1.0, -1.0, 0.0], dtype=np.float64)
    covariance = np.array(
        [
            [2.0, 0.3, 0.2],
            [0.3, 1.5, 0.1],
            [0.2, 0.1, 0.8],
        ],
        dtype=np.float64,
    )
    return feature_matrix, contrast, covariance


def test_gaussian_selected_projection_path_preserves_conditioned_coordinates() -> None:
    feature_matrix, contrast, covariance = _gaussian_path_inputs()
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        covariance,
        projection_dimension=2,
    )

    np.testing.assert_allclose(path.projection @ path.projection.T, np.eye(2), atol=1e-12)
    np.testing.assert_allclose(path.reconstruct(path.observed_radius), feature_matrix, atol=1e-12)
    assert path.observed_statistic == pytest.approx(path.observed_radius**2)

    candidate_radius = 0.75
    candidate = path.reconstruct(candidate_radius)
    candidate_nuisance = candidate - np.outer(
        contrast,
        contrast @ candidate,
    ) / float(contrast @ contrast)
    np.testing.assert_allclose(candidate_nuisance, path.nuisance_matrix, atol=1e-12)

    candidate_contrast = contrast @ candidate / path.contrast_norm
    candidate_z = np.linalg.solve(path.covariance_factor, candidate_contrast)
    candidate_projected = path.projection @ candidate_z
    candidate_orthogonal = candidate_z - path.projection.T @ candidate_projected
    np.testing.assert_allclose(
        candidate_projected,
        candidate_radius * path.projected_direction,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        candidate_orthogonal,
        path.orthogonal_whitened_contrast,
        atol=1e-12,
    )

    with pytest.raises(ValueError, match="non-negative"):
        path.reconstruct(-0.1)


def test_gaussian_projection_is_derived_only_from_contrast_orthogonal_nuisance() -> None:
    feature_matrix, contrast, covariance = _gaussian_path_inputs()
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        covariance,
        projection_dimension=2,
    )
    covariance_factor = np.linalg.cholesky(covariance)
    contrast_norm = float(np.linalg.norm(contrast))
    altered_feature_matrix = feature_matrix + np.outer(
        contrast / contrast_norm,
        covariance_factor @ np.array([2.0, -1.0, 0.5]),
    )

    altered_path = build_gaussian_selected_projection_path(
        altered_feature_matrix,
        contrast,
        covariance,
        projection_dimension=2,
    )

    np.testing.assert_allclose(altered_path.nuisance_matrix, path.nuisance_matrix, atol=1e-12)
    np.testing.assert_allclose(
        altered_path.projection.T @ altered_path.projection,
        path.projection.T @ path.projection,
        atol=1e-12,
    )


def test_gaussian_selected_projection_path_rejects_undefined_conditioning() -> None:
    feature_matrix, contrast, covariance = _gaussian_path_inputs()

    with pytest.raises(ValueError, match="non-zero contrast"):
        build_gaussian_selected_projection_path(
            feature_matrix,
            np.zeros_like(contrast),
            covariance,
            projection_dimension=2,
        )

    with pytest.raises(ValueError, match="positive-definite"):
        build_gaussian_selected_projection_path(
            feature_matrix,
            contrast,
            np.diag([1.0, 1.0, 0.0]),
            projection_dimension=2,
        )

    contrast_only_matrix = np.outer(contrast, np.array([1.0, -0.5, 2.0]))
    with pytest.raises(ValueError, match="nuisance PCA"):
        build_gaussian_selected_projection_path(
            contrast_only_matrix,
            contrast,
            covariance,
            projection_dimension=2,
        )


@pytest.mark.parametrize("contrast_scale", [1e-300, 1e300])
def test_gaussian_selected_projection_path_is_stable_under_extreme_contrast_rescaling(
    contrast_scale: float,
) -> None:
    feature_matrix, contrast, covariance = _gaussian_path_inputs()
    reference_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        covariance,
        projection_dimension=2,
    )

    scaled_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast * contrast_scale,
        covariance,
        projection_dimension=2,
    )

    np.testing.assert_allclose(scaled_path.nuisance_matrix, reference_path.nuisance_matrix)
    np.testing.assert_allclose(
        scaled_path.projection.T @ scaled_path.projection,
        reference_path.projection.T @ reference_path.projection,
    )
    assert scaled_path.observed_radius == pytest.approx(reference_path.observed_radius)


def test_selected_gaussian_tail_integrates_a_disconnected_radial_region() -> None:
    result = compute_selected_gaussian_tail(
        observed_radius=3.0,
        projection_dimension=2,
        selection_intervals=(
            RadialSelectionInterval(0.0, 1.0),
            RadialSelectionInterval(2.0, np.inf),
        ),
    )
    expected_selection_probability = (1.0 - np.exp(-0.5)) + np.exp(-2.0)
    expected_selected_tail_probability = np.exp(-4.5)

    assert result.statistic == 9.0
    assert result.degrees_of_freedom == 2
    assert result.selection_probability == pytest.approx(expected_selection_probability)
    assert result.selected_tail_probability == pytest.approx(
        expected_selected_tail_probability
    )
    assert result.p_value == pytest.approx(
        expected_selected_tail_probability / expected_selection_probability
    )


def test_selected_gaussian_tail_reduces_to_unconditional_chi_tail() -> None:
    result = compute_selected_gaussian_tail(
        observed_radius=2.0,
        projection_dimension=2,
        selection_intervals=(RadialSelectionInterval(0.0, np.inf),),
    )

    assert result.selection_probability == 1.0
    assert result.selected_tail_probability == pytest.approx(np.exp(-2.0))
    assert result.p_value == pytest.approx(np.exp(-2.0))


def test_selected_gaussian_tail_requires_observed_radius_inside_selection_event() -> None:
    with pytest.raises(ValueError, match="outside the declared selection event"):
        compute_selected_gaussian_tail(
            observed_radius=1.5,
            projection_dimension=2,
            selection_intervals=(
                RadialSelectionInterval(0.0, 1.0),
                RadialSelectionInterval(2.0, np.inf),
            ),
        )


@pytest.mark.parametrize("dimension", [1, 2, 4])
def test_selected_tail_survives_absolute_probability_underflow(dimension: int) -> None:
    lower, observed = 40.0, 40.01
    result = compute_selected_gaussian_tail(
        observed_radius=observed,
        projection_dimension=dimension,
        selection_intervals=(RadialSelectionInterval(lower, np.inf),),
    )
    if dimension == 1:
        expected = np.exp(log_ndtr(-observed) - log_ndtr(-lower))
    else:
        expected = np.exp(-0.5 * (observed - lower) * (observed + lower))
        if dimension == 4:
            expected *= (1.0 + observed**2 / 2.0) / (1.0 + lower**2 / 2.0)
    assert result.p_value == pytest.approx(expected, rel=1e-11)


def test_selected_tail_handles_disconnected_underflowed_intervals() -> None:
    result = compute_selected_gaussian_tail(
        observed_radius=40.015,
        projection_dimension=2,
        selection_intervals=(
            RadialSelectionInterval(40.0, 40.01),
            RadialSelectionInterval(40.012, 40.02),
        ),
    )

    def scaled_mass(lower, upper):
        return np.exp(-0.5 * (lower - 40.0) * (lower + 40.0)) * (
            -np.expm1(-0.5 * (upper - lower) * (upper + lower))
        )

    expected = scaled_mass(40.015, 40.02) / (
        scaled_mass(40.0, 40.01) + scaled_mass(40.012, 40.02)
    )
    assert result.p_value == pytest.approx(expected, rel=1e-11)


@pytest.mark.parametrize("fraction", [0.0, 0.75, 1.0])
def test_selected_tail_handles_underflowed_lower_tail_and_endpoints(fraction: float) -> None:
    result = compute_selected_gaussian_tail(
        observed_radius=fraction * 1e-90,
        projection_dimension=4,
        selection_intervals=(RadialSelectionInterval(0.0, 1e-90),),
    )
    # In this interval exp(-r^2/2) equals one to floating-point precision;
    # integrating the remaining r^3 density gives this independent ratio.
    assert result.p_value == pytest.approx(1.0 - fraction**4, abs=1e-12)
