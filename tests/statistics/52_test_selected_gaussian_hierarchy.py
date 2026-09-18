"""Exact hierarchy-region contracts for the Gaussian selected radial path."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_hierarchy import (
    construct_squared_euclidean_average_linkage_region,
    replay_squared_euclidean_average_linkage,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    build_gaussian_selected_projection_path,
    compute_selected_gaussian_tail,
)


@pytest.mark.parametrize("epsilon", [1e-6, 1e-8, 1e-10])
def test_nearly_coincident_leaves_preserve_selected_interval_and_tail(epsilon) -> None:
    a = 1.0 - 0.2 * epsilon
    matrix = np.array([[a], [-a], [-1.0], [-1.0 + epsilon]])
    path = build_gaussian_selected_projection_path(
        matrix, np.array([1.0, -1.0, 0.0, 0.0]), np.ones((1, 1)),
        projection_dimension=1,
    )
    region = construct_squared_euclidean_average_linkage_region(path)
    assert len(region.selection_intervals) == 1
    interval = region.selection_intervals[0]
    expected_lower = np.sqrt(2.0) * (1.0 - epsilon / 2.0)
    expected_upper = np.sqrt(2.0) * (1.0 + epsilon)
    assert interval.lower == pytest.approx(expected_lower, abs=1e-15, rel=0)
    assert interval.upper == pytest.approx(expected_upper, abs=1e-15, rel=0)
    for radius in np.linspace(interval.lower, interval.upper, 11)[1:-1]:
        direct = linkage(pdist(path.reconstruct(radius), "sqeuclidean"), method="average")
        assert set(direct[0, :2]) == {1, 2}
        assert set(direct[1, :2]) == {3, 4}
    outside_radius = expected_upper + epsilon / 2.0
    assert replay_squared_euclidean_average_linkage(path, radius=outside_radius) != (
        region.observed_merge_signature
    )
    result = compute_selected_gaussian_tail(
        observed_radius=path.observed_radius, projection_dimension=1,
        selection_intervals=region.selection_intervals,
    )
    assert result.p_value == pytest.approx(0.8, abs=5e-6)


def _three_leaf_radial_path():
    root_two = float(np.sqrt(2.0))
    feature_matrix = np.array(
        [
            [1.0 / root_two],
            [-1.0 / root_two],
            [3.0],
        ],
        dtype=np.float64,
    )
    return build_gaussian_selected_projection_path(
        feature_matrix,
        np.array([1.0, -1.0, 0.0], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
        projection_dimension=1,
    )


def test_squared_euclidean_average_linkage_replays_the_full_merge_signature() -> None:
    path = _three_leaf_radial_path()
    observed_signature = (
        ((0,), (1,)),
        ((0, 1), (2,)),
    )

    assert path.observed_radius == pytest.approx(1.0)
    assert replay_squared_euclidean_average_linkage(
        path,
        radius=path.observed_radius,
    ) == observed_signature
    assert replay_squared_euclidean_average_linkage(path, radius=2.0)[0] == (
        (0,),
        (2,),
    )


def test_squared_euclidean_average_linkage_constructs_exact_observed_region() -> None:
    path = _three_leaf_radial_path()

    region = construct_squared_euclidean_average_linkage_region(path)

    assert region.observed_merge_signature == (
        ((0,), (1,)),
        ((0, 1), (2,)),
    )
    assert region.constraint_count == 2
    np.testing.assert_allclose(
        region.nonnegative_boundary_radii,
        np.array([np.sqrt(2.0), 3.0 * np.sqrt(2.0)]),
        rtol=1e-12,
        atol=1e-12,
    )
    assert len(region.selection_intervals) == 1
    assert region.selection_intervals[0].lower == 0.0
    assert region.selection_intervals[0].upper == pytest.approx(np.sqrt(2.0))

    selected_tail = compute_selected_gaussian_tail(
        observed_radius=path.observed_radius,
        projection_dimension=path.projection_dimension,
        selection_intervals=region.selection_intervals,
    )
    expected_p_value = (
        math.erf(1.0) - math.erf(1.0 / np.sqrt(2.0))
    ) / math.erf(1.0)
    assert selected_tail.p_value == pytest.approx(expected_p_value)


def test_two_leaf_hierarchy_has_the_full_nonnegative_radial_region() -> None:
    feature_matrix = np.array([[1.0], [3.0]], dtype=np.float64)
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        np.array([1.0, 0.0], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
        projection_dimension=1,
    )

    region = construct_squared_euclidean_average_linkage_region(path)

    assert region.observed_merge_signature == (((0,), (1,)),)
    assert region.constraint_count == 0
    assert region.nonnegative_boundary_radii == ()
    assert len(region.selection_intervals) == 1
    assert region.selection_intervals[0].lower == 0.0
    assert region.selection_intervals[0].upper == np.inf


def test_hierarchy_region_preserves_a_narrow_unselected_gap() -> None:
    root_two = float(np.sqrt(2.0))
    observed_radius = 3.0
    feature_matrix = np.array(
        [
            [observed_radius / root_two],
            [-observed_radius / root_two],
            [-5.0],
            [-5.0 + 1e-7],
        ],
        dtype=np.float64,
    )
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        np.array([1.0, -1.0, 0.0, 0.0], dtype=np.float64),
        np.ones((1, 1), dtype=np.float64),
        projection_dimension=1,
    )
    region = construct_squared_euclidean_average_linkage_region(path)
    zero_competitor_radius = 5.0 * root_two
    reselected_radius = 8.0

    assert replay_squared_euclidean_average_linkage(
        path,
        radius=zero_competitor_radius,
    ) != region.observed_merge_signature
    assert not any(
        interval.lower <= zero_competitor_radius <= interval.upper
        for interval in region.selection_intervals
    )
    assert replay_squared_euclidean_average_linkage(
        path,
        radius=reselected_radius,
    ) == region.observed_merge_signature
    assert any(
        interval.lower <= reselected_radius <= interval.upper
        for interval in region.selection_intervals
    )
    assert len(region.selection_intervals) == 2
