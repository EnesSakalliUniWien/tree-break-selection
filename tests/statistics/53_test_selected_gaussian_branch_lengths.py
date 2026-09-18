"""Branch-length compatibility checks for the selected Gaussian prototype."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from tree_break_selection.hierarchy_analysis.statistics.contrast_covariance import (
    build_contrast_covariance,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_hierarchy import (
    construct_squared_euclidean_average_linkage_region,
    replay_squared_euclidean_average_linkage,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    build_gaussian_selected_projection_path,
    compute_selected_gaussian_tail,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.tree_estimator import (
    INTERNAL_DISTRIBUTION_BRANCH_LENGTH_STATE,
    compute_spectral_decomposition,
)
from tree_break_selection.tree.construction.hierarchical import tree_from_linkage
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.optimized_branch_lengths import (
    BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
    BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC,
    BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN,
    fit_fixed_topology_nnls_branch_lengths,
)


def _four_leaf_inputs() -> tuple[np.ndarray, np.ndarray]:
    feature_matrix = np.array(
        [
            [-3.0, 0.0],
            [-2.0, 0.0],
            [2.0, 1.0],
            [4.0, 1.0],
        ],
        dtype=np.float64,
    )
    contrast = np.array([1.0, -1.0, 0.0, 0.0], dtype=np.float64)
    return feature_matrix, contrast


def _observed_linkage(feature_matrix: np.ndarray) -> np.ndarray:
    return linkage(
        pdist(feature_matrix, metric="sqeuclidean"),
        method="average",
    )


def _tree_merge_signature(tree) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    members_by_node: dict[object, tuple[int, ...]] = {
        node_id: (int(tree.nodes[node_id]["label"]),)
        for node_id in tree.nodes
        if bool(tree.nodes[node_id]["is_leaf"])
    }
    signature: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    internal_nodes = sorted(
        (node_id for node_id in tree.nodes if not bool(tree.nodes[node_id]["is_leaf"])),
        key=lambda node_id: int(str(node_id).removeprefix("N")),
    )
    for node_id in internal_nodes:
        children = list(tree.successors(node_id))
        left_members, right_members = sorted(members_by_node[child] for child in children)
        signature.append((left_members, right_members))
        members_by_node[node_id] = tuple(sorted((*left_members, *right_members)))
    return tuple(signature)


def _branch_length_state_projection(
    feature_matrix: np.ndarray,
    *,
    branch_length_method: str,
) -> tuple[np.ndarray, np.ndarray]:
    leaf_labels = [str(index) for index in range(len(feature_matrix))]
    columns = ["x", "y"]
    leaf_data = pd.DataFrame(
        feature_matrix,
        index=leaf_labels,
        columns=columns,
    )
    tree = tree_from_linkage(_observed_linkage(feature_matrix), leaf_names=leaf_labels)
    if branch_length_method == BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS:
        optimization = fit_fixed_topology_nnls_branch_lengths(
            tree,
            leaf_data,
            target_metric=BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN,
            pair_sample_size=None,
            solver_tolerance=1e-10,
        )
        assert optimization.status == "ok"
        assert optimization.applied_to_tree is True
    else:
        assert branch_length_method == BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC

    feature_space = continuous_feature_space_from_columns(columns)
    tree.populate_node_divergences(leaf_data, feature_space=feature_space)
    spectral = compute_spectral_decomposition(
        tree,
        leaf_data,
        minimum_projection_dimension=1,
        feature_space=feature_space,
        include_internal_barycenters=True,
        internal_distribution_mode=INTERNAL_DISTRIBUTION_BRANCH_LENGTH_STATE,
        projection_basis_dimension=1,
    )
    projection = spectral.principal_component_projections_by_node[tree.root()]
    branch_lengths = np.asarray(
        [attributes["branch_length"] for *_edge, attributes in tree.edges(data=True)],
        dtype=np.float64,
    )
    return projection.T @ projection, branch_lengths


def test_fixed_branch_time_covariance_rescales_the_exact_radial_region() -> None:
    root_two = float(np.sqrt(2.0))
    feature_matrix = np.array(
        [
            [1.0 / root_two],
            [-1.0 / root_two],
            [3.0],
        ],
        dtype=np.float64,
    )
    contrast = np.array([1.0, -1.0, 0.0], dtype=np.float64)
    branch_time_multiplier = 1.0 + 0.4 / 0.2
    feature_space = continuous_feature_space_from_columns(["x"])
    covariance_inputs = (
        np.array([0.25], dtype=np.float64),
        np.array([0.75], dtype=np.float64),
        2.0,
        2.0,
    )
    base_covariance = build_contrast_covariance(
        *covariance_inputs,
        comparison="sibling",
        feature_space=feature_space,
        continuous_covariance_by_block={"continuous": np.ones(1, dtype=np.float64)},
        ridge=0.0,
    ).covariance_blocks[0]
    adjusted_covariance = build_contrast_covariance(
        *covariance_inputs,
        comparison="sibling",
        feature_space=feature_space,
        continuous_covariance_by_block={"continuous": np.ones(1, dtype=np.float64)},
        tree_time=0.4,
        tree_time_normalizer=0.2,
        ridge=0.0,
    ).covariance_blocks[0]
    np.testing.assert_array_equal(base_covariance, np.ones((1, 1), dtype=np.float64))
    np.testing.assert_array_equal(
        adjusted_covariance,
        np.array([[branch_time_multiplier]], dtype=np.float64),
    )
    base_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        base_covariance,
        projection_dimension=1,
    )
    adjusted_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        adjusted_covariance,
        projection_dimension=1,
    )

    base_region = construct_squared_euclidean_average_linkage_region(base_path)
    adjusted_region = construct_squared_euclidean_average_linkage_region(adjusted_path)
    radius_scale = float(np.sqrt(branch_time_multiplier))

    assert adjusted_path.observed_radius == pytest.approx(
        base_path.observed_radius / radius_scale
    )
    for base_radius in (0.0, 0.5, 1.0, np.sqrt(2.0)):
        np.testing.assert_allclose(
            adjusted_path.reconstruct(base_radius / radius_scale),
            base_path.reconstruct(base_radius),
            atol=1e-12,
        )
    assert len(base_region.selection_intervals) == 1
    assert len(adjusted_region.selection_intervals) == 1
    assert adjusted_region.selection_intervals[0].lower == 0.0
    assert adjusted_region.selection_intervals[0].upper == pytest.approx(
        base_region.selection_intervals[0].upper / radius_scale
    )

    adjusted_tail = compute_selected_gaussian_tail(
        observed_radius=adjusted_path.observed_radius,
        projection_dimension=adjusted_path.projection_dimension,
        selection_intervals=adjusted_region.selection_intervals,
    )
    expected_p_value = (
        math.erf(1.0 / np.sqrt(3.0)) - math.erf(1.0 / np.sqrt(6.0))
    ) / math.erf(1.0 / np.sqrt(3.0))
    assert adjusted_tail.p_value == pytest.approx(expected_p_value)


def test_scalar_covariance_rescaling_does_not_create_a_remote_hierarchy_interval() -> None:
    feature_matrix = np.array(
        [
            [-0.2344109572800355, 0.05870205710974445, 0.11991411011599169],
            [1.6903709879284714, -0.270289915672194, -0.7388767629392597],
            [0.04821780309684464, -1.2351328394153571, 0.32433872604675495],
            [0.4904554073932952, 1.0448685552435102, -0.5772782481357334],
            [0.06717005068541126, 1.0840957934332742, -1.5913222500584427],
        ],
        dtype=np.float64,
    )
    contrast = np.array([1.0, -1.0, 0.5, -0.5, 0.0], dtype=np.float64)
    covariance = np.array(
        [
            [7.257206881857945, -4.246344088677516, 1.293347555646921],
            [-4.246344088677516, 3.950276408170196, -1.2846652987669274],
            [1.293347555646921, -1.2846652987669274, 1.4121793432384928],
        ],
        dtype=np.float64,
    )
    multiplier = 1.456028978217311
    base_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        covariance,
        projection_dimension=2,
    )
    adjusted_path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        multiplier * covariance,
        projection_dimension=2,
    )

    base_region = construct_squared_euclidean_average_linkage_region(base_path)
    adjusted_region = construct_squared_euclidean_average_linkage_region(adjusted_path)
    scale = float(np.sqrt(multiplier))

    assert len(base_region.selection_intervals) == 1
    assert len(adjusted_region.selection_intervals) == 1
    assert adjusted_region.selection_intervals[0].lower == pytest.approx(
        base_region.selection_intervals[0].lower / scale,
    )
    assert adjusted_region.selection_intervals[0].upper == pytest.approx(
        base_region.selection_intervals[0].upper / scale,
    )


def test_fixed_topology_nnls_branch_lengths_do_not_change_the_selected_merge_signature() -> None:
    feature_matrix, contrast = _four_leaf_inputs()
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        np.eye(2, dtype=np.float64),
        projection_dimension=1,
    )
    region = construct_squared_euclidean_average_linkage_region(path)
    linkage_matrix = _observed_linkage(feature_matrix)
    leaf_labels = [str(index) for index in range(len(feature_matrix))]
    tree = tree_from_linkage(linkage_matrix, leaf_names=leaf_labels)
    initial_signature = _tree_merge_signature(tree)
    initial_branch_lengths = np.asarray(
        [attributes["branch_length"] for *_edge, attributes in tree.edges(data=True)],
        dtype=np.float64,
    )
    leaf_data = pd.DataFrame(
        feature_matrix,
        index=leaf_labels,
        columns=["x", "y"],
    )

    optimization = fit_fixed_topology_nnls_branch_lengths(
        tree,
        leaf_data,
        target_metric=BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN,
        pair_sample_size=None,
        solver_tolerance=1e-10,
    )
    adjusted_branch_lengths = np.asarray(
        [attributes["branch_length"] for *_edge, attributes in tree.edges(data=True)],
        dtype=np.float64,
    )

    assert initial_signature == region.observed_merge_signature
    assert _tree_merge_signature(tree) == region.observed_merge_signature
    assert optimization.status == "ok"
    assert optimization.applied_to_tree is True
    assert not np.allclose(adjusted_branch_lengths, initial_branch_lengths)
    assert np.isfinite(adjusted_branch_lengths).all()
    assert np.all(adjusted_branch_lengths >= 0.0)


@pytest.mark.parametrize(
    "branch_length_method",
    [
        BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC,
        BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
    ],
)
def test_branch_length_internal_projection_changes_inside_one_hierarchy_region(
    branch_length_method: str,
) -> None:
    feature_matrix, contrast = _four_leaf_inputs()
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        np.eye(2, dtype=np.float64),
        projection_dimension=1,
    )
    region = construct_squared_euclidean_average_linkage_region(path)
    first_radius = 0.25
    second_radius = 1.0

    assert all(
        any(interval.lower < radius < interval.upper for interval in region.selection_intervals)
        for radius in (first_radius, second_radius)
    )
    assert replay_squared_euclidean_average_linkage(
        path,
        radius=first_radius,
    ) == region.observed_merge_signature
    assert replay_squared_euclidean_average_linkage(
        path,
        radius=second_radius,
    ) == region.observed_merge_signature

    first_projection, first_branch_lengths = _branch_length_state_projection(
        path.reconstruct(first_radius),
        branch_length_method=branch_length_method,
    )
    second_projection, second_branch_lengths = _branch_length_state_projection(
        path.reconstruct(second_radius),
        branch_length_method=branch_length_method,
    )

    assert not np.allclose(first_branch_lengths, second_branch_lengths, atol=1e-6)
    assert not np.allclose(first_projection, second_projection, atol=1e-6)
