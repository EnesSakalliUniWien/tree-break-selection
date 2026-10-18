"""The child-parent edge test must hold its chi-square level on null data."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence import (
    annotate_child_parent_divergence,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.spectral_context import (
    compute_child_parent_spectral_context,
    split_edge_test_rows,
)
from tree_break_selection.tree.construction.build import build_tree
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def test_edge_test_rejects_near_alpha_on_iid_bernoulli_null() -> None:
    """Edge p-values on iid Bernoulli leaves must reject at about their nominal rate.

    The tree is built from independent noise, so any excess rejection comes from the
    edge test itself, not from selecting the tree on the tested data.
    """
    rng = np.random.default_rng(0)
    alpha = 0.05
    n_samples, n_features = 120, 40
    rejected = tested = 0
    for _ in range(3):
        data = pd.DataFrame(
            (rng.random((n_samples, n_features)) < 0.5).astype(float),
            index=[f"S{i}" for i in range(n_samples)],
            columns=[f"F{j}" for j in range(n_features)],
        )
        tree = build_tree(
            data, pdist(rng.normal(size=(n_samples, 5))), linkage_method="average"
        ).tree
        tree.populate_node_distributions(data)
        annotated, _ = annotate_child_parent_divergence(
            tree, tree.annotations_df, leaf_data=data, significance_level_alpha=alpha
        )
        has_test = annotated["Child_Parent_Divergence_df"] > 0
        p_values = annotated.loc[has_test, "Child_Parent_Divergence_P_Value"]
        rejected += int((p_values < alpha).sum())
        tested += len(p_values)
        # Both halves of a parent with five or more leaves are large enough to test its edges.
        large_parent_children = [
            child for parent, child in tree.edges if tree.nodes[parent]["leaf_count"] >= 5
        ]
        assert has_test[large_parent_children].all()

    assert alpha / 5 <= rejected / tested <= 2 * alpha


def test_edge_p_values_do_not_depend_on_leaf_row_order() -> None:
    """Shuffling leaf_data rows must leave every edge p-value unchanged."""
    rng = np.random.default_rng(1)
    n_samples, n_features = 60, 20
    group_probabilities = rng.uniform(0.2, 0.8, size=(3, n_features))
    group_of_sample = np.arange(n_samples) % 3
    data = pd.DataFrame(
        (rng.random((n_samples, n_features)) < group_probabilities[group_of_sample]).astype(float),
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"F{j}" for j in range(n_features)],
    )
    tree = build_tree(data, pdist(data.to_numpy()), linkage_method="average").tree
    tree.populate_node_distributions(data)

    def edge_p_values(leaf_data: pd.DataFrame) -> pd.Series:
        annotated, _ = annotate_child_parent_divergence(
            tree, tree.annotations_df, leaf_data=leaf_data, significance_level_alpha=0.05
        )
        return annotated["Child_Parent_Divergence_P_Value"]

    ordered = edge_p_values(data)
    shuffled = edge_p_values(data.sample(frac=1.0, random_state=7))

    pd.testing.assert_series_equal(shuffled.reindex(ordered.index), ordered, rtol=1e-9, atol=1e-12)


def test_edge_basis_sees_between_child_signal_on_separated_continuous_mixture() -> None:
    """The root's edge basis must not be whitened by the between-group spread it should detect."""
    rng = np.random.default_rng(3)
    n_per_group, n_groups, n_features = 40, 4, 24
    shift = np.zeros((n_groups, n_features))
    for group in range(n_groups):
        shift[group, 6 * group : 6 * group + 6] = 4.0
    group_of_sample = np.repeat(np.arange(n_groups), n_per_group)
    data = pd.DataFrame(
        shift[group_of_sample] + rng.normal(size=(n_groups * n_per_group, n_features)),
        index=[f"S{i}" for i in range(n_groups * n_per_group)],
        columns=[f"F{j}" for j in range(n_features)],
    )
    feature_space = continuous_feature_space_from_columns(tuple(data.columns))
    tree = build_tree(data, pdist(data.to_numpy()), linkage_method="average").tree
    tree.populate_node_distributions(data, feature_space=feature_space)
    root = tree.root()

    _, basis_rows_by_parent = split_edge_test_rows(tree, data)
    edge_basis = compute_child_parent_spectral_context(
        tree, data, feature_space=feature_space, basis_row_indices_by_node=basis_rows_by_parent
    )
    annotated, _ = annotate_child_parent_divergence(
        tree,
        tree.annotations_df,
        leaf_data=data,
        feature_space=feature_space,
        significance_level_alpha=0.05,
    )

    assert edge_basis.raw_mp_signal_counts_by_node[root] >= 1
    root_edge_p_values = annotated.loc[
        list(tree.successors(root)), "Child_Parent_Divergence_P_Value"
    ]
    assert (root_edge_p_values < 0.05).all()


def test_edge_test_does_not_over_reject_on_iid_gaussian_null() -> None:
    """Edge p-values on iid continuous leaves must not reject above twice their level.

    Every tested edge counts. There is no lower bound: parents too small for a
    within-child covariance fall back to parent-total whitening, which is conservative.
    """
    rng = np.random.default_rng(2)
    alpha = 0.05
    n_samples, n_features = 240, 5
    rejected = tested = 0
    for _ in range(40):
        data = pd.DataFrame(
            rng.normal(size=(n_samples, n_features)),
            index=[f"S{i}" for i in range(n_samples)],
            columns=[f"F{j}" for j in range(n_features)],
        )
        feature_space = continuous_feature_space_from_columns(tuple(data.columns))
        tree = build_tree(data, pdist(rng.normal(size=(n_samples, 5))), linkage_method="ward").tree
        tree.populate_node_distributions(data, feature_space=feature_space)
        annotated, _ = annotate_child_parent_divergence(
            tree,
            tree.annotations_df,
            leaf_data=data,
            feature_space=feature_space,
            significance_level_alpha=alpha,
        )
        has_test = annotated["Child_Parent_Divergence_df"] > 0
        p_values = annotated.loc[has_test, "Child_Parent_Divergence_P_Value"]
        rejected += int((p_values < alpha).sum())
        tested += len(p_values)

    assert tested > 0
    assert rejected / tested <= 2 * alpha


def test_edge_test_detects_separated_clusters_in_a_node_with_fewer_leaves_than_features() -> None:
    """A node with fewer leaves than features must still test its split.

    Whitening the leaves by their own singular covariance puts them on a simplex, where
    the basis half and the test half are exactly orthogonal and the statistic is 0.
    """
    rng = np.random.default_rng(3)
    n_per_cluster, n_features = 15, 40
    values = rng.normal(size=(2 * n_per_cluster, n_features))
    values[n_per_cluster:] += 3.0
    data = pd.DataFrame(
        values,
        index=[f"S{i}" for i in range(2 * n_per_cluster)],
        columns=[f"F{j}" for j in range(n_features)],
    )
    feature_space = continuous_feature_space_from_columns(tuple(data.columns))
    tree = build_tree(data, pdist(values), linkage_method="average").tree
    tree.populate_node_distributions(data, feature_space=feature_space)

    annotated, _ = annotate_child_parent_divergence(
        tree,
        tree.annotations_df,
        leaf_data=data,
        feature_space=feature_space,
        significance_level_alpha=0.05,
    )

    root = next(node for node in tree.nodes if tree.in_degree(node) == 0)
    root_children = list(tree.successors(root))
    assert sorted(tree.nodes[child]["leaf_count"] for child in root_children) == [15, 15]
    assert (annotated.loc[root_children, "Child_Parent_Divergence_P_Value"] < 1e-3).all()
