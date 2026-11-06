"""Variance and fixed-topology controls for the attribution diagnostic."""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.diffusion_attribution import (
    covariance_for_node,
    full_space_sibling_records,
    root_covariance_diagnostics,
)
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.poset_tree import PosetTree
from tree_break_selection.tree.topology import tree_root


def test_identity_full_space_statistic_uses_sibling_mean_variance_and_preserves_tree():
    data = pd.DataFrame([[-2.0, 0.0], [-1.0, 1.0], [1.0, 2.0], [3.0, -1.0]], columns=["A", "B"])
    data.index = [f"S{i}" for i in range(len(data))]
    space = continuous_feature_space_from_columns(data.columns)
    tree = build_tree(data, pdist(data.values), linkage_method="average").tree
    before = tuple(tree.edges())
    lengths = [tree.edges[e]["branch_length"] for e in before]
    populated, records = full_space_sibling_records(tree, data, space, "known_identity")
    root = tree_root(populated)
    record = next(record for record in records if record.parent == root)
    left, right = (populated.nodes[node] for node in (record.left, record.right))
    variance = 1 / left["leaf_count"] + 1 / right["leaf_count"]
    expected = np.sum((left["distribution"] - right["distribution"]) ** 2) / variance
    np.testing.assert_allclose(record.stat, expected, rtol=1e-10)
    np.testing.assert_allclose(record.p_value, chi2.sf(expected, 2), rtol=1e-10)
    assert record.degrees_of_freedom == 2
    assert tuple(tree.edges()) == before
    assert [tree.edges[e]["branch_length"] for e in before] == lengths
    assert "distribution" not in tree.nodes[root]
    assert all(populated.edges[e]["branch_length"] == 0 for e in before)


def test_independent_gaussian_fixed_membership_identity_control():
    rng = np.random.default_rng(19)
    # Conditional on a fixed split, independent Gaussian means have this exact covariance.
    n_left, n_right, p = 3, 77, 24
    difference = rng.normal(size=(50000, p)) * np.sqrt(1 / n_left + 1 / n_right)
    statistics = np.sum(difference**2, axis=1) / (1 / n_left + 1 / n_right)
    assert abs(statistics.mean() - p) < 0.15
    assert 0.008 < np.mean(chi2.sf(statistics, p) < 0.01) < 0.012


@pytest.mark.parametrize(
    "n_left,source,representation",
    [
        (7, "parent_total_fallback", "dense"),
        (9, "pooled_within_child", "diagonal"),
        (25, "pooled_within_child", "dense"),
    ],
)
def test_actual_guarded_covariance_regimes(n_left, source, representation):
    data = pd.DataFrame(
        np.random.default_rng(8).normal(size=(80, 24)), index=[f"S{i}" for i in range(80)]
    )
    space = continuous_feature_space_from_columns(data.columns)
    tree = PosetTree()
    tree.add_node("root", is_leaf=False)
    for name, rows in (("left", range(n_left)), ("right", range(n_left, 80))):
        tree.add_node(name, is_leaf=False)
        tree.add_edge("root", name, branch_length=0)
        for row in rows:
            tree.add_node(data.index[row], is_leaf=True, label=data.index[row])
            tree.add_edge(name, data.index[row], branch_length=0)
    populated, _ = full_space_sibling_records(tree, data, space, "production_guarded")
    meta = root_covariance_diagnostics(populated, "root", space, "production_guarded")
    assert meta["covariance_source"] == source
    assert meta["covariance_representation"] == representation
    actual = next(
        iter(covariance_for_node(populated, "root", space, "production_guarded").values())
    )
    if source == "parent_total_fallback":
        expected = np.cov(data.values, rowvar=False)
    else:
        groups = (data.values[:n_left], data.values[n_left:])
        expected = (
            sum(
                (values - values.mean(axis=0)).T @ (values - values.mean(axis=0))
                for values in groups
            )
            / 78
        )
        if representation == "diagonal":
            expected = np.diag(expected)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
