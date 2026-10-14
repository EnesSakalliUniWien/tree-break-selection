from __future__ import annotations

import numpy as np
from tree_break_selection.hierarchy_analysis.statistics.distributional_action import (
    edge_distributional_action_summary,
    split_distributional_action_summary,
)


def test_binary_split_action_matches_child_parent_variance_identity() -> None:
    left_mean = np.array([0.0, 2.0])
    right_mean = np.array([4.0, 2.0])
    left_count = 2
    right_count = 6

    action = split_distributional_action_summary(
        child_means=[left_mean, right_mean],
        child_masses=[left_count, right_count],
    ).action

    parent_mean = (left_count * left_mean + right_count * right_mean) / (left_count + right_count)
    child_parent_action = left_count * np.sum((left_mean - parent_mean) ** 2)
    child_parent_action += right_count * np.sum((right_mean - parent_mean) ** 2)
    pairwise_action = (left_count * right_count / (left_count + right_count)) * np.sum(
        (left_mean - right_mean) ** 2
    )

    assert np.isclose(action, child_parent_action)
    assert np.isclose(action, pairwise_action)


def test_binary_split_edges_are_mass_tied_to_one_sibling_contrast() -> None:
    summary = split_distributional_action_summary(
        child_means=[np.array([0.0]), np.array([10.0])],
        child_masses=[1, 99],
    )
    left_edge_action, right_edge_action = summary.child_edge_actions

    assert np.allclose(summary.parent_mean, np.array([9.9]))
    assert np.isclose(left_edge_action, 98.01)
    assert np.isclose(right_edge_action, 0.99)
    assert np.isclose(summary.action, 99.0)
    assert np.isclose(left_edge_action / summary.action, 0.99)
    assert np.isclose(right_edge_action / summary.action, 0.01)


def test_split_action_decomposes_total_inertia_recursively_across_topology() -> None:
    leaf_values = np.array([[0.0], [2.0], [10.0], [14.0]])
    left = split_distributional_action_summary(
        child_means=[leaf_values[0], leaf_values[1]],
        child_masses=[1, 1],
    )
    right = split_distributional_action_summary(
        child_means=[leaf_values[2], leaf_values[3]],
        child_masses=[1, 1],
    )
    root = split_distributional_action_summary(
        child_means=[left.parent_mean, right.parent_mean],
        child_masses=[left.parent_mass, right.parent_mass],
    )

    root_mean = leaf_values.mean(axis=0)
    total_inertia = float(np.sum((leaf_values - root_mean) ** 2))

    assert np.isclose(left.action, 2.0)
    assert np.isclose(right.action, 8.0)
    assert np.isclose(root.action, 121.0)
    assert np.isclose(total_inertia, left.action + right.action + root.action)


def test_edge_action_changes_when_distribution_moves_with_same_branch_length() -> None:
    parent_mean = np.array([0.0, 0.0])

    near_action = edge_distributional_action_summary(
        parent_mean, np.array([1.0, 0.0]), parent_leaf_count=10, child_leaf_count=5
    ).action
    far_action = edge_distributional_action_summary(
        parent_mean, np.array([2.0, 0.0]), parent_leaf_count=10, child_leaf_count=5
    ).action

    assert far_action == 4.0 * near_action


def test_edge_action_summary_exposes_distribution_masses() -> None:
    summary = edge_distributional_action_summary(
        parent_mean=np.array([0.0, 0.0]),
        child_mean=np.array([3.0, 4.0]),
        parent_leaf_count=20,
        child_leaf_count=5,
    )

    assert summary.parent_mass == 20.0
    assert summary.child_mass == 5.0
    assert summary.child_parent_mass_fraction == 0.25
    assert summary.squared_displacement == 25.0
    assert summary.action == 125.0


def test_split_summary_exposes_child_and_parent_masses() -> None:
    summary = split_distributional_action_summary(
        child_means=[np.array([0.0, 0.0]), np.array([4.0, 0.0])],
        child_masses=[2, 6],
    )

    assert summary.child_masses == (2.0, 6.0)
    assert summary.parent_mass == 8.0
    assert summary.child_parent_mass_fractions == (0.25, 0.75)
    assert summary.squared_child_parent_displacements == (9.0, 1.0)
    assert np.isclose(summary.action, 24.0)
