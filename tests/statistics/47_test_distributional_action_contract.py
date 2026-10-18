from __future__ import annotations

import numpy as np

from tree_break_selection.hierarchy_analysis.statistics.distributional_action import (
    edge_distributional_action_summary,
)


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
