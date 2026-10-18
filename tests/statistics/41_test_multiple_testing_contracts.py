from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from tree_break_selection.hierarchy_analysis.statistics.multiple_testing import (
    benjamini_hochberg_correction,
)
from tree_break_selection.hierarchy_analysis.statistics.multiple_testing.tree_bh import (
    apply_tree_bh_correction,
)


def _two_child_tree() -> nx.DiGraph:
    tree = nx.DiGraph()
    tree.add_edge("root", "A")
    tree.add_edge("root", "B")
    return tree


def test_bh_correction_rejects_nonfinite_p_values() -> None:
    with pytest.raises(ValueError, match="finite before correction"):
        benjamini_hochberg_correction(np.array([0.01, np.nan]), alpha=0.05)


def test_bh_correction_rejects_out_of_range_p_values() -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        benjamini_hochberg_correction(np.array([0.01, 1.2]), alpha=0.05)


def test_tree_bh_requires_p_values_to_align_with_child_ids() -> None:
    with pytest.raises(ValueError, match="align one-to-one with child_ids"):
        apply_tree_bh_correction(
            _two_child_tree(),
            np.array([0.01]),
            ["A", "B"],
            alpha=0.05,
        )


def test_tree_bh_rejects_nonfinite_p_values() -> None:
    with pytest.raises(ValueError, match="finite before correction"):
        apply_tree_bh_correction(
            _two_child_tree(),
            np.array([0.01, np.nan]),
            ["A", "B"],
            alpha=0.05,
        )


def test_tree_bh_rejects_invalid_alpha() -> None:
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        apply_tree_bh_correction(
            _two_child_tree(),
            np.array([0.01, 0.02]),
            ["A", "B"],
            alpha=0.0,
        )


def test_tree_bh_requires_child_ids_to_be_non_root_edges() -> None:
    with pytest.raises(ValueError, match="exactly one parent"):
        apply_tree_bh_correction(
            _two_child_tree(),
            np.array([0.01]),
            ["root"],
            alpha=0.05,
        )


def test_tree_bh_requires_child_ids_to_exist_in_tree() -> None:
    with pytest.raises(ValueError, match="present in the tree"):
        apply_tree_bh_correction(
            _two_child_tree(),
            np.array([0.01]),
            ["missing"],
            alpha=0.05,
        )


def test_flat_bh_matches_hand_computed_adjustment():
    rejected, adjusted = benjamini_hochberg_correction(
        np.array([0.001, 0.01, 0.03, 0.05, 0.1]), alpha=0.05
    )

    np.testing.assert_allclose(adjusted, [0.005, 0.025, 0.05, 0.0625, 0.1])
    assert rejected.tolist() == [True, True, True, False, False]


def test_tree_bh_scales_child_alpha_by_parent_family_rejections():
    tree = nx.DiGraph()
    tree.add_edges_from(
        [("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2")]
    )
    child_ids = ["A", "B", "A1", "A2", "B1", "B2"]
    p_values = np.array([0.001, 0.5, 0.012, 0.03, 1e-9, 1e-9])

    result = apply_tree_bh_correction(tree, p_values, child_ids, alpha=0.05)

    # Root family: BH-adjusted [0.002, 0.5] at alpha 0.05 rejects only A.
    # A's family is tested at 0.05 * 1/2; BH-adjusted [0.024, 0.03] rejects only A1.
    # B was not rejected, so its children are never tested despite tiny p-values.
    assert result.child_parent_edge_null_rejected_by_tree_bh.tolist() == [
        True,
        False,
        True,
        False,
        False,
        False,
    ]
    assert result.child_parent_edge_tested_by_tree_bh.tolist() == [
        True,
        True,
        True,
        True,
        False,
        False,
    ]
    np.testing.assert_allclose(
        result.child_parent_edge_corrected_p_values_by_tree_bh,
        [0.002, 0.5, 0.024, 0.03, 1.0, 1.0],
    )
