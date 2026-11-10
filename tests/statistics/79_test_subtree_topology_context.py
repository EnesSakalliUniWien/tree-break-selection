"""Exact labeled witnesses for subtree context changes."""

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    topology_from_linkage,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_context import (
    compare_node_contexts,
)


def node_row(reference, alternative, node):
    return next(row for row in compare_node_contexts(reference, alternative) if row.node == node)


def test_absent_node_and_both_absent_children_can_overlap():
    row = node_row((((0, 1), (2, 3)), (4, 5)), (((0, 2), (1, 4)), (3, 5)), (0, 1, 2, 3))
    assert row.case_1 and row.case_2
    assert not row.case_3 and not row.case_4
    assert not row.node_exists and row.original_parent_exists
    assert row.original_children_exist == (False, False)
    assert row.both_children_absent and row.at_least_one_child_absent
    assert row.alternative_parent is None and row.alternative_children is None
    assert not row.direct_parent_preserved and not row.direct_children_preserved


def test_node_with_same_children_but_original_parent_absent():
    reference = ((((0, 1), (2, 3)), 4), 5)
    alternative = ((((0, 1), (2, 3)), 5), 4)
    row = node_row(reference, alternative, (0, 1, 2, 3))
    assert row.node_exists and row.case_3
    assert not row.case_1 and not row.case_2 and not row.case_4
    assert row.original_parent == (0, 1, 2, 3, 4)
    assert row.alternative_parent == (0, 1, 2, 3, 5)
    assert not row.direct_parent_preserved and row.direct_children_preserved


def test_six_leaf_witness_preserves_node_and_loses_parent_and_both_children():
    reference = ((((0, 1), (2, 3)), 4), 5)
    alternative = ((((0, 2), (1, 3)), 5), 4)
    row = node_row(reference, alternative, (0, 1, 2, 3))
    assert row.node_exists and not row.case_1
    assert row.case_2 and row.case_3 and row.case_4
    assert row.original_children == ((0, 1), (2, 3))
    assert row.original_children_exist == (False, False)
    assert row.alternative_children == ((0, 2), (1, 3))
    assert row.alternative_parent == (0, 1, 2, 3, 5)
    assert not row.direct_parent_preserved and not row.direct_children_preserved


def test_six_leaf_case_four_is_realized_by_distinct_covariance_eigendirections():
    x = [0, 1, 10, 11, 30, 100]
    y = [0, 10, 1, 11, 1000, -12751 / 224]
    data = np.column_stack([x, y])
    covariance = np.cov(data, rowvar=False)
    assert covariance[0, 1] == pytest.approx(0.0, abs=1e-12)
    assert covariance[1, 1] > 100 * covariance[0, 0]
    reference = topology_from_linkage(linkage(data[:, :1], method="ward"))
    alternative = topology_from_linkage(linkage(data[:, 1:], method="ward"))
    row = node_row(reference, alternative, (0, 1, 2, 3))
    assert row.case_4
    assert row.original_children == ((0, 1), (2, 3))
    assert row.alternative_children == ((0, 2), (1, 3))


def test_parent_and_one_child_can_exist_without_the_original_adjacencies():
    reference = (((0, 1), (2, 3)), 4)
    alternative = ((((0, 1), 2), 3), 4)
    row = node_row(reference, alternative, (0, 1))
    assert row.node_exists and row.original_parent_exists
    assert row.alternative_parent == (0, 1, 2) and not row.direct_parent_preserved
    assert row.original_children_exist == (True, True) and row.direct_children_preserved

    row = node_row(reference, alternative, (0, 1, 2, 3))
    assert row.original_children_exist == (True, False)
    assert row.at_least_one_child_absent and not row.both_children_absent
    assert not row.case_2 and not row.direct_children_preserved
    assert row.alternative_children == ((0, 1, 2), (3,))


def test_root_and_child_reversal_equal_tree():
    a = ((0, 1), (2, 3))
    reversed_a = ((3, 2), (1, 0))
    rows = compare_node_contexts(a, reversed_a)
    assert rows == compare_node_contexts(reversed_a, a)
    assert {row.node for row in rows} == {(0, 1), (2, 3), (0, 1, 2, 3)}
    assert all(row.node_exists and row.direct_children_preserved for row in rows)
    assert all(not (row.case_1 or row.case_2 or row.case_3 or row.case_4) for row in rows)
    root = node_row(a, reversed_a, (0, 1, 2, 3))
    assert root.original_parent is root.original_parent_exists is None
    assert root.alternative_parent is root.direct_parent_preserved is None
    assert all(row.direct_parent_preserved for row in rows if row.original_parent is not None)


def test_root_can_lose_children_but_has_no_missing_parent_case():
    root = node_row(((0, 1), (2, 3)), ((0, 2), (1, 3)), (0, 1, 2, 3))
    assert root.case_2 and root.node_exists
    assert not root.case_3 and not root.case_4


@pytest.mark.parametrize(
    "invalid", [(0, 1, 2), ((0, 1), 0), ((0, 1), 2.0), ((False, 1), 2), ((0, 1), np.bool_(True))]
)
def test_invalid_tree_fails_explicitly(invalid):
    with pytest.raises(ValueError):
        compare_node_contexts(((0, 1), 2), invalid)


def test_different_leaf_universe_fails_and_integer_labels_are_normalized():
    with pytest.raises(ValueError, match="same"):
        compare_node_contexts(((0, 1), 2), ((0, 1), 3))
    rows = compare_node_contexts(((np.int64(0), 1), 2), ((0, 1), 2))
    assert all(type(leaf) is int for row in rows for leaf in row.node)
