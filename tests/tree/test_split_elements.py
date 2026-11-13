import networkx as nx
import numpy as np
import pytest

from tree_break_selection.tree.construction import tree_from_linkage
from tree_break_selection.tree.elements.partition import Partition
from tree_break_selection.tree.elements.partition_set import PartitionSet
from tree_break_selection.tree.elements.split_semantics import (
    same_rooted_split,
    same_unrooted_split_within_pivot,
    unrooted_difference_within_pivot,
)
from tree_break_selection.tree.spr import subtree_prune_regraft


def test_same_rooted_split_requires_exact_bitmask():
    encoding = {"A": 0, "B": 1, "C": 2, "D": 3}

    assert same_rooted_split(Partition((0, 1), encoding), Partition((0, 1), encoding))
    assert not same_rooted_split(Partition((0, 1), encoding), Partition((2, 3), encoding))


def test_same_unrooted_split_uses_pivot_scope():
    encoding = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    pivot = Partition((0, 1, 2, 3), encoding)

    assert same_unrooted_split_within_pivot(
        Partition((0, 1), encoding),
        Partition((2, 3), encoding),
        pivot,
        max_complement_ratio=2,
    )


def test_unrooted_difference_preserves_unbalanced_rooted_clades():
    encoding = {name: index for index, name in enumerate("ABCDEFGH")}
    pivot = Partition(tuple(range(8)), encoding)
    large = Partition((0, 1, 2, 3, 4, 5), encoding)
    small = Partition((6, 7), encoding)

    difference = unrooted_difference_within_pivot(
        PartitionSet({large}, encoding=encoding),
        PartitionSet({small}, encoding=encoding),
        pivot,
        max_complement_ratio=2,
    )

    assert large in difference


def test_partition_membership_rejects_matching_bitmask_with_different_encoding():
    partition_set = PartitionSet({Partition((0,), {"A": 0})})

    with pytest.raises(ValueError, match="different encoding"):
        Partition((0,), {"B": 0}) in partition_set


def test_partition_membership_accepts_matching_bitmask_with_same_encoding():
    encoding = {"A": 0}
    partition_set = PartitionSet({Partition((0,), encoding)})

    assert Partition((0,), {"A": 0}) in partition_set


def test_partition_membership_miss_skips_encoding_validation():
    partition_set = PartitionSet({Partition((0,), {"A": 0, "B": 1})})

    assert Partition((1,), {"X": 0, "Y": 1}) not in partition_set


def _linkage_tree(linkage):
    return tree_from_linkage(np.array(linkage, dtype=float), leaf_names=["a", "b", "c", "d"])


def test_poset_tree_nodes_are_identified_by_their_leaf_sets():
    # Same clades {a, b} and {c, d}, merged in a different order: the node ids differ.
    first = _linkage_tree([[0, 1, 1.0, 2], [2, 3, 2.0, 2], [4, 5, 3.0, 4]])
    second = _linkage_tree([[2, 3, 1.0, 2], [0, 1, 2.0, 2], [4, 5, 3.0, 4]])
    encoding = {"a": 0, "b": 1, "c": 2, "d": 3}
    first.initialize_split_indices(encoding)
    second.initialize_split_indices(encoding)

    ab = Partition((0, 1), encoding)
    assert first.nodes["L2"]["split_indices"] == Partition((2,), encoding)
    assert first.nodes[first.root()]["split_indices"].taxa == frozenset("abcd")
    assert first.find_node_by_split(ab) == "N4"
    assert second.find_node_by_split(ab) == "N5"
    assert first.to_splits() == second.to_splits()
    assert len(first.to_splits()) == 3 and len(first.to_splits(with_leaves=True)) == 7


def test_poset_tree_reports_a_leaf_set_that_another_tree_lacks():
    balanced = _linkage_tree([[0, 1, 1.0, 2], [2, 3, 2.0, 2], [4, 5, 3.0, 4]])
    caterpillar = _linkage_tree([[0, 1, 1.0, 2], [4, 2, 2.0, 3], [5, 3, 3.0, 4]])

    missing = balanced.to_splits() - caterpillar.to_splits()
    assert [split.taxa for split in missing] == [frozenset("cd")]
    assert caterpillar.find_node_by_split((2, 3)) is None
    assert caterpillar.find_node_by_split((0, 1, 2)) == "N5"


def test_poset_tree_rejects_a_split_with_another_encoding():
    tree = _linkage_tree([[0, 1, 1.0, 2], [2, 3, 2.0, 2], [4, 5, 3.0, 4]])
    with pytest.raises(ValueError, match="different encoding"):
        tree.find_node_by_split(Partition((0, 1), {"a": 3, "b": 2, "c": 1, "d": 0}))


def test_spr_changes_only_the_leaf_sets_on_the_path():
    # ((a, b), (c, (d, e))): move a onto the edge above d -> (b, (c, ((a, d), e))).
    tree = tree_from_linkage(
        np.array([[0, 1, 1.0, 2], [3, 4, 2.0, 2], [2, 6, 3.0, 3], [5, 7, 4.0, 5]]),
        leaf_names=list("abcde"),
    )
    before = tree.to_splits()
    assert tree.find_node_by_split((0, 1)) == "N5"

    subtree_prune_regraft(tree, "L0", "L3")

    after = tree.to_splits()
    assert {split.taxa for split in before - after} == {
        frozenset("ab"),
        frozenset("de"),
        frozenset("cde"),
    }
    assert {split.taxa for split in after - before} == {
        frozenset("ad"),
        frozenset("ade"),
        frozenset("acde"),
    }
    # The node id N5 stays, its identity does not: it held {a, b} and now holds {a, d}.
    assert tree.find_node_by_split((0, 3)) == "N5"
    assert tree.find_node_by_split((0, 1)) is None
    assert sorted(tree.successors(tree.root())) == ["L1", "N7"]
    assert nx.is_arborescence(tree) and tree.number_of_edges() == 8


def test_spr_rejects_a_regraft_inside_the_pruned_subtree():
    tree = _linkage_tree([[0, 1, 1.0, 2], [2, 3, 2.0, 2], [4, 5, 3.0, 4]])
    with pytest.raises(ValueError, match="onto itself"):
        subtree_prune_regraft(tree, "N4", "L0")
