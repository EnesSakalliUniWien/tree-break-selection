"""Contract for the peel mapping that the edge gate, Tree-BH and the sibling frontier share."""

from __future__ import annotations

import networkx as nx

from tree_break_selection.tree.topology import peel_continuations


def _tree() -> nx.DiGraph:
    # root peels the leaf x; C peels the 2-leaf group S; D splits 3 | 3 leaves.
    return nx.DiGraph(
        [("root", "x"), ("root", "C"), ("C", "S"), ("C", "D"), ("S", "s1"), ("S", "s2")]
        + [("D", "E"), ("D", "F")]
        + [("E", f"e{i}") for i in range(3)]
        + [("F", f"f{i}") for i in range(3)]
    )


def test_without_a_limit_only_singleton_peels_are_mapped() -> None:
    assert peel_continuations(_tree()) == {"root": "C"}


def test_a_limit_adds_nodes_whose_smaller_child_is_below_it() -> None:
    # S has 2 leaves (< 3) and D has 6: C continues in D. D's children both have 3 leaves.
    assert peel_continuations(_tree(), small_peel_leaf_limit=3) == {"root": "C", "C": "D"}


def test_a_node_whose_larger_child_is_below_the_limit_is_not_a_small_peel() -> None:
    # With limit 7 no child of C or D reaches the limit, so only the singleton peel remains.
    assert peel_continuations(_tree(), small_peel_leaf_limit=7) == {"root": "C"}
