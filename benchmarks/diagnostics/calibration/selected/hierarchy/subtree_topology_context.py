"""Compare node and adjacency identities across labeled binary topologies."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .subtree_topology_geometry import Tree

type Clade = tuple[int, ...]
type Children = tuple[Clade, Clade]


@dataclass(frozen=True)
class NodeContext:
    """Context of one reference internal node in the alternative tree."""

    node: Clade
    original_parent: Clade | None
    original_children: Children
    node_exists: bool
    original_parent_exists: bool | None
    original_children_exist: tuple[bool, bool]
    both_children_absent: bool
    at_least_one_child_absent: bool
    alternative_parent: Clade | None
    alternative_children: Children | None
    direct_parent_preserved: bool | None
    direct_children_preserved: bool
    case_1: bool
    case_2: bool
    case_3: bool
    case_4: bool


def _index_tree(tree: Tree) -> tuple[Clade, dict[Clade, Clade | None], dict[Clade, Children]]:
    parents: dict[Clade, Clade | None] = {}
    children: dict[Clade, Children] = {}
    leaves: set[int] = set()

    def visit(node: Tree) -> Clade:
        if isinstance(node, Integral) and not isinstance(node, (bool, np.bool_)):
            leaf = int(node)
            if leaf in leaves:
                raise ValueError("Tree leaf labels must be unique")
            leaves.add(leaf)
            return (leaf,)
        if not isinstance(node, tuple) or len(node) != 2:
            raise ValueError("Expected a binary nested tuple with integer leaf labels")
        left, right = sorted(map(visit, node))
        clade = tuple(sorted(left + right))
        children[clade] = left, right
        parents[left] = parents[right] = clade
        return clade

    root = visit(tree)
    parents[root] = None
    return root, parents, children


def compare_node_contexts(reference: Tree, alternative: Tree) -> list[NodeContext]:
    """Compare exact descendant leaf sets for every reference internal node.

    Cases can overlap. Case 1 means the node is absent. Case 2 means both
    original child clades are absent, independently of node presence. Case 3
    means the node remains and its original parent clade is absent. Case 4
    means case 3 with both original children absent. The root's parent and
    direct-parent fields are not applicable. A clade can remain elsewhere
    in the alternative without preserving its direct adjacency.
    """
    root, parents, children = _index_tree(reference)
    other_root, other_parents, other_children = _index_tree(alternative)
    if root != other_root:
        raise ValueError("Reference and alternative must have the same labeled leaf set")

    rows = []
    for node, original_children in sorted(
        children.items(), key=lambda item: (len(item[0]), item[0])
    ):
        parent = parents[node]
        exists = node in other_parents
        parent_exists = None if parent is None else parent in other_parents
        left_exists = original_children[0] in other_parents
        right_exists = original_children[1] in other_parents
        both_absent = not (left_exists or right_exists)
        alternative_parent = other_parents.get(node)
        alternative_children = other_children.get(node)
        case_3 = exists and parent_exists is False
        rows.append(
            NodeContext(
                node=node,
                original_parent=parent,
                original_children=original_children,
                node_exists=exists,
                original_parent_exists=parent_exists,
                original_children_exist=(left_exists, right_exists),
                both_children_absent=both_absent,
                at_least_one_child_absent=not (left_exists and right_exists),
                alternative_parent=alternative_parent,
                alternative_children=alternative_children,
                direct_parent_preserved=None if parent is None else alternative_parent == parent,
                direct_children_preserved=alternative_children == original_children,
                case_1=not exists,
                case_2=both_absent,
                case_3=case_3,
                case_4=case_3 and both_absent,
            )
        )
    return rows
