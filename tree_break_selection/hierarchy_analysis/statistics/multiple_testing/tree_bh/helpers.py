"""Helper utilities for the Tree-BH hierarchical multiple-testing procedure."""

from __future__ import annotations

from collections import defaultdict

import networkx as nx


def find_root_nodes(tree: nx.DiGraph) -> list[str]:
    """Return root nodes of the directed tree."""
    return [str(node_id) for node_id in tree.nodes() if tree.in_degree(node_id) == 0]


def group_child_indices_by_parent(
    tree: nx.DiGraph,
    child_ids: list[str],
) -> dict[str, list[int]]:
    """Map each parent node to indices of its child hypotheses."""
    missing_child_ids = [child_id for child_id in child_ids if child_id not in tree]
    if missing_child_ids:
        raise ValueError(
            f"Tree-BH child_ids must all be present in the tree; missing {missing_child_ids}."
        )
    sibling_group_indices: dict[str, list[int]] = defaultdict(list)

    for index, child_id in enumerate(child_ids):
        predecessors = list(tree.predecessors(child_id))
        if len(predecessors) != 1:
            raise ValueError(
                f"Tree-BH child {child_id!r} must have exactly one parent; "
                f"found {len(predecessors)}."
            )
        parent_id = str(predecessors[0])
        sibling_group_indices[parent_id].append(index)

    return dict(sibling_group_indices)


__all__ = [
    "find_root_nodes",
    "group_child_indices_by_parent",
]
