"""Strict topology helpers for rooted ``PosetTree`` instances."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Iterator

import networkx as nx


def tree_root(tree: nx.DiGraph) -> object:
    """Return the cached ``graph["root"]`` or the single node without a parent."""
    if "root" in tree.graph:
        return tree.graph["root"]
    roots = [node for node, degree in tree.in_degree() if degree == 0]
    if len(roots) != 1:
        raise ValueError(f"Expected one root, got {roots!r}.")
    return roots[0]


def bottom_up_nodes(tree: nx.DiGraph) -> Iterator[object]:
    """Yield nodes in bottom-up order (leaves first, root last).

    This is the single canonical traversal for all bottom-up aggregation
    passes across the codebase — descendant sets, distribution propagation,
    split-flag precomputation, and spectral index collection should all
    call this instead of inlining their own topological sort.

    Equivalent to ``reversed(list(nx.topological_sort(tree)))`` and to
    ``nx.dfs_postorder_nodes(tree, source=root)`` on a rooted tree, but
    expressed as a single source of truth.

    Parameters
    ----------
    tree
        Directed acyclic graph representing the hierarchy.

    Yields
    ------
    object
        Node identifiers, leaves first, root last.
    """
    ordered = nx.lexicographical_topological_sort(
        tree,
        key=lambda node: (type(node).__qualname__, repr(node)),
    )
    return reversed(list(ordered))


def compute_node_depths(tree: nx.DiGraph) -> dict[object, int]:
    """Compute depth of each node from the root via BFS.

    Parameters
    ----------
    tree
        Directed acyclic graph representing the hierarchy.

    Returns
    -------
    Dict[str, int]
        Mapping from node_id to depth (root = 0).

    Raises
    ------
    ValueError
        If the tree has no root node (all nodes have parents).
    """
    roots = sorted(
        (node for node, degree in tree.in_degree() if degree == 0),
        key=lambda node: (type(node).__qualname__, repr(node)),
    )
    if not roots:
        raise ValueError("Tree has no root node (all nodes have parents)")

    depths: dict[object, int] = {root: 0 for root in roots}
    queue = deque(roots)
    while queue:
        node = queue.popleft()
        for child in tree.successors(node):
            if child not in depths:
                depths[child] = depths[node] + 1
                queue.append(child)

    return depths


def is_leaf(tree: nx.DiGraph, node_id: object) -> bool:
    """Return the explicit leaf flag for a node."""
    return bool(tree.nodes[node_id]["is_leaf"])


def get_leaf_nodes(tree: nx.DiGraph, node: object | None = None) -> list[object]:
    """Collect leaf node ids globally or inside one subtree."""
    if node is None:
        candidates = list(tree.nodes)
    elif is_leaf(tree, node):
        candidates = [node]
    else:
        candidates = list(nx.descendants(tree, node))
    return [candidate for candidate in candidates if is_leaf(tree, candidate)]


def get_leaf_values(
    tree: nx.DiGraph,
    *,
    node: object | None = None,
    use_labels: bool = True,
    sort: bool = True,
) -> list[object]:
    """Collect leaf labels or ids globally or inside one subtree."""
    leaf_nodes = get_leaf_nodes(tree, node)
    values = [
        tree.nodes[leaf_node]["label"] if use_labels else leaf_node for leaf_node in leaf_nodes
    ]
    return sorted(values) if sort else values


def compute_descendant_leaf_sets(
    tree: nx.DiGraph,
    *,
    use_labels: bool = True,
) -> dict[object, frozenset]:
    """Map each node to the set of descendant leaf labels or ids."""
    descendant_sets: dict[object, frozenset] = {}
    for node in bottom_up_nodes(tree):
        if is_leaf(tree, node):
            leaf_value = tree.nodes[node]["label"] if use_labels else node
            descendant_sets[node] = frozenset([leaf_value])
            continue

        child_sets = [descendant_sets[child] for child in tree.successors(node)]
        descendant_sets[node] = frozenset.union(*child_sets)
    return descendant_sets


def _single_parent(tree: nx.DiGraph, node_id: object) -> object:
    parents = list(tree.predecessors(node_id))
    if len(parents) != 1:
        raise ValueError(f"Expected node {node_id!r} to have exactly one parent.")
    return parents[0]


def lowest_common_ancestor(
    tree: nx.DiGraph,
    node_a: object,
    node_b: object,
    *,
    depths: dict[object, int] | None = None,
) -> object:
    """Find the lowest common ancestor of two nodes in a rooted tree."""
    if node_a == node_b:
        return node_a

    node_depths = depths if depths is not None else compute_node_depths(tree)
    try:
        depth_a = node_depths[node_a]
        depth_b = node_depths[node_b]
    except KeyError as exc:
        raise ValueError(
            f"Cannot compute LCA for nodes outside the rooted tree: {node_a!r}, {node_b!r}."
        ) from exc

    current_a, current_b = node_a, node_b
    for _ in range(max(depth_a - depth_b, 0)):
        current_a = _single_parent(tree, current_a)
    for _ in range(max(depth_b - depth_a, 0)):
        current_b = _single_parent(tree, current_b)

    while current_a != current_b:
        current_a = _single_parent(tree, current_a)
        current_b = _single_parent(tree, current_b)

    return current_a


def lowest_common_ancestor_for_set(
    tree: nx.DiGraph,
    nodes: Iterable[object],
) -> object:
    """Find the lowest common ancestor for a non-empty node collection."""
    node_iterator = iter(nodes)
    try:
        lca = next(node_iterator)
    except StopIteration as exc:
        raise ValueError("Cannot find LCA for an empty node set.") from exc

    root = tree.graph["root"]
    depths = compute_node_depths(tree)
    for node in node_iterator:
        lca = lowest_common_ancestor(tree, lca, node, depths=depths)
        if lca == root:
            return root
    return lca


__all__ = [
    "bottom_up_nodes",
    "compute_descendant_leaf_sets",
    "compute_node_depths",
    "get_leaf_nodes",
    "get_leaf_values",
    "is_leaf",
    "lowest_common_ancestor",
    "lowest_common_ancestor_for_set",
    "tree_root",
]
