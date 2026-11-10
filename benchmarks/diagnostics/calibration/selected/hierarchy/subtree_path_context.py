"""Paths and disappearing neighborhoods between exact shared labeled clades."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .subtree_topology_context import Clade, _index_tree
from .subtree_topology_geometry import Tree


@dataclass(frozen=True)
class PathContext:
    start: Clade
    end: Clade
    relation: str
    reference_path: tuple[Clade, ...]
    alternative_path: tuple[Clade, ...]
    reference_lca: Clade
    alternative_lca: Clade
    missing_intermediates: tuple[Clade, ...]
    surviving_intermediates_off_path: tuple[Clade, ...]
    inserted_intermediates: tuple[Clade, ...]
    reference_neighborhood: dict[str, tuple[Clade, ...]]
    alternative_neighborhood: dict[str, tuple[Clade, ...]]
    missing_neighborhood: tuple[Clade, ...]
    surviving_neighborhood_displaced: tuple[Clade, ...]


def _ancestors(parents):
    result = {}
    for node in parents:
        chain = []
        current = node
        while current is not None:
            chain.append(current)
            current = parents[current]
        result[node] = tuple(chain)
    return result


def _path(start, end, ancestors):
    left, right = ancestors[start], ancestors[end]
    right_set = set(right)
    i = next(i for i, node in enumerate(left) if node in right_set)
    lca = left[i]
    j = right.index(lca)
    return left[: i + 1] + right[:j][::-1], lca


def _neighborhood(path, parents, children):
    path_set = set(path)
    roles = {
        "parents": {parents[node] for node in path} - {None},
        "children": {child for node in path for child in children.get(node, ())},
        "siblings": {
            child for node in path for child in children.get(parents[node], ()) if child != node
        },
    }
    return {role: tuple(sorted(nodes - path_set)) for role, nodes in roles.items()}


def compare_surviving_paths(reference: Tree, alternative: Tree):
    """Enumerate each unordered pair of shared clades, including singletons.

    Ancestor paths are oriented toward the descendant. Cross-branch paths use
    sorted endpoint order. Neighborhoods are parents, children and siblings of
    path vertices, excluding the path itself; role sets can overlap. A surviving
    vertex moved off the alternative path is distinct from a missing clade.
    """
    root, parents, children = _index_tree(reference)
    other_root, other_parents, other_children = _index_tree(alternative)
    if root != other_root:
        raise ValueError("Reference and alternative must have the same labeled leaf set")
    ancestors, other_ancestors = _ancestors(parents), _ancestors(other_parents)
    common = sorted(parents.keys() & other_parents.keys())
    for start, end in combinations(common, 2):
        if start in ancestors[end]:
            relation = "ancestor_descendant"
        elif end in ancestors[start]:
            start, end = end, start
            relation = "ancestor_descendant"
        else:
            relation = "cross_branch"
        path, lca = _path(start, end, ancestors)
        other_path, other_lca = _path(start, end, other_ancestors)
        original_neighbors = _neighborhood(path, parents, children)
        new_neighbors = _neighborhood(other_path, other_parents, other_children)
        old_set = set().union(*original_neighbors.values())
        new_set = set().union(*new_neighbors.values())
        yield PathContext(
            start=start,
            end=end,
            relation=relation,
            reference_path=path,
            alternative_path=other_path,
            reference_lca=lca,
            alternative_lca=other_lca,
            missing_intermediates=tuple(n for n in path[1:-1] if n not in other_parents),
            surviving_intermediates_off_path=tuple(
                n for n in path[1:-1] if n in other_parents and n not in other_path
            ),
            inserted_intermediates=tuple(n for n in other_path[1:-1] if n not in path),
            reference_neighborhood=original_neighbors,
            alternative_neighborhood=new_neighbors,
            missing_neighborhood=tuple(sorted(old_set - other_parents.keys())),
            surviving_neighborhood_displaced=tuple(
                sorted((old_set & other_parents.keys()) - new_set)
            ),
        )
