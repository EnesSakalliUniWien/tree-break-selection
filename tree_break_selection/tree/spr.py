"""Subtree prune and regraft (spr) on a rooted binary ``PosetTree``."""

from __future__ import annotations

import networkx as nx

from tree_break_selection.tree.poset_tree import PosetTree
from tree_break_selection.tree.topology import tree_root


def _parent(tree: PosetTree, node: object) -> object:
    (parent,) = tree.predecessors(node)
    return parent


def _join(tree: PosetTree, parent: object, child: object, *lengths: float | None) -> None:
    if any(length is None for length in lengths):
        tree.add_edge(parent, child)
    else:
        tree.add_edge(parent, child, branch_length=sum(lengths))


def subtree_prune_regraft(tree: PosetTree, prune: object, regraft: object) -> None:
    """Move the subtree at ``prune`` onto the edge above ``regraft``, in place.

    The parent of ``prune`` keeps its node id and becomes the parent of ``prune`` and
    ``regraft``; the old sibling of ``prune`` takes the parent's place. Only nodes on the path
    between the two places change their leaf set; every other node keeps its
    ``split_indices``. Node distributions and annotations are stale afterwards.
    """
    root = tree_root(tree)
    if prune == root or regraft == root:
        raise ValueError("Cannot prune the root or regraft above it.")
    parent = _parent(tree, prune)
    if regraft in (prune, parent) or regraft in nx.descendants(tree, prune):
        raise ValueError("Cannot regraft a subtree onto itself or onto its own parent edge.")
    children = list(tree.successors(parent))
    if len(children) != 2:
        raise ValueError("spr needs a binary tree.")
    sibling = children[0] if children[1] == prune else children[1]
    if regraft == sibling:
        return

    sibling_length = tree.edges[parent, sibling].get("branch_length")
    tree.remove_edge(parent, sibling)
    if parent == root:
        tree.graph["root"] = sibling
    else:
        grandparent = _parent(tree, parent)
        parent_length = tree.edges[grandparent, parent].get("branch_length")
        tree.remove_edge(grandparent, parent)
        _join(tree, grandparent, sibling, parent_length, sibling_length)

    target = _parent(tree, regraft)
    length = tree.edges[target, regraft].get("branch_length")
    half = None if length is None else length / 2
    tree.remove_edge(target, regraft)
    # shortcut: the regraft edge is halved; refit branch lengths before testing edges.
    _join(tree, target, parent, half)
    _join(tree, parent, regraft, half)

    if "split_index" in tree.graph:
        tree.initialize_split_indices(tree.graph["taxa_encoding"])


__all__ = ["subtree_prune_regraft"]
