---
title: PosetTree
type: entity
status: reviewed
updated: 2026-11-11
sources:
  - tree_break_selection/tree/README.md
  - tree_break_selection/tree/poset_tree.py
  - tree_break_selection/tree/construction/build.py
  - tree_break_selection/tree/elements/partition.py
  - tree_break_selection/tree/spr.py
tags:
  - entity
  - tree
---

# PosetTree

## Summary

`PosetTree` is the central directed tree structure for the Tree-Break Selection pipeline. It
stores parent-child topology, node metadata, descendant sets, distributions,
and helpers for decomposition and sample cluster assignments.

## Details

The tree package README describes `PosetTree` as a NetworkX `DiGraph`
subclass. Construction is separate: the explicit construction interface builds
linkage topologies and returns a `PosetTree` plus
construction evidence. `PosetTree` itself owns rooted hierarchy operations,
distribution population, leaf counts, and decomposition access; it no longer
contains pass-through constructor aliases.

Node identity is set based since 2026-11-11. The node ids `L<i>` and `N<j>`
are labels that follow the merge order, so the same cells get another id in
another tree. `initialize_split_indices` gives every node its leaf set as
`split_indices`, a `Partition` (one bit per leaf under a shared encoding)
taken over from BranchArchitect together with `PartitionSet` and the rooted
and unrooted split comparison. `to_splits` returns the leaf sets of a tree,
`find_node_by_split` the node holding a leaf set, or `None`. The four node
cases of [[subtree-descriptions-and-statistical-tests]] are lookups of leaf
sets in another tree's `to_splits`. `subtree_prune_regraft` moves a subtree
in place; only the nodes on the path change their leaf set. The split
indices are computed on request; one bitmask of 48,879 cells is about 6 KB.

`PosetTree` is the implementation entity that carries candidate hierarchy
structure into [[tree-decomposition]].

## Evidence

- `tree_break_selection/tree/README.md` maps the public methods and related
  topology and distribution helpers.
- `tree_break_selection/tree/poset_tree.py` contains the class
  implementation.

## Links

- [[tree-break-selection]]
- [[tree-decomposition]]
- [[top-down-traversal]]
- [[spr-alternatives-best-flat-cut]]
- [[spr-refit]]

## Open Questions

- Should wiki coverage split distribution population and topology helpers into
  separate entity pages as the wiki grows?
