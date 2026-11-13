---
title: Subtree Descriptions Traced to the Statistical Tests
type: analysis
status: draft
updated: 2026-11-11
sources:
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_context.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_spectral_alignment.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_geometry.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_path_context.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-real-20261111.md
tags:
  - subtree
  - hypotheses
  - production
  - sibling-gate
---

# Subtree Descriptions Traced to the Statistical Tests

## Summary

The four structural cases, the alternative partitions, the distance fit and
the inherited and rebuilt geometries of 9 November are traced here to the
sibling tests, including the opt-in Gaussian best-flat-cut reference of 11
November. One relation is usable: the reference depends only on a node's cells
and its split, so one reference covers every alternative child partition of a
surviving node in the same invariant directions. Three relations are gaps: the
tree's splits use inherited geometry while the reference uses rebuilt
geometry, statistics along a path use different covariances and directions,
and no reference covers the choice of a node's cells.

## Details

### What a sibling test is made of

A sibling test at node C is fixed by its cells, its split into two children
and the directions in which the statistic is measured. Its decision also
depends on reachability from the root and on the BH family. The native
statistic adds the branch factor g, a ridge, an adaptive rank and an empirical
inflation, which depend on the parent and on global branch fits. The reference
computes `rho` from C's cells and split alone, in `invariant_directions` of
C's cells.

### The four cases

Cases compare a node of one tree with an alternative tree
(`compare_node_contexts`); counts are from 16,128 dependent comparisons.

- Case 1, node does not exist (15,316): no test of that node exists in the
  alternative. Most tested nodes exist in one tree only, so their tests exist
  through the tree's choice of cells. No reference built so far covers that
  choice.
- Case 2, both children do not exist (5,271): where C remains, the cells are
  the same and the split differs. Native statistics of the two splits are not
  comparable. Under the reference both splits have the same directions and the
  same null, and both are bounded by the best flat cut.
- Case 3, node exists and its parent does not (401): the incoming edge test,
  the branch factor and reachability change; all 19 comparisons with the same
  parent and children still changed raw T and g. The reference's `rho` and
  p-value at C are unchanged; its BH family and traversal still change.
- Case 4, node exists, parent and children do not (102): only the cells
  survive. Under the reference the null is the same and `rho` differs.

### The three production measurements

- Alternative partitions each subtree preserves: `spectral_subspaces`,
  `spectral_geometry` and `ward_topology_family` give them, exactly up to 6 to
  10 leaves and through a bounded candidate selector up to 256 leaves, at 14.10
  and 24.88 times the baseline runtime. Null datasets also give 63
  alternatives each, so a count of alternatives is not evidence. An
  alternative that is a split of the same node can be tested with `rho` in the
  node's full invariant directions against the same reference.
- Distance fit: `fit_additive_distances` reports RMSE and the largest
  residual for a supplied topology. It measures representation; a better
  distance fit did not repair calibration. It belongs beside the test, not in
  it.
- Inherited parent geometry and rebuilt child geometry: `restrict_parent_basis`
  and `compare_parent_child` compare them; restriction lost rank in 239 of 320
  spans and 316 of 320 induced and refitted operators differ. The tree is built
  once from the root's diffusion distances, so every split is inherited
  geometry. The reference takes its directions from the subtree's own
  covariance in feature space, so it is rebuilt geometry. Each test therefore
  asks about a split chosen in one geometry with directions from another.
  [[split-vs-best-flat-cut]] measures this per split: the median ARI between
  the tree split and the best flat cut is 0.02 to 0.08 on Markov/Ward.

### Where the 11 November results fit

- The impure cancer_pc splits and the 27 merged line pairs whose true-label
  split passes are instances of the inherited and rebuilt mismatch.
- A subtree without substructure has no invariant directions; with the
  covariance identity of parent and children, directions arise from
  between-group means. The refusals at 157 to 256 cells in cancer_pc are this
  case.
- Weak directions matter: PC2 with a median variance share of 0.2885 kept a
  crossing factor in 20 of 20 datasets and leading-rank truncation reduced
  recovery from 20 of 20 to 6 of 20. The reference's limit k up to 10 refuses
  every real-data split with 1,000 cells or more; restricting to the 10 leading
  directions instead would repeat the truncation.
- Along a path every node uses its own pooled covariance and its own
  directions; summed edge scores and the endpoint score differ (42 against
  4.5). Pass-through combines such decisions without a derivation and decided
  the disconnected cancer_pc results.

### Status of the 10 November plan

A1 (bound) was replaced by the simulated Gaussian reference; A3 exists as the
opt-in. A2 (define k), B (separated spaces and topology alternatives), C
(neighborhood with one common Sigma) and D (fixtures and frozen panels) were
not started. The plan's text is not recorded in the wiki.

## Evidence

[[subtree-context-and-eigenvalue-statistics]],
[[subtree-spectral-alignment-and-crossing-clusters]],
[[subtree-paths-and-missing-neighborhoods]] and
[[subtree-clustering-production-gates]] hold the 9 November descriptions and
counts. [[best-flat-cut-reference-benchmark]] and
[[best-flat-cut-reference-real]] hold the 11 November results. The relation of
the reference to cases 2 to 4 follows from its definition in
`best_flat_cut_reference.py`; it was not measured on the saved fits.

## Links

- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[subtree-paths-and-missing-neighborhoods]]
- [[subtree-clustering-production-gates]]
- [[subtree-distance-topology-families]]
- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-benchmark]]
- [[best-flat-cut-reference-real]]
- [[split-vs-best-flat-cut]]
- [[best-flat-cut-divisive]]
- [[spr-alternatives-best-flat-cut]]
- [[spr-subspaces-equidistances]]
- [[spr-eigenvalue-blocks]]
- [[spr-refit]]

## Open Questions

- [[best-flat-cut-divisive]] separates the lines (ARI 0.874, 16 lines in a
  merge) and splits 88 lines internally. Are those internal splits real
  structure or non-Gaussian shape, and which reference covers them?
- Which reference covers the choice of a node's cells?
- Can the search be made reliable above k = 10, so that weak directions are
  kept without refusing large subtrees?
- Does pass-through have a derivation under one common covariance per
  neighborhood?
