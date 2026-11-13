---
title: Prototype of k with a Margin at the Edge, also for More Features than Cells
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/k-margin-prototype-20261112.md
  - raw/inbox/choir-cell-lines/k-margin-prototype-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/k-margin-choir-20261112.md
  - raw/inbox/choir-cell-lines/k-margin-choir-scripts-20261112.py.txt
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/split_count_rule_out.py
tags:
  - invariant-directions
  - marchenko-pastur
  - calibration
  - benchmark
  - prototype
---

# Prototype of k with a Margin at the Edge, also for More Features than Cells

## Summary

A prototype rule for k, outside production: eigenvalue j divided by the
median of the remaining eigenvalues is compared with the simulated 99% point
of the same ratio for one spherical Gaussian population with the same cells
and features. For an independent sample its level is about 1% instead of 16%
to 23%, and it finds directions with more features than cells. Plugged into
the opt-in `sibling_best_flat_cut_reference`, the continuous suite goes from
mean ARI 0.6950 to 0.8880 while the null cases stay at 1 cluster. In subtrees
chosen by the tree it still finds directions in 40% to 75% of nodes. On the
CHOIR cells it fails: ARI 0.8219 to 0.1285, because the cells of one cell
line are far from spherical.

## Key Points

- Level, one Gaussian population, 2,000 draws with k >= 1: prototype 10 to
  20 at seven sizes including (30, 40) and (40, 2000); current rule 291 to
  472 where it applies and 0 with more features than cells.
- Power with more features than cells: 15 and 15 cells in 40 features, `rho`
  5: 972 of 1,000; 20 and 20 cells in 2,000 features, `rho` 50: 1,000 of
  1,000; current rule 0.
- Cost of the margin with more cells than features (64 and 256 cells, 30
  features): `rho` 0.6 gives 568 of 1,000 against 901; `rho` 1.0 gives 997.
- Four balanced groups: k = 3 in 990 of 1,000; two groups: k >= 2 in 14 of
  1,000 against 57 for the current rule.
- Opt-in with the prototype rule, 14 cases: gauss_clear_medium 3 to 4
  clusters, gauss_moderate_3c 2 to 3, gauss_dense_signal_highd 1 to 4
  (default tree) and 1 to 3 (Markov/Ward); every other case unchanged; no
  case above its true groups.
- Continuous suite (9 cases of 11 November), mean ARI: 0.6950 to 0.8880
  (default tree) and 0.8544 (Markov/Ward); today's default test 0.8802 and
  0.8525. `cancer_pc` unchanged at 0.7221.
- `no_invariant_directions`: 58 to 44 of 151 calls. One-population grid: 1
  cluster in 80 of 80 runs; smallest reference p-value 0.0185 against 1.0
  with the current rule.
- Selected subtrees of one Gaussian population, nodes with k >= 1:
  symmetric/average 67.4% to 40.2%, symmetric/ward 94.8% to 72.9%,
  markov/average 83.9% to 56.2%, markov/ward 93.5% to 75.0%.
- Remaining benchmark losses are refusals by the reference (dim_diffuse,
  connected `cancer_pc`), not missing directions.
- CHOIR cells, Markov/Ward tree, same 1,838 calls: 230 clusters to 64, ARI
  0.8219 to 0.1285, lines recovered 116 to 36, accepted by BH 134 to 9. k is
  higher at 1,815 nodes (1 to 28 against 3 to 15); 417 calls have k above 15
  and are refused, among them all 103 nodes with at least 1,000 cells.
- Cause: within a single cell line (175 lines with at least 100 cells) the
  largest eigenvalue is 7.0 to 94.1 times the median eigenvalue, against a
  spherical 99% point of 1.30 to 2.93. k of one line's cells: 5 to 12 with
  the current rule, 6 to 27 with the prototype rule; no line has k = 0.
- So a count of eigenvalues above a spherical edge does not separate one
  population from several on real cells, under either rule; the synthetic
  benchmark has spherical noise by construction and hid this.
- Limits: the null is one spherical Gaussian population; 2,000 simulated
  draws give a realized level of 0.5% to 1.75%; selection by the tree is not
  addressed; one real dataset and tree; no alternative definition of k
  tried; the harness ran once before the rule was checked and was rerun
  with equal results.

## Evidence

`raw/inbox/choir-cell-lines/k-margin-prototype-20261112.md` holds the rule,
the tables of level, power and count, the independent spot check, the
integration with the opt-in, the selected subtrees, the four predictions, the
deviations and the limits; the protocol, the rule, the scripts and their
outputs are in
`raw/inbox/choir-cell-lines/k-margin-prototype-scripts-20261112.py.txt`.
The run on the CHOIR cells and its cause are in
`raw/inbox/choir-cell-lines/k-margin-choir-20261112.md` with protocol,
scripts and outputs in
`raw/inbox/choir-cell-lines/k-margin-choir-scripts-20261112.py.txt`.

## Links

- [[line-non-sphericity]]
- [[non-heuristical-statistic-plan]]
- [[marchenko-pastur-multiplicity-literature]]
- [[best-flat-cut-reference-benchmark]]
- [[gaussian-best-flat-cut-reference]]
- [[selected-tree-null-over-splitting]]
- [[subtree-spectrum-topology]]
