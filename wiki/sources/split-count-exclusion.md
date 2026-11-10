---
title: Split Counts with Exclusions from Subspaces and Equidistant Cells
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/split-count-exclusion-20261108.md
  - raw/inbox/choir-cell-lines/split-count-exclusion-scripts-20261108.py.txt
tags:
  - calibration
  - selective-inference
  - sibling-gate
  - tree-choice
  - spectral
---

# Split Counts with Exclusions from Subspaces and Equidistant Cells

## Summary

Counting all possible splits of a subtree, then ruling a split out as a real
separation unless it lies in the subtree's invariant directions and leaves at most
5% of cells equidistant (`rule_out`), is the first arm that meets the frozen rule:
no-group false passes at or below 2% for both large and small splits, with every
simulated real group kept. Shrinking the count instead (`shrink`) makes the test
looser and fails.

## Key Points

- Markov/Ward no-group splits >= 40 cells pass 72-73% today, 0.0-0.9% with
  rule_out; 10-39 cells 1-3% today, 0.7-1.6% with rule_out (count_all alone 85-93%).
- False partitions of 300: Markov/Ward 298/296 to 2/25 (sd 0.25/0.42); default tree
  2/50/259 to 0/0/0.
- Markov/Ward on real groups: ARI 0.03-0.32 today to 0.87-1.00, with 4-5 clusters.
- Default-tree real-group ARI unchanged (density variation 0.651 to 0.667).
- Separability does most of the excluding; under about 40 cells the tie share is
  uninformative (whitened 30-dimensional metric).
- Limits: error control measured, not guaranteed; no split inside a subtree of
  <= 31 cells can be accepted; 25/300 false partitions remain at sd 0.42 on
  Markov/Ward; real CHOIR data not yet run end to end.

## Production opt-in

`sibling_split_count_rule_out` (default off) applies `rule_out` in production:
`tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/split_count_rule_out.py`,
wired through `annotate_sibling_divergence`, `run_gate_annotation_pipeline` (projected-Wald
gate only, continuous features required), `run_tbs_on_distance` and benchmark dispatch.
Each tested split's reason is written to `Sibling_Test_Method`
(`split_count_rule_out:kept`, `:no_invariant_directions`, `:not_separable`,
`:equidistant_cells`). Tests: `tests/statistics/70_test_split_count_rule_out.py`,
`tests/integration/65_test_split_count_rule_out_runner.py`. Measured with 30 test
coordinates only: in a 10-coordinate check one population still gave 6 clusters,
because subtrees above 11 cells keep invariant directions there.

## Evidence

`raw/inbox/choir-cell-lines/split-count-exclusion-20261108.md` holds the arms,
tables, exclusion shares and limits; protocol, code and output are in
`raw/inbox/choir-cell-lines/split-count-exclusion-scripts-20261108.py.txt`.

## Links

- [[split-count-correction]]
- [[subspace-distance-ties]]
- [[subtree-spectrum-topology]]
- [[exact-gate-sibling-test]]
