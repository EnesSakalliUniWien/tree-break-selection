---
title: Split-Count Rule-Out on the Continuous Benchmark Suite
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/rule-out-continuous-benchmark-20261108.md
  - raw/inbox/choir-cell-lines/rule-out-continuous-benchmark-scripts-20261108.py.txt
tags:
  - benchmark
  - calibration
  - sibling-gate
---

# Split-Count Rule-Out on the Continuous Benchmark Suite

## Summary

On the 9 continuous benchmark cases the opt-in rule-out removes over-splitting
where today's test over-splits but loses real groups wherever the subtree to be
split has no more cells than features plus one; mean ARI drops from 0.85-0.88
(today) to 0.70.

## Key Points

- Only 9 of 126 benchmark cases are continuous; the rule-out does not apply to the
  binary and categorical cases.
- Gains: dim_consolidated 4 groups, 72 features: 5-6 clusters to 4, ARI 0.82-0.91 to
  1.000; 24 features on Markov/Ward 0.920 to 1.000.
- Losses: 60 cells in 40 features (4 to 3 clusters), 45 in 30 (3 to 2), 40 in 20,000
  (4 to 1), 180 in 136 (already failing today, now 1 cluster).
- Null and outlier cases unchanged; both trees give identical rule-out results.
- The module now measures equidistant cells only for separable splits, so wide data
  need no p x p matrix for refused splits.

## Evidence

`raw/inbox/choir-cell-lines/rule-out-continuous-benchmark-20261108.md` holds the setup
and per-case table; protocol, runner and result rows are in
`raw/inbox/choir-cell-lines/rule-out-continuous-benchmark-scripts-20261108.py.txt`.

## Links

- [[split-count-exclusion]]
- [[split-count-exclusion-real]]
