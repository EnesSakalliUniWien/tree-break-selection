---
title: The Gradient Behind the Split-Count Rule-Out
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/rule-out-gradient-20261108.md
  - raw/inbox/choir-cell-lines/rule-out-gradient-scripts-20261108.py.txt
tags:
  - calibration
  - spectral
  - sibling-gate
  - power
---

# The Gradient Behind the Split-Count Rule-Out

## Summary

The rule-out's acceptance of a real split is governed by three terms: the group
direction must rise above the noise band (needs separation that shrinks as cells per
feature grow), the split must beat the count of all splits (fixed separation, any
size), and the two sides must not overlap at all (needs separation that grows with
the number of cells). The third term binds: the required separation is about
2 Phi^-1(1 - 1/n) noise SDs, 4.7 at 100 cells and 8.2 at 48,879.

## Key Points

- Benchmark gradient: recovery rises with cells of a two-group subtree per feature
  (failures all at <= p + 1 cells, the implementation edge).
- Invariant direction found in 28% to 99% of runs as the transition variable
  log(pi(1-pi) delta^2) - 0.5 log(p/n) goes from below 0 to above 1.
- Count of all splits: pass 0% at pi(1-pi) delta^2 = 1, 66% at 2.25, 95% at 4.
- Separable share follows Phi(delta/2)^n (p = 10, delta = 4: 63% at 20 cells, 5% at
  160); kept falls with n at fixed separation.
- Logistic fit: log(n/p) adds pseudo-R2 0.43 to 0.74 with a negative sign.
- The zero-overlap demand is also what removes noise splits.

## Evidence

`raw/inbox/choir-cell-lines/rule-out-gradient-20261108.md` holds the grid,
relationships and reading; protocol, deviation, code and outputs are in
`raw/inbox/choir-cell-lines/rule-out-gradient-scripts-20261108.py.txt`.

## Links

- [[rule-out-continuous-benchmark]]
- [[split-count-exclusion]]
- [[split-count-exclusion-real]]
