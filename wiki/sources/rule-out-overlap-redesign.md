---
title: Redesigning the Zero-Overlap Rule
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/rule-out-overlap-redesign-20261108.md
  - raw/inbox/choir-cell-lines/rule-out-overlap-redesign-scripts-20261108.py.txt
tags:
  - calibration
  - sibling-gate
  - power
  - tree-choice
---

# Redesigning the Zero-Overlap Rule

## Summary

Replacing the rule-out's zero-overlap demand with an allowed overlap fraction (the
overlap of the weakest separation the count of all splits accepts) stops the
tightening with cell count on clean Gaussian groups, but end to end it lets more of
the tree's noise splits through. Tree-chosen noise splits misplace only 10-20% of
cells, about as little as weak real groups, so no such fraction separates them.
Production keeps the zero-overlap rule.

## Key Points

- Grid: with the new overlap and tie rules, kept rises with n at fixed separation
  (p = 10, delta = 3: 25% to 85%); random noise splits kept 0 of 1,080.
- End to end, Markov/Ward false partitions of 300 rise from 2/25/15 (zero overlap) to
  13/42/36 (new overlap) and 17/46/141 (new overlap and ties); density-variation ARI
  falls from 0.87 to 0.53 and 0.26.
- All arms meet the frozen per-split rule, which does not check false partitions per
  dataset or Markov/Ward ARI; the protocol's replacement clause was not applied.
- Tree-chosen noise splits with an invariant direction: median misplaced 13-20%,
  15-41% at or below the 12% allowance, fully separable 0-12%.
- Correction: `allowed_overlap` had its Bayes-error terms swapped (right at 50/50, 7.3%
  instead of 1.3% at 10/90). Fixed, false partitions are 13/17/28 instead of 13/42/36;
  the conclusion stands ([[single-population-split-bound]]).

## Evidence

`raw/inbox/choir-cell-lines/rule-out-overlap-redesign-20261108.md` holds the rules,
grid, end-to-end tables and reading; protocol, deviations, code and outputs are in
`raw/inbox/choir-cell-lines/rule-out-overlap-redesign-scripts-20261108.py.txt`.

## Links

- [[rule-out-gradient]]
- [[split-count-exclusion]]
