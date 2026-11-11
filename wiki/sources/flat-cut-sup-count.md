---
title: The Tail of One Population's Best Cut
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/root-validation-calibration-20261109.md
  - raw/inbox/choir-cell-lines/flat-cut-sup-count-20261109.md
  - raw/inbox/choir-cell-lines/flat-cut-sup-count-scripts-20261109.py.txt
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
tags:
  - calibration
  - sibling-gate
  - tree-choice
---

# The Tail of One Population's Best Cut

## Summary

The count of all flat cuts in the ceiling rule can be replaced by the tail of one
population's best cut: an even Gaussian field over the sphere of directions, whose
maximum is approximated by an expected-Euler-characteristic tail. The original
Gaussian simulations were close to nominal and tighter than the count. The later
[[three-tests-construction-audit]] finds finite-sample miscalibration and an unmet
C2 regularity assumption; the envelope is not a proved survival bound.
For ceiling 3 the field's smoothness depends on
the least favourable shape, which the ceiling does not fix; with the k = 3 value it is
no tighter below about 1,000 cells. None of the 27 hidden-boundary splits passes.

## Key Points

- The tree's split is never above the best flat cut. Its exact survival law under
  a specified iid population would control a selected split of that sample;
  the implemented approximation and selected descendants need separate justification.
- Frozen version (rho scale) failed its null check (up to 22% at 1%); amended to the R
  scale: fresh seeds, k = 1-3, 0-1.3% at 1%. A sign-swing bug for k >= 4 (p = 0 at the
  ceiling) was fixed with a non-increasing envelope before use.
- Gaussian null, k = 6-9: needed rho 2.23-2.95 vs null 99th percentile 2.18-2.77 (count:
  2.52-3.21). Ceiling 3, uniform-product null: 3.71-5.37 vs 3.34-4.06 (count 3.82-4.83).
- Ceiling 3, k = 9: needed rho 16.2 at 100 cells, 6.0 at 300, 4.2 at 1,000, 3.5 at 5,000.
- At the merged real clusters' sizes even an exact count needs rho near 4; the tree's
  splits there reach 1.6-4.5, so blurred tree splits, not the count, hide the boundary.
- Not done in the original study: end-to-end runs. The later
  [[three-tests-construction-audit]] includes full-cell C fits; a least-favourable
  law above k = 3 remains unresolved, and uniform marginals in every direction
  are impossible there.
- In the follow-up all-coordinate n=32,k=20 diagnostic, C clips to 1 for
  every attainable statistic: the envelope already exceeds one at the upper
  endpoint R=1. Its zero rejections there are a power failure. B rejects
  818/1,000 searched Gaussian coordinate-threshold cuts at nominal .01 in the
  same diagnostic. These are isolated tail checks, not production traversal.

## Evidence

`raw/inbox/choir-cell-lines/flat-cut-sup-count-20261109.md` holds the derivation,
checks, application and literature; protocol, deviations, code and outputs are in
`raw/inbox/choir-cell-lines/flat-cut-sup-count-scripts-20261109.py.txt`.

## Links

- [[root-validation-generalisation-evidence]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[single-population-split-bound]]
- [[ceiling-rule-real-and-benchmark]]
- [[ceiling-rule-literature-comparison]]
- [[three-tests-construction-audit]]
- [[gaussian-best-flat-cut-reference]]
