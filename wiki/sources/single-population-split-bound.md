---
title: The Largest Split One Population Allows
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/single-population-split-bound-20261108.md
  - raw/inbox/choir-cell-lines/single-population-split-bound-scripts-20261108.py.txt
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
tags:
  - calibration
  - sibling-gate
  - power
  - tree-choice
---

# The Largest Split One Population Allows

## Summary

The overlap of a split does not decide whether it is real: a flat cut through one
population has zero overlap. What decides is the split's between/within ratio along its
own discriminant, rho = pi(1-pi) D^2, against the most a cut of one population can
produce. For one Gaussian population that ceiling is g(pi) = r/(1-r),
r = phi(Phi^-1(pi))^2/(pi(1-pi)) (1.75 balanced); it is 3 for cuts at the mode of any
unimodal population and reaches 8 for long flat tails. With a count of flat cuts for
the tree's choice, the allowed overlap grows with cells toward 9.3% (balanced, Gaussian
ceiling) instead of shrinking.

## Key Points

- Proof sketch and checks: the best of all splits is a flat cut (140/140 exhaustive
  cases); a flat Gaussian cut gives g(pi) (simulated within 0.2%).
- Finite-n rule: p = min{1, sum_{i<=k} binom(n-1, i) P(chi2_k(n bound(pi)) >= n rho)} on
  the k invariant directions; the original Gaussian median-cut check gave 0.95% at
  1%. This is not a general guarantee: [[three-tests-construction-audit]] finds
  an incorrect uniform-boundary variance and 40/1,000 rejections after the count
  on the sphere-3 boundary at n=100, nominal 0.01.
- Allowed overlap, balanced, Gaussian ceiling: 4.5% at 100 cells, 7.2% at 1,000, 8.9%
  at 48,879; ceiling 3: 1.8% to 4.0%.
- End to end on 935 datasets, Markov/Ward false partitions of 300 (sd 0.25/0.42/1.0): zero
  overlap 2/25/15, Gaussian ceiling 0/0/5, ceiling 3 0/0/0; ARI with ceiling 3 is
  1.00 for all four group types (density variation 0.87 before); the Gaussian ceiling
  over-splits skewed groups (density ARI 0.57).
- The earlier `allowed_overlap` formula had swapped terms (7.3% instead of 1.3% at a
  10/90 split); fixed, the overlap rule still gives 13/17/28 false partitions.
- Real CHOIR data and the continuous benchmark: see [[ceiling-rule-real-and-benchmark]]
  (Markov/Ward lines recovered 80 to 119; benchmark unchanged). Not run: long-tailed
  populations. No production change.

## Evidence

`raw/inbox/choir-cell-lines/single-population-split-bound-20261108.md` holds the
derivation, checks, tables and reading; protocol, deviations, code and outputs are in
`raw/inbox/choir-cell-lines/single-population-split-bound-scripts-20261108.py.txt`.
Literature looked up: Liu et al. 2008 (doi:10.1198/016214508000000454) and Cover 1965
(doi:10.1109/PGEC.1965.264137) via Crossref; the Khintchine representation via
arXiv:2208.08813; the original 1938 paper and Cover's count formula were not checked.

## Links

- [[rule-out-overlap-redesign]]
- [[rule-out-gradient]]
- [[split-count-exclusion]]
- [[three-tests-construction-audit]]
- [[gaussian-best-flat-cut-reference]]
