---
title: Marchenko-Pastur, Eigenvalue Multiplicity and the Definition of k
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/marchenko-pastur-multiplicity-literature-20261112.md
  - raw/inbox/choir-cell-lines/marchenko-pastur-multiplicity-scripts-20261112.py.txt
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/split_count_rule_out.py
  - tree_break_selection/eigenvalue_groups.py
tags:
  - literature
  - marchenko-pastur
  - multiplicity
  - invariant-directions
---

# Marchenko-Pastur, Eigenvalue Multiplicity and the Definition of k

## Summary

Looked-up literature and four numerical checks on k, the number of eigenvalues
of a node's covariance above the edge of `invariant_directions`. The edge is a
limit, not a quantile: one independent Gaussian sample gets a direction in 16%
to 23% of draws. A direction exists in the limit only when the population
`rho` exceeds sqrt(p/n). The rule "no direction when n <= p + 1" is not
required by the theory. Geometrical and algebraic multiplicity are equal for
every covariance matrix, and sample eigenvalues do not tie.

## Key Points

- Baik and Silverstein 2006 (read): eigenvalues of i.i.d. samples fill
  [(1 - sqrt(c))^2, (1 + sqrt(c))^2], c = p/n; a population eigenvalue alpha
  shows only above 1 + sqrt(c), at alpha + c alpha / (alpha - 1); the same
  holds for c > 1; a population eigenvalue of multiplicity m moves m sample
  eigenvalues to that one limit.
- Johnstone 2001 (abstract): the largest eigenvalue fluctuates with a
  Tracy-Widom law of order 1. Check C1: k >= 1 in 16.0% to 23.2% of 2,000
  independent draws of one Gaussian population, at five sizes.
- Gavish and Donoho 2014 (abstract): with unknown noise the threshold is a
  multiple of the median singular value; this is the device behind the median
  in the edge, stated for a low-rank matrix in noise.
- Own derivation with check C2: two groups add the population eigenvalue
  1 + pi(1 - pi) d^2, so a direction exists only for population `rho` above
  sqrt(p/n). At sqrt(p/n) = 0.306: `rho` 0.15 gives k >= 1 in 19.7%, 0.6 in
  89.9%, 1.44 in 100%.
- Check C3, 30 cells in 40 features: the current rule returns no direction
  in 1,000 of 1,000 draws; a prototype on the nonzero eigenvalues finds one in
  75.4% (`rho` 2.3) and 100% (`rho` 5 and 50), and in 21.5% for one population.
- A covariance matrix is real symmetric, so geometrical and algebraic
  multiplicity agree; the Markov operator is similar to a symmetric matrix.
- Okamoto 1973 (abstract): eigenvalues of a sample with a density are distinct
  with probability one. Bai and Yao 2008 (read): the m sample eigenvalues of a
  multiple spike fluctuate jointly like a Gaussian random matrix. Check C4:
  0 numerical ties in 1,000 draws; relative spread 0.115 at 400 cells and
  0.038 at 4,000.
- The exact ties of 12 November (29 tied eigenvalues) come from whitening
  rows with their own covariance, not from a population multiplicity.
- Limits: Marcenko and Pastur 1967 and Anderson 1963 verified as records
  only; Tracy-Widom thresholds in Kritchman and Nadler 2009 unverified; no
  source on subsets selected by clustering; checks use identity noise.

## Evidence

`raw/inbox/choir-cell-lines/marchenko-pastur-multiplicity-literature-20261112.md`
lists the eleven looked-up sources with the DOI or link used and what each
retrieved text confirms, the three meanings of multiplicity in the project,
the table of the four checks and the reading for k; protocol, script and
output are in
`raw/inbox/choir-cell-lines/marchenko-pastur-multiplicity-scripts-20261112.py.txt`.

## Links

- [[non-heuristical-statistic-plan]]
- [[local-marchenko-pastur-rule]]
- [[subtree-spectrum-topology]]
- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-benchmark]]
- [[cancer-pc-edge-test]]
- [[spr-eigenvalue-blocks]]
