---
title: Independent Gaussian Production Benchmark
type: source
status: deprecated
updated: 2026-10-18
sources:
  - reports/production_pvalue_replacement_20260928/README.md
  - reports/production_pvalue_replacement_20260928/protocol.md
  - reports/production_pvalue_replacement_20260928/evaluation/manifest.json
  - reports/production_pvalue_replacement_20260928/evaluation/summary.csv
  - reports/production_pvalue_replacement_20260928/verification.json
tags:
  - calibration
  - projection
  - independent-data
---

# Independent Gaussian Production Benchmark

## Summary

An explicit production fit/test API separates discovery topology/PCA from
independent Gaussian measurements of the same leaves. It requires externally
known common covariance and independent rows. Conditional energy p-values
replace nominal adaptive-rank tails before both edge and sibling corrections;
no empirical inflation is fitted. Ordinary single-matrix calls lack these
inputs and retain their existing behavior.

The API was removed on 2026-10-18: nothing outside tests and the report scripts called it. The report evidence stays, but its scripts no longer run against current code.

## Key Points

- Fit consumes only discovery measurements and known covariance. Test consumes
  frozen discovery context plus a separate matrix with matching leaf/feature
  labels. Observation IDs must be disjoint; copied matrices are rejected.
  These checks do not prove independence or make estimated covariance known.
- The reference integrates chi-square candidate-radius survival over directional
  shares conditional on the retained complete eigenvalue group. It generalizes
  the two-dimensional angular formula using conditional Beta shares. Projection
  dimension remains distinct from the mixture's reference-law identity.
- Existing tree construction, eigendecomposition, numerical eigenvalue groups,
  Tree-BH, sibling BH and TreeDecomposition are reused. Under known common
  covariance, binary child-parent contrasts equal the sibling contrast up to
  sign. No branch-length stochastic variance or same-data inflation is added.
- Evaluation uses 900 independent measurement pairs across three Gaussian
  profiles and eight tree methods. There are 9,267 completed runs per independent
  reference, with 1,533 unsupported tree cases per reference. Existing
  centroid/median construction rejects nonmonotone merge heights. No unexpected
  runtime errors occurred. A separate current single-matrix comparison completes
  2,700 runs and changes data quantity/covariance as well as the equation.
- High-dimensional average-linkage root rejections at .01 fall from 5/250 to
  1/250 on matched independent contrasts. High-dimensional neighbor joining
  remains 8/250; this pointwise excess remains visible. A post-evaluation Holm
  diagnostic over 96 planned cells yields minimum adjusted p .3864. Rare-tail
  precision is limited and absence of rejection is not proof of calibration.
- Equation changes alter cluster counts in 32 matched primary-alpha outcomes.
  Signal ARI sometimes decreases; high-dimensional recovery remains weak.
  Large gains against single-matrix methods cannot be assigned to the equation
  alone because independent data and known covariance also change.

## Evidence

146 focused tests pass, including 24 new equation/data-flow tests. Twelve
independent adaptive-quadrature checks agree within 7.11e-15. Source hashes,
10,800 paired variant rows, 175,677 node rows, root joins and p-value ranges are
verified. Report CSVs retain unsupported denominators and pointwise intervals.
The earlier six-family baseline completes 1,960 of 2,380 runs and explicitly
skips 420, without runtime or audit errors. It does not validate the new Gaussian
contract on Bernoulli/categorical inputs or binary diffusion.

## Links

- [[distribution-equation-fitting-comparison]]
- [[independent-projection-vandermonde-control]]
- [[lower-dimension-degrees-of-freedom-audit]]
- [[open-mathematical-questions]]
