---
title: Naming and Claims Review
type: source
status: reviewed
updated: 2026-10-17
sources:
  - reports/naming_claims_review_20260928/README.md
  - reports/naming_claims_review_20260928/fixes.md
  - README.md
  - tree_break_selection/hierarchy_analysis/tree_decomposition.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/space_separation/leading_axes.py
  - tree_break_selection/tree/optimized_branch_lengths.py
  - tree_break_selection/space_separation/diffusion.py
  - applications/scrna/plots/selected_nnls_fit_summary.py
tags:
  - terminology
  - review
  - calibration
---

# Naming and Claims Review

## Summary

Three reviews covered statistical names, geometry/solver names and
documentation claims. Seven actionable mismatches were verified against
implementations and callers. The subsequent user-authorized
corrections addressed all seven through clarified descriptions, plot
labels and metadata. Numerical formulas and existing API identifiers remain
unchanged; the original review is explicitly marked as the pre-fix snapshot.

## Key Points

- Fixed-subspace descriptions now separate fixed projection/rank from
  same-data covariance estimation and tree selection.
- README and TreeDecomposition now describe the active support-eligible
  plug-in rule and BH without claiming a validated-reference prerequisite.
- The passthrough report now describes empirical-inflation plug-in values and
  explicitly leaves selected-tail validity unproven.
- Retained invariant/equivariant identifiers are documented as operational
  labels; visible scRNA labels describe leading and following PCA axes.
- NNLS output now records the target definition, normalization and divisor.
  The plotting caption describes fitted branch lengths rather than time.
- Adaptive diffusion metadata now separates the pydiffmap kernel from the
  repository symmetric-coordinate computation and SciPy eigensolver.
- Properly qualified internal-admissibility, exact-F, retained-rank, numerical
  tie-group and NNLS terms should remain. Historical evidence was not treated
  as current merely because it remains in the repository.

## Evidence

The report gives seven numbered findings with source locations, observed
behavior, consequences, confidence and suggested clarifications. Review scope
includes existing dirty/untracked work at HEAD
`f34f59cacd09b19c34fafa3a4a6e783dd3973f72`. This is a targeted source review,
not exhaustive coverage or a new calibration experiment. The follow-up records
71 passing focused tests, numerical equality in two NNLS snapshot comparisons
and one diffusion comparison, and compatibility of retained identifiers.
Documentation and metadata corrections do not prove statistical calibration.

## Links

- [[eigenvector-calibration-construction-review]]
- [[empirical-null-calibration-reference-law-contract]]
- [[independent-diffusion-nnls-benchmark]]
