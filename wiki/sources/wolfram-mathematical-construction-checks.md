---
title: Wolfram Mathematical Construction Checks
type: source
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/wolfram-mathematical-construction-checks-20261102.md
tags:
  - mathematics
  - wolfram
  - verification
---

# Wolfram Mathematical Construction Checks

## Summary

Captured Wolfram Language evaluator inputs and outputs for the Tree-Break
Selection mathematical glossary and construction review on 2026-11-02.

## Key Points

- Evaluation 1 verifies barycentric contrasts, nested variance cancellation,
  standardized edge/sibling identities, and scalar split action.
- A fixed unit projection of a correlated unit-variance vector has variance
  \(1+\rho\); orthonormal projection alone does not ensure a chi-square law.
- A truncated chi-square counterexample demonstrates the distinction between
  selected and unselected tails.
- Evaluation 3 verifies a categorical covariance interior example, the
  branch-factor ratio, adaptive-projection algebra, weighted quadratic moments,
  and equal means with unequal variances for two chi-square constructions.
- Evaluation 2 contains an unevaluated scaled-distribution expression; evaluation
  3 uses `TransformedDistribution` to resolve it. Its moment check also derives
  results from the MGF rather than returning assigned expressions.
- These are symbolic checks under stated assumptions, not TBS execution,
  empirical calibration, or a selected-inference theorem.

## Evidence

- `raw/inbox/wolfram-mathematical-construction-checks-20261102.md` preserves
  all three evaluator calls and returned text without normalizing away the
  intermediate unevaluated expression.
- The initial context-tool internal error is recorded; the successful evaluator
  calls supplied the computation evidence.

## Links

- [[mathematical-glossary]]
- [[mathematical-construction-wolfram-review]]
