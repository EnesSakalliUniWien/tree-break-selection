---
title: Selective and Rebuilding Replacement Comparison
type: source
status: reviewed
updated: 2026-11-04
sources:
  - reports/selective_rebuild_comparison_20261104/README.md
  - reports/selective_rebuild_comparison_20261104/protocol.md
  - reports/selective_rebuild_comparison_20261104/results.csv
  - reports/selective_rebuild_comparison_20261104/summary.csv
  - reports/selective_rebuild_comparison_20261104/manifest.json
  - reports/selective_rebuild_comparison_20261104/verification.json
  - reports/selective_rebuild_comparison_20261104/chain_null_control.json
  - raw/inbox/selective-rebuild-primary-sources-20261104.md
tags:
  - selective-inference
  - diffusion
  - calibration
---

# Selective and Rebuilding Replacement Comparison

## Summary

A paired root-only known-I Gaussian panel evaluates conditional inference
and whole-procedure rebuilding on the native diffusion/average tree. Changing
the conditional sampler recovers some strong signals without changing the
observed tree or full-history event. This supports a development candidate,
not a calibrated cancer-data replacement.

## Key Points

- At 0.05, ordinary full-history sampling detects zero of six strong cases,
  a random-position reversible-chain rank detects three, and conditioning the
  chain only on the root split detects five. Both chain versions reject zero
  of 12 null cases and one of six moderate cases. This panel is too small to
  establish empirical error control.
- Importance estimates produce the same strong-case counts but reject one
  null case. They are finite-simulation approximations; their strong-case
  root effective sample sizes range from 12.67 to 61.93 out of 999 proposals.
- Full-history matches imply root-split matches. A root-split event fixes the
  selected two children, not every merge or the entire descendant test family.
- Rebuilding known-I null datasets detects six strong cases and rejects one
  null case. Fitting sample covariance for the reference detects none of the
  strong cases. This is a different hypothesis: one Gaussian population,
  rather than equality of the two selected groups' average means.
- A reversible chain can give a valid conservative random-anchor rank under
  an exact event and Gaussian model without having mixed. Poor movement still
  reduces detection. Native numerical event membership is not exact-real
  certification. Settings were fixed before the panel.
- The 2026-10-12 Chen–Gao unknown-covariance preprint is relevant further work,
  but its common-mean null is stronger and its Gaussian/rank assumptions do
  not resolve the current cancer-data model.

## Evidence

The fixed panel contains 12 null, six moderate and six strong nine-by-three
datasets. All 192 candidate statistics completed. Native replay baseline:
16 tests. Verification recomputes the saved ranks, checks 144 chain points
and 112 source hashes. Additional finite-state, truncated-chi and 8,000
known-law chain draws check the numerical construction. The six strong
input hashes match the earlier 999-proposal zero-acceptance examples.

This is a diagnostic adaptation, not an SHC package reproduction. No production
statistic was replaced. Cancer covariance and preprocessing, numerical event
certification, descendant hypotheses and final-partition error control remain
open. The report retains simulation budgets, paired-dependence limits and the
distinction between rank tests and importance approximations.

## Links

- [[root-conditional-sampler]]
- [[conditional-diffusion-monte-carlo-evidence]]
- [[statistical-construction-audit]]
- [[choir-cell-lines-recovery-plan]]
- [[analytic-selected-root-laws]]
