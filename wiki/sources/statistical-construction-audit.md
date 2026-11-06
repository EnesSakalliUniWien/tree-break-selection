---
title: Statistical Construction Audit
type: source
status: reviewed
updated: 2026-11-04
sources:
  - reports/statistical_construction_audit_20261104/README.md
  - reports/statistical_construction_audit_20261104/mp_audit.md
  - reports/statistical_construction_audit_20261104/reference_law_audit.md
  - reports/statistical_construction_audit_20261104/knowledge_audit.md
  - reports/statistical_construction_audit_20261104/mp_results.json
  - reports/statistical_construction_audit_20261104/reference_results.json
  - reports/statistical_construction_audit_20261104/downstream_results.json
  - reports/statistical_construction_audit_20261104/verification.json
  - raw/inbox/tbs-mp-primary-sources-20261104.md
  - raw/inbox/tbs-reference-law-primary-sources-20261104.md
  - raw/inbox/tbs-knowledge-primary-sources-20261104.md
  - raw/inbox/tbs-multiplicity-primary-sources-20261104.md
tags:
  - statistics
  - spectral
  - calibration
  - audit
---

# Statistical Construction Audit

## Summary

Three Sol-agent audits and coordinator boundary checks trace the statistical
construction from row ownership through covariance, Marchenko–Pastur, PCA,
degrees of freedom, calibration, multiple testing and traversal. Correct rank
bookkeeping does not establish a chi-square law for the selected statistic.
No production inference was replaced.

## Key Points

- MP spike count, nominal projection dimension, retained rank and reference
  df are distinct. At energy fraction 0.90 the tested contrast can replace a
  positive nominal rank. Fitted full whitening on the same 60-by-30 Gaussian
  rows produces zero MP spikes but nominal rank 30 through numerical ties.
- The edge PCA uses assigned basis halves, while its contrast denominator can
  use full-data covariance. Sibling PCA uses all leaves and optional internal
  rows. Leaf-only PCA and splitting after selecting the tree do not establish
  independence of the statistic and its nuisance estimates.
- Fixed Gaussian controls identify chi-square, scaled-Beta and Hotelling-F
  laws under different assumptions. With the production guard satisfied at
  8+8 rows and six features, full means give 6,540/20,000 nominal chi-square
  rejections at 0.05 versus 1,027 under the restricted exact F reference.
  That F reference does not apply to the alternating edge means with the same
  full-data denominator. Native quadratic forms match their formulas.
- Adaptive energy rank alone gives 1,562/20,000 rejections at 0.05 with known
  covariance and a fixed eight-axis basis. A dimension change alone cannot
  repair covariance, rank and tree-selection dependence.
- Calibration roles can include mean-different records. Group ESS describes
  support concentration, while scale estimation remains record-weighted.
  Tree-BH arithmetic, per-depth sibling frontiers and live traversal are
  separate objects; published error guarantees require further premises.
- Current wiki synthesis was corrected for the leaf-only independence claim,
  spectral contexts, nominal/final ranks, zero-basis edge sentinel and exact-K
  oracle scope. Manuscript discrepancies remain explicitly listed. Historical
  numerical evidence was preserved.

## Evidence

The three agent reports give file/line traces, formulas, primary sources and
claim-verdict tables. Independent review verified the restricted Gaussian
derivations and caught a production-covariance-guard mismatch in an initial
mathematical control; an additional native 8+8 replay resolved that gap.
The MP controls use 3,000 draws per shape; reference controls use 20,000 per
configuration. Focused test batches passed 37, 92 and 33 tests, with overlap
between the first two batches. These checks establish their stated mechanisms,
not full-pipeline calibration.

Two software-level gaps are recorded without a production change: inconsistent
augmented zero-rank threshold-row metadata, and missing orthonormality
validation for arbitrary direct kernel callers. Full cancer-data covariance,
selected descendant laws, numerical event certification and final-partition
error control remain unverified.

## Links

- [[local-marchenko-pastur-rule]]
- [[selected-pca-projected-wald-validation]]
- [[mathematical-glossary]]
- [[method-proof-web]]
- [[empirical-null-calibration-reference-law-contract]]
- [[cancer-pc-statistics-diagnosis]]
