---
title: Selected PCA Projected-Wald Validation
type: analysis
status: reviewed
updated: 2026-11-04
sources:
  - benchmarks/validation/manifests/selected_pca_projected_wald_validation_manifest.json
  - raw/assets/selected-pca-projected-wald-validation/20260601-selected-pca-projected-wald-calibration.json
  - raw/assets/selected-pca-projected-wald-validation/20260601-selected-pca-projected-wald-calibration.csv
  - benchmarks/results/validation/selected_pca_projected_wald_calibration_check_20260605.json
  - benchmarks/results/validation/selected_pca_projected_wald_calibration_check_20260605.csv
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/tree_estimator.py
  - reports/statistical_construction_audit_20261104/knowledge_audit.md
  - reports/statistical_construction_audit_20261104/reference_law_audit.md
  - reports/cancer_pc_statistics_debug_20261104/README.md
tags:
  - method
  - projection
  - validation
---

# Selected PCA Projected-Wald Validation

## Summary

The three recorded fixed-membership Gaussian settings have rejection rates
near the nominal level with leaf-only PCA, and much larger rejection rates
when deterministic child-mean rows are included in the spectral matrix.
These are bounded empirical findings, not an independence theorem for
leaf-only PCA. The locked run used commit
`492c8520809e6cfad9e4853e91c89a887eac5a21`, seed `20260601`, 1000 replicates
per setting, and an empty recorded worktree status.

This is not full-pipeline evidence. It measures one local question: whether
the chi-square projected-Wald reference remains calibrated when PCA rows and the
MP dimension are selected from the same fixed-membership Gaussian null context.
It does not validate hierarchy construction, tree-selected sibling pairs,
sibling FDR, traversal, empirical-null inflation, categorical blocks, or
real-data misspecification.

A 2026-06-05 rerun on the then-current worktree reproduces the locked
pattern: leaf-only rejection rates at nominal `0.05` are `0.046`, `0.048`, and
`0.059`, while the corresponding child-mean-row settings reject at `0.671`,
`0.992`, and `1.000`.

## Details

The diagnostic simulates two sibling samples from the same Gaussian
distribution. For each replicate, it computes the production
empirical-Gaussian sibling Wald contrast, maps the local rows into the
null-whitened tangent coordinates used by the spectral backend, selects PCA rows
and projection dimension through the production MP path, and evaluates the
production projected-Wald kernel.

The paired validation grid compares leaf-only local spectral rows against the
same setting with child-mean rows appended as deterministic internal rows:

| setting | row mode | rejection rate at 0.05 | 95% Wilson interval | mean p-value | KS p-value |
| --- | --- | ---: | ---: | ---: | ---: |
| 8-dimensional identity Gaussian, 80x80 | leaves | 0.047 | 0.0355--0.0619 | 0.5023 | 0.9673 |
| 8-dimensional identity Gaussian, 80x80 | child means | 0.648 | 0.6179--0.6770 | 0.0625 | 0 |
| 16-dimensional AR(1) Gaussian, 60x60 | leaves | 0.045 | 0.0338--0.0597 | 0.5003 | 0.6877 |
| 16-dimensional AR(1) Gaussian, 60x60 | child means | 0.988 | 0.9791--0.9931 | 0.0036 | 0 |
| 32-dimensional ill-conditioned Gaussian, 50x50 | leaves | 0.050 | 0.0381--0.0653 | 0.4861 | 0.1987 |
| 32-dimensional ill-conditioned Gaussian, 50x50 | child means | 1.000 | 0.9962--1.0000 | 0.0000053 | 0 |

All six settings used mean projection dimension 2 and mean raw MP signal count
0. The anti-conservative behavior is therefore not caused by the MP rule
counting many false spikes in this run. It arises even when the selected
dimension is only the configured spectral floor. The leaf-only rows are
compatible with the fixed-subspace approximation in this diagnostic; the
child-mean rows are not.

This means Q12 and Q13 should be kept separate. Q12's fixed-projection
chi-square statement is valid only under fixed or conditionally independent
valid projections. Q13's selected-PCA question is not globally negative:
leaf-only fixed-membership Gaussian selected PCA is calibrated in this
diagnostic. The failure is the leakage mode where deterministic summaries of
the tested sibling means enter the spectral row set, plus the broader
selected-tree setting where the hierarchy and focal sibling are selected using
the same evidence.

Leaf-only PCA is not generally independent of the sibling contrast, even
with fixed groups. The total centered scatter equals pooled within-group
scatter plus \(n_Ln_R/(n_L+n_R)\,\Delta\Delta^\top\), where
\(\Delta=\bar X_L-\bar X_R\). The tested contrast therefore enters the
PCA input. The observed near-nominal rates apply to the sizes and covariance
conditions in the table. In the later two-leaf control, a fitted rank-one
direction can retain all 30 coordinates' contrast energy. Appending child
means introduces another form of contrast reuse, but removing those rows
alone does not establish the fixed-subspace law. See
[[statistical-construction-audit]] for the covariance and rank controls.

## Evidence

- `raw/assets/selected-pca-projected-wald-validation/20260601-selected-pca-projected-wald-calibration.json`
  records the locked JSON report, including commit, command, seed, worktree
  status, grid, confidence intervals, p-value uniformity summaries, and
  limitations.
- `raw/assets/selected-pca-projected-wald-validation/20260601-selected-pca-projected-wald-calibration.csv`
  records the row-level summary table.
- `benchmarks/results/validation/selected_pca_projected_wald_calibration_check_20260605.csv`
  records the current-worktree rerun used to check that the locked conclusion
  still reproduces.
- `benchmarks/validation/statistics/selected_pca_projected_wald_calibration.py` defined
  the scaffold and validation contract (retired 2026-10-17).
- `tree_break_selection/hierarchy_analysis/statistics/projection/spectral/tree_estimator.py`
  implements separate contexts: the edge basis uses its assigned leaf rows;
  the full-descendant sibling context uses leaves by default and can append
  internal distributions when explicitly configured. The archived child-mean
  comparison remains evidence for its particular construction.
- `benchmarks/validation/statistics/selected_pca_projected_wald_calibration.py` was retired on 2026-10-17 (last present at commit 5c554bc5).

## Links

- [[open-mathematical-questions]]
- [[local-marchenko-pastur-rule]]
- [[projected-wald-statistic]]
- [[statistical-construction-audit]]

## Open Questions

The edge basis excludes internal rows, while the sibling context retains an
opt-in internal-row path. Neither leaf-only sibling PCA nor splitting rows
after selecting a tree establishes a calibrated reference. Covariance fitting,
contrast-energy rank, tree selection and sibling calibration support remain
separate questions. High-dimensional and high-cardinality contexts can have no
strict empirical-null calibration records even when the raw sibling evidence is
strong. That requires either a validated external/null calibration object or an
explicit unsupported status; it should not be patched by reintroducing internal
summary rows into the inferential basis.
