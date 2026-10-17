---
title: Passthrough P-Value and Normalisation Audit
type: source
status: reviewed
updated: 2026-10-17
sources:
  - reports/passthrough_pvalue_normalisation_20260928/README.md
  - reports/passthrough_pvalue_normalisation_20260928/verification.json
  - reports/passthrough_pvalue_normalisation_20260928/node_reach_pvalue_summary.csv
  - reports/passthrough_pvalue_normalisation_20260928/paired_partition_changes.csv
  - reports/passthrough_pvalue_normalisation_20260928/upstream_null_root_factor_summary.csv
  - reports/passthrough_pvalue_normalisation_20260928/normalisation_audit.md
  - reports/passthrough_pvalue_normalisation_20260928/pvalue_audit.md
  - reports/passthrough_pvalue_normalisation_20260928/production_probe.json
  - reports/passthrough_pvalue_normalisation_20260928/production_verification.json
  - tree_break_selection/hierarchy_analysis/tree_decomposition.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py
tags:
  - calibration
  - traversal
  - normalisation
---

# Passthrough P-Value and Normalisation Audit

## Summary

Passthrough changes which precomputed sibling gates can affect the final
partition. The inspected production path computes p-values and edge-reachable
frontier BH before the passthrough boolean is applied. A changed reached-only
p-value distribution does not imply changed values for individual node tests.

## Key Points

- Reanalysis of the existing 10,000 matched Gaussian datasets finds identical
  raw and independently reconstructed BH p-values for all 30,000 node pairs,
  while 3,971 datasets change partition count and 7,942 extra internal-node
  occurrences are reached. Raw equality is by construction in this benchmark.
- Under global null, reached false-split counts rise from 12/1000 and 15/1000
  to 20/1000 and 23/1000. For signal only below A, A recovery rises from
  11/1000 and 12/1000 to 1000/1000 in both covariance profiles.
- Raw p-values can change through child membership, sample sizes, covariance,
  whitening, projection basis/rank and calibration. Adjusted p-values also
  depend on the correction family. The inspected sibling BH function does
  not accept the traversal passthrough boolean.
- The same matched global-null datasets give 33/1000 and 37/1000 root
  rejections under known-covariance same-data PCA/MP, versus 12/1000 and
  15/1000 under full coordinates. Same-data average-linkage root selection
  gives 765/1000 and 762/1000. These are upstream changes, not passthrough
  effects; PCA basis/rank and selected-root membership effects are bundled.
- P-values describe a node hypothesis under its reference law; they are not
  eigenvalue-cluster distances or final-clustering quality measures.
- A real-data Hamming/average-linkage probe on 626 observations and 456 binary
  features leaves all 1,251 annotation rows unchanged, including 345 calibrated
  sibling p-values and their BH corrections. Passthrough changes 1 cluster to
  109, with 182 reached internal nodes, 35 reached leaves, 32 significant splits
  and 76 pass-through decisions. No population truth is available for these
  real clusters; this probe does not test a diffusion tree.
- The GO workflow normalizes cosine/diffusion geometry for tree construction
  while passing original typed features to TBS. Bernoulli, categorical and
  continuous covariance models require their own coordinate contracts; no
  additional normalization follows merely from passthrough.
- Custom diffusion uses symmetric-operator coordinates rather than directly
  using paired random-walk left/right coordinates. Correlated Bernoulli
  features, ridge placement across whitening stages, parent-total covariance
  and same-data tree/projection selection remain calibration questions.

## Evidence

The report preserves per-node reach summaries, paired partition-count changes,
paired upstream rejection changes, and input/script hashes. Its source audit
distinguishes annotation construction from traversal. The known-tree benchmark
uses known Gaussian covariance and deliberately open edge gates; it does not
establish error control for the full inferred-tree pipeline. Reanalysis adds
no independent simulation datasets.

Two independent audits completed the source review and production probe. A
separate check verified the recorded per-node equalities, traversal counts, input hash and
seven current source hashes independently. Production hashes are post-run
records, not a complete pre-run snapshot. Targeted root baseline: 22 tests;
production probe baseline: six overlapping tests. Report scripts pass Ruff,
local report links resolve, and wiki lint passes 204 pages. No production
source was changed by this audit.

## Links

- [[known-branch-skew-diagnostic]]
- [[independent-projection-vandermonde-control]]
- [[eigenspace-benchmark-calibration-reassessment]]
- [[diffusion-nnls-calibration-open-points]]
- [[current-adaptive-diffusion-subspace-tree-pipeline]]
