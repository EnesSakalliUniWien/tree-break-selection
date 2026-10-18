---
title: Naming Consistency Audit
type: source
status: reviewed
updated: 2026-10-17
sources:
  - reports/naming_audit_20261017/README.md
  - raw/inbox/naming-audit-reference-notes-20261017.md
  - tree_break_selection/tree/poset_tree.py
  - tree_break_selection/tree/optimized_branch_lengths.py
  - tree_break_selection/space_separation/adaptive_cosine.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py
  - benchmarks/shared/metrics.py
  - benchmarks/validation/statistics/traversal_sibling_fdr_null.py
  - benchmarks/cloud/aws_traversal_sibling_fdr_null.py
  - applications/scrna/plots/selected_nnls_report.py
  - reports/eigenvalue_block_comparison_20260926/README.md
  - reports/eigenspace_reassessment_20260928/spectral_blocks/README.md
tags:
  - terminology
  - review
  - aliases
---

# Naming Consistency Audit

## Summary

The 2026-10-17 audit checked active names against implementations, callers,
and official statistical definitions. The verified corrections were
implemented with no compatibility aliases for retired names.

## Key Points

- `PosetTree.populate_node_distributions` names distribution preparation;
  it replaces `populate_node_divergences`, which computed no divergence.
- Generic coordinate APIs use `Adjustment`, covering both FDR and FWER
  methods. The redundant BH-only forwarding wrapper is removed.
- Unweighted spectral input uses `weighting="raw"`; actual binary data
  labels keep their original meaning.
- NNLS uses `mean_squared_standardized_euclidean` to expose its
  feature-count average. Numerical target calculations are unchanged.
- NNLS plots distinguish fitted length from time and disclose radial
  normalization. The metric entry point is `calculate_clustering_metrics`.
- Diagnostic schema v2 names the observed event `any_false_rejection_rate`
  and uses `empirical_control_observed` for finite below-alpha results.
  Merging rejects older or mixed schemas.
- Historical artifacts retain their identifiers. Current strict spectral
  readers reject retired weighting labels instead of silently discarding rows.
- Two hash-pinned report helpers keep `binary`: their recorded working tree
  matches no commit, so they stop at the source-hash check on later source.

## Evidence

The report records eight resolved findings, API/CLI/schema migrations,
focused checks, numerical comparisons, and scope limits. Official R
documentation distinguishes Holm/Bonferroni family-wise error control from
BH/BY FDR; SciPy defines the unaveraged standardized Euclidean distance.
Local captured reference notes preserve that source context. The full
`make check` gate passed with 895 tests and one existing scikit-learn warning;
the quick-start pipeline also completed. Numerical comparisons, migration
tests, and scope limits are recorded in the report and [[log]].

## Links

- [[naming-and-claims-review]]
- [[dead-code-alias-duplication-audit]]
- [[poset-tree]]
- [[projected-wald-statistic]]
