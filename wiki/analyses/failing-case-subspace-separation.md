---
title: Failing-Case Subspace Separation
type: analysis
status: draft
updated: 2026-10-18
sources:
  - reports/failing_case_subspace_analysis_20261018/README.md
  - reports/failing_case_subspace_analysis_20261018/subspace_failure_analysis.py
  - reports/failing_case_subspace_analysis_20261018/pair_subspace_check.py
  - reports/failing_case_subspace_analysis_20261018/failing_case_subspace_summary.csv
  - reports/failing_case_subspace_analysis_20261018/failing_case_top_splits.csv
  - reports/failing_case_subspace_analysis_20261018/pair_subspace_heldout.csv
  - reports/failing_case_subspace_analysis_20261018/tbs_rerun_failing_cases.csv
  - reports/failing_case_subspace_analysis_20261018/bug_check/bug_harness.py
  - reports/failing_case_subspace_analysis_20261018/bug_check/bug_checks.csv
  - reports/failing_case_subspace_analysis_20261018/bug_check/bug_verdicts.csv
  - reports/failing_case_subspace_analysis_20261018/bug_check/edge_basis_experiment.py
  - reports/failing_case_subspace_analysis_20261018/bug_check/edge_basis_experiment.csv
  - reports/failing_case_subspace_analysis_20261018/bug_check/mahalanobis_separation_check.csv
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py
  - tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py
  - tree_break_selection/tree/distributions.py
  - tree_break_selection/tree/continuous_distance.py
  - benchmarks/shared/runners/tbs_support.py
  - benchmarks/shared/util/method_execution.py
tags:
  - analysis
  - benchmark
  - subspace
  - distance
  - gates
---

# Failing-Case Subspace Separation

## Summary

Of the 49 full-suite cases where canonical `tbs` fails (`unsupported` or ARI
below 0.5), only 10 fail because of the geometry of the space the tree is built
in. In 21 cases the TBS tree already holds every true cluster as a clade and the
gate fails; 11 cases have no location signal that any distance tree recovers;
5 spread their separation over many weak features; and 2 are null cases.

The geometric failures take three forms. Four `dim_consolidated` cases separate
in one shared subspace of 12 or 13 variables that the full-space distance drowns
in noise features. Four `dim_diffuse` cases separate each cluster pair in its
own small subspace, so no single tree serves all pairs. Two cases depend on the
metric or linkage.

A per-case bug check found no fault in data, distances, trees, labels or
scoring, but three method defects. The edge test fits its projection basis to
the leaves it tests, so it rejects nearly every edge, even on null data; the
continuous edge statistic cannot exceed n_parent - 1; and `mahalanobis_time`
removes mean separation. These defects are the cause in 17 cases, are present
in 22 more where the data would still limit recovery, and possibly contribute
in 2.

## Details

The analysis uses the distance TBS actually used per case (Hamming for binary,
median-binary and one-hot data; precomputed `standardized_euclidean`,
`mahalanobis_time` or `sbm_shifted_modularity` otherwise) and `build_tree` with
average linkage. It measures near-equidistance, the separating features of each
true cluster pair (Welch t, BH q < 0.01), oracle cuts at the true K on the TBS
tree, on the separating subspace and under other metrics, held-out pair
subspaces, and the features that drive the top TBS splits.

### Groups

| Group | Cases | Fails at |
| --- | ---: | --- |
| Tree holds the clusters | 21 | gate |
| One shared subspace drowned by noise | 4 | full-space distance |
| Pair-specific subspaces | 4 | any single tree |
| Metric or linkage choice | 2 | TBS distance or linkage |
| Separation spread over many weak features | 5 | tree in part, then gate |
| No location signal in the space used | 11 | every distance tree |
| Null cases (K = 1) | 2 | gate |

### Gate failures on a correct tree

All 21 trees cut at the true K with ARI 1.0 (0.987 once), and the ten
strongest features of every root split are true separating features. Twelve
cases are `unsupported` because the child-parent edge test rejects every tested
edge, which leaves no null-like sibling record to calibrate; the null cases
`gauss_null_small` and `gauss_null_large` show the same pattern on pure noise.
This matches the same-data edge selection in [[oracle-gate-path-diagnostic]].
The three `binary_perfect` cases are `unsupported` because every cluster is
copies of one row, so all nearest neighbours tie at distance 0 and no record has
positive degrees of freedom. Two continuous Gaussian cases stop at a root edge
test with p = 0.027 and 0.048 against edge alpha 0.001 although the clusters are
4.7 to 8.0 times farther apart between than within. Four cases get wrong sibling
decisions.

### Subspace drowning

In `dim_consolidated_4c_24f`, `_72f`, `_272f` and
`dim_consolidated_4c_72f_continuous`, the between/within distance ratio is 1.26,
1.08, 1.02 and 1.08 over all features but 1.60 to 2.07 on the separating
features. The TBS oracle cut falls from 0.50 to about 0 as noise features are
added, while held-out pair subspaces separate the pairs at ARI 0.86 to 1.0. At
72 and 272 features none of the root split's ten strongest features separates
any true pair: the root splits along noise. A cosine distance recovers the
continuous case completely.

### Pair-specific subspaces

In the four `dim_diffuse_6c` cases each pair separates on a median of 0 to 2
held-out features and the feature sets of different pairs barely overlap
(Jaccard 0.10 to 0.19). A pair's own features separate it better (held-out ARI
0.12 to 0.38) than the union of all pairs' features or all features (at most
0.20). What separates one pair is noise for the others, so one tree cannot
order all six clusters.

### Equidistance and metric

Exact equidistance appears when the signal is in the covariance:
`mp_spike_above_bbp_continuous` and `cont_lowrank_pggn_shrinkage` have a
distance coefficient of variation of 0.003 and 0.000 under `mahalanobis_time`,
and their trees peel single points. The metric is decisive in three cases:
cosine for `dim_consolidated_4c_72f_continuous`, Ward for `cat_mod_3cat_4c`
(0.96 against 0.69), and Euclidean for `phylo_brownian_null_16taxa`, whose
`mahalanobis_time` root splits one point from 159 and blocks the whole tree.

### Bug check

`bug_check/bug_harness.py` reruns each case through `run_single_method_once`
and checks every stage; all 49 cases pass. Three defects remain:

- **Edge basis fitted to the tested leaves.** The projected Wald kernel compares
  the projection of the child-parent difference onto the parent's own leading
  components with chi-square on k degrees of freedom, a law its docstring
  states only for a fixed basis. On pure-null data the production basis
  rejects 98 to 100% of edges on the data tree and 73 to 80% on a random tree;
  a basis from an independent null replicate rejects 0 to 3%. This refines the
  edge-selection explanation in [[oracle-gate-path-diagnostic]]: the basis,
  not the tree selection, drives the rejections. 26 cases reject every edge and
  return `unsupported`.
- **Continuous statistic capped at n_parent - 1.** Whitening by the parent's own
  covariance bounds the statistic; with df = p the root of four continuous
  cases can never reject. Above a memory limit only per-feature variances are
  stored and the first defect applies instead.
- **`mahalanobis_time` removes separation.** Whitening by the total covariance
  gives ARI 0 on two clusters shifted 3 standard deviations on every feature,
  where Euclidean gives 1.0.

The 8 cases without a defect are `binary_perfect_*` (identical rows; designed
fail-closed), the two Gaussian outlier cases, `dim_consolidated_4c_72f_continuous`
and the two SBM cases.

## Evidence

- `reports/failing_case_subspace_analysis_20261018/README.md` gives the
  method, every group's case list and numbers, the limits and the reproduction
  commands; a second run reproduced all tables byte for byte.
- `failing_case_subspace_summary.csv` holds the per-case measurements and
  diagnosis; `failing_case_top_splits.csv` the top four TBS splits per case;
  `pair_subspace_heldout.csv` the held-out pair checks;
  `tbs_rerun_failing_cases.csv` the rerun on the current code, including the
  `unsupported` evidence counts.

## Links

- [[oracle-gate-path-diagnostic]]
- [[tree-construction-method-map]]
- [[selected-hierarchy-selection-geometry]]
- [[traceable-mathematical-benchmark-suite]]
- [[open-mathematical-questions]]

## Open Questions

- Which independent basis should the edge test use (a sample split, a basis
  that leaves out the child, or a fixed design basis), and does it restore
  calibration support without opening false splits on the null cases?
- Should the continuous edge test whiten by a within-child covariance (the
  existing guarded policy) so its statistic is no longer capped?
- Should `mahalanobis_time` whiten by a within-cluster or shrunk diagonal
  covariance, or be replaced for cluster recovery?
- Would a tree per selected subspace, or a subspace-weighted distance, recover
  the `dim_consolidated` and `dim_diffuse` families without using labels?
- Covariance-only and frequency-shape signals (group 6) need a non-location
  representation before any distance tree can see them.
