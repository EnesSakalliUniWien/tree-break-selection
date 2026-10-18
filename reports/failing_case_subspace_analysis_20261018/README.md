# Failing-case subspace analysis (2026-10-18)

## Question

Canonical `tbs` fails on 49 of the 122 full-suite cases: in the 2026-10-16 run
it returned `unsupported` or an ARI below 0.5. Do these cases fail because the
separation lives in a subspace, because different cluster pairs separate in
different subspaces, because points are near-equidistant in the space the tree
uses, or because the distance metric does not suit the data? And which tree
splits show differences in which dimensions?

## Method

All measurements use the exact distance TBS used for each case: Hamming for
binary, median-binary and one-hot data, the precomputed `standardized_euclidean`
or `mahalanobis_time` distance for continuous data, and the precomputed
`sbm_shifted_modularity` distance for graphs. Trees are built with `build_tree`
and average linkage, as in TBS.

- **Equidistance:** coefficient of variation of the pairwise distances, the
  share of points whose nearest neighbour is tied, the share of nearest
  neighbours with the same true label, and mean between-cluster over mean
  within-cluster distance.
- **Subspaces:** for every pair of true clusters, the features whose Welch
  t-test passes BH at q < 0.01; one-hot columns `F<v>_c<k>` are also counted per
  variable `F<v>`. Overlap is the mean Jaccard index between the feature sets
  of different pairs.
- **Topology:** ARI of the TBS tree cut at the true K (oracle cut), and the
  share of true clusters that are a clade (Jaccard >= 0.9). The same oracle cut
  on a tree built from the separating features only, and on trees with other
  metrics (Hamming, Jaccard, cosine; or Euclidean, standardized Euclidean,
  cosine, correlation, city-block) and Ward linkage.
- **Held-out subspaces:** the separating features above are chosen with the
  true labels, so they flatter the subspace. `pair_subspace_check.py` selects
  each pair's features on half of its rows and scores a 2-cluster cut on the
  other half, using all features, the union of all pairs' features, or only
  that pair's features.
- **Splits:** for the TBS root and the next three largest splits, the child
  sizes, the true-label mix of each child, the ten features with the largest
  Welch |t| between the children, and the share of those that separate some
  true pair.
- **Gate:** the root diagnosis and edge counts from the 2026-10-16 audit, and
  the `unsupported` evidence counts from the rerun.

## Results

| Group | Cases | What fails |
| --- | ---: | --- |
| 1. Tree holds the clusters | 21 | the gate |
| 2. One shared subspace drowned by noise features | 4 | the full-space distance |
| 3. Each pair separates in its own small subspace | 4 | any single tree |
| 4. Metric or linkage choice | 2 | the TBS distance or linkage |
| 5. Separation spread over many weak features | 5 | partly the tree, then the gate |
| 6. No location signal in the space used | 11 | every distance tree |
| 7. Null cases (K = 1) | 2 | the gate |

**1. The tree holds the clusters; the gate does not split (21).** The oracle
cut of the TBS tree has ARI 1.0 (0.987 for `cat_mod_4cat_6c`), every true
cluster is a clade, and the root split's ten strongest features are all true
separating features. The failure is the gate:

- 12 cases return `unsupported` because the child-parent edge test rejects
  every tested edge, so no null-like sibling record is left to calibrate the
  empirical null: `binary_balanced_low_noise__2`, `binary_many_features`,
  `binary_noise_feat_30i_500n`, `binary_noise_feat_80i_400n`,
  `phylo_divergent_8taxa`, `phylo_dna_16taxa_low_mut`, `phylo_protein_12taxa`,
  `phylo_large_32taxa`, `phylo_large_64taxa`, `gauss_dense_signal_highd`,
  `gauss_dense_signal_highd_continuous`, `gauss_extreme_noise_many`. This is
  the same-data edge selection documented in the wiki's oracle gate-path
  analysis (99% edge rejection on pure Bernoulli noise).
- 3 cases (`binary_perfect_2c`, `_4c`, `_8c`) return `unsupported` because each
  cluster is copies of one row: every point has tied nearest neighbours at
  distance 0, the subtrees have zero variance, and no record has positive
  degrees of freedom.
- 2 continuous cases stop at the root because the root edge test is not
  significant at edge alpha 0.001 although the clusters are far apart
  (between/within distance 8.0 and 4.7): `gauss_clear_medium_continuous`
  (p = 0.027) and `gauss_moderate_3c_continuous` (p = 0.048).
- 4 cases get wrong sibling decisions: `gauss_single_outlier_4c_continuous`
  (root rejected, P = 0.016, 1 of 5 clusters), `gauss_outlier_cluster_4c_continuous`
  (2 of 5), `cat_mod_4cat_6c` (3 of 6) and `binary_low_noise_2c` (15 for 2;
  92.5% of its points have tied nearest neighbours).

**2. One shared subspace drowned by noise features (4).** In
`dim_consolidated_4c_24f`, `_72f`, `_272f` (median-binary) and
`dim_consolidated_4c_72f_continuous`, 12 or 13 variables separate the clusters
and the rest are noise. Over all features the points are close to equidistant
(between/within distance 1.26, 1.08, 1.02, 1.08); on the separating features it
is 1.65, 1.67, 1.60, 2.07. The oracle cut of the TBS tree falls with the number
of noise features (0.50, 0.013, -0.009; 0.0 for the continuous case), while the
held-out pair subspaces separate the pairs (ARI 0.89, 0.89, 0.86, 1.0, against
0.66, 0.03, 0.02, 0.07 for all features). At 24 features the root split is
still driven by separating features; at 72 and 272 features none of its ten
strongest features is a separating feature, so the root splits along noise. In
the continuous case the standardized-Euclidean root splits 3 against 157 and
then peels single points; a cosine distance on all features gives oracle ARI 1.0.

**3. Each pair separates in its own small subspace (4).** In
`dim_diffuse_6c_36f`, `_136f`, `_536f` and `dim_diffuse_6c_136f_continuous`
each of the 15 pairs is separated by a median of 0 to 2 held-out features, up to
7 pairs have no separating feature even when selected on all rows, and the feature sets of different
pairs barely overlap (Jaccard 0.10 to 0.19). Held out, a pair's own features
separate it better (mean ARI 0.12 to 0.38) than the union of all pairs'
features (0.02 to 0.20) or all features (-0.02 to 0.09). No single distance
puts all six clusters in one tree: what separates one pair is noise for the
others. The signal is also thin, so even the pair subspaces stay weak.

**4. Metric or linkage choice (2).** In `cat_mod_3cat_4c` all four clusters are
clades, but the Hamming average-linkage cut at K = 4 reaches ARI 0.69, while
Ward on Euclidean reaches 0.96. In `phylo_brownian_null_16taxa` all 16 taxa
are clades, but the `mahalanobis_time` root splits one point against 159; that
root edge test is not significant (p = 0.14), so every node below is blocked
and TBS returns one cluster. A Euclidean average-linkage tree cuts at ARI 1.0.

**5. Separation spread over many weak features (5).** `bar_binary_balanced_4c`,
`bar_binary_unbalanced_95_5`, `traversal_deep_signal_under_same_parent`,
`cat_highcard_10cat_4c` and `cat_highcard_20cat_4c` separate best with all
features (held-out ARI 0.80 to 0.91) and no better on any selected subspace: each
feature carries little, and many rare one-hot columns carry the rest. The TBS
tree reaches an oracle ARI of 0.55 to 0.88, but the clusters are mostly not
clean clades (clade share 0 to 0.5). Four of the five are `unsupported` with
every edge rejected; `bar_binary_balanced_4c` stops at a rejected root
(P = 0.022).

**6. No location signal in the space used (11).** No metric, linkage or
subspace tried here recovers these clusters (best oracle ARI 0.33, on the
label-selected subspace of `overlap_heavy_4c_med_feat`):

- `mp_spike_above_bbp_continuous` and `cont_lowrank_pggn_shrinkage`: the
  clusters differ in covariance, not in mean. Under `mahalanobis_time` all
  pairwise distances are equal to within a coefficient of variation of 0.003
  and 0.000, so the tree peels one point at a time.
- `cat_simplex_face_rare_20cat` and `cat_dirichlet_overdispersed_20cat`: the
  clusters differ in the shape of rare-category frequencies; 4 and 5 of 6 pairs
  have no separating column, and 37% of points have tied nearest neighbours.
- `sbm_moderate` and `sbm_hard`: no adjacency column separates any pair.
- `overlap_heavy_4c_small_feat`, `overlap_heavy_4c_med_feat`,
  `overlap_heavy_8c_large_feat` and `overlap_extreme_4c`: between/within
  distance 1.00 to 1.05; the designed overlap leaves no separable structure.
- `gauss_sparse_signal_highd_noise`: 40 samples and 20,000 features; no feature
  passes for any pair.

**7. Null cases (2).** `gauss_null_small` and `gauss_null_large` have one true
cluster. Every edge test rejects, so TBS returns `unsupported` instead of one
cluster.

## Answers

- **Subspaces:** yes for 8 cases. Four have one shared subspace that the
  full-space distance drowns (group 2); four have pair-specific subspaces that
  no single tree can serve (group 3).
- **Equidistance:** yes in three forms. Near-equidistance over all features in
  groups 2 and 3 (between/within 1.00 to 1.08; 1.26 at 24 features); exact equidistance under
  `mahalanobis_time` when the signal is in the covariance (group 6); and
  equidistance inside clusters made of identical rows, which leaves no degrees
  of freedom to calibrate (`binary_perfect_*`).
- **Metric:** decisive in 3 cases: cosine for
  `dim_consolidated_4c_72f_continuous`, Ward for `cat_mod_3cat_4c`, Euclidean
  for `phylo_brownian_null_16taxa`. Elsewhere the TBS metric is as good as any
  alternative tried.
- **Splits and dimensions:** in group 1 every root split is driven by true
  separating features; in group 2 at 72 and more features the root is driven
  only by noise features; in the covariance and rare-category cases of group 6
  and the continuous cases of groups 2 and 3 the root splits off 1 to 3 points.
- **The largest group is not geometric.** In 21 cases (23 with the null cases)
  the tree already holds the clusters and the gate fails, mostly because the
  same-data edge test rejects every edge and leaves nothing to calibrate.

## Bug check

Could the failures come from a bug rather than the data or the method? The
scripts in `bug_check/` answer this per case.

**No bug in data, distances, trees, labels or scoring.** `bug_harness.py` runs
canonical `tbs` on all 49 cases through `run_single_method_once`, as the
benchmark does, and captures what the dispatcher receives and returns. Every
case passes every check:

- the labels have the metadata's K; the data are finite; sample names are unique;
- the TBS distance equals an independent recomputation (Hamming exactly,
  standardized Euclidean to 6e-14, `mahalanobis_time` to 4e-15);
- the tree's leaves are exactly the samples, and every split is binary;
- each returned cluster is exactly the leaf set of one decomposition cluster
  root, and the ARI recomputes to the reported value;
- the `unsupported` support count recomputes from the annotations.

**Three method defects, each confirmed by an experiment.**

1. **The edge test fits its projection basis to the leaves it tests.** The
   child-parent statistic projects the whitened child-parent difference onto
   the top principal components of the parent's own leaves and compares it with
   chi-square on k degrees of freedom. The kernel's docstring says that law holds
   for a fixed basis. A basis fitted to the same leaves puts the difference into
   its first components: for a parent with two leaves, the first component *is*
   the difference. On the per-feature-variance path that makes the statistic
   exactly p (2000.00 on 1 df for 58 edges of `cont_lowrank_pggn_shrinkage`). On
   pure-null data (`edge_basis_experiment.py`, all 118 edges, no ancestor
   blocking):

   | Basis | Tree built from the data | Random tree |
   | --- | ---: | ---: |
   | Fitted to the same leaves (production) | 98-100% rejected | 73-80% rejected |
   | Fitted to an independent null replicate | 0-3% rejected | 0-2% rejected |

   The statistic averages 17 to 19 times its degrees of freedom with the
   production basis and 0.35 to 0.56 times with an independent one. The random
   tree shows that tree selection is not the main cause. As a result, 26 cases
   reject every edge, leave no null-like record to calibrate and return
   `unsupported`, both null cases included; 8 more reject at least 95% of
   internal edges.
2. **The continuous edge statistic is capped at n_parent - 1.** When full
   covariances are stored, the difference is whitened by the parent's own sample
   covariance, which contains the split under test, so the statistic cannot
   exceed n_parent - 1. Root statistics sit at the cap (58.98 for 59, 43.98 for
   44, 119.00 for 119, 159.7 for 160, 164.4 for 165). With df = p, the 0.999
   cutoff lies above the cap in `gauss_clear_medium_continuous` (73.4 > 59),
   `gauss_moderate_3c_continuous` (59.7 > 44), `mp_spike_above_bbp_continuous`
   (172.4 > 119) and `dim_diffuse_6c_136f_continuous` (192.7 > 179): their root
   can never reject. Above a memory limit only per-feature variances are stored
   and defect 1 applies instead, so the edge test's behaviour switches with a
   storage threshold.
3. **`mahalanobis_time` removes mean separation.** It whitens by the covariance
   of all samples, including the between-cluster spread, with 5% shrinkage. In
   `mahalanobis_separation_check.py` two clusters shifted 3 standard deviations
   on every feature (n = 120) give Euclidean ARI 1.0 for every p and
   `mahalanobis_time` ARI 0.00 for p = 40 to 300. It is the distance of
   `mp_spike_above_bbp_continuous`, `cont_lowrank_pggn_shrinkage` and
   `phylo_brownian_null_16taxa`.

**Per-case verdicts** (`bug_verdicts.csv`):

| Verdict | Cases |
| --- | ---: |
| A defect is the cause: the tree already holds the clusters | 15 |
| A defect is the cause: null data should give one cluster | 2 |
| A defect is present; the data or geometry would still limit recovery | 22 |
| Edge over-rejection possibly contributes (80-90% of internal edges) | 2 |
| No defect implicated | 8 |

The 17 caused cases are the 12 all-edges-rejected cases with a correct tree,
`cat_mod_4cat_6c`, the two capped continuous Gaussian cases and the two null
cases. The 8 without a defect are the three `binary_perfect` cases (identical
rows leave zero variance; failing closed is the designed outcome), the two
Gaussian outlier cases and `dim_consolidated_4c_72f_continuous` (their root edge
can reject), and `sbm_hard` and `sbm_moderate`.

**Side findings.** `tests/integration/60_test_benchmark_methods_smoke.py` calls
`benchmark_cluster_algorithm` without an output folder, so each `make check`
writes `case_1_tbs_diffusion_*` audit files into `benchmarks/results/audit/`.
The two reruns for this report did the same and overwrote 20
`case_*_tbs_stats.csv` files there (gitignored; the 2026-10-16 run folder is
untouched). `bug_harness.py` calls the runner directly and writes nothing there.

## Limits

- The Welch test sees mean differences only; covariance and frequency-shape
  differences (group 6) have no separating feature by construction.
- Oracle cuts use the true K, so they bound what a tree can give; they are not
  TBS results.
- The audit columns come from the 2026-10-16 run; `unsupported` cases have no
  audit file, so their gate evidence comes from the rerun's evidence counts.

## Reproduce

Run from the repository root with the local
`benchmarks/results/run_20261016_163716Z_full/` run present:

```bash
uv run python reports/failing_case_subspace_analysis_20261018/rerun_failing_cases.py
uv run python reports/failing_case_subspace_analysis_20261018/subspace_failure_analysis.py
uv run python reports/failing_case_subspace_analysis_20261018/pair_subspace_check.py
```

For the bug check, run `bug_check/bug_harness.py`, then `bug_check/bug_verdicts.py`,
`bug_check/edge_basis_experiment.py` and `bug_check/mahalanobis_separation_check.py`
the same way; each writes its CSV next to itself.

The three analysis scripts write `tbs_rerun_failing_cases.csv`,
`failing_case_subspace_summary.csv` with `failing_case_top_splits.csv`, and
`pair_subspace_heldout.csv`. A second run reproduced all three tables byte for
byte.
