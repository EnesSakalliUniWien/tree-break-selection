# Split-count rule-out on the continuous benchmark suite

Run on 2026-11-08 (user choice: one comparison run, benchmark defaults unchanged). Local folder:
`results/analyses/rule_out_continuous_benchmark_20261108` (git-ignored); protocol, runner and the result rows in
`rule-out-continuous-benchmark-scripts-20261108.py.txt`.

## Setup

The 9 cases of `get_test_cases_by_suite("continuous")` (the other 117 of 126 benchmark cases are binary or
categorical, where the rule-out does not apply) with `tbs_family_diffusion` and
`tbs_family_diffusion_markov_ward`, each in its registered parameter cell and the same cell with
`sibling_split_count_rule_out=True`, through the shared `benchmark_cluster_algorithm` (default alphas). All 36
runs finished with status ok. Before the run the module was changed to measure the equidistant share only for
splits that pass the separability check (otherwise the 20,000-feature case would need a 20,000 x 20,000
matrix); unit and end-to-end tests pass after the change.

## Results (found clusters, ARI)

| Case | Cells | Features | True groups | Default today | Markov/Ward today | Rule-out (both trees) |
| --- | ---: | ---: | ---: | --- | --- | --- |
| dim_consolidated_4c_24f | 160 | 24 | 4 | 4, 1.000 | 5, 0.920 | 4, 1.000 |
| dim_consolidated_4c_72f | 160 | 72 | 4 | 5, 0.913 | 6, 0.824 | 4, 1.000 |
| dim_diffuse_6c_136f | 180 | 136 | 6 | 5, 0.017 | 4, 0.006 | 1, 0.000 |
| gauss_clear_medium | 60 | 40 | 4 | 4, 1.000 | 4, 1.000 | 3, 0.704 |
| gauss_dense_signal_highd | 40 | 20,000 | 4 | 4, 1.000 | 4, 0.932 | 1, 0.000 |
| gauss_moderate_3c | 45 | 30 | 3 | 3, 1.000 | 3, 1.000 | 2, 0.560 |
| gauss_null_large | 200 | 80 | 1 | 1, 1.000 | 1, 1.000 | 1, 1.000 |
| gauss_outlier_cluster_4c | 166 | 40 | 5 | 5, 1.000 | 5, 1.000 | 5, 1.000 |
| gauss_single_outlier_4c | 161 | 40 | 5 | 4, 0.992 | 4, 0.992 | 4, 0.992 |

Mean ARI: default 0.880 today, 0.695 rule-out; Markov/Ward 0.852 today, 0.695 rule-out. Cases with the exact
number of clusters: 6 and 5 today, 4 with rule-out.

## Reading

The rule-out removes the extra splits where today's test over-splits (dim_consolidated 72 features: 5 and 6
clusters to 4, ARI to 1.000) and keeps the null case at one cluster, but it loses real groups where the
subtree to be split has p + 1 cells or fewer (60 cells in 40 features, 45 in 30, 40 in 20,000): such a subtree
has no invariant directions, so its split is refused. The protocol expected at most the first splits to pass;
that was too pessimistic, because the limit applies to the subtree being split (an 80-cell subtree in 72
features can still be split into two 40-cell groups). On these small, wide benchmark cases the rule-out is
worse than today on average.
