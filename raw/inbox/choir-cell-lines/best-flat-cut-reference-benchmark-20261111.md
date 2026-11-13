# Gaussian best-flat-cut reference on the continuous and cancer_pc benchmark suites

Run on 2026-11-11, Europe/Vienna. Local folder:
`results/analyses/best_flat_cut_reference_benchmark_20261111` (git-ignored); protocol, runner, result rows and
the per-split reasons are in `best-flat-cut-reference-benchmark-scripts-20261111.py.txt`. The opt-in
`sibling_best_flat_cut_reference` was uncommitted working-tree code on top of `bc044e07`. No production code
was changed by this run and nothing was promoted.

## Setup

The 9 cases of `get_test_cases_by_suite("continuous")` and the 4 cases of
`get_test_cases_by_suite("cancer_pc")` with `tbs_family_diffusion` and `tbs_family_diffusion_markov_ward`, each
in its registered parameter cell ("today") and the same cell with `sibling_best_flat_cut_reference=True`
("reference"), through the shared `benchmark_cluster_algorithm` (default alphas, plots off). All 52 runs
finished with status ok. The protocol was written before the run and names the expected limit: a subtree with
no invariant directions is not split. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

A second pass (`reasons.py`) replayed the 26 reference runs one case and method at a time and recorded every
call of `best_flat_cut_reference` (subtree size, reason, k, `rho`, p-value). Its found clusters and ARI equal
the first run in 26 of 26.

The "today" cells reproduce earlier runs: the continuous cells equal the 8 November run row by row, and the
default tree's mean ARI on cancer_pc is 0.8218 as on 10 November.

## Results (found clusters, ARI)

| Case | Cells | Features | True groups | Default today | Default reference | Markov/Ward today | Markov/Ward reference |
| --- | ---: | ---: | ---: | --- | --- | --- | --- |
| gauss_clear_medium | 60 | 40 | 4 | 4, 1.000 | 3, 0.704 | 4, 1.000 | 3, 0.704 |
| gauss_moderate_3c | 45 | 30 | 3 | 3, 1.000 | 2, 0.560 | 3, 1.000 | 2, 0.560 |
| gauss_null_large | 200 | 80 | 1 | 1, 1.000 | 1, 1.000 | 1, 1.000 | 1, 1.000 |
| gauss_dense_signal_highd | 40 | 20,000 | 4 | 4, 1.000 | 1, 0.000 | 4, 0.932 | 1, 0.000 |
| dim_consolidated_4c_24f | 160 | 24 | 4 | 4, 1.000 | 4, 1.000 | 5, 0.920 | 4, 1.000 |
| dim_consolidated_4c_72f | 160 | 72 | 4 | 5, 0.913 | 4, 1.000 | 6, 0.824 | 4, 1.000 |
| dim_diffuse_6c_136f | 180 | 136 | 6 | 5, 0.017 | 1, 0.000 | 4, 0.006 | 1, 0.000 |
| gauss_single_outlier_4c | 161 | 40 | 5 | 4, 0.992 | 4, 0.992 | 4, 0.992 | 4, 0.992 |
| gauss_outlier_cluster_4c | 166 | 40 | 5 | 5, 1.000 | 5, 1.000 | 5, 1.000 | 5, 1.000 |
| cancer_pc_hub_disconnected | 496 | 30 | 6 | 6, 0.998 | 6, 0.998 | 10, 0.446 | 6, 0.998 |
| cancer_pc_hub_connected | 496 | 30 | 6 | 5, 0.593 | 4, 0.446 | 13, 0.265 | 4, 0.446 |
| cancer_pc_nohub_disconnected | 496 | 30 | 6 | 6, 0.998 | 6, 0.998 | 16, 0.268 | 6, 0.998 |
| cancer_pc_nohub_connected | 496 | 30 | 6 | 6, 0.697 | 4, 0.446 | 15, 0.234 | 4, 0.446 |

| Suite | Tree | Mean ARI today | Mean ARI reference | Exact count today | Exact count reference |
| --- | --- | ---: | ---: | ---: | ---: |
| continuous (9) | default | 0.8802 | 0.6950 | 6 | 4 |
| continuous (9) | Markov/Ward | 0.8525 | 0.6950 | 5 | 4 |
| cancer_pc (4) | default | 0.8218 | 0.7221 | 3 | 2 |
| cancer_pc (4) | Markov/Ward | 0.3033 | 0.7221 | 0 | 2 |

On the continuous suite the reference gives the same found clusters and ARI as the split-count rule-out
(8 November) and the ceiling-3 rule (9 November) in every case. With the reference both trees give the same
found clusters and the same ARI to three decimals in all 13 cases.

## Per-split reasons (reference cell, 151 calls)

93 calls were `tested` and 58 were `no_invariant_directions`; `too_few_cells` and
`too_many_invariant_directions` did not occur. Of the 93 tested p-values, 56 were 0.0005 (the lowest possible
value), 36 were 1.0 (`rho` below the stored upper 10% of the null) and one was 0.056 (41 cells, k = 2,
`rho` = 3.42, cancer_pc_nohub_connected on Markov/Ward). A first version of this file said that every tested
p-value was 0.0005 or 1.0; that sentence was written from the eight largest subtrees per run and was wrong.

- Losses in the three small, wide Gaussian cases come from `no_invariant_directions`: the last real split has
  30 cells in 40 features (gauss_clear_medium), 30 cells in 30 features (gauss_moderate_3c) and 40 cells in
  20,000 features (gauss_dense_signal_highd, root). The splits above them pass with `rho` of 549 to 1,069.
- dim_diffuse_6c_136f: the root is tested with k = 3 and `rho` = 1.61 (p = 1.0); all 6 lower subtrees have
  117 cells or fewer in 136 features and no invariant directions.
- dim_consolidated: every real split passes (`rho` 12.3 to 23.3, k = 2 to 5); the extra splits of today's test
  sit in subtrees of 40 cells in 72 features (no invariant directions) or get p = 1.0.
- gauss_null_large: the root is tested with k = 1, `rho` = 0.03, p = 1.0.
- cancer_pc, disconnected: real splits pass with `rho` 2.64 to 55.8; the Markov/Ward tree's extra splits are
  refused, as `no_invariant_directions` at 157 to 256 cells in 30 features (no eigenvalue above the noise band,
  not the n <= p + 1 limit) or as tested with `rho` below 0.5.
- cancer_pc, connected: two splits pass (`rho` 2.64 and 5.36); the remaining tested splits have k = 1 and `rho`
  of 0.15 to 1.03 (p = 1.0), so two true groups are not separated. Today's default tree finds 5 and 6 clusters
  there (ARI 0.593 and 0.697).

## Do the refused splits follow a true group boundary? (cancer_pc, follow-up on the same day)

`boundaries.py` replayed the 8 cancer_pc reference runs (found clusters and ARI equal the first run in 8 of 8)
and recorded, for every call, the true-group counts on both sides and, for each true group inside the subtree,
the `rho` of "this group against the rest" in the same invariant directions and in all 30 features.

The generator puts groups 0 to 3 (256, 64, 32, 16 cells) at 0, 6, 12, 18 on PC1 and groups 4 and 5 (64 cells
each) at -offset and +offset on PC2, with unit noise; offset is 20 in the disconnected and 3 in the connected
cases.

Connected cases, subtree holding groups 0, 4 and 5 (382 or 384 cells, k = 1). The table's 99% point at the
grid n = 300, k = 1 is 2.151 and its lowest stored value (90% point) is 1.971.

| Case, tree | Tree split (true-group counts) | `rho` of the tree split | Group 4 against rest, invariant / all 30 features | Group 5 against rest, invariant / all 30 features |
| --- | --- | ---: | --- | --- |
| hub, default | 46 of group 5 + 7 of group 0 against 248 + 63 + 18 | 0.76 | 0.81 / 0.92 | 0.88 / 0.98 |
| hub, Markov/Ward | 43 of group 5 + 4 of group 0 against 252 + 64 + 21 | 0.70 | 0.82 / 0.93 | 0.88 / 0.98 |
| nohub, default | 42 of group 5 + 6 of group 0 against 250 + 64 + 22 | 0.64 | 0.82 / 0.93 | 0.87 / 0.97 |
| nohub, Markov/Ward | 52 of group 5 + 9 of group 0 against 247 + 64 + 12 | 0.93 | 0.82 / 0.93 | 0.87 / 0.97 |

- The tree split does not follow the true boundary: it isolates 42 to 52 of the 64 cells of group 5 together
  with 4 to 9 cells of group 0. The next split isolates 38 to 49 of the 64 cells of group 4 (`rho` 0.95 to 1.03).
- The exact true-group split would be refused as well: `rho` is 0.81 to 0.88 in the invariant direction and
  0.92 to 0.98 in all 30 features, against a 99% point of 2.151. The invariant direction keeps most of the
  all-feature value, so the directions are not the cause at this level.
- Arithmetic for the cause: with `rho = pi(1-pi) D^2`, a group of 64 next to one of 256 (pi = 0.2) at offset 3
  has population `rho` = 0.16 x 9 = 1.44; a second group on the opposite side enlarges the pooled within
  variance and lowers it further. The median cut of one Gaussian population has population `rho` =
  2 / (pi - 2) = 1.752. The 99% point 2.151 needs an offset of 3.67 for sizes 64 and 256 and 2.93 for equal
  sizes.
- Disconnected cases: on the Markov/Ward tree the root split is exactly group 4 (64 cells, offset 20) against
  the other 432 cells and is refused (`rho` = 1.51, p = 1.0), because the rest holds group 5 at the opposite
  offset. On the default tree the same happens one level down for group 5 against the rest (`rho` = 1.54,
  p = 1.0). The next split (the other offset group against the rest) passes with `rho` of 49.4 to 55.8, and
  the pass-through traversal then separates both groups, so the final partition is correct (ARI 0.998).
- Deeper subtrees of the connected cases (41 to 96 cells): the true-group `rho` is 0.00 to 0.11 in the
  invariant directions and 0.88 to 9.14 in all 30 features. With so few cells per feature the all-feature value
  is inflated by noise, so this does not show that the directions miss the difference.
- Density along the offset axis (calculation, not a simulation): for 256 cells at 0 and 64 at 3 the mixture
  density has a valley at 2.49 whose height is 0.998 of the neighbouring mode; for the connected layout (256 at
  0, 64 at -3, 64 at 3) the two valleys are at 0.999, the skewness is 0 and the kurtosis is 3.000. At offset 2.5
  or 2 the two-group density has one mode. For equal sizes at offset 3 the valley is at 0.642. So the
  connected layout has no gap in its density; a measure of the gap between the two sides has almost nothing to
  detect there, while it would detect the disconnected layout.

Answer: both. The refused tree splits are impure, and the pure true-group splits would be refused too. A group
is accepted only if its `rho` against everything else in the subtree exceeds what one Gaussian population's
best flat cut reaches; an unequal-sized group at offset 3, or a group whose complement holds another distant
group, stays below that.

## Reading

- The reference removes over-splitting: dim_consolidated to the exact 4 clusters, and the Markov/Ward tree on
  cancer_pc from 10 to 16 clusters to 6, 4, 6, 4 (mean ARI 0.303 to 0.722).
- It loses real groups in two ways. (1) A subtree with no invariant directions gets p = 1; this is the limit
  the protocol named and is identical to the rule-out and the ceiling-3 rule, because all three use
  `invariant_directions`. (2) In the connected cancer_pc cases the split is tested but its `rho` inside the
  one invariant direction is below the null's upper 10%.
- The losses are located in the choice of directions, not in the reference: no accepted or refused split
  changed relative to the rule-out on the continuous suite.

## Limits

- One run, one seed per case, 13 synthetic cases; no real data.
- The true-group check below covers the cancer_pc cases only and compares one group against the rest; it
  does not search other splits of a subtree.
- The other 113 benchmark cases are binary or categorical and were not run; the reference is derived for
  Gaussian coordinates.
- Only one tested split had a p-value between 0.0005 and 1.0, so this run barely exercises intermediate
  p-values or the rounding of n to the table grid.
