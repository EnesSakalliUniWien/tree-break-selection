# Prototype of k with a margin at the edge, also for more features than cells

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/k_margin_prototype_20261112` (git-ignored);
the protocol, the deviations of both parts, the rule, the scripts and their outputs are in
`k-margin-prototype-scripts-20261112.py.txt`. Follow-up of `marchenko-pastur-multiplicity-literature-20261112.md`
and step A2 of `plan-non-heuristical-statistic-20261110.md`. Working tree on top of `bc044e07` with the
uncommitted opt-in reference. The rule exists in the analysis folder only; no production code was changed.
Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## The rule (`k_rule.py`, frozen in the protocol before any run)

For n cells in p features: lambda_1 >= ... >= lambda_r are the nonzero eigenvalues of the covariance
(r = min(p, n - 1); from the n x n Gram matrix when p > n - 1). Step j = 1, 2, ...:
T_j = lambda_j / median(lambda_j, ..., lambda_r) is compared with the 99% point of T_1 for one spherical Gaussian
population of n cells in p - j + 1 features (2,000 simulated draws from a Wishart matrix, seed from (n, p - j + 1),
cached). k is the number of steps passed before the first failure. No noise variance is estimated. The returned
rows are the k leading eigenvectors, the same contract as `invariant_directions`.

The work was split in two parts against this frozen interface, built separately and checked before use: Part A
(the rule and synthetic checks) and Part B (a harness that plugs a rule into `best_flat_cut_reference` at run
time, validated with the current rule). The integration was run after both were checked.

## Part A: synthetic checks (Gaussian data, identity noise; current rule beside prototype rule)

Self-checks passed: the 99% point is reproducible; the Wishart route equals drawing cells (two-sample
Kolmogorov-Smirnov p 0.33, 0.67, 0.44 at (64, 30), (30, 40), (40, 2000)); the rows are orthonormal leading
eigenvectors; k is unchanged by a shift and by scaling.

Level, one Gaussian population, 2,000 draws, draws with k >= 1:

| Cells, features | Current rule | Prototype rule |
| --- | ---: | ---: |
| 64, 30 | 472 | 20 |
| 128, 30 | 350 | 18 |
| 256, 30 | 311 | 17 |
| 800, 30 | 291 | 19 |
| 200, 80 | 433 | 10 |
| 30, 40 | 0 | 11 |
| 40, 2000 | 0 | 19 |

Power, 1,000 draws, draws with k >= 1:

| Setting | Population `rho` | Current rule | Prototype rule |
| --- | ---: | ---: | ---: |
| 64 and 256 cells, 30 features (sqrt(p/n) = 0.306) | 0.15 | 194 | 16 |
| | 0.306 | 358 | 44 |
| | 0.6 | 901 | 568 |
| | 1.0 | 1,000 | 997 |
| | 1.44 | 1,000 | 1,000 |
| 15 and 15 cells, 40 features | 2.3 | 0 | 234 |
| | 5 | 0 | 972 |
| | 50 | 0 | 1,000 |
| 20 and 20 cells, 2,000 features | 10 | 0 | 263 |
| | 50 | 0 | 1,000 |
| | 200 | 0 | 1,000 |

Count of directions, 1,000 draws: four balanced groups (three true directions), 400 cells, 30 features: k = 3 in
997 (current) and 990 (prototype); two groups, 64 and 256 cells, `rho` 1.44: k >= 2 in 57 (current) and 14
(prototype). One population with variances 4, 2, 1, ..., 1: k = 2 in 970 (current) and 943 (prototype), as the
rule's stated limit says. Cost: 0.04 to 0.10 s for a new pair of cells and features, about 1 ms afterwards, also
at (48879, 30) and (40, 20000).

All seven predictions of the protocol for Part A held.

Independent spot check (`verify_rule.py`, other seed, own draws of cells): own draws above the rule's 99% point
in 47, 23 and 45 of 4,000 at (128, 30), (30, 40), (500, 12); k >= 1 in 35, 27, 21 and 24 of 2,000 at (100, 30),
(300, 30), (45, 30), (20, 60); the first direction equals the current rule's up to sign in 200 of 200; with 15
and 15 cells in 40 features the direction lies on the separating feature (component 0.976 to 0.987).

## Part B: the harness, validated with the current rule

With the current rule the harness reproduces 11 November in 26 of 26 rows (found clusters and ARI), with 58
`no_invariant_directions` among 151 calls, and gives 1 cluster in 80 of 80 runs of the one-population grid.

## Integration: the opt-in reference with the prototype rule

14 cases (10 continuous, 4 `cancer_pc`), both registered methods, `sibling_best_flat_cut_reference=True`.
Found clusters and ARI, current rule to prototype rule; cases not listed are unchanged in both methods:

| Case (true groups) | Default tree | Markov/Ward tree |
| --- | --- | --- |
| gauss_clear_medium (4) | 3, 0.704 to 4, 1.000 | 3, 0.704 to 4, 1.000 |
| gauss_moderate_3c (3) | 2, 0.560 to 3, 1.000 | 2, 0.560 to 3, 1.000 |
| gauss_dense_signal_highd (4) | 1, 0.000 to 4, 1.000 | 1, 0.000 to 3, 0.698 |

| | Current rule | Prototype rule |
| --- | ---: | ---: |
| Continuous, 9 cases of 11 November, mean ARI, default tree | 0.6950 | 0.8880 |
| Continuous, 9 cases, mean ARI, Markov/Ward tree | 0.6950 | 0.8544 |
| Continuous, 9 cases, exact number of clusters, default / Markov/Ward | 4 / 4 | 7 / 6 |
| cancer_pc, 4 cases, mean ARI, both trees | 0.7221 | 0.7221 |
| `no_invariant_directions` among the 151 calls of the 13 old cases | 58 | 44 |
| `no_invariant_directions` among the 165 calls of all 14 cases | 66 | 54 |
| Null cases `gauss_null_large` and `gauss_null_256n_30f` | 1 cluster | 1 cluster |
| One-population grid (4 sizes x 5 seeds x 4 cells), runs with 1 cluster | 80 of 80 | 80 of 80 |
| Smallest reference p-value among the tested calls of that grid | 1.0 (279 calls) | 0.0185 (222 calls) |

For comparison, today's default test on the same 9 continuous cases: mean ARI 0.8802 and 0.8525 (11 November).

Remaining losses with the prototype rule are refusals by the reference, not missing directions: dim_diffuse_6c_136f
(root tested with k = 3, `rho` 1.61, p = 1.0), the connected `cancer_pc` cases (`rho` 0.76, p = 1.0),
gauss_single_outlier_4c (the single cell), and one split of gauss_dense_signal_highd on the Markov/Ward tree
(20 cells split 9 and 11, `rho` 4.19, p = 0.014).

Selected subtrees of one Gaussian population (496 cells, 30 features, seeds 401 to 405, nodes with at least 40
cells), nodes with k >= 1:

| Operator/linkage | Nodes | Current rule | Prototype rule | Largest k, current / prototype |
| --- | ---: | ---: | ---: | --- |
| symmetric/average | 132 | 89 (67.4%) | 53 (40.2%) | 3 / 5 |
| symmetric/ward | 96 | 91 (94.8%) | 70 (72.9%) | 4 / 6 |
| markov/average | 112 | 94 (83.9%) | 63 (56.2%) | 4 / 4 |
| markov/ward | 92 | 86 (93.5%) | 69 (75.0%) | 4 / 4 |

Predictions of the protocol: P1 held (three cases find more clusters, no case finds more than its true groups);
P2 held (null cases and grid at 1 cluster); P3 held (58 to 44); P4 held (share with k >= 1 above 0.01 for both
rules and lower for the prototype rule in all four cells).

## Reading

- The margin and the Gram route do what they were built for: the level for an independent sample is about 1%
  instead of 16% to 23%, directions are found with more features than cells, and the opt-in recovers the three
  cases it lost for that reason while both null cases and the grid stay at 1 cluster.
- The margin costs power near the threshold with more cells than features (`rho` 0.6: 568 against 901 of
  1,000).
- The margin does not address selection: in subtrees the tree chose from one Gaussian population, 40% to 75%
  of nodes still have a direction against a nominal 1%. The directions of a subtree are shaped by the tree's
  choice of cells; that part of step A2 is open.
- With the prototype rule the reference is asked more often in the one-population grid's small subtrees and
  its smallest p-value is 0.0185, where the current rule's is 1.0. No run split, but the distance to the level
  0.01 is smaller.

## Deviations

Part A: choices the protocol left open (Bartlett decomposition for the Wishart draw; "within Monte Carlo error"
defined as Kolmogorov-Smirnov p above 0.001; "near sqrt(p/n)" judged at `rho` 0.306); no assert or prediction
changed. Part B: the harness ran once with the prototype rule before the rule was checked, because the rule file
already existed; those outputs were set aside (`premature/`) and the integration was rerun after the check. The
rerun equals them in found clusters, ARI, reasons, k and p-values; `rho` differs by at most 2.8e-12 relative.

## Not checked

The realized level depends on the simulated 99% point with 2,000 draws (0.5% to 1.75% in the 11 sizes
checked). At step 2 and later with more features than cells the observed median runs over fewer eigenvalues
than the null's; not measured apart (k >= 2 in 0 of 4,000 null draws there). Rank-deficient or duplicate cells,
`alpha` other than 0.01 and `max_directions` were not exercised. One population with unequal variances has
directions by definition. Identity noise in all synthetic checks; the CHOIR cells, the reference's own level
with this k, and the directions being chosen from the tested cells were not examined. One seed per benchmark
case.
