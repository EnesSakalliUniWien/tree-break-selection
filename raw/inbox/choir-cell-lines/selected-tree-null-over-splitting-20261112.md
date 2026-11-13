# Systematic debugging of the accepted extra splits: the tree is built from the tested cells

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/cancer_pc_systematic_debugging_20261112`
(git-ignored); protocol with its amendment, the three scripts and their complete outputs are in
`selected-tree-null-over-splitting-scripts-20261112.py.txt`. Follow-up of `cancer-pc-over-splits-20261112.md`,
`cancer-pc-operator-linkage-20261112.md`, `cancer-pc-edge-test-20261112.md` and `cancer-pc-edge-k-20261112.md`.
Working tree on top of `bc044e07` with the uncommitted opt-in reference. No production code was changed and no
fix was made. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Symptom

Today's test splits one true group (one Gaussian population) on the `cancer_pc` cases with operator "markov"
and linkage "ward": 13 accepted splits with all cells from one true group, 0 in the other three cells.

## Phase 1: root cause investigation (no new run)

- Reproduction: deterministic; the four earlier folders reproduce found clusters and ARI.
- Backward trace: accepted split <- `Sibling_BH_Different` <- sibling p-value present <- edge significant for a
  child <- edge reached by Tree-BH <- edge raw p <= 0.001 <- statistic in k directions <- k from the whitening
  of the basis rows and `extend_rank_through_ties` <- split chosen by the tree from the same cells.
- Recent changes: `extend_rank_through_ties` entered with `bb96e2b5` (2026-10-14); the basis half, the test half
  and the within-child whitening of the basis rows entered with `9c5cd805` (2026-10-18). That commit's message
  says a basis fitted to the tested leaves rejected nearly every null edge.

## Phase 2: pattern analysis (no new run)

- Working reference: `tests/statistics/58_test_edge_null_calibration.py` (5 tests pass today). Its two null
  tests build the tree from independent noise, `pdist(rng.normal(size=(n_samples, 5)))`; the docstring says
  the excess rejection then does not come "from selecting the tree on the tested data".
- No test builds the tree from the tested data. The docstring of `annotate_child_parent_divergence` says the
  row split "does not make child membership independent when the tree was selected using these data" and that
  the chi-square reference "is not a calibrated selected-tree reference".
- The registered null case `gauss_null_large_continuous` has 200 cells in 80 features, seed 401.
- Difference between working and broken: where the tree comes from.

## Phase 3: one hypothesis per experiment, one variable changed

### E1. H1: the edge test over-rejects because the tree is built from the tested data

Setup of the existing test `test_edge_test_does_not_over_reject_on_iid_gaussian_null` (240 cells, one Gaussian
population, alpha 0.05, 40 datasets, seed 2, same data in both arms). One variable: the distance the tree is
built from. Raw edge p-values.

| Features | Linkage | Tree from | All edges, p < 0.05 | Parent of at least 40 cells, p < 0.05 | Parent of at least 40 cells, p <= 0.001 |
| ---: | --- | --- | ---: | ---: | ---: |
| 5 | ward | independent noise | 106 of 7,672 (1.4%) | 56 of 682 (8.2%) | 4 of 682 (0.6%) |
| 5 | ward | the data | 1,326 of 7,684 (17.3%) | 678 of 698 (97.1%) | 658 of 698 (94.3%) |
| 5 | average | independent noise | 96 of 7,990 (1.2%) | 54 of 1,138 (4.7%) | 0 of 1,138 |
| 5 | average | the data | 1,340 of 7,936 (16.9%) | 810 of 1,164 (69.6%) | 600 of 1,164 (51.5%) |
| 30 | ward | independent noise | 218 of 7,710 (2.8%) | 76 of 684 (11.1%) | 6 of 684 (0.9%) |
| 30 | ward | the data | 2,032 of 7,782 (26.1%) | 394 of 664 (59.3%) | 368 of 664 (55.4%) |
| 30 | average | independent noise | 208 of 7,994 (2.6%) | 50 of 1,144 (4.4%) | 4 of 1,144 (0.3%) |
| 30 | average | the data | 1,202 of 7,930 (15.2%) | 186 of 2,292 (8.1%) | 172 of 2,292 (7.5%) |

H1 is supported in all four settings: with the tree from independent noise the test's bound of 0.10 holds on
all edges (1.2% to 2.8%); with the tree from the data the share is 15.2% to 26.1%, and 7.5% to 94.3% at 0.001
for parents of at least 40 cells. Self-check: the arm of the existing test (5 features, ward, independent)
gives 1.4%, within the test's bound.

Side observation: with independent noise, ward and 30 features, parents of at least 40 cells reject at 11.1%
at 0.05 and 0.9% at 0.001, above the level; this is the fixed-split excess with the full covariance already
seen in `cancer-pc-edge-test-20261112.md`.

### E2. H2: the tie extension is necessary for that over-rejection

Tree from the data, ward. One variable: `extend_rank_through_ties` replaced by no extension while the edge
test's basis is computed.

| Features | Tie extension | k at parents of at least 40 cells | All edges, p < 0.05 | Parent of at least 40 cells, p < 0.05 | p <= 0.001 |
| ---: | --- | --- | ---: | ---: | ---: |
| 5 | on | 5 at 604, 2 at 94 | 1,326 of 7,684 | 678 of 698 (97.1%) | 658 of 698 (94.3%) |
| 5 | off | 2 at 698 | 1,326 of 7,684 | 678 of 698 (97.1%) | 658 of 698 (94.3%) |
| 30 | on | 2 to 30 | 2,032 of 7,782 (26.1%) | 394 of 664 (59.3%) | 368 of 664 (55.4%) |
| 30 | off | 2 at 664 | 1,978 of 7,782 (25.4%) | 376 of 664 (56.6%) | 350 of 664 (52.7%) |

H2 is not supported. Without the extension k is 2 everywhere and the rejections stay (identical counts at 5
features, 59.3% to 56.6% at 30 features). The k rule decides which nodes reject, not whether the test is
valid: the first basis direction is the between-child direction of the basis half, which is the tree's own
split.

### E3. H3: the symptom needs none of the cancer_pc structure

Today's complete method on one Gaussian population (`blobs_continuous`, 1 cluster, the registered null case
with other sizes), seeds 401 to 405, default alphas. Found clusters per seed:

| Cells, features | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | --- | --- | --- | --- |
| 200, 80 (registered) | 1, 3, 1, 1, 1 | 1, 2, 1, 1, 1 | 1, 1, 1, 1, 1 | 1, 1, 1, 1, 2 |
| 200, 30 | 7, 2, 6, 6, 1 | 5, 5, 5, 4, 6 | 4, 1, 1, 5, 1 | 4, 3, 6, 4, 3 |
| 256, 30 | 10, 5, 8, 2, 9 | 3, 6, 3, 6, 11 | 1, 5, 8, 1, 1 | 4, 7, 4, 4, 7 |
| 496, 30 | 1, 1, 1, 1, 1 | 12, 8, 15, 11, 10 | 9, 1, 11, 1, 7 | 14, 12, 10, 12, 12 |

All 80 runs finished ok. More than 1 cluster: 3 of 20 runs at 80 features and 46 of 60 at 30 features
(symmetric/average 9 of 15, symmetric/ward 15 of 15, markov/average 7 of 15, markov/ward 15 of 15).

H3 is supported, and the failure is wider than the symptom: with 30 features one Gaussian population is split
in all four cells, including the default cell at 200 and 256 cells. The registered null case (200 cells, 80
features, seed 401) gives 1 cluster in all four cells and does not show it.

## Root cause

The edge test and the sibling test are valid for a split that does not depend on the tested cells, and the
repository's null tests check exactly that. In the method the tree is built from the tested cells, so every
split it offers is a selected split. Operator, linkage, the whitening of the basis rows, k and Tree-BH decide
where the selected splits pass; none of them makes the tests valid for a selected split.

## Phase 4: verification of an existing candidate (amendment written before the run, no new code)

E4. The E3 grid with the opt-in `sibling_best_flat_cut_reference=True` (the split is compared with the Gaussian
null of the best flat cut; the edge gate stays as today): 1 cluster in 80 of 80 runs, as predicted.

No failing test was added to the test suite and no default was changed. Known cost of the opt-in on the
benchmark (`best-flat-cut-reference-benchmark-20261111.md`): real groups are lost where a subtree has no
invariant directions and in the connected `cancer_pc` cases.

## Follow-up on the same day: benchmark null case added

`gauss_null_256n_30f_continuous` was added to `benchmarks/shared/cases/gaussian.py` (`blobs_continuous`, 256
cells, 30 features, 1 cluster, `cluster_std` 1.0, seed 401, no discretized baseline); the continuous suite has
10 cases. `tests/core/test_benchmark_case_data_contract.py`, `benchmarks/README.md` and `CHANGELOG.md` follow.
158 tests pass in `tests/core` and the other test files that read the case lists; ruff passes.

The registered case, run by name (`new_case.py`):

| Method | Today: found clusters, ARI | With `sibling_best_flat_cut_reference=True` |
| --- | --- | --- |
| `tbs_family_diffusion` | 10, 0.0 | 1, 1.0 |
| `tbs_family_diffusion_markov_ward` | 4, 0.0 | 1, 1.0 |

The case fails for today's default test by design; no default was changed.

## Not checked

The sibling test alone on a tree from independent noise. More than 5 seeds, other sizes, correlated noise,
real data. Whether the edge gate is still needed when the opt-in is on. E4 covers one Gaussian population
only.
