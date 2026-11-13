# The Gaussian best-flat-cut reference extended to 11 to 15 invariant directions

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/best_flat_cut_reference_k15_20261111`
(git-ignored); protocol (frozen before any run), the recorded deviations of part 2, the check scripts, their
output and the real-data table are in `best-flat-cut-reference-k15-scripts-20261111.py.txt`. Production change of
the opt-in `sibling_best_flat_cut_reference` only (uncommitted working-tree code on top of `bc044e07`); no
default changed. Part 1 and part 2 were written separately in parallel; both were checked afterwards with
separate comparisons before the real-data runs.

## Part 1: the change

- `scripts/build_best_flat_cut_reference.py`: `K_GRID` is 1 to 15; `--k-min` computes only the cells with k at or
  above it and copies the others from the existing table (same n grid, B and seed required). Search, seeds and
  B = 1,999 are unchanged. For k above 10 the cells with n below 50 are stored as NaN (at n = 15 the whitening
  fails for k = 15 and `rho` is 3,000 to 1e10 for k = 11 and 13).
- `best_flat_cut_reference.py`: `MAX_DIRECTIONS` = 15, `MIN_CELLS_ABOVE_10_DIRECTIONS` = 50 (k above 10 with
  fewer cells fails closed as `too_few_cells`; `reference_p_value` raises instead of reading a NaN cell); the
  docstring states what is measured for k above 10.
- Tests: `tests/statistics/82_test_best_flat_cut_reference.py` (table shape and NaN cells; 11 to 15 directions
  are tested; too few cells above 10 directions; more than 15 directions). The three reference test files give
  13 passed (11 before); `tests/statistics` 534 passed; ruff passes.
- Checks: the cells for k = 1 to 10 are bit-identical to the table before the change; (n = 15, k = 10) and
  (n = 50, k = 5) recomputed from their seeds equal the stored values; the build of the 40 new cells took
  1,216 s on 8 workers.

99% points of the new cells (k = 10 is the old table):

| n | k=10 | k=11 | k=12 | k=13 | k=14 | k=15 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 9.678 | 11.089 | 12.449 | 13.955 | 15.000 | 17.339 |
| 75 | 6.477 | 6.994 | 7.386 | 8.123 | 8.622 | 9.322 |
| 100 | 5.317 | 5.648 | 5.893 | 6.235 | 6.591 | 6.824 |
| 150 | 4.186 | 4.333 | 4.544 | 4.716 | 4.946 | 5.123 |
| 200 | 3.688 | 3.797 | 3.942 | 4.095 | 4.171 | 4.343 |
| 300 | 3.160 | 3.235 | 3.379 | 3.462 | 3.531 | 3.603 |
| 400 | 2.901 | 2.978 | 3.041 | 3.120 | 3.179 | 3.256 |
| 800 | 2.484 | 2.531 | 2.566 | 2.608 | 2.623 | 2.671 |

They rise with k and fall with n in every cell. At n = 800 they differ from the run-time values of the divisive
trial (399 datasets: 2.533, 2.562, 2.600, 2.659, 2.682) by at most 0.036 (k = 14).

## Part 2: validity under one Gaussian population (k = 11, 13, 15)

The stored null for k above 10 is the null of the searched value, because the search does not reliably find the
best flat cut there. Measured (exact binomial 95% intervals; rejection = `rho` above the 99% point):

| Cell (n, k) | Searched value, fresh data | Ward root split | `rho` of the Ward split / 99% point: median, maximum |
| --- | --- | --- | --- |
| 200, 11 | 10/1,000 | 0/1,000 | 0.299, 0.498 |
| 200, 13 | 6/1,000 | 0/1,000 | 0.283, 0.457 |
| 200, 15 | 6/1,000 | 0/1,000 | 0.271, 0.416 |
| 800, 11 | 10/1,000 | 0/1,000 | 0.262, 0.386 |
| 800, 13 | 9/1,000 | 0/1,000 | 0.240, 0.346 |
| 800, 15 | 5/1,000 | 0/1,000 | 0.227, 0.313 |

- With k chosen by `invariant_directions` in 30 features (variance 25 on k coordinates, n = 400 and 2,000,
  1,000 datasets per cell): the intended k in 6,000 of 6,000 datasets and 0 rejections of the Ward root split
  in all six cells (largest ratio to the 99% point 0.456). In all, 0 rejections in 12,000 null datasets.
- Stronger searches (30,000 directions; re-sorts and single-cell moves from the Ward split; the same from 20
  random half-splits): the best of them beats the table's search in 191, 300 and 350 of 600 datasets at n = 200
  and in 172, 187 and 208 of 300 at n = 800, by a median of 0.6 to 2.6% and at most 16.8% in `rho`. It exceeds
  the stored 99% point in 25 of 2,700 datasets (0.93%) against 20 of 2,700 for the searched value; the largest
  cell is (800, 15) with 6/300 = 2.0% [0.74%, 4.30%] against 3/300.
- A separate spot check with other seeds agrees: 0 of 300 at (200, 13), largest `rho` 1.771 against the 99%
  point 4.095; through the whole function at n = 400 with 12 strong directions, 0 of 200.

Part 2 deviations (recorded before the summary): moves started from the table's result add nothing, so three
other searches were run; the searched value through the chosen k was run on 200 datasets per cell; the same
quantity as the whole function was computed directly because the table was being rebuilt.

## Part 3: the real CHOIR cells (Markov/Ward tree, harness of spr-refit)

| | k <= 10, unchanged tree | k <= 15, unchanged tree | k <= 10, tree after sprs | k <= 15, tree after sprs |
| --- | ---: | ---: | ---: | ---: |
| Found clusters | 138 | 230 | 257 | 316 |
| ARI against the 190 cell lines | 0.2975 | 0.8219 | 0.3425 | 0.7376 |
| NMI | 0.8512 | 0.9358 | 0.8846 | 0.9316 |
| Lines recovered at 90/90 | 76 | 116 | 82 | 107 |
| Lines split only | 9 | 21 | 47 | 58 |
| Lines merged only | 91 | 38 | 51 | 19 |
| Lines both | 14 | 15 | 10 | 6 |
| Largest cluster | 7,678 | 1,990 | 6,346 | 1,451 |
| Clusters at least 90% pure | 98 | 174 | 202 | 251 |
| Clusters under 20 cells | 10 | 25 | 30 | 39 |
| Splits at p <= 0.01 | 57 | 149 | 184 | 275 |
| Splits accepted by BH | 31 | 134 | 104 | 210 |
| Lines that are the majority of no cluster | 75 | 23 | 41 | 11 |
| ARI after joining clusters with the same majority line | 0.3008 | 0.8372 | 0.4397 | 0.8914 |

- Unchanged tree: all 190 splits with 11 to 15 invariant directions are tested; 92 have p <= 0.01, all at the
  lowest possible 0.0005, and BH accepts 88 (50 of 103 with 1,000 leaves or more, 38 of 87 with 160 to 999).
  Among the 1,285 splits with 1 to 10 directions the same 57 have p <= 0.01 and BH now accepts 46 (31 before);
  7 splits at p 0.0005 are not accepted (17 before).
- For comparison on the same tree: ceiling 3 gave 217 clusters, ARI 0.8329, 119 lines recovered, largest
  cluster 1,990; the divisive trial 323 clusters, ARI 0.874.
- Tree after the sprs (built under the k <= 10 reference): fewer merged lines (19 against 38) and fewer lines
  without a cluster (11 against 23), more split lines (58 against 21) and a lower ARI (0.7376 against 0.8219).

## Reading

- The limit of 10 directions was the cause of the poor real-data result of the opt-in: with 11 to 15 directions
  tested, the ARI on the unchanged tree goes from 0.30 to 0.82, about the value of ceiling 3.
- Under one Gaussian population the tree's split stays far below the reference at k = 11 to 15, so the test is
  conservative in the way the pipeline uses it. The stored null is slightly below the null of the true best flat
  cut; the measured excess is within about one percentage point except in one cell.
- The sprs trade merged lines for split lines. With the extended reference they lower the ARI.

## Limits

One Gaussian population, three k and two n per part, Ward root split only (not deeper splits); the strongest
search is still a search, so the true best flat cut's null may be higher than measured; k is chosen from the
same data; one real dataset and one tree; the lowest possible p is 1/2,000; no benchmark panel and no
complete-fit null calibration were run; the default tree was not rerun.
