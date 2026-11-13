# A fit that splits each set of cells at its best flat cut (CHOIR cell lines)

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/best_flat_cut_divisive_20261111`
(git-ignored); protocol (written before the run), scripts and summary output are in
`best-flat-cut-divisive-scripts-20261111.py.txt`. Trial, descriptive; no production change.

## Procedure

Not Tree-Break Selection: no diffusion distances, no linkage tree, no edge gate, no pass-through, no BH step.
Start with all 48,879 cells in 30 PCs. For a set of cells: k = number of `invariant_directions` of these
cells; stop if k = 0 or fewer than 15 cells; find the best flat cut in those directions with the table's
search; p-value under one Gaussian population at (n, k) with n rounded down to the table grid; split and
repeat on both sides if p <= 0.01. For k <= 10 the p-value comes from the stored table; for k > 10 from 399
null datasets drawn at run time and searched with the same search. The searched value has the same null
distribution for every covariance (the whitened Gaussian sample does), so the run-time null is a valid
reference for the searched value at fixed (n, k) also where the search gives only a lower bound.

Self-check on simulated data (`selfcheck.py`): three separated Gaussian groups of 300, 150 and 150 cells in 30
features give 3 clusters with ARI 1.0; one Gaussian population of 600 cells gives 1 cluster.

Run-time nulls, 99% point: n = 800: 2.533, 2.562, 2.600, 2.659, 2.682 for k = 11 to 15 (the table has 2.48
at k = 10); n = 400: 2.954, 3.059, 3.081 for k = 11 to 13; n = 300: 3.246, 3.292 for k = 11, 12; n = 200:
3.777 for k = 11. The whole fit took 282 s with 10 workers; the guard of 6,000 tested sets was not reached.

## Result (190 cell lines)

| Fit | Clusters | ARI | NMI | Lines 90/90 | ARI vs CHOIR_1 | Clusters >= 90% pure | Clusters < 20 cells | Largest cluster | Lines split only / merged only / both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Divisive best flat cut | 323 | 0.874 | 0.941 | 93 | 0.617 | 264 | 44 | 1,990 | 81 / 8 / 8 |
| Markov/Ward, ceiling 3 (9 Nov) | 217 | 0.833 | 0.940 | 119 | 0.538 | 164 | 25 | 1,990 | 16 / 33 / 22 |
| Markov/Ward, rule-out | 664 | 0.765 | 0.901 | 80 | 0.514 | 494 | 253 | 1,990 | 72 / 11 / 27 |
| Markov/Ward, reference on the tree | 138 | 0.297 | 0.851 | 76 | 0.163 | 98 | 10 | 7,678 | 9 / 91 / 14 |
| Markov/Ward, today (8 Nov) | 954 | 0.549 | 0.863 | 41 | 0.445 | 761 | 374 | 1,451 | not recorded |

## Tested sets

| Cells | Sets | No invariant directions | Tested | Split | Median k (largest) | Median `rho` | Best flat cut with fewer than 5 cells on one side |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| at most 14 | 33 | 33 | 0 | 0 | | | |
| 15 to 31 | 37 | 37 | 0 | 0 | | | |
| 32 to 159 | 176 | 0 | 176 | 43 | 7 (10) | 4.13 | 21 |
| 160 to 999 | 320 | 0 | 320 | 201 | 10 (13) | 4.25 | 4 |
| 1,000 and more | 79 | 0 | 79 | 78 | 13 (15) | 7.00 | 0 |

322 cuts were accepted: 169 with k above 10; 290 at the lowest possible p-value (0.0005 or 0.0025); 9 with
fewer than 5 cells on the smaller side and 44 with fewer than 20. 253 tested sets were not split, 36 of them
with p between 0.01 and 0.05. The first 13 cuts (sets of 10,000 cells or more) have `rho` of 3.88 to 8.65 at
k of 13 to 15 and a smaller side of 250 to 12,291 cells.

## Where the errors are (added after the endpoints)

- Lines by the number of clusters that hold at least 5% of the line: 102 lines in one cluster, 65 in two, 23
  in three or more. Line size does not separate them (median 230, 221 and 213 cells).
- Every one of the 190 lines is the majority of at least one cluster; under ceiling 3 on Markov/Ward, 22 lines
  are the majority of no cluster.
- Joining clusters with the same majority line (uses the truth, so it is a description, not a method) gives
  ARI 0.957 for the divisive fit and 0.841 for ceiling 3.
- ARI between the divisive fit and ceiling 3 on Markov/Ward: 0.755.

## Reading

- The divisive fit separates the cell lines from each other almost completely: 16 lines are involved in a
  merge (8 merged only, 8 both) against 55 under ceiling 3 and 105 under the reference on the tree. Its
  remaining error is splitting inside lines: 88 lines are spread over two or more clusters.
- It has the highest ARI against the cell lines (0.874) and against CHOIR_1 (0.617) of all fits recorded on
  these data, with fewer lines recovered at 90/90 than ceiling 3 (93 against 119) because of the splits
  inside lines.
- This agrees with `split-vs-best-flat-cut-20261111.md`: the rebuilt geometry holds the line boundaries the
  tree hides, and the Gaussian reference rejects inside single lines more often than one Gaussian population
  allows.
- Expectation check: fewer merged lines, more than 190 clusters, lines split and clusters under 20 cells all
  held. ARI above ceiling 3 was not predicted.

## Limits

- Per-test level 0.01 with no family-wise or FDR step; the number of false splits over the fit is not
  controlled, and no complete-fit null calibration was run beyond the one simulated population of the
  self-check.
- k is chosen from the same cells, and the sides of an accepted cut are not Gaussian samples; the reference
  covers neither.
- Whether the splits inside lines are real structure or non-Gaussian shape cannot be told from these data.
- One dataset; cell lines are the only truth. Not run: the benchmark cases, other values of the level, a
  comparison of distance fit with the tree, any integration with the tree or with pass-through.
