# Gaussian reference for the best flat cut: table, small n, search, selection of k

Date: 2026-11-11, Europe/Vienna. Scope: simulation checks of a proposed reference
for the sibling test, run on 10 and 11 November. No production statistical code
was changed and nothing was promoted.

## Finding

For one Gaussian population the statistic `rho = pi(1-pi) D^2` of the best flat
cut inside k directions does not depend on shift or invertible linear maps, so
its distribution depends only on (n, k). A table of its 99% point per (n, k) is
therefore a candidate reference, and because the best flat cut is the best of
all splits, the same reference covers every other split in those k directions.

Four limits of that table were checked:

1. n below 50: the table is finite and ordered from n = 15; n = 10 is unusable
   for k of 3 and above; at n = k + 1 `rho` is infinite in some datasets.
2. Search: the search that built the first table (fixed grid of directions,
   re-sorts from the single best one) misses the best flat cut for k of 4 and
   above, and its found `rho` is then not invariant. Re-sorts from the top 50
   directions plus single-point moves repair this up to k = 6.
3. More than 30 features, k chosen from the same data: no excess at n = 400; a
   small excess at n = 100. With n - 1 features or more, `invariant_directions`
   returns no directions and the test is skipped.
4. k above 6: the ordering continues. The 99% points are settled for k = 8 and
   10 with a stronger search; for k = 15 and 20 they are lower bounds only.

## Reproducibility and scope

Seed 20261110, one stream id per job kind; all scripts are deterministic and are
archived with the protocol in
`gaussian-best-flat-cut-reference-scripts-20261111.py.txt`. Numbers are in
`gaussian-best-flat-cut-reference-evidence-20261111.json`.

The first run of the small-n table and of the search comparison at n = 50, 100,
200 was stopped before it finished, and its per-dataset files were not retained.
Its printed summaries are kept verbatim in the evidence file. A same-seed re-run
of 10 small-n cells and 3 search cells reproduced the printed 99% points and
counts exactly; the other 20 small-n cells and 9 search cells rest on the printed
summary alone. The table for n = 50 to 800, the capped run and the enumerator
check were executed completely and their output is in the evidence file.

Baseline: wiki-lint passed on 272 pages. Environment: Python 3.12.11, uv 0.9.5;
NumPy and SciPy versions are recorded in the evidence file. No dependency
installation, commit, push or persistent service.

Everything assumes one Gaussian population. Nothing here tests non-Gaussian
populations, a complete tree fit, selected descendants or the final partition.

## Searches

- S0: all thresholds on a fixed grid of directions (k = 1 exact; k = 2, 720
  angles; k = 3, 2,000 points on the sphere; k of 4 and above, 3,000 random
  directions), then re-sorts from the single best direction.
- S1: the same grid, re-sorts from the top 50 directions, then single-point moves.
- S2 and S3: 30,000 and 300,000 directions, top 50, moves.
- S4: 30,000 directions, top 500, moves.
- S5: 200 random labelings, moves. It does not use the direction grid.
- EX: exact enumeration of every flat cut (hyperplanes through k of the n points
  with all assignments of the points on the plane).

Checks of EX: the number of distinct cuts equals the Cover count in 8 of 8 cases
(n = 8 to 16, k = 2 to 6); every enumerated cut is separable (129 of 129 at
n = 10, k = 3); on 100 datasets with n = 10 to 14 and k = 2 to 6 the maximum over
all 2^(n-1) - 1 splits equals the maximum over flat cuts to 1.1e-15.

## Table for n = 50 to 800 (S0, 4,000 datasets per cell)

99% point of the best flat cut's `rho`:

| k \ n | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|
| 1 | 2.987 | 2.515 | 2.247 | 2.086 | 1.989 |
| 2 | 3.623 | 2.856 | 2.493 | 2.219 | 2.066 |
| 3 | 4.298 | 3.149 | 2.640 | 2.332 | 2.136 |
| 4 | 4.769 | 3.469 | 2.788 | 2.416 | 2.200 |
| 5 | 5.452 | 3.787 | 2.932 | 2.503 | 2.250 |
| 6 | 6.173 | 4.028 | 3.102 | 2.586 | 2.302 |

The point rises with k and falls with n in every cell. A fresh S1 table at
n = 100 and n = 400 (3,000 datasets) agrees within Monte Carlo error: 2.491,
2.863, 3.171, 3.423, 3.740, 4.013, 4.280, 4.610 for k = 1 to 8 at n = 100 and
2.093, 2.222, 2.316, 2.425, 2.521, 2.600, 2.674, 2.752 at n = 400.

First selection test (k chosen by `invariant_directions` from the same data, 10
or 30 features, 2,000 datasets per cell, rejections at the table's 99% point for
n = 100 and n = 400): isotropic with 10 features 5 and 5; isotropic with 30
features 7 and 5; one weak direction (2.5) 28 and 24; two strong directions
(16, 4) 16 and 24; geometric 100 to 1 over 10 features 23 and 22. Pooled
159/20,000 (95% interval 0.0068 to 0.0093).

## n below 50 (first run, EX, 4,000 datasets per cell)

| k \ n | 10 | 15 | 20 | 30 | 40 |
|---|---|---|---|---|---|
| 1 | 7.744 | 5.359 | 4.355 | 3.521 | 3.207 |
| 2 | 16.914 | 8.039 | 6.319 | 4.773 | 3.907 |
| 3 | 36.484 | 13.018 | 8.311 | 5.882 | 4.748 |
| 4 | 106.848 | 19.592 | 11.762 | 7.211 | 5.692 |
| 5 | 376.018 | 32.662 | 15.984 | 8.691 | 6.522 |
| 6 | 3481.894 | 62.040 | 21.628 | 10.471 | 7.685 |

- Bootstrap standard errors: at n = 10 from 0.48 (k = 1) to 437 (k = 6); at
  n = 15 at most 2.8; at n = 20 at most 0.50; at n = 40 at most 0.10.
- S0 equals EX in every fully enumerated cell (0 of 4,000 datasets below, apart
  from one difference of 3.5e-9 at n = 10, k = 6).
- Three cells were too costly to enumerate fully and use S0: n = 30, k = 6;
  n = 40, k = 5; n = 40, k = 6. EX on their first 100, 100 and 20 datasets found
  no miss. At n = 40, k = 6 the stronger searches beat S0 in 6 of 1,000 datasets,
  by at most 0.234.
- Degenerate sizes, 200 datasets each: at n = k + 1 `rho` is infinite in 27, 12,
  6, 3 and 1 datasets for k = 2 to 6 and undefined for k = 1; at n = k + 2 it is
  finite with a 99% point between 1.7e3 and 1.0e9.

## Search for k of 4 and above (first run, 2,000 datasets per cell)

Datasets where S1 found a larger `rho` than S0, with the largest difference:

| k | n = 50 | n = 100 | n = 200 |
|---|---|---|---|
| 3 | 0 | 0 | 1 (0.000) |
| 4 | 0 | 24 (0.033) | 289 (0.041) |
| 5 | 1 (0.045) | 311 (0.117) | 727 (0.111) |
| 6 | 67 (0.394) | 634 (0.269) | 1,075 (0.201) |

- 99% point under S0 and S1: equal to three decimals except n = 100, k = 5
  (3.766 and 3.772), n = 200, k = 5 (2.940 and 2.944) and n = 200, k = 6 (3.067
  and 3.069). All differences are inside the bootstrap standard error.
- S1 against S2 on 200 random datasets per cell: S1 lower in at most 1, by at
  most 0.003. S2 against S3 on 60 random datasets per cell: 0 lower in all cells.
- Invariance, 200 datasets per cell, same data after a random invertible linear
  map and shift: the S0 value changed by more than 1e-6 in 0 datasets for k = 3,
  in up to 42 for k = 4, 73 for k = 5 and 138 for k = 6 (largest change 0.294);
  the S1 value changed in at most 2 (largest change 0.014).
- n = 400, k = 4, 5, 6 (capped run, 60 datasets per cell): S1 below the best of
  S1, S2, S4, S5 in 0, 1 and 5 datasets, by at most 0.002.

Consequence: with S0 the table's 99% point is nearly unchanged, but the found
value is not the best flat cut, so "any other split's `rho` is at most the best
flat cut's `rho`" does not hold for the found value. With S1 it holds up to
k = 6 in all but a few datasets per hundred.

## More than 30 features (capped run, S1, 2,000 datasets per cell)

`invariant_directions` returns no directions when n <= p + 1. Probe under
isotropic noise: 0 directions at n = 100 with 98, 99, 100, 200 and 500 features
and at n = 400 with 398, 399 and 500 features.

Rejections at the fresh S1 table's 99% point (20 would be 0.01):

| Covariance | n=100, p=50 | n=100, p=80 | n=400, p=50 | n=400, p=100 | n=400, p=200 |
|---|---|---|---|---|---|
| isotropic | 8 | 4 | 4 | 4 | 4 |
| one weak direction (2.5) | 31 | 16 | 15 | 18 | 14 |
| two strong directions (16, 4) | 21 | 23 | 20 | 19 | 21 |
| geometric 100 to 1 over 10, rest 1 | 28 | 33 | 17 | 20 | 19 |

- Isotropic: no direction is chosen in 1,483 to 1,675 of 2,000 datasets; pooled
  24/10,000; among datasets with at least one direction 24/2,066 (0.0075 to 0.0172).
- Three structured covariances at n = 400: 163/18,000 (0.0077 to 0.0105).
- Three structured covariances at n = 100: 152/12,000 (0.0107 to 0.0148). Two
  single cells exclude 0.01: 31/2,000 (0.0106 to 0.0219) and 33/2,000 (0.0114 to
  0.0231). The grouping by n was chosen after the cells were seen.
- The geometric covariance chose k = 9 in 245, 51 and 10 datasets at n = 400;
  the table stops at k = 8, so these 306 datasets were not tested and count as
  not rejected.

## k above 6 (capped run, S1, 1,000 datasets per cell)

99% point under S1:

| k \ n | 100 | 200 | 400 |
|---|---|---|---|
| 8 | 4.566 | 3.415 | 2.733 |
| 10 | 5.276 | 3.692 | 2.882 |
| 15 | 6.879 | 4.361 | 3.232 |
| 20 | 8.511 | 4.919 | 3.492 |

Search ladder on the first 60 datasets of each cell; "below" counts datasets
where a search is below the best of S1, S2, S4 and S5:

| k | n | S1 below (mean shortfall) | S4 below (largest) | S4 equals S5 |
|---|---|---|---|---|
| 8 | 100 | 2 (0.003) | 0 | 57 |
| 8 | 200 | 5 (0.002) | 0 | 53 |
| 8 | 400 | 16 (0.002) | 0 | 54 |
| 10 | 100 | 8 (0.004) | 0 | 54 |
| 10 | 200 | 17 (0.014) | 0 | 49 |
| 10 | 400 | 26 (0.010) | 2 (0.006) | 39 |
| 15 | 100 | 42 (0.286) | 1 (0.016) | 35 |
| 15 | 200 | 37 (0.077) | 4 (0.032) | 16 |
| 15 | 400 | 50 (0.061) | 8 (0.015) | 11 |
| 20 | 100 | 56 (0.759) | 8 (1.198) | 12 |
| 20 | 200 | 50 (0.275) | 15 (0.240) | 4 |
| 20 | 400 | 54 (0.146) | 16 (0.155) | 0 |

For k = 8 and 10, S4 is not beaten by more than 0.006 and the S1 table is short
by at most 0.014 on average. For k = 15 the S1 table is a lower bound, short by
0.06 to 0.29 on average. For k = 20 no search of the ladder is adequate: S4 is
beaten in 8 to 16 of 60 datasets and two independent search families agree in at
most 12 of 60.

## Limits

- One Gaussian population only. Earlier checks in
  `three-tests-audit-20261109.md` show that the variance of `rho` and the largest
  population `rho` depend on the population.
- The invariance argument covers fixed k. For k chosen from the same data there
  are only the measurements above, on five covariance shapes.
- Not run: n = 50 and 800 for k above 6; k = 9 and k = 11 to 14; the cause of the
  excess at n = 100; any change of production code.
