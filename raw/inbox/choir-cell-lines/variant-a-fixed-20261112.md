# Variant A of the draft definition of k with the findings of the check fixed

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/variant_a_fixed_20261112` (git-ignored); the
protocol with its amendment, the two versions of the rule, the builder of the exact null, the run scripts and the
comparison outputs are in `variant-a-fixed-scripts-20261112.py.txt`. Follow-up of `variant-a-check-20261112.md`.
The rule exists in the analysis folder only; no production code was changed. Environment: Python 3.12.11, NumPy
2.3.4, SciPy 1.16.3, uv 0.9.5.

## Looked-up sources (retrieved 2026-11-12)

1. N. Meinshausen, L. Meier, P. Buehlmann, "p-Values for High-Dimensional Regression", Journal of the American
   Statistical Association 104(488), 1671-1681, 2009. https://doi.org/10.1198/jasa.2009.tm08647 (record); text
   read at https://arxiv.org/abs/0811.2177, page 5. Confirmed: one random split gives a "p-value lottery"; over
   B splits Q(gamma) = min{1, gamma-quantile of the p-values divided by gamma} is a p-value for any fixed gamma
   (asymptotically correct in the paper's setting); for gamma = 0.5 it is twice the median of the p-values.
2. V. Vovk, R. Wang, "Combining p-values via averaging", https://arxiv.org/abs/1212.4966 (abstract). Confirmed:
   twice the arithmetic mean of p-values is valid without any assumption on their dependence.

## The two findings and the fixes

1. A feature constant in the basis half raised an uncaught `LinAlgError`. Fix: the direction uses the
   pseudo-inverse of the pooled within covariance; a division without a usable direction gives p = 1. On the two
   constructed inputs of the check there is no error (p = 1.0 both).
2. The halves alternated in the order of the rows. Fix: each node is divided at random 20 times (random number
   generator seeded by the cells of the node, the cells of the left side and the index of the division); the
   node's p-value is twice the median of the 20 p-values (source 1).

Second version: these two fixes, p-values from the production table (1,999 draws per cell). Third version
(amendment, written after the second version's CHOIR result and before any run of the third): the same, but the
p-value of a division comes from an exact null of the best cut of one direction, 99,999 Gaussian datasets per n
of the table's grid.

## Why the third version: the smallest possible p-value

Second version on the CHOIR cells, data order: 151 calls have the smallest possible p-value, 0.001 = 2 x 0.0005,
and Benjamini-Hochberg accepts 97 of them. Of the 54 splits variant A accepted and the second version does not,
46 have p = 0.001. The table cannot give a p-value below 1/2,000, and doubling it is too much for the BH step,
which works per depth level of the tree. With the exact null the smallest node p-value is 0.00002: 138 calls
have it and 138 are accepted.

Exact 99% points of the best cut of one direction against the production table at k = 1: n = 15: 5.343 (table
4.888); n = 100: 2.531 (2.439); n = 150: 2.356 (2.389); n = 400: 2.091 (2.073); n = 800: 1.983 (1.984).

## Results

Benchmark, 14 cases, both registered methods: found clusters and ARI equal variant A in 28 of 28 rows for the
second and the third version; both null cases 1 cluster; one-population grid 1 cluster in 80 of 80 runs
(smallest p-value 1.0 and 0.188; variant A 0.062).

Level, one Gaussian population, 800 fits ((256, 30) and (496, 30), identity and the unequal variances of the
CHOIR root, 50 seeds, 4 operator and linkage cells):

| | Fits with more than 1 cluster | Tested calls with p <= 0.01 | Smallest p |
| --- | ---: | ---: | ---: |
| Variant A (one division in row order) | 8 of 800 | 20 of 7,249 | 0.001 |
| Second version | 0 of 800 | 0 of 7,249 | 0.0215 |
| Third version | 0 of 800 | 0 of 7,249 | 0.0332 |

CHOIR cells, Markov/Ward tree, five orders of the rows (data order and seeds 1, 2, 3, 777):

| | Found clusters | ARI | Lines recovered | Lines merged only | Accepted by BH |
| --- | --- | --- | --- | --- | --- |
| Current rule (no dependence on the order) | 230 | 0.8219 | 116 | 38 | 134 |
| Variant A | 196 to 225 | 0.8043 to 0.8397 | 114 to 125 | 27 to 45 | 139 to 156 |
| Second version | 133 to 135 | 0.5371 to 0.5374 | 89 to 91 | 79 to 82 | 96 to 97 |
| Third version | 206 to 215 | 0.8349 to 0.8411 | 126 to 127 | 29 to 31 | 150 to 155 |

Third version: per order 206, 211, 215, 213, 210 clusters and ARI 0.8359, 0.8401, 0.8411, 0.8404, 0.8349; the
partitions of two orders have an ARI of 0.9859 to 0.9989; 149 nodes are accepted in every order and 9 in only
some; 3 or 4 accepted splits per order are at nodes with at least 90% of the cells from one line; the largest
cluster has 1,990 cells in every order; one run takes 143 s (current rule 126 s).

## Predictions

| Prediction | Result |
| --- | --- |
| G0 no error on the constructed inputs | held |
| G1 second version: benchmark as variant A, nulls at 1 cluster | held (28 of 28, 80 of 80) |
| G2 second version: at most 6 of 800 fits with more than 1 cluster | held (0) |
| G3 second version: CHOIR spread at most 10 clusters and 0.01 ARI | held (2 clusters, 0.0002) |
| G4 second version: CHOIR ARI at least 0.80 in every order | not held (0.5371 to 0.5374) |
| H1 third version: benchmark as variant A, nulls at 1 cluster | held (28 of 28, 80 of 80) |
| H2 third version: at most 6 of 800 fits with more than 1 cluster | held (0) |
| H3 third version: CHOIR spread at most 10 clusters and 0.01 ARI | held (9 clusters, 0.0062) |
| H4 third version: CHOIR ARI at least 0.80 in every order | held (0.8349 at the lowest) |

## Reading

- Both findings of the check are fixed in the third version: no error, and the spread over the order of the
  rows is 9 clusters and 0.0062 ARI instead of 29 clusters and 0.035.
- The loss of the second version was a loss of resolution of the p-value, not of evidence: with an exact null of
  99,999 draws the same rule recovers it.
- The third version is above the current rule on the CHOIR cells in all five orders (ARI at least 0.8349
  against 0.8219, 126 to 127 lines recovered against 116), equals variant A on the benchmark, and rejects in 0
  of 800 fits of one Gaussian population.
- It uses k = 1, no count of variance directions, and a null that is exact for one direction.

## Deviations

The level runs of the second version did not start at first (two arguments were passed as one) and were
restarted unchanged. The protocol states 6 fits with more than 1 cluster for variant A on the 800 fits; the
recount from the recorded files gives 8, so the bound of G2 and H2 was stricter than intended. The third version was not in the protocol; it was added by an amendment after the second
version's CHOIR result, with its cause measured and four predictions stated before its runs.

## Not checked

The third version's level beyond these 800 fits and the 80 grid runs (the 1,960 fits and the direct calls at
every node were measured for variant A only). Populations that are not Gaussian. The default tree on the CHOIR
cells; other real datasets. More than five orders. The exact null exists for k = 1 only; n above 800 uses the
n = 800 null. The direction still comes from the tree's split, so the level is measured, not derived; twice the
median is a p-value in the source's setting, where each split's p-value is valid.
