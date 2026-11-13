# The level of variant A of the draft definition of k on one-population data with the method's own trees

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/k_split_direction_level_20261112` (git-ignored);
the protocol, the deviations, the scripts, the summary and an independent spot check are in
`k-split-direction-level-scripts-20261112.py.txt`. Follow-up of `k-split-direction-runs-20261112.md`. The draft
exists in the analysis folder only; no production code was changed. Environment: Python 3.12.11, NumPy 2.3.4,
SciPy 1.16.3, uv 0.9.5.

## Question

Variant A of the draft takes its one direction from the tree's split on the basis half, so its level is not
derived. How often does it reject on data of one Gaussian population when the tree is the method's own?

## Setup

One Gaussian population through the registered generator (`blobs_continuous`, 1 cluster), method
`tbs_family_diffusion` with operator symmetric or markov and linkage average or ward, opt-in on, variant A
plugged in at run time. Identity covariance at (cells, features) = (60, 30), (100, 30), (200, 30), (256, 30),
(496, 30), (1000, 30), (200, 80), (60, 200); unequal variances with correlated features (the 30 eigenvalues of
the CHOIR root, 53.09 to 5.12, and a fixed rotation) at (256, 30), (496, 30), (1000, 30). 50 seeds per setting,
20 for 1,000 cells; 1,960 fits. The protocol with four predictions was written before any run; the measurement
was built separately from it and checked before use.

Measured: the found clusters of every fit; the calls made during the fits; and variant A and the current rule
applied directly to every internal node with at least 30 cells of every fitted tree.

## Results

Fits with more than 1 cluster:

| | symmetric/average | symmetric/ward | markov/average | markov/ward | All |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fits with more than 1 cluster | 2 of 490 | 6 of 490 | 2 of 490 | 2 of 490 | 12 of 1,960 |

By covariance: identity 5 of 1,480 fits, unequal variances 7 of 480. The 12 fits find 3, 3, 6, 6, 6, 7, 7, 8,
8, 8, 9 and 9 clusters. In 11 of them exactly one call has p <= 0.01 and in one fit three calls; these nodes
have 29 to 65 cells, `rho` 3.64 to 7.87 and p 0.001 to 0.008. The other clusters come from the nodes above the
accepted node, which the traversal passes through.

Called nodes (12,302 tested calls):

| | Tested calls | With p <= 0.01 | Smallest p |
| --- | ---: | ---: | ---: |
| symmetric/average | 1,627 | 3 | 0.0025 |
| symmetric/ward | 3,983 | 12 | 0.002 |
| markov/average | 2,577 | 4 | 0.002 |
| markov/ward | 4,115 | 15 | 0.001 |
| All | 12,302 | 34 (0.28%) | 0.001 |

By node size: 28 of 2,405 calls below 60 cells, 6 of 7,094 with 60 to 199 cells, 0 of 2,803 with 200 and more.

Every node with at least 30 cells (32,038 nodes):

| | Variant A: tested, with p <= 0.01 | Current rule: tested, with p <= 0.01 |
| --- | --- | --- |
| symmetric/average | 8,490, 22 | 6,208, 4 |
| symmetric/ward | 7,236, 25 | 5,793, 6 |
| markov/average | 8,434, 23 | 6,210, 1 |
| markov/ward | 7,322, 39 | 5,924, 14 |
| All | 31,482, 109 (0.35%) | 24,135, 25 (0.10%) |

Variant A by the number of test cells: 15 to 19 test cells 60 of 6,983 (0.86%); 20 to 29: 40 of 8,125 (0.49%);
30 to 49: 8 of 6,406 (0.12%); 50 to 99: 1 of 4,856 (0.02%); 100 and more: 0 of 5,112. By covariance: identity
65 of 19,093 tested nodes (0.34%), unequal variances 44 of 12,389 (0.36%).

## Predictions

| Prediction | Result |
| --- | --- |
| V1 more than 1 cluster in at most 19 of the 1,960 fits | held (12) |
| V2 called nodes: p <= 0.01 in at most 2% of the tested calls | held (34 of 12,302, 0.28%) |
| V3 every node: variant A p <= 0.01 in at most 2% of the tested nodes | held (109 of 31,482, 0.35%) |
| V4 markov/ward has the largest share under variant A | held (0.53% against 0.26%, 0.35%, 0.27%) |

## Independent spot check (`verify.py`)

Inside the recorded files a called node and the direct call on the same node have the same p-value in 12,043 of
12,043 cases. An own run of (256, 30) identity, seeds 1 to 10, reproduces the found clusters and the p-values of
the tested calls in 40 of 40 fits.

## Reading

- With one Gaussian population variant A stays at or below the level 0.01 in every class: 0.35% of all tested
  nodes, 0.86% at the smallest test halves (15 to 19 test cells), none at 100 test cells and more.
- The level does not depend on the covariance of the population (0.34% and 0.36%).
- Variant A is tested at more nodes than the current rule (31,482 against 24,135) and rejects more often (0.35%
  against 0.10%); both are below 0.01.
- The complete method finds more than 1 cluster in 12 of 1,960 fits (0.61%). When it does, it finds 3 to 9
  clusters: one accepted node of 29 to 65 cells is enough, because the nodes above it are passed through.
- All false acceptances are at small nodes. A larger smallest number of test cells would remove most of them;
  this was not tried.

## Deviations

- The protocol's self-check "largest eigenvalue of the data's covariance below 5 for identity settings" fails
  for (60, 200) in all 50 seeds (7.28 to 8.48). The check was written for more cells than features and does
  not apply to 60 cells in 200 features; the setting's results are kept and counted.
- One fit, (1000, 30) with unequal variances, seed 1, symmetric/ward, stopped with a branch-length solver that
  did not converge. It is counted among the 1,960 fits as not more than 1 cluster, with no calls and no nodes.
  A catch for this error was added to the script after the first run of that setting stopped.
- The left side of a node in the direct calls is the first child in the tree's order; the draft is symmetric in
  the two sides, and the spot check confirms equal p-values.

## Not checked

Why the one fit did not converge. Populations that are not Gaussian. A larger smallest number of test cells.
The complete method with the current rule on the same 1,960 datasets (only its direct calls were recorded).
Real data beyond the CHOIR run of the same day.
