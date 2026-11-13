# The draft definition of k on the benchmark cases and the CHOIR cells

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/k_split_direction_runs_20261112` (git-ignored);
the protocol, the run-time plug-in, the scripts and their outputs are in
`k-split-direction-runs-scripts-20261112.py.txt`. Follow-up of `k-split-direction-draft-20261112.md`. Working tree
on top of `bc044e07` with the uncommitted opt-in reference. The draft exists in the analysis folder only; no
production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Setup

Draft: k = 1; the node's cells are divided into a basis half and a test half; one direction from the basis half
(variant A from the tree's split, variant B from a 2-means clustering of the basis half); `rho` of the tree's
split of the test half along it against the reference at (test cells, 1); p = 1 with fewer than 15 test cells.

At run time `best_flat_cut_reference` inside the opt-in's module is replaced by the draft, and each side's cells
are passed in the order of their rows in the data, so the halves alternate in data order and not in the tree's
order. The protocol with five predictions was written before the runs.

## Benchmark: 10 continuous and 4 cancer_pc cases, both registered methods, opt-in on

Found clusters, with the true groups in brackets; cases not listed are unchanged under all four rules:

| Case | Current rule | Prototype rule | Draft A | Draft B |
| --- | --- | --- | --- | --- |
| gauss_clear_medium (4) | 3 | 4 | 4 | 4 |
| gauss_moderate_3c (3) | 2 | 3 | 3 | 3 |
| gauss_dense_signal_highd (4), default tree / Markov/Ward | 1 / 1 | 4 / 3 | 2 / 2 | 1 / 1 |

| | Current rule | Prototype rule | Draft A | Draft B |
| --- | ---: | ---: | ---: | ---: |
| Continuous, 9 cases of 11 November, mean ARI, default tree | 0.6950 | 0.8880 | 0.8302 | 0.7768 |
| Continuous, 9 cases, mean ARI, Markov/Ward tree | 0.6950 | 0.8544 | 0.8302 | 0.7768 |
| cancer_pc, 4 cases, mean ARI | 0.7221 | 0.7221 | 0.7221 | 0.7221 |
| Cases with the exact number of clusters, of 14, default tree | 7 | 10 | 9 | 9 |
| Cases with more clusters than true groups | 0 | 0 | 0 | 0 |
| One-population grid, runs with 1 cluster | 80 of 80 | 80 of 80 | 80 of 80 | 80 of 80 |
| Smallest p-value among the tested calls of that grid | 1.0 | 0.0185 | 0.062 | 0.0735 |

Both null cases give 1 cluster under both variants. Of the 165 calls on the 14 cases, 155 are tested and 10 fail
closed with fewer than 15 test cells (current rule: 99 tested, 66 `no_invariant_directions`).

- gauss_dense_signal_highd stops at 2 clusters with variant A: the root (40 cells, 20 test cells) passes with
  `rho` 22.8; the two splits of 20 cells have 10 test cells and fail closed.
- Variant B refuses splits of nodes with more than two groups: in gauss_clear_medium the splits of 60 and 45
  cells have `rho` 0.38 and 0.31 along the 2-means direction (282 and 1,213 along variant A's direction);
  pass-through still reaches the 4 clusters. In gauss_dense_signal_highd the root has `rho` 0.39 and the case
  stays at 1 cluster.
- Remaining losses are refusals by the reference in both variants: dim_diffuse_6c_136f (root `rho` 1.61 with
  variant A), the connected cancer_pc cases, the single cell of gauss_single_outlier_4c.

## CHOIR cells: Markov/Ward tree, 48,879 cells, the same 1,838 calls

| | Current rule | Draft A | Draft B |
| --- | ---: | ---: | ---: |
| Found clusters | 230 | 218 | 153 |
| ARI against the cell lines | 0.8219 | 0.8311 | 0.5830 |
| NMI | 0.9358 | 0.9390 | 0.9051 |
| Lines recovered | 116 | 120 | 98 |
| Lines split only / merged only / both | 21 / 38 / 15 | 19 / 35 / 16 | 9 / 67 / 16 |
| Largest cluster | 1,990 | 1,990 | 3,124 |
| Clusters below 20 cells | 25 | 17 | 8 |
| Calls tested | 1,475 | 1,486 | 1,486 |
| Calls with p <= 0.01 | 149 | 170 | 81 |
| Accepted by BH | 134 | 151 | 55 |
| Accepted by BH at nodes with at least 1,000 cells (103 calls) | 50 | 68 | 8 |
| Run time | 126 s | 103 s | 105 s |

ARI between the assignments of the current rule and draft A: 0.9817. Node by node against the current rule:
draft A accepts 124 of the same splits, 27 further ones and not 10 of the current rule's; draft B accepts 53 of
the same, 2 further ones and not 81.

## Predictions

| Prediction | Result |
| --- | --- |
| R1 null cases and the one-population grid at 1 cluster, both variants | held (80 of 80 for each) |
| R2 gauss_clear_medium 4, gauss_moderate_3c 3, gauss_dense_signal_highd 2 clusters | held for variant A; variant B gives 4, 3 and 1 |
| R3 no benchmark case with more clusters than true groups | held |
| R4 CHOIR, variant A: ARI at least 0.75 | held (0.8311) |
| R5 CHOIR, variant A: more than 230 clusters | not held (218) |

## Reading

- Variant A, with k = 1 and no count of variance directions, gives on the CHOIR cells the result of the current
  rule or slightly better (ARI 0.8311 against 0.8219, 120 against 116 lines recovered), and on the continuous
  benchmark suite a mean ARI of 0.8302 against 0.6950, with the null cases and the one-population grid at 1
  cluster.
- Against the prototype rule it loses gauss_dense_signal_highd in part (2 of 4 groups), because a node needs 15
  test cells, that is about 30 cells.
- Variant B, whose level is exact, loses power wherever a node holds more than two groups, because its direction
  does not follow the tree's split: ARI 0.5830 on the CHOIR cells.
- So the working variant is the one whose level is measured and not derived. Measured so far: at most 13 of
  1,000 simulated datasets of one Gaussian population, 80 of 80 grid runs at 1 cluster, smallest p-value 0.062.

## Not checked

The level of variant A beyond these simulations; the default tree on the CHOIR cells; other real datasets; more
than one direction; halves formed differently; whether the 27 further splits accepted by variant A on the CHOIR
cells are between or inside cell lines. One seed per benchmark case.
