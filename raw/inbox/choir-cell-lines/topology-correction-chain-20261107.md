# Tree topology, the estimated correction and the chain test on the same no-group data

Analysis on 2026-11-07, code at `cbc317b2` plus the staged working tree. Local
folder: `results/analyses/topology_correction_20261107` (git-ignored); scripts
captured in `topology-correction-chain-scripts-20261107.py.txt`. Protocol frozen
before results; deviations listed at the end.

## Data and methods

Count nulls with one true group, 400 cells, 3,000 genes, NB dispersion 8, 300
replicates each: `low` (log-depth sd 0.25), `high` (sd 1.0), `matched_median`
(sd 0.42, the median within-line spread of the real cell lines). Frozen
preparation (log-normalised counts, 2,000 HVG, 30 PCs, kernel, symmetric and
Markov diffusion coordinates). Trees: `symmetric_average` (default method) and
`markov_ward`, production alphas (edge 1e-3, sibling 1e-2). Every split in these
data is a no-difference split.

## Q1 Topology on no-group data

| Tree | Depth sd | Root smaller side (median) | Root peels <= 5 cells | Peel share, nodes >= 40 | False partitions /300 | Median clusters |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Default (symmetric/average) | 0.25 | 2 | 88% | 0.32 | 2 | 1 |
| Default | 0.42 | 5.5 | 50% | 0.14 | 50 | 1 |
| Default | 1.0 | 125.5 | 0.7% | 0.10 | 259 | 5 |
| Markov/Ward | 0.25 | 138.5 | 0% | 0 | 298 | 8 |
| Markov/Ward | 0.42 | 140 | 0% | 0 | 296 | 8 |
| Markov/Ward | 1.0 | 161 | 0% | 0 | 300 | 8 |

The default tree peels when the table has no structure and divides when depth
puts structure into the table. Ward divides pure noise into balanced halves.

Root gates (share of datasets): the default tree's child-parent (edge) gate
opens at the root in 1%, 18% and 87% (sd 0.25, 0.42, 1.0) and the root splits
in 0%, 6.7% and 84.7%. Markov/Ward opens it in 100% and splits in 95-99%.

## Q2 Applied versus needed correction

Needed factor at 0.01 = 99th percentile over tested no-difference splits of
`T / (s * chi2.ppf(0.99, df))`; applied = today's empirical-null inflation
factor. Score = `T / (s * df)`, about 1 for one population with no selection.

| Tree | Depth sd | Root score (median) | Needed factor | Applied factor (median) | Adjusted p <= 0.01 |
| --- | --- | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 9.9 | 7.8 | 2.6 | 96.7% |
| Markov/Ward | 0.42 | 8.7 | 7.7 | 2.6 | 95.0% |
| Markov/Ward | 1.0 | 17.6 | 22.1 | 2.2 | 99.0% |
| Default | 0.42 (54 roots tested) | 3.9 | 4.7 | 2.6 | 37.0% |
| Default | 1.0 (261 tested) | 10.1 | 12.2 | 1.9 | 97.3% |

By parent size (Markov/Ward, sd 0.25): 10-39 cells applied 9.2 vs needed 9.4,
1.0% below 0.01; 40-159 applied 2.5 vs needed 7.0, 61%; 160-400 applied 2.6 vs
needed 7.6, 94%. The estimate is close for small splits and three times too
small for large ones. The sd 1.0 rows mix selection with real depth structure
in the table.

Code reading (not measured): splits enter the reference set with weight equal
to their child-parent edge p-value, so splits whose sides clearly differ from
the whole carry almost no weight
(`child_parent_edge_metadata.estimate_sibling_null_weight_from_child_parent_edges`).
The stored annotations do not record the reference splits' scores, so this
cannot be checked from the saved output.

## Q3 Chain test at the default tree's root

`selected_diffusion_root_chain_test`, 30 PCs, 199 steps, sigma^2 = mean PC
variance (fixed plug-in rule; the test's contract of known spherical Gaussian
noise does not hold for count data).

| Depth sd | Rejected at 0.01 | p <= 0.05 | Stalled | Median accepted moves /199 | Median sec | Production root splits |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.25 | 0/300 | 1.0% | 6.0% | 7 | 9.7 | 0 |
| 0.42 | 0/300 | 0.3% | 4.3% | 6 | 17.5 | 20 |
| 1.0 | 0/300 | 0.7% | 19.7% | 2 | 9.3 | 254 |

The replayed root equals the production root in all 900 datasets.

Added power check (deviation): 400-cell random subsets of the signal datasets
(4 true groups), 10 replicates each. Production splits every root and recovers
the groups (mean ARI 1.000, 0.993, 0.993). The chain test rejects 2/10
(balanced), 1/10 (imbalanced), 1/10 (singleton); median p 0.02; median accepted
moves 72, 24.5, 33.5. At 1,000 cells a 9-step chain did not finish within
10 minutes.

## Deviations

- Q3 restarted after 197 tests with 8 instead of 4 workers; seeded by replicate,
  so results are unchanged.
- The signal power check was added after the Q3 results.
