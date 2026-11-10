# Rule-out test on the real CHOIR cell-line data, end to end

Run on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/split_count_exclusion_real_20261108` (git-ignored); protocol, runner, summary and the metric
lines in `split-count-exclusion-real-scripts-20261108.py.txt`. Protocol frozen before the run; no deviations.

## Setup

Both registered methods as in the 2026-11-06 full-data runs (`choir_cell_lines_tbs_lbfgsb_20261106`): 48,879
cells, 190 lines, 30 PCs, saved diffusion coordinates, registry parameters, edge alpha 1e-3, sibling alpha 1e-2.
The sibling p-value of each tested split is swapped by the unchanged harness of
`split-count-exclusion-20261108.md`: today (production), count_all (union bound over 2^(n-1) - 1 splits) and
rule_out (count_all, and p = 1 unless the split is separable in its invariant directions with at most 5% of
cells equidistant). Leaf rows of tested splits are found on demand (scale only). Replay: today reproduces the
2026-11-06 results exactly (75 and 954 clusters, ARI 0.0275 and 0.549, 27 and 41 lines); cell counts and
p-values of tested splits match (90/90 default, 1,838/1,838 Markov/Ward). Runtime 109-163 s per fit.

## Results

| Method | Arm | Clusters | ARI (lines) | NMI | Lines 90/90 | ARI (published CHOIR_1) | Clusters >= 90% one line | Clusters < 20 cells | Lines with >= 90% in one cluster | Largest cluster |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Default | today | 75 | 0.028 | 0.548 | 27 | 0.014 | 48 | 20 | 177 | 28,611 |
| Default | count_all | 82 | 0.027 | 0.549 | 29 | 0.014 | 54 | 23 | 174 | 28,611 |
| Default | rule_out | 69 | 0.018 | 0.497 | 29 | 0.009 | 49 | 12 | 174 | 32,680 |
| Markov/Ward | today | 954 | 0.549 | 0.863 | 41 | 0.445 | 761 | 374 | 56 | 1,451 |
| Markov/Ward | count_all | 1,605 | 0.402 | 0.833 | 24 | 0.334 | 1,330 | 940 | 28 | 1,451 |
| Markov/Ward | rule_out | 664 | 0.765 | 0.901 | 80 | 0.514 | 494 | 253 | 91 | 1,990 |

Lines 90/90: a cluster holds >= 90% of the line and the line is >= 90% of that cluster.

Tested splits under rule_out (kept as real / accepted after the traversal):

| Method | <= 31 cells | 32-159 | 160-999 | >= 1,000 |
| --- | --- | --- | --- | --- |
| Default | 8 tested, 0% / 0% | 12, 67% / 50% | 21, 43% / 33% | 49, 25% / 18% |
| Markov/Ward | 363, 0% / 0% | 867, 25% / 17% | 505, 26% / 14% | 103, 36% / 34% |

Separability decides most exclusions (separable 0-31% for 32-999 cells on Markov/Ward); the tie share leaves
47% (32-159) and 35% (160-999) of Markov/Ward splits.

## Reading

On Markov/Ward, rule_out roughly doubles the lines recovered (41 to 80), raises ARI from 0.55 to 0.77 and
removes 290 clusters, but 664 clusters remain for 190 lines: lines are still split (median 2 clusters per
line) and 253 clusters have fewer than 20 cells. Whether the remaining within-line splits are biology or false
cannot be read from these data, which have no known no-group parts. The default tree does not change: its giant
cluster (28,611-32,680 cells) comes from the topology, which no test of a split can open. count_all alone is
worse than today on real data (1,605 clusters, 24 lines) because small splits pass, as on simulated data.
