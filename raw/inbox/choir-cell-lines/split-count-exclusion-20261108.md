# Counting all possible splits with exclusions from subspaces and equidistant cells

Analysis on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/split_count_exclusion_20261108` (git-ignored); protocol, code and summary output in
`split-count-exclusion-scripts-20261108.py.txt`. Protocol frozen before the run; no deviations. Self-checks:
replay 11,488/11,488, cell counts 11,488/11,488.

## Question and arms

The user asked to count all possible splits of a subtree's cells and exclude possibilities based on the numbers
of different subspaces and equidistant relationships, and to test which reading of "exclude" makes more sense.
p_split = chi2.sf(T, df) of the production sibling statistic, no inflation factor. Separable: a flat cut in the
subtree's invariant directions (eigenvectors above its noise band) puts the two children on different sides;
k = 0 (including every subtree of <= 31 cells, p = 30) is not separable. Equidistant cells: margin |m| < 0.1 in
the test metric; tie share = their fraction.

| Arm | p-value |
| --- | --- |
| today | production |
| count_all | min(1, (2^(n-1) - 1) p_split) |
| rule_out | count_all, and p = 1 unless separable with tie share <= 0.05 |
| shrink | count of flat cuts in the invariant directions (if separable) or of all splits, over n minus the equidistant cells, times 2 |
| both | shrink, and p = 1 unless separable with tie share <= 0.05 |

The 0.05 and 0.1 cut-offs are conventions fixed in the protocol before any result.

## Results

Tested no-group splits with p <= 0.01 (Markov/Ward, sd 0.25 / 0.42 / 1.0):

| Arm | 10-39 cells | >= 40 cells | False partitions of 300 | Median clusters |
| --- | --- | --- | --- | --- |
| today | 1.0% / 3.1% / 9.3% | 71.8% / 73.1% / 84.8% | 298 / 296 / 300 | 8 |
| count_all | 92.9% / 85.4% / 90.1% | 14.3% / 10.6% / 37.6% | 279 / 246 / 297 | 7 |
| rule_out | 0.7% / 1.6% / 2.3% | 0.0% / 0.9% / 0.4% | 2 / 25 / 15 | 1 |
| shrink | 92.5% / 87.0% / 88.4% | 20.9% / 17.0% / 46.4% | 294 / 265 / 299 | 7-8 |
| both | 1.4% / 3.1% / 2.9% | 1.8% / 3.1% / 0.8% | 48 / 77 / 25 | 1 |

Default tree: false partitions 2 / 50 / 259 (today) to 0 / 0 / 0 (rule_out), 0 / 3 / 2 (both); large splits at
sd 0.42 54.6% to 0.0% (rule_out).

Real groups, mean ARI (median clusters):

| Tree | Data | today | count_all | rule_out | shrink | both |
| --- | --- | --- | --- | --- | --- | --- |
| Markov/Ward | balanced (4) | 0.315 (20) | 0.582 | 1.000 (4) | 0.484 | 0.989 (4) |
| Markov/Ward | imbalanced (4) | 0.079 (25.5) | 0.137 | 1.000 (4) | 0.110 | 0.813 (4) |
| Markov/Ward | singleton (4 + 1 cell) | 0.092 (22) | 0.157 | 0.879 (4) | 0.134 | 0.595 (8) |
| Markov/Ward | density variation (4) | 0.025 (81) | 0.043 | 0.867 (5) | 0.032 | 0.700 (6) |
| Default | balanced / imbalanced / singleton | 0.800 / 1.000 / 0.880 | same | same | same | same |
| Default | density variation | 0.651 (13) | 0.652 | 0.667 (7) | 0.670 | 0.685 (8) |

Decision rule (frozen): rule_out fixes large splits, fixes small splits and keeps real groups. count_all and
shrink fail both calibration parts; both fails the large-split part at sd 0.42 (3.1%) and the "within 0.02"
real-group part only because the default density-variation ARI rose by 0.034.

## What does the excluding

Share of today's tested no-group splits that rule_out keeps as real (Markov/Ward): 3.2% (10-39 cells), 4.1%
(40-159), 0% (160-400). Separable in the invariant directions: 3.2% / 4.9% / 0%. Tie share <= 0.05: 100% /
27% / 0.2%. For subtrees under about 40 cells the whitened test metric (30 directions) leaves no cell
equidistant, so there only separability decides.

## Limits

- rule_out uses the observed split (separability, equidistant cells), so its error control is measured, not
  guaranteed by construction.
- Subtrees of <= 31 cells have no invariant directions here (p = 30), so rule_out can never accept a split
  inside them, real or not.
- Markov/Ward no-group data at sd 0.42 still split 25 of 300 datasets (8%); at sd 1.0, 15 (sequencing depth puts
  structure into the table there).
- Not run: rule_out with today's factor or without the count, and the real CHOIR data end to end.

## Reading

Excluding possibilities as real separations works; shrinking the number of possibilities does not. The
count of all splits handles the largest noise splits; requiring that a split lie in the subtree's invariant
directions and leave almost no cell equidistant removes nearly all remaining noise splits at every size,
including the small ones the count let through, while every simulated real group is kept and Markov/Ward
recovers the true groups (ARI 0.87-1.00).
