# Correction for the tree's choice from split counts and invariant directions

Analysis on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/split_count_correction_20261108` (git-ignored); protocol (with one amendment), deviations,
scripts and summary output in `split-count-correction-scripts-20261108.py.txt`.

## Question and arms

User-confirmed idea: build the sibling test's correction for the tree's choice from (1) how many splits are
possible for the known leaf set of a subtree and (2) how many lie in its invariant directions (eigenvalue above
the subtree's noise band) versus its fuzzy directions. p_split = chi2.sf(T, df), production statistic, no
inflation factor.

- `count_all`: min(1, (2^(n-1) - 1) p_split), a union bound over every split of the leaf set.
- `count_invariant`: if the split is linearly separable in the k invariant directions (linear program),
  min(1, 2 M_inv p_split) with M_inv = sum_{i=0}^{k} binom(n-1, i) hyperplane splits (Cover 1965, from
  secondary sources); otherwise min(1, 2 (2^(n-1) - 1) p_split).

Amendment before the full run: production's per-subtree Marchenko-Pastur count is 0 for every subtree here
(each subtree is whitened by its own covariance first), so invariant directions were taken from the subtree's
covariance in the shared coordinates, with the noise level read from the eigenvalue bulk (median matched to
the Marchenko-Pastur median). Replicate 0 check: no-group root 0 directions, four-group root 3.

Data: the 900 no-group and 35 real-group datasets of 2026-11-08, both trees, production alphas. Replay
11,488/11,488.

## Results: tested no-group splits with p <= 0.01

| Tree | sd | Size | today | count_all | count_invariant |
| --- | --- | --- | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 10-39 | 1.0% | 92.9% | 92.5% |
| Markov/Ward | 0.25 | 40-159 | 61.0% | 21.0% | 23.5% |
| Markov/Ward | 0.25 | 160-400 | 94.3% | 0.4% | 0.4% |
| Markov/Ward | 0.42 | 10-39 | 3.2% | 85.4% | 87.0% |
| Markov/Ward | 0.42 | 40-159 | 63.3% | 16.1% | 19.9% |
| Markov/Ward | 0.42 | 160-400 | 91.5% | 0.2% | 0.2% |
| Markov/Ward | 0.25 / 0.42 | >= 40 | 71.8% / 73.1% | 14.3% / 10.6% | 16.0% / 13.1% |
| Default | 0.42 | >= 40 (163) | 54.6% | 1.2% | 3.1% |

False partitions of 300: Markov/Ward 298 / 296 (today) to 279 / 246 (`count_all`); default tree 2 to 2
(sd 0.25), 50 to 12 (sd 0.42), 259 to 219 (sd 1.0). Real groups: default-tree mean ARI unchanged (0.800,
1.000, 0.880; density variation 0.651 to 0.652 / 0.670); Markov/Ward ARI rises (balanced 0.315 to 0.582).

Decision rule (frozen): no arm fixes large splits (Markov/Ward >= 40 cells pass 11-16%); both keep real groups.

## Invariant directions

Share of tested no-group splits that are separable in their subtree's invariant directions: 0-9% (median k
0-7 for 40-159-cell subtrees, 0-1 for 160-400). Noise splits almost never follow invariant directions, so
`count_invariant` falls back to the all-splits count; on real-group data the real splits are separable in 1-3
invariant directions (replicate 0 check).

## Best split of one population (descriptive)

For tested no-group splits whose smaller side holds >= 25% of the cells, observed T against (2/pi) n + df - 1
(the share of variance a split at the mean of a normal captures along its axis is 2/pi; arithmetic, not a
cited result): median ratio 1.1-1.3 for 160-400 cells at sd 0.25/0.42 (Markov/Ward), 1.9-2.2 for 40-159, 3.7-4.1
for 10-39.

## Reading

Counting every split of the known leaf set is enough for the largest subtrees: the count grows as 2^n while a
noise split's score grows about like (2/pi) n, so the bound wins for large n. It fails for small subtrees,
where the score is already about four times what a best split of one population gives, so the fixed-split
chi-square law itself does not hold there and no count can repair it. Today's estimated factor handles the
small subtrees (1% at 10-39 cells) and fails the large ones; the count handles the large ones and fails the
small. Invariant directions add nothing on no-group data because noise splits do not follow them. Not run: a
rule combining the two.
