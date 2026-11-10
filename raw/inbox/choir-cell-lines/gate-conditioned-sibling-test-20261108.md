# Sibling test conditioned on the edge gate having opened

Analysis on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/gate_conditioned_sibling_test_20261108` (git-ignored); protocol, deviations and scripts
captured in `gate-conditioned-sibling-test-scripts-20261108.py.txt`. Protocol frozen before the full run; no
deviations.

## Question

Does the sibling test of large splits become calibrated on no-group data, without losing real groups, when
its p-value is computed given that the split's edge gate opened?

## Why the gate can be written as a cutoff on the sibling score

- For continuous data the parent mean is exactly the pooled mean of its two sides (checked: maximum absolute
  difference 0.0), so both child-parent contrasts are fixed multiples of the sibling contrast
  d = mean(A) - mean(B).
- Scaling d by r (moving the two sides apart rigidly, parent mean fixed) multiplies both edge statistics and
  the sibling statistic by r^2, with the whitening covariance, PCA basis, branch lengths and the
  energy-fraction projection dimension held fixed, as today's chi-square law already holds them.
- Tree-BH tests the children of each parent as one family at a level set by the ancestors. Ancestor edges,
  deeper edges and family levels do not change under the move, so the gate is a function of s = r^2 alone,
  and it is monotone. A copy of Tree-BH reproduced production rejections in 5,300/5,300 fits.
- s* = the value of s where the gate flips, by bisection. Tested splits: gate open at s = 1 in
  45,952/45,952 evaluations; references with own tested edges: closed in 12,296/12,296.
- Every tested split of 40 or more cells has sibling df 30 = number of features, so T is the full squared
  length of the whitened contrast. Given the contrast's direction, T ~ c * chi2(30) truncated to
  T >= s* T_obs. This is exact under the fixed-quantity assumptions above. For splits with df < 30 (mostly
  under 40 cells) it is an approximation.

## Arms

Production code runs end to end; only each tested split's p-value is swapped.

| Arm | p-value |
| --- | --- |
| today | production (replay 11,488/11,488) |
| gate | sf(T) / sf(s* T), no correction for the tree's choice |
| gate_c | sf(T/c) / sf(s* T/c), today's c |
| gate_both | as gate_c, c from a truncated likelihood (own-tested references observed below their gate); without truncation it reproduces today's c (11,488/11,488) |

Data: the 900 no-group count datasets (400 cells; sd 0.25 `low`, 0.42 `matched_median`, 1.0 `high`; 300
each) and 35 real-group datasets (balanced, imbalanced, singleton: 1,000 cells, 4 groups, 10 each; density
variation: 3,000 cells, 4 groups, 5). Both trees, production alphas.

## E1 Tested no-group splits of 40 or more cells with p <= 0.01 (target about 1%)

| Tree | sd | today | gate | gate_c | gate_both |
| --- | --- | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 71.8% | 87.5% | 63.5% | 63.5% |
| Markov/Ward | 0.42 | 73.1% | 90.2% | 63.4% | 63.3% |
| Markov/Ward | 1.0 | 84.8% | 91.4% | 76.2% | 76.0% |
| Default | 0.42 (163 tested) | 54.6% | 73.0% | 32.5% | 32.5% |
| Default | 1.0 | 91.6% | 90.3% | 79.2% | 78.4% |

Markov/Ward sd 0.25 by size: 10-39 cells 1.0% (today) to 0.0% (gate_c); 40-159 cells 61.0% to 49.4%;
160-400 cells 94.3% to 92.8%. Median factor: today 2.5, gate_both 2.5 (the truncated likelihood moved it by
less than 0.5%).

## E2 False partitions of 300

| Tree | sd | today | gate | gate_c | gate_both |
| --- | --- | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 298 | 300 | 297 | 297 |
| Markov/Ward | 0.42 | 296 | 300 | 293 | 293 |
| Markov/Ward | 1.0 | 300 | 300 | 300 | 300 |
| Default | 0.25 | 2 | 5 | 0 | 0 |
| Default | 0.42 | 50 | 68 | 34 | 34 |
| Default | 1.0 | 259 | 261 | 259 | 258 |

## E3 Real groups, mean ARI

| Tree | Signal | today | gate | gate_c | gate_both |
| --- | --- | ---: | ---: | ---: | ---: |
| Default | balanced | 0.800 | 0.800 | 0.800 | 0.800 |
| Default | imbalanced | 1.000 | 1.000 | 1.000 | 1.000 |
| Default | singleton | 0.880 | 0.880 | 0.880 | 0.880 |
| Default | density variation | 0.651 | 0.663 | 0.651 | 0.651 |
| Markov/Ward | balanced | 0.315 | 0.276 | 0.327 | 0.327 |
| Markov/Ward | imbalanced | 0.079 | 0.073 | 0.092 | 0.092 |

On real splits the gate sits far below the observed score (s* about 1e-4), so conditioning leaves their
p-values unchanged.

## E4 Where the gate sits on no-group tested splits (median s*)

| Tree | sd | 10-39 | 40-159 | 160-400 |
| --- | --- | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 0.83 | 0.64 | 0.41 |
| Markov/Ward | 0.42 | 0.83 | 0.62 | 0.46 |
| Default | 0.42 | 0.69 | 0.84 | 0.74 |

Median score T/df of these splits: 13.5 to 14.2 (10-39 cells), 4.8 to 6.4 (40-159), 7.8 to 13.6 (160-400).

## Decision rule (fixed before results)

- Fixes large splits (Markov/Ward, sd 0.25 and 0.42, splits >= 40 cells pass <= 2%): no arm.
- Keeps real groups (default tree, ARI within 0.02 of today on every signal model): all arms.

## Reading

Conditioning on the gate is exact here and removes part of the excess: with today's c it cuts large-split
passes from 72-73% to 63% on Markov/Ward and from 55% to 33% on the default tree at sd 0.42, and false
partitions on the default tree from 2 to 0 (sd 0.25) and 50 to 34 (sd 0.42), without changing real-group
results. Large no-difference splits still score 5 to 14 times the one-population value while the gate opens
at about half the observed score, so most of the excess comes from the tree's choice of the split, not from
the gate. Without a correction for that choice (arm gate) the test is worse than today. The truncated
likelihood barely changes c, so the gate is not why the reference splits score low. The reference splits
that drive c are mostly edge-blocked and carry no truncation of their own (2026-11-08 depth/weight/paths
analysis).

## Correction (2026-11-08, after the results)

The derivation above is wrong in one step. The edge test does not use all of a parent's cells:
`split_edge_test_rows` alternates each child's cells between a test half and a basis half, fits the edge PCA
basis on the basis half and computes the child-parent contrast on the test half only. The sibling test uses
all cells. So the edge statistics are not fixed multiples of the sibling statistic, and s* * T is not where
the gate cuts the sibling score. The "parent = pooled mean" check above used all-cell means, not the
test-half means.

Check (`halves_check.py`, 92 tested Markov/Ward splits of >= 40 cells, `low` replicates 0-9): the test-half
contrast points almost the same way as the full contrast (whitened cosine median 0.95, range 0.90-0.98) and
less so as the basis-half contrast (0.80). Holding the basis half fixed and scaling the test-half contrast
until the gate flips, the gate cuts the sibling score at a median 0.72 of the observed score (range 0.54 to
0.99), where the analysis used s* = 0.47 (0.28 to 0.98).

Consequences:

- The arms `gate`, `gate_c` and `gate_both` are approximations that put the cut too low, so they correct too
  little. The tables above stand as results of that approximation, not of an exact conditional test.
- An exact version must condition on the basis half and on the test half's direction. Then the gate is a
  cutoff on the test half's length (which scales exactly), and the sibling score is a quadratic function of
  that length with a fixed basis-half part, not a scaled chi-square(30).
- Expected direction, not measured: an exact version removes more of the excess than the approximation.
  For a typical no-group large split (T about 300, c about 2.5, df 30) a cut at 0.72 T instead of 0.47 T
  still leaves p far below 0.01, so the reading that the tree's choice carries most of the excess is
  unchanged.
