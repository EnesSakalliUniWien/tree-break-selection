# The sprs under the third version of the draft definition of k, CHOIR tree

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/spr_third_version_20261112` (git-ignored);
protocol, deviations, both scripts and their outputs are in `spr-third-version-scripts-20261112.py.txt`.
Follow-up of `variant-a-fixed-20261112.md` and `spr-refit-20261111.md`. No production code was changed.
Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

On 11 November sprs were applied at 166 refused splits whose best flat cut passed under the rule of that day
(4,371 sprs, `tree_spr.npz`). With that rule the unchanged tree was better than the tree after the sprs. The
third version of the draft takes its one direction from the split itself, so an spr changes the direction too.
What do the sprs do under the third version?

## Setup

48,879 CHOIR cells, 30 PCs, Markov/Ward tree, third version (20 random divisions per node, twice the median
p-value, exact null for one direction), data order of the rows. The protocol with three predictions was written
before the run.

## T1: the fit on the tree after the sprs

| | Found clusters | ARI | Lines recovered | Split only | Merged only | Both | Largest cluster | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Third version, unchanged tree | 206 | 0.8359 | 127 | 16 | 30 | 17 | 1,990 | 150 |
| Third version, tree after the sprs | 258 | 0.7691 | 127 | 42 | 13 | 8 | 1,451 | 193 |
| Rule of 11 November (k up to 15), unchanged tree | 230 | 0.8219 | 116 | 21 | 38 | 15 | 1,990 | 134 |
| Rule of 11 November (k up to 15), tree after the sprs | 316 | 0.7376 | 107 | 58 | 19 | 6 | 1,451 | 210 |

Under the third version the sprs lower the ARI from 0.8359 to 0.7691. They separate merged lines (merged only 30
to 13) and split lines internally (split only 16 to 42); the number of recovered lines stays at 127.

## T2: the 166 nodes where sprs were applied

110 of the 166 candidates have the cells of a node of the unchanged tree; the other 56 are nodes of the tree as
it was already changed by sprs above them. Third version, nodes with p <= 0.01:

| Candidates | All 166, after the sprs | Of those that are nodes of the unchanged tree | Before the sprs | After the sprs |
| --- | ---: | ---: | ---: | ---: |
| All | 77 | 110 | 2 | 60 |
| Cells from one cell line | 12 of 58 | 38 | 0 | 10 |
| Cells from several cell lines | 65 of 108 | 72 | 2 | 50 |

The smaller side of the split grows from 28 to 61 cells at the median. The directly computed p-value equals the
recorded one of the fit at the 86 candidates that were called and tested in the unchanged tree.

## T3: room for sprs under the third version itself

For every called node of the unchanged tree with at least 30 cells that the third version refuses, the same 20
divisions and directions with `rho` of the best cut of the test cells along the direction, twice the median:

| Refused nodes | Nodes | Best cut along the direction has p <= 0.01 |
| --- | ---: | ---: |
| All | 1,302 | 32 (2.5%) |
| At least 90% of the cells from one line | 1,099 | 10 |
| Several lines | 203 | 22 |
| 30 to 99 cells | 618 | 8 |
| 100 to 999 cells | 652 | 21 |
| 1,000 cells and more | 32 | 3 |

## Predictions

| Prediction | Result |
| --- | --- |
| U1 the fit on the tree after the sprs has a lower ARI than the unchanged tree | held (0.7691 against 0.8359) |
| U2 at the 166 nodes more splits have p <= 0.01 after the sprs than before | held (2 to 60 at the 110 comparable nodes) |
| U3 fewer than 10% of the refused nodes have a passing best cut along their direction | held (2.5%) |

## Reading

- The third version follows the sprs: a split moved to the best flat cut of 11 November gets a new direction and
  passes at 60 of 110 nodes where it was refused before.
- That is not a gain as a whole. The sprs of 11 November were applied wherever a best flat cut passed, including
  inside one cell line; there 10 of 38 refused splits pass afterwards, and lines are split internally.
- The sprs do one useful thing: they separate lines the tree had merged (merged only 30 to 13).
- Under the third version itself there is little left for sprs: at 32 of 1,302 refused nodes a better split
  exists along the node's direction, 22 of them at nodes with several lines.

## Deviations

After the sprs a candidate's cells sit at the node recorded as `holder_id`, so "after" is read there. "Before"
exists only for the 110 candidates that are nodes of the unchanged tree. Both were found when the script's own
checks stopped it.

## Not checked

Sprs derived under the third version, for example only at the 22 refused nodes with several lines. Other orders
of the rows for the tree after the sprs. The default tree. T3 looks only along the direction of the tree's
split, so a separation in another direction is not counted.
