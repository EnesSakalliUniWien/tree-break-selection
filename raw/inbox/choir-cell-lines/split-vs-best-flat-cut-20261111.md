# How far each tested split is from the best flat cut in the subtree's own invariant directions

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/split_vs_best_flat_cut_20261111`
(git-ignored); protocol (written before the run), scripts and summary output are in
`split-vs-best-flat-cut-scripts-20261111.py.txt`. Descriptive; no production change.

## Question

The tree is built once from the root's diffusion distances (inherited geometry). The Gaussian best-flat-cut
reference measures `rho` in `invariant_directions` of the subtree's own cells (rebuilt geometry). How far is
the split the tree offers from the best flat cut in those directions, and which of the two respects the cell
lines?

## Setup

48,879 CHOIR cells, 30 PCs, both registered methods refitted with `sibling_best_flat_cut_reference=True` to
capture the cells of every tested split. Replay: 138 and 53 clusters, 1,838 and 90 tested splits, as in
`best-flat-cut-reference-real-20261111.md`. For each of the 1,475 (Markov/Ward) and 82 (default tree) tested
splits with at least one invariant direction, in those directions:

- tree split: its `rho` (equal to `split_rho` in every case);
- nearest flat cut: start at the tree split, re-sort along its mean difference and move single cells until
  stable;
- best flat cut: the table's search (grid, top 50 re-sorts, single-cell moves); a lower bound for k above 10;
- agreement: adjusted Rand index (ARI) between two partitions of the subtree's cells;
- reference p-values where n >= 15 and k <= 10;
- line cut share: share of the subtree's cells on the minority side of their own cell line.

## Results by tree and subtree size (medians unless counts)

| | M/W 32-159 | M/W 160-999 | M/W 1,000+ | Default 32-159 | Default 160-999 | Default 1,000+ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Splits | 867 | 505 | 103 | 12 | 21 | 49 |
| k | 6 | 9 | 14 | 6 | 11 | 14 |
| `rho`, tree split | 0.97 | 1.00 | 3.13 | 2.49 | 0.73 | 0.11 |
| `rho`, nearest flat cut | 3.19 | 2.65 | 4.97 | 5.32 | 6.35 | 3.99 |
| `rho`, best flat cut | 4.19 | 3.33 | 7.20 | 6.17 | 7.94 | 8.15 |
| ARI, tree split and nearest flat cut | 0.305 | 0.303 | 0.769 | 0.626 | 0.476 | 0.020 |
| ARI, tree split and best flat cut | 0.036 | 0.020 | 0.077 | 0.232 | 0.100 | -0.009 |
| Cells moved to reach the nearest flat cut | 20.6% | 21.0% | 4.1% | 3.1% | 6.8% | 16.2% |
| Tree split is itself a flat cut (count) | 24 | 19 | 4 | 6 | 4 | 0 |
| Best flat cut equals the tree split (count) | 10 | 12 | 2 | 4 | 1 | 0 |
| ARI with best flat cut at least 0.9 (count) | 14 | 52 | 12 | 4 | 5 | 3 |
| Best flat cut has fewer than 5 cells on one side (count) | 38 | 6 | 0 | 0 | 0 | 0 |

The 37 splits accepted by BH in the two fits have a median ARI of 0.968 with the best flat cut; 10 of them
are flat cuts themselves. For the 250 splits with k above 10 the median `rho` is 2.05 for the tree split and
at least 8.14 for the best flat cut, with median ARI 0.034.

## Reference p-values, tree split against best flat cut (n >= 15, k <= 10, p <= 0.01)

| Tree, size | Splits | Both | Best flat cut only | Neither | Tree split only |
| --- | ---: | ---: | ---: | ---: | ---: |
| Markov/Ward, 32-159 | 867 | 15 | 101 | 751 | 0 |
| Markov/Ward, 160-999 | 418 | 42 | 120 | 256 | 0 |
| Default, 32-159 | 12 | 3 | 4 | 5 | 0 |
| Default, 160-999 | 10 | 3 | 6 | 1 | 0 |

"Tree split only" is 0 by construction, because no split has a larger `rho` than the best flat cut. Of the 231
subtrees where only the best flat cut passes, 20 have a best flat cut with fewer than 5 cells on one side.

## Cell lines (Markov/Ward, n >= 15, k <= 10; added after the endpoints)

Subtrees holding one cell line only (765):

| Size | Subtrees | Best flat cut p <= 0.01 | Best flat cut p <= 0.05 | Tree split p <= 0.01 | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32-159 | 570 | 38 | 58 | 1 | 1 |
| 160-999 | 195 | 48 | 60 | 2 | 0 |

Subtrees holding two or more cell lines (520):

| Size | Subtrees | Best flat cut p <= 0.01 | Tree split p <= 0.01 | Significant best flat cut cuts no line | Line cut share at most 5% | Mean line cut share, best / tree | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| 32-159 | 297 | 78 | 14 | 38 | 62 | 0.056 / 0.136 | 4 |
| 160-999 | 223 | 114 | 40 | 42 | 86 | 0.074 / 0.123 | 26 |

In 138 subtrees with two or more lines the best flat cut passes and the tree split does not; there the line
cut share is at most 5% for the best flat cut in 95 and for the tree split in 44, and the best flat cut's
share is lower in 98.

Over all tested splits with two or more lines in the subtree (not restricted to significant ones) the best
flat cut does not cut fewer lines on average: mean line cut share tree / best 0.175 / 0.241 (32-159, 297
splits), 0.112 / 0.171 (160-999, 299), 0.010 / 0.008 (1,000 and more, 98).

## Reading

- The split the tree offers is, in most subtrees, unrelated to the best flat cut in the subtree's own
  invariant directions: median ARI 0.02 to 0.08 on Markov/Ward, and the tree split is itself a flat cut in 47
  of 1,475 cases. The tree's split has a median `rho` of about 1 in subtrees below 1,000 cells, the best flat
  cut 3.3 to 4.2.
- Where the tree split and the best flat cut agree, the split is accepted: the accepted splits have a median
  ARI of 0.968 with the best flat cut.
- In subtrees with several cell lines the rebuilt geometry holds a significant split far more often than the
  tree offers one: 192 of 520 against 54 of 520. Of the 192 significant best flat cuts, 148 cut at most 5% of
  the cells from their line and 80 cut no line. This measures the hidden boundaries of the 9 November merge
  check on every tested split.
- The Gaussian null does not hold inside single cell lines: the best flat cut is significant at 0.01 in 38 of
  570 (6.7%) and 48 of 195 (24.6%) single-line subtrees, where about 1% would be expected if a line were one
  Gaussian population. Whether this is non-Gaussian shape or real structure inside a line cannot be told from
  these data. The tree's own split is significant in 3 of the 765.
- Default tree: at 1,000 cells and more the tree split has a median `rho` of 0.11 and ARI -0.009 with the
  best flat cut; its large splits are unrelated to the feature-space structure of their cells.
- Expectation check: low agreement in small subtrees held. Higher agreement in large subtrees held only for
  the nearest flat cut (ARI 0.769, 4.1% of cells moved on Markov/Ward), not for the best flat cut (0.077).

## Limits

- One dataset, two trees; nested subtrees, so the splits are dependent and the counts are not independent
  trials.
- Cell lines are the only truth; structure inside a line counts as a cut line.
- The best flat cut is found by search; for k above 10 it is a lower bound. Splits with k above 10 (250) and
  subtrees without invariant directions (371) have no reference p-value here.
- The best flat cut's p-value uses the reference of the best flat cut, so it needs no further correction
  inside a subtree; across subtrees no BH step was applied.
- Not run: fits that use the best flat cut as the split; other alternative partitions (second-best cuts,
  splits in single spectral blocks); the distance fit of the alternatives.
