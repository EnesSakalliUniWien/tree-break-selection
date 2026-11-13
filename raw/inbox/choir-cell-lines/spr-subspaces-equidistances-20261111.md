# How the sprs relate to the subspaces and the equidistances (CHOIR cell lines)

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/spr_subspaces_equidistances_20261111`
(git-ignored); protocol (written before the run), scripts, self-check output and summary output are in
`spr-subspaces-equidistances-scripts-20261111.py.txt`. Descriptive; no gate or default was changed and nothing
was applied to the fit. Same data, saved Markov/Ward tree and 1,285 tested splits (32 cells or more, 1 to 10
invariant directions) as `spr-alternatives-best-flat-cut-20261111.md`; `rho` of the tree split, of the greedy
sprs and of the best flat cut equal that run in all 1,285 splits.

## Trace

Node P with children A (a cells) and B (b cells); a clade C (c cells) below A; centroids mu in a Euclidean space
S. The spr of C to B's side changes the within-children sum of squares by the transfer cost

    dW_S(C) = c b / (b + c) |mu_C - mu_B|^2  -  c a / (a - c) |mu_C - mu_A|^2 .

1. In P's whitened invariant directions the between/total ratio of the split changes by exactly -dW / n
   (self-check: largest deviation 5.2e-17 over 236 clades). An spr raises `rho` if and only if dW < 0.
2. dW = 0 is equidistance of the clade: the spr leaves `rho` unchanged and `rho` cannot rank the two
   alternatives. For one cell this is the cell margin of 8 November, m = (distance to the other centroid -
   distance to its own) / distance between the centroids; tied means |m| < 0.1 (convention).
3. Every clade has a margin in every space. The tree was built in the tree's space (30 Markov diffusion
   coordinates, Ward). The sprs that raise `rho` are the clades that the invariant directions put nearer the
   other side.
4. Ward joined A's children A1 and A2 at the lowest merge cost. The tie ratio (cost of joining a grandchild with
   the other child / cost of the merge made) measures how close the tree was to the spr of that grandchild.

Spaces: tree; pcs (30 PCs); invariant (P's invariant directions, whitened: the space of `rho`); fuzzy (PCs with
the invariant directions projected out).

## Results

Cells, pooled over the splits (tied / wrong side):

| Class | Splits | tree | pcs | invariant | fuzzy | Side differs, tree against invariant |
| --- | ---: | --- | --- | --- | --- | ---: |
| All | 1,285 | 3.2% / 2.4% | 11.8% / 12.0% | 12.3% / 10.3% | 25.7% / 32.2% | 10.7% |
| Tree split and best flat cut pass | 57 | 0.8% / 0.8% | 2.3% / 1.8% | 2.1% / 1.2% | 28.5% / 43.7% | 1.5% |
| Only best passes, greedy sprs pass | 122 | 3.1% / 2.7% | 9.3% / 13.7% | 12.1% / 9.3% | 25.8% / 33.6% | 9.8% |
| Only best passes, greedy sprs fail | 99 | 3.4% / 2.5% | 13.2% / 14.2% | 14.3% / 11.9% | 25.9% / 31.3% | 11.9% |
| Neither passes | 1,007 | 3.5% / 2.5% | 13.2% / 12.5% | 13.3% / 11.4% | 25.3% / 30.6% | 11.8% |

Clades whose spr raises `rho` against all other clades (all splits):

| Clades | Group | n | Tied in the tree's space | Wrong side in the tree's space | Tied in the fuzzy space |
| --- | --- | ---: | ---: | ---: | ---: |
| Single cells | spr raises `rho` | 24,766 | 10.5% | 8.7% | 21.0% |
| Single cells | other | 163,560 | 2.1% | 1.4% | 26.4% |
| 10 cells or more | spr raises `rho` | 4,485 | 4.9% | 2.7% | 3.5% |
| 10 cells or more | other | 23,109 | 1.0% | 0.5% | 11.5% |

Cells moved by the sprs against the cells that stay (all splits; tied / wrong side):

| Sprs | Cells moved | tree: moved | tree: stay | pcs: moved | pcs: stay |
| --- | ---: | --- | --- | --- | --- |
| Greedy | 50,054 | 5.7% / 4.6% | 2.3% / 1.6% | 23.0% / 30.8% | 7.7% / 5.1% |
| Targeted (to the best flat cut) | 62,372 | 4.6% / 3.6% | 2.5% / 1.8% | 17.6% / 21.5% | 8.9% / 7.3% |

Grandchild sprs by Ward tie ratio in the tree's space (5,052 grandchildren; median ratio 2.96; 1 below 1):

| Tie ratio | Grandchildren | Spr raises `rho` | Split passes after this one spr (failing splits) |
| --- | ---: | ---: | --- |
| below 1.1 | 175 | 43.4% | 2 of 168 |
| 1.1 to 1.5 | 687 | 37.8% | 8 of 669 |
| 1.5 to 3 | 1,689 | 26.8% | 7 of 1,635 |
| 3 and more | 2,501 | 12.8% | 6 of 2,356 |

All 23 passes after one grandchild spr lie in the 122 splits that the greedy sprs rescue.

Sprs derived without the best flat cut and without a search, in one pass: "tied sprs" move every maximal clade
that is tied in the tree's space and has dW < 0 in the invariant directions; "untied sprs" move every maximal
clade with dW < 0 that is not tied in the tree's space.

| Split | Passing at 0.01 | Median `rho` | Cells moved (all splits) | Median line cut share |
| --- | ---: | ---: | ---: | ---: |
| Tree split | 57 | 0.95 | 0 | 0.268 |
| Tied sprs | 61 | 1.12 | 6,218 | 0.282 |
| Untied sprs | 74 | 1.12 | 81,779 | 0.331 |
| Greedy sprs | 179 | 3.00 | 50,054 | (see the spr capture) |
| Best flat cut | 278 | 3.87 | 62,372 | (see the spr capture) |

Tied sprs exist in 837 splits (median 1 clade, 0.018 of the cells) and never lower `rho`. Untied sprs lower
`rho` below the tree split in 489 splits, because all clades move at once.

## Reading

- The relationship is exact in the algebra: an spr's change of `rho` is the clade's transfer cost in the
  invariant directions, and equidistance is the spr that changes nothing.
- Equidistance in the tree's space is related to the sprs but is a small part of them. Cells whose spr raises
  `rho` are tied in the tree's space five times as often as other cells (10.5% against 2.1%), yet at least 81%
  of them are neither tied nor on the wrong side there. The tree placed them clearly; the invariant
  directions place them on the other side.
- The larger part is the change of space. The tree's space and the invariant directions put 10.7% of the cells
  on different sides (1.5% where tree split and best flat cut both pass). The moved cells are tied or on the
  wrong side in the PCs (23.0% / 30.8% against 7.7% / 5.1% for the cells that stay).
- Near-tied merges mark the sprs in direction: 43.4% of near-tied grandchild sprs raise `rho` against 12.8% at
  a ratio of 3 and more. One grandchild spr rarely makes a split pass (23 of 4,828).
- Sprs derived from equidistance in the tree's space alone add 4 passing splits (61 against 57). The fuzzy
  space carries no side information (25.7% tied, 32.2% on the wrong side), as on 8 November.

## Deviation

Protocol item 3 names "cells that differ between the tree split and the best flat cut". The side assignment
with fewer differing cells was used; the earlier spr run chose the assignment with fewer clades.

## Limits

One dataset and one tree; dependent splits; the tie cuts 0.1 and 1.1 are conventions, not tests; centroid
distances only, so equidistance to a curved group is not seen; the saved diffusion coordinates are taken as
the tree's space; eigenvalue blocks inside the invariant directions were not separated; per-split p-values
without the BH step; splits with more than 10 invariant directions are not covered; nothing was applied to
the fit.
