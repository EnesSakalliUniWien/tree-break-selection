# Alternative sprs on the tree against the best flat cut (CHOIR cell lines)

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/spr_alternatives_best_flat_cut_20261111`
(git-ignored); protocol (written before the run, amended once after the self-check and before the CHOIR run),
scripts, self-check output and summary output are in `spr-alternatives-best-flat-cut-scripts-20261111.py.txt`.
Descriptive; no gate or default was changed. The sprs were not applied to the fit.

## Setup

48,879 CHOIR cells, 30 PCs, `tbs_family_diffusion_markov_ward` with `sibling_best_flat_cut_reference=True`
(uncommitted working-tree code on top of `bc044e07`). The fit was repeated to save the complete tree; replay:
138 clusters, 1,838 tested splits, as before. The default tree was not run. 1,475 tested splits have 32 cells or
more; 1,285 of them have 1 to 10 invariant directions and are measured; the other 190 have 11 to 15 (168 to
48,879 cells; every split with 1,000 cells or more is among them) and are outside the stored reference.

Identity: every node is its leaf set (`Partition`, taken over from BranchArchitect into
`tree_break_selection/tree/elements` on the same day). `subtree_prune_regraft` in
`tree_break_selection/tree/spr.py` is the spr on a `PosetTree`.

## What an spr does to a tested split

Node P with children A and B. An spr prunes a clade C below one child and regrafts it on the other side. P keeps
its cells; its split becomes (A minus C | B plus C). P's invariant directions, k, n, best flat cut and reference
null do not change, so every spr alternative is one more split of the same cells and its `rho` cannot exceed the
best flat cut's `rho`. Only the children change their cells.

- Greedy sprs: move the clade (leaves included) that raises `rho` most; repeat until no spr raises it.
- Targeted sprs: the fewest clades of the fitted tree that must change sides so that P's split equals the best
  flat cut (the maximal clades lying completely on the wrong side).

## Self-check (simulated, 120 cells, tree built from noise coordinates)

The bound held and replaying every move with `subtree_prune_regraft`, finding the pruned clade by its leaf set,
gave the same child leaf sets as the bookkeeping in 3 of 3 datasets. With two groups at offset 6 the greedy
sprs stopped after 15 and 16 moves at `rho` 1.78; the best flat cut has `rho` 9.2 to 9.9 (ARI of the greedy
result with the groups 0.02). Greedy sprs end at a flat cut near the tree split.

## Results (1,285 splits)

Medians of `rho`, and counts passing the reference at 0.01:

| Subtree size | Splits | Tree split | First spr | Greedy sprs | Best flat cut | Pass: tree | Pass: greedy | Pass: best |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32-159 | 867 | 0.97 | 1.60 | 3.19 | 4.19 | 15 | 77 | 116 |
| 160-999 | 418 | 0.92 | 1.49 | 2.61 | 3.22 | 42 | 102 | 162 |

| Verdict of the reference | Splits | Pass after greedy sprs | Median greedy sprs | Median ARI greedy/best | Median targeted sprs | Median cells moved (targeted) | Line cut share tree / greedy / best |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Tree split and best flat cut pass | 57 | 57 | 2 | 1.000 | 2 | 0.010 of n | 0.008 / 0.000 / 0.000 |
| Only the best flat cut passes | 221 | 122 | 11 | 0.743 | 16 | 0.336 of n | 0.201 / 0.204 / 0.096 |
| Neither passes | 1,007 | 0 | 11 | 0.071 | 20 | 0.400 of n | 0.289 / 0.378 / 0.415 |

- The bound: none of the 1,007 splits whose best flat cut is refused passes after greedy sprs. In 4 of 1,285
  splits (k = 9 or 10, n = 171 to 261) the greedy sprs reached a higher `rho` than the table's search (3.18
  against 2.92, 3.07 against 3.06, 2.94 against 2.92, 5.09 against 4.52). The bound is for the true best flat
  cut, so the search missed it there; no verdict changed in these 4.
- Greedy sprs end exactly at the best flat cut (ARI above 0.999) in 260 of 1,285 splits. No split reached the
  cap of 1,000 moves.
- Of the 221 splits where only the best flat cut passes, greedy sprs make the tree's split pass in 122. Sprs
  until the first pass: median 3, quartiles 1 to 8, largest 195; 1 spr in 38 (median size of that clade 0.092
  of n), 2 to 5 in 46, 6 to 20 in 29, more in 9. Cells moved until the first pass: median 0.267 of n. In 96 of
  the 122 the greedy end equals the best flat cut. The other 99 stop at median `rho` 2.56 against 3.99 for the
  best flat cut (median ARI 0.027).
- Cell lines in the 122 rescued splits: median line cut share 0.144 for the tree split, 0.038 at the first
  pass, 0.021 at the greedy end and 0.021 for the best flat cut. With several lines in the subtree: 94 of 138
  rescued, 66 of them cut at most 5% of the cells from their line at the first pass (74 with the best flat
  cut). With one line: 28 of 83 rescued, 2 at most 5%.
- Targeted sprs: single cells are half or more of the moved clades (median share 0.50 to 0.71) and the largest
  moved clade holds 0.03 to 0.06 of n. For the 138 only-best splits with several lines: 5 sprs or fewer in 54,
  20 or fewer in 98, median 9.
- Moved clades of 10 cells or more with a majority line of at least 90%: greedy 883 of 953 (all splits) and
  175 of 229 (only-best); targeted 466 of 500 and 180 of 206. Such clades hold 39,255 of 67,992 cells moved by
  greedy sprs and 11,615 of 66,589 moved by targeted sprs.
- Reaction below the split (122 rescued splits, 244 children before and after the greedy sprs): the best flat
  cut of a child passes in 92 of 170 testable children before and in 25 of 180 after. No child has a passing
  best flat cut in 46 splits before and in 101 after. k is unchanged in 94 of 244 children; children without
  invariant directions: 73 before, 64 after.

## Reading

- At the split itself the best flat cut does not react to sprs: same cells, same directions, same reference. It
  is the upper bound of every spr alternative, so its verdict decides whether any spr can make the split pass.
- Below the split it reacts: after the sprs the children hold different cells, and a passing best flat cut is
  left in far fewer of them (25 against 92). The structure the best flat cut saw in the children was moved up
  into the parent's split.
- The tree holds the best flat cut's sides as a few subtrees in the 57 agreeing splits and in part of the
  splits with several lines (54 of 138 need 5 sprs or fewer). Elsewhere most moves are single cells.

## Limits

One dataset and one tree; dependent splits; per-split p-values without the BH step; the sprs were not applied to
the fit, so the edge gate, pass-through and the complete partition were not measured; greedy order, so the spr
counts are upper bounds of the fewest needed for a pass; the best flat cut is found by search (missed in at
least 4 splits); splits with more than 10 invariant directions, which include every split with 1,000 cells or
more, are not covered; cell lines are the only truth; the rescued splits were selected by the same data that
gives their p-value, and the reference does not cover that selection.
