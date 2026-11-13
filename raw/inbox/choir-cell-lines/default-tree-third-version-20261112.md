# The third version of the draft definition of k on the default tree of the CHOIR cells

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/default_tree_third_version_20261112`
(git-ignored); protocol, scripts and outputs are in `default-tree-third-version-scripts-20261112.py.txt`.
Follow-up of `variant-a-fixed-20261112.md`. No production code was changed. Environment: Python 3.12.11, NumPy
2.3.4, SciPy 1.16.3, uv 0.9.5.

## Setup

Default tree (`tbs_family_diffusion`: symmetric operator, average linkage) of the 48,879 CHOIR cells, rebuilt as
a parent array from the linkage saved on 6 November. Check of the conversion: the Markov/Ward tree rebuilt the
same way from its saved linkage equals the saved Markov/Ward tree. Fits through the harness of
`spr_refit_20261111` with the opt-in reference on: the current rule (table, k up to 15) and the third version in
the data order of the rows and two random orders. The protocol with three predictions was written before the run.

## Result

| | Found clusters | ARI | NMI | Lines recovered | Split only | Merged only | Both | Largest cluster | Calls | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Current rule | 64 | 0.0275 | 0.5495 | 29 | 4 | 145 | 12 | 28,611 | 90 | 16 |
| Third version, data order | 63 | 0.0275 | 0.5493 | 29 | 4 | 145 | 12 | 28,611 | 90 | 18 |
| Third version, random order 1 | 63 | 0.0275 | 0.5493 | 29 | 4 | 145 | 12 | 28,611 | 90 | 18 |
| Third version, random order 2 | 63 | 0.0275 | 0.5493 | 29 | 4 | 145 | 12 | 28,611 | 90 | 18 |

For comparison on the Markov/Ward tree: 1,838 calls, current rule 230 clusters and ARI 0.8219, third version 206
to 215 clusters and ARI 0.8349 to 0.8411. Today's default test on the default tree: 75 clusters, ARI 0.0275.

Predictions: X1 held (29 lines recovered under both rules); X2 held (ARI 0.0275 under both); X3 held (63 clusters
in all three orders).

## Why nothing changes (`giant.txt`)

- The default tree is a chain. From the root down, the first 14 calls split 1, 151, 170, 185, 5, 196, 226, 208,
  250, 249, 267, 278, 288 and 293 cells from a remainder of 46,000 to 48,877 cells.
- Such a split has a small `rho` because one side holds less than 1% of the cells: under the current rule all
  14 have p = 1.0 (k 13 or 14); under the third version 2 of them pass (250 against 47,486 and 267 against
  46,970, p = 0.00002) and 12 have p = 1.0. The traversal passes through the refused ones.
- The sibling test is called at only 90 nodes of the whole tree (1,838 on the Markov/Ward tree).
- The largest cluster has 28,611 cells from 144 cell lines. The node above it (32,679 cells, 4,069 against
  28,610) is accepted. The next split below takes 4 cells from 28,606, and the sibling test is not called there
  or at any of the next 14 steps down the larger side, although those steps separate 2,525, 917, 11,579 and
  other cells.

## Reading

- On the default tree the rule for k and the sibling statistic are not the limit: both rules give the same
  clusters. The limits are the topology (a chain of small groups peeled from one large remainder; 42 of the 190
  lines are clean subtrees) and the gate before the sibling test, which stops at a peel of 4 cells above 28,606
  cells.
- A split that peels one small group from a very large remainder is refused by `rho` under both rules even when
  the small group may be a whole cell line (151 cells against 48,726, p = 1.0).

## Follow-up on the same day: the edge test's values at the peel that stops the tree (`edge_values.py`, `edge_values.txt`)

The fit with the third version was run again with the edge and sibling columns of every node saved (63 clusters,
ARI 0.0275 as before). Edge alpha 0.001.

- Of the 97,756 edges of the tree, Tree-BH tests 382 and finds 171 significant; 6,006 edges have a raw p-value
  of at most 0.001.
- N97671 (32,679 cells, 4,069 against 28,610): both edges significant (statistics 28,558 and 36,084, k = 5), the
  sibling test accepts.
- N97667, the next node (4 against 28,606): both edges have statistic 44.4 with k = 30, raw p 0.044, tested and
  not significant. The gate closes here and the sibling test is skipped.
- Everything below is not tested because this ancestor is closed, although the edges there are far beyond the
  threshold: 2,525 against 26,081 has edge statistics 15,117 and 22,769 (k = 6); 129 against 25,952 has 1,044
  (k = 19, raw p 1.8e-209); 917 against 25,035 has 11,255 and 11,194; 173 against 24,862 has 2,350; 279 against
  24,580 has 2,202.
- All nodes with at least 1,000 cells (155): 106 have no significant edge. 5 of them were tested and failed, and
  their smaller side has 1, 1, 4, 3 and 14 cells. The other 101 were never tested because an ancestor is closed;
  94 of these have a raw p-value of at most 0.001 on at least one edge.
- The branch multiplier at these edges is 1.00 to 1.83.

Reading: the default tree is stopped by peels of very few cells. An edge to 4 cells cannot be significant at
0.001 with 30 directions, and Tree-BH then closes the whole subtree below it. The code lets a peel of one cell
pass this gate; peels of 3, 4 or 14 cells do not pass.

## Not checked

Which cell lines the peeled groups are. What the result would be if peels of a few cells passed the gate like
a peel of one cell. The tree is the one of 6 November, not rebuilt today. Sprs
on this tree.
