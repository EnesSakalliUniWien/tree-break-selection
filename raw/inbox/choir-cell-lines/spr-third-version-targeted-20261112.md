# Sprs only at refused nodes where the third version sees a better split, CHOIR tree

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/spr_third_version_targeted_20261112`
(git-ignored); protocol, both scripts and their outputs are in
`spr-third-version-targeted-scripts-20261112.py.txt`. Follow-up of `spr-third-version-20261112.md`. No production
code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

Under the third version of the draft definition of k the best cut along a node's direction passes at 32 of
1,302 refused nodes of the unchanged Markov/Ward tree, 22 with cells from several cell lines and 10 from one
line. Do sprs at these nodes only help, where the sprs of 11 November at 166 nodes did not?

## Setup

Set (a): the 22 nodes with several lines, chosen with the labels, so a diagnostic and not a rule. Set (b): all
32, without labels. Target split of a node: the third version's direction from all cells of the node and the
tree's split, the best cut of all cells along it, its smaller side; that side is made a child of the node by
sprs with the code of 11 November. Nodes are processed from the largest to the smallest; a node whose cells an
earlier spr changed is skipped. Refit with the third version, data order. The protocol with three predictions
was written before the run.

## The sprs

| | Set (a) | Set (b) |
| --- | ---: | ---: |
| Nodes chosen | 22 | 32 |
| Processed | 17 | 25 |
| Skipped because an earlier spr changed their cells | 5 | 7 |
| Sprs | 677 | 757 |
| Sprs per node, median (largest) | 11 (348) | 11 (348) |
| Share of cells whose line is cut by the split, median, before to after | 0.065 to 0.047 | 0.069 to 0.054 |
| Nodes where that share falls / rises | 9 / 3 | 12 / 6 |
| New split has p <= 0.01 under the third version | 10 of 17 | 14 of 25 |

In set (b) the new split passes at 4 of the 8 processed nodes from one line and at 10 of the 17 with several
lines. Both trees are rooted binary trees on the same leaves.

## The refit, third version

| Tree | Found clusters | ARI | Lines recovered | Split only | Merged only | Both | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Unchanged | 206 | 0.8359 | 127 | 16 | 30 | 17 | 150 |
| Set (a): sprs at 17 nodes with several lines | 213 | 0.8266 | 126 | 17 | 28 | 19 | 154 |
| Set (b): sprs at 25 nodes, no labels | 218 | 0.8257 | 126 | 18 | 27 | 19 | 156 |
| All sprs of 11 November (166 nodes) | 258 | 0.7691 | 127 | 42 | 13 | 8 | 193 |

The unchanged tree's ARI over five orders of the rows is 0.8349 to 0.8411; both trees after the targeted sprs
are below that range.

## Predictions

| Prediction | Result |
| --- | --- |
| W1 set (a): ARI at least 0.8359 and fewer than 30 lines merged only | not held (ARI 0.8266; merged only 28) |
| W2 set (b): ARI at least 0.8259 | not held (0.8257) |
| W3 the new split has p <= 0.01 at more than half of the processed nodes | held (10 of 17, 14 of 25) |

## Reading

- Sprs at the nodes where the third version itself sees a better split do not help either, not even when the
  nodes are chosen with the labels: 2 to 3 merged lines are separated, 1 to 2 more are split, and the ARI falls
  by about 0.01.
- Together with the run on the sprs of 11 November: on this tree and with this statistic the sprs do not improve
  the result. The alternative splits they produce pass the test more often than they follow the cell lines.

## Not checked

Other orders of the rows. The default tree, where the topology holds far fewer lines as subtrees. A target
split in another direction than that of the tree's split. Sprs that move whole clades between distant parts of
the tree instead of rearranging one node.
