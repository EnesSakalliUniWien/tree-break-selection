# The spr measure per eigenvalue block of the invariant directions (CHOIR cell lines)

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/spr_eigenvalue_blocks_20261111` (git-ignored);
protocol (written before the run), scripts, self-check output and summary output are in
`spr-eigenvalue-blocks-scripts-20261111.py.txt`. Descriptive; no gate or default was changed and nothing was
applied to the fit. Same data, saved Markov/Ward tree and 1,285 tested splits as
`spr-alternatives-best-flat-cut-20261111.md`; `rho` of the tree split and of the best flat cut and the verdict
of the greedy sprs equal that run in the 30 splits compared before the full run.

## Blocks and trace

Blocks are the numerical eigenvalue groups of `tied_eigenvalue_neighbors`. The 1,285 nodes have 9,250 invariant
directions (median 7 per node) and no numerical tie, so every invariant direction is its own block. 781 of the
7,965 neighbouring pairs have an eigenvalue ratio below 1.1 (in 618 nodes); single directions inside such a
pair are not identified by the data. No statistical test of equal eigenvalues exists in the code.

The invariant directions are uncorrelated in the node, so whitening scales each one by its own standard
deviation and squared distances add over blocks. For every split r2 = sum of r2_block, and for every spr
dW = sum of dW_block (self-check: differences 0.0 and 8.5e-14). The self-check on four simulated groups on a
2 x 2 grid finds two blocks whose own cuts pass and cross (ARI -0.002, each cut equal to one grid axis); under
one Gaussian population 2 of 401 blocks in 200 datasets pass at 0.01 / k.

## Results

Share of the first block in r2, and the block with the largest share:

| Class | Splits | First block's share: tree split | greedy sprs | best flat cut | Tree split and best flat cut have the same largest block |
| --- | ---: | ---: | ---: | ---: | ---: |
| All | 1,285 | 0.27 | 0.43 | 0.27 | 396 |
| Tree split and best flat cut pass | 57 | 0.89 | 0.88 | 0.88 | 52 |
| Only best passes, greedy sprs pass | 122 | 0.42 | 0.73 | 0.73 | 65 |
| Only best passes, greedy sprs fail | 99 | 0.23 | 0.31 | 0.32 | 19 |
| Neither passes | 1,007 | 0.25 | 0.36 | 0.20 | 260 |

The block alone (one direction, reference at (n, 1); passing at 0.01 / k, in brackets at 0.01):

| Block | Blocks | Median `rho`: tree split | own best cut | Tree split passes | Own cut passes | Cells tied | Cells on the wrong side |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1,285 | 0.12 | 1.95 | 39 (47) | 139 (217) | 2.8% | 27.4% |
| 2 | 1,285 | 0.04 | 1.79 | 1 (1) | 25 (60) | 2.3% | 36.6% |
| 3 | 1,285 | 0.02 | 1.71 | 0 (0) | 6 (19) | 1.8% | 38.9% |
| later | 5,395 | 0.01 | 1.63 | 0 (0) | 2 (28) | 1.4% | 41.7% |

Greedy sprs inside one block end at that block's own cut (passing counts 217, 60, 18 of 19, 28 at 0.01).

Nodes by the number of blocks whose own cut passes at 0.01 / k:

| Class | Splits | 0 blocks | 1 block | 2 blocks | 3 or more |
| --- | ---: | ---: | ---: | ---: | ---: |
| All | 1,285 | 1,123 | 153 | 8 | 1 |
| Tree split and best flat cut pass | 57 | 10 | 42 | 4 | 1 |
| Only best passes, greedy sprs pass | 122 | 56 | 64 | 2 | 0 |
| Only best passes, greedy sprs fail | 99 | 78 | 19 | 2 | 0 |
| Neither passes | 1,007 | 979 | 28 | 0 | 0 |

- 172 own cuts pass in 162 nodes. Their median ARI is 0.854 with the joint best flat cut and 0.242 with the
  tree split; median line cut share 0.067; 9 isolate fewer than 5 cells.
- Cell lines: with several lines in the subtree 131 cuts pass in 122 of 520 nodes, and 71 of them cut at most
  5% of the cells from their line. With one line 41 cuts pass in 40 of 765 nodes.
- 28 nodes whose joint best flat cut is refused have one block whose own cut passes at 0.01 / k (135 at 0.01).
  Bonferroni over a node's blocks allows about 1% of nodes under one Gaussian population per node, about 10 of
  1,007.
- Crossing: 11 pairs of passing blocks in 9 nodes; all 11 have an ARI below 0.1 (median 0.004). Where a node
  holds two passing blocks, the two cuts are unrelated partitions, which one split of a tree cannot hold.
- Sprs per block: of the 47,423 clades whose spr raises `rho`, 55.4% of the blocks agree (dW_block < 0); all
  blocks agree in 0.6%; the first block agrees in 61.3% (75.3% for clades of 10 cells or more, 84.5% in the
  rescued class).
- Equidistance per block: tied cells fall from 2.8% (block 1) to 1.4% (later blocks); cells on the wrong side
  rise from 27.4% to 41.7%. Where tree split and best flat cut both pass, block 1 has 4.8% of the cells on the
  wrong side and the other blocks 38.8% to 46.0%.

## Reading

- A split that passes lives in one block: the first block holds 0.89 of r2 where tree split and best flat cut
  agree. The tree split of a failing node spreads over the blocks (first block 0.23 to 0.42).
- The sprs that rescue a split move it into the first block (share 0.42 to 0.73, the best flat cut's share).
- In the other blocks the tree split carries no side information: about 40% of the cells are on the wrong
  side, near the 50% of an unrelated partition. So an spr that helps in one block is opposed by about half of
  the others, and the joint transfer cost is a sum of terms with both signs.
- Nodes with two passing blocks are rare (9 of 1,285) and their cuts are different partitions (ARI below 0.1).
  Corrected after the inspection below: a low ARI does not mean the two cuts exclude each other. In about half
  of the pairs they are two cell-line boundaries that a tree holds as two successive splits.

## Inspection of the 9 nodes with two passing blocks (follow-up, same day)

`inspect_crossing.py` cross-tabulates the cell lines against the two cuts (output `crossing.txt` in the scripts
file). A line is "whole" on a side with at least 90% of its cells in the node there; otherwise it is cut.

| Node | Cells | Lines in the node | First cut | Second cut | Kind |
| --- | ---: | --- | --- | --- | --- |
| N97405 | 478 | HCC1419 219, MDAMB361 135, BT474 122 (breast) | block 1: BT474 against HCC1419 (MDAMB361 115 to 20) | block 2: MDAMB361 against the other two, no line cut | two line boundaries |
| N97405 | 478 | same | blocks 1 and 2 | block 3: cuts HCC1419 140 to 79 and MDAMB361 89 to 46 | through the lines |
| N97237 | 163 | COLO792 57, SKMEL30 34, COLO741 26, SKMEL5 23, A2058 16 (skin) | block 1: COLO792 against A2058 and COLO741 | block 2: SKMEL5 against the rest, no line cut | two line boundaries |
| N96842 | 76 | TE10 32, SNU46 30, KYSE520 6 | block 1: KYSE520 against the rest | block 3: SNU46 against TE10 | two line boundaries |
| N96808 | 43 | SNU46 29, KYSE520 6, single cells of 6 lines | block 1: KYSE520 against the rest | block 2: 7 single cells of 5 lines against the rest | two line boundaries |
| N97713 | 368 | OVK18 178, SBC5 102, COV434 88 | block 1: OVK18 against the other two, no line cut | block 2: COV434 against SBC5, cuts OVK18 125 to 53 | line boundary plus a cut line |
| N97705 | 523 | NCIH2077 445, TOV112D 78 | block 1: the two lines, no line cut | block 2: cuts NCIH2077 298 to 147 and TOV112D 31 to 47 | through the lines |
| N94004 | 84 | VMCUB1 73, BICR56 11 | block 1: the two lines, no line cut | block 2: cuts VMCUB1 39 to 34 and BICR56 6 to 5 | through the lines |
| N96475 | 114 | SNU1214 101, TE14 10 | block 1: the two lines, no line cut | block 3: cuts SNU1214 43 to 58 and TE14 6 to 4 | through the lines |
| N93303 | 127 | SNU1196 127 (one line) | block 1: cuts the line 43 to 84 | block 6: cuts the line 75 to 52 | both through one line |

- Two line boundaries (4 nodes): each cut isolates a different line, and the lines fall into three of the four
  combinations of sides. Both cuts are real separations and a tree holds them as two successive splits. They are
  not alternatives to each other.
- Through the lines (4 nodes, counting N97405 by its third block): the first cut is the line boundary; the second
  cut divides the cells of each line in similar proportions, so all four combinations of sides are filled. This
  is a second partition of the same cells that the cell lines do not explain and that no tree holding the line
  boundary can also hold as one node.
- N97713 is mixed; N93303 has two unrelated cuts inside one line.
- In the fit, the nodes N97705, N97405, N97713 and N97237 are split (sibling p 0.0005, BH different); N94004 has
  sibling p 0.0005 without BH different; the other four are not split (p 1.0, 1.0, 0.0445, 0.021).

Not checked: what the second cut through the lines corresponds to (for example a state shared by cells of
different lines); whether the fit separates the second line boundary further down the tree.

## N97405 as leaf-set splits (follow-up, same day)

`node_splits.py` builds the subtree of N97405 as a `PosetTree`, gives every node its leaf set
(`initialize_split_indices`), writes the tree split and the three block cuts as `Partition`s and looks each
side up with `find_node_by_split`; a side that is not a node is written as the union of the largest tree nodes
inside it (`PartitionSet.maximals_under`). Output `splits_N97405.txt` in the scripts file; one row per leaf in
the local file `leaves_N97405.csv`.

The tree below N97405 (478 leaves: HCC1419 219, MDAMB361 135, BT474 122, EFM192A 2), with the fit's verdicts:

| Node | Leaves | Cell lines | Sibling test of its split |
| --- | ---: | --- | --- |
| N97405 | 478 | all | p 0.0005, BH different |
| - N96707 | 220 | HCC1419 219, EFM192A 1 | p 1.0 (not split) |
| - N97129 | 258 | MDAMB361 135, BT474 122, EFM192A 1 | p 0.0005, BH different |
| - - N96510 | 136 | BT474 122, MDAMB361 13, EFM192A 1 | p 1.0 (not split) |
| - - N96138 | 122 | MDAMB361 122 | p 1.0 (not split) |

So the fit ends with three sets here: N96707, N96510 and N96138. 15 of the 478 leaves are not with their line
(13 MDAMB361 leaves inside N96510 and the 2 EFM192A leaves).

How the tree holds each line: HCC1419 is the union of 6 nodes (largest N95969 with 177 leaves), MDAMB361 of 10
nodes (N96138 with 122, then 4, 2 and 7 single leaves), BT474 of 22 nodes (largest N94337 with 44). Moving the 9
MDAMB361 nodes outside N96138 (13 leaves) by sprs would make MDAMB361 one node.

The block cuts as leaf sets:

| Cut | Low side | High side | A node of the tree? |
| --- | --- | --- | --- |
| Block 1 | 249: BT474 122, MDAMB361 115, HCC1419 11, EFM192A 1 | 229: HCC1419 208, MDAMB361 20, EFM192A 1 | no; unions of 21 and 23 nodes |
| Block 2 | 143: MDAMB361 126, HCC1419 16, BT474 1 | 335: HCC1419 203, BT474 121, MDAMB361 9, EFM192A 2 | no; unions of 25 and 54 nodes |
| Block 3 | 343: HCC1419 140, BT474 112, MDAMB361 89, EFM192A 2 | 135: HCC1419 79, MDAMB361 46, BT474 10 | no; unions of 96 and 66 nodes, most of them single leaves |

- Block 1's cut is the tree split of N97405 with 31 leaves on the other side (11 HCC1419, 20 MDAMB361). Block
  2's cut is the tree split of the child N97129 seen from N97405, with 26 leaves on the other side. The tree's
  own two splits are purer than the two block cuts.
- Exact compatibility (one of the four intersections of the sides empty) fails for every pair, by 11 to 63
  leaves. The reading "two line boundaries that a tree holds" is true for the lines and approximate for the
  leaves.
- Block 3's cut is not held by the tree anywhere: its sides are unions of 96 and 66 nodes.
- This answers one open point for this node: the fit does separate the second line boundary, one level down.

## All 9 nodes as leaf-set splits (follow-up, same day)

`node_splits.py` was run for the other 8 nodes (outputs `splits_<node>.txt` in the scripts file; one row per
leaf in the local files `leaves_<node>.csv`). "Compatible" means one of the four intersections of the sides is
empty, so one tree can hold both splits.

| Node | Tree split (leaf sets) | Sibling test in the fit | Block cuts against the tree |
| --- | --- | --- | --- |
| N97705 | N97553: NCIH2077 443; N97587: TOV112D 78, NCIH2077 2 | p 0.0005, BH different | block 1 (424 / 99): compatible, moves 19 NCIH2077 leaves to the TOV112D side; block 2: not compatible, unions of 196 and 150 nodes |
| N97405 | N96707: HCC1419 219, EFM192A 1; N97129: MDAMB361 135, BT474 122, EFM192A 1 | p 0.0005, BH different | none compatible (11 to 63 leaves); block 3: unions of 96 and 66 nodes |
| N97713 | N97204: OVK18 178; N97616: SBC5 102, COV434 88 | p 0.0005, BH different | block 1 is exactly the tree split; block 2: not compatible, unions of 56 and 53 nodes |
| N97237 | N96797: COLO792 57, SKMEL30 33 and 8 more; N96890: COLO741 26, SKMEL5 21, A2058 15 and 3 more | p 0.0005, BH different | block 1: not compatible by 10 leaves; block 2 (SKMEL5 21 / rest): compatible, a union of 5 nodes |
| N93303 | N88150: SNU1196 60; N90669: SNU1196 67 | p 1.0 | neither compatible; unions of 14 to 44 nodes |
| N96475 | N94813: TE14 10, SNU1214 6, BICR16 3; N95374: SNU1214 95 | p 1.0 | block 1 (TE14 10 / rest): compatible, TE14 is a union of 2 nodes; block 3: not compatible, unions of 34 and 38 nodes |
| N94004 | N90482: BICR56 11; N92777: VMCUB1 73 | p 0.0005, not BH different | block 1 is exactly the tree split; block 2: not compatible, unions of 12 and 16 nodes |
| N96842 | N96516: TE10 32, SNU46 1; N96808: SNU46 29, KYSE520 6 and 8 more | p 0.0445 | block 1 (KYSE520 6, BICR16 1 / rest) and block 3 (SNU46 28 and 3 more / rest): both compatible with the tree split |
| N96808 | N96070: 7 leaves of 5 lines; N96190: SNU46 29, KYSE520 6, BICR16 1 | p 0.021 | block 2 is exactly the tree split; block 1 (KYSE520 6, BICR16 1 / rest): compatible, a union of 2 nodes |

Counts over the 19 block cuts of the 9 nodes:

- 3 are exactly the node's tree split (N97713, N94004, N96808).
- 6 are compatible with the tree split without being a node (unions of 2 to 16 tree nodes on the smaller side).
- 10 are not compatible with the tree split: the 7 cuts through the lines or inside one line (unions of 12 to
  196 tree nodes per side) and 3 cuts that follow a line boundary less purely than the tree (N97405 blocks 1
  and 2, N97237 block 1).

Lines that are one tree node here: OVK18 (N97204), the 73 VMCUB1 leaves (N92777), BICR56 (N90482), SNU1196's
127 leaves (N93303). NCIH2077 is N97553 plus 2 single leaves.

Reading: where a block cut follows a cell-line boundary, the tree already has that boundary as a node or as a
union of a few nodes, and in 4 nodes the tree's split is at least as pure as the block cut. Where the tree's
split is refused although a line boundary exists (N96475: TE14 against SNU1214; N96842 and N96808: KYSE520
against SNU46), the boundary is a union of 2 tree nodes, so 1 or 2 sprs would make it a node. The cuts through
the lines are spread over tens to hundreds of tree nodes and are not sprs of the tree in any useful sense.

## Sprs applied in N96475, N96842 and N96808 and the splits re-tested (follow-up, same day)

Follow-up protocol written before the run (in the scripts file). `apply_sprs.py` builds the node's subtree as a
`PosetTree` with set based identity, finds the two tree nodes that make up the smaller side of block 1's cut
(`PartitionSet.maximals_under`), joins them with `subtree_prune_regraft` and regrafts the joined node on the
edge above the node. The sibling reference is `best_flat_cut_reference` (per split; the fit is not rerun).

| Node | Leaves | Tree split before | k | `rho`, p before | Sprs | Split after | `rho`, p after | Leaf sets changed |
| --- | ---: | --- | ---: | --- | ---: | --- | --- | --- |
| N96475 | 114 | SNU1214 95 against TE14 10, SNU1214 6, BICR16 3 | 9 | 2.221, 1.0 | 2 | TE14 10 against SNU1214 101, BICR16 3 | 7.651, 0.0005 | 3 of 114 |
| N96842 | 76 | SNU46 29, KYSE520 6 and 8 more against TE10 32, SNU46 1 | 8 | 5.174, 0.0445 | 2 | KYSE520 6, BICR16 1 against the other 69 | 10.305, 0.0005 | 4 of 76 |
| N96808 | 43 | SNU46 29, KYSE520 6, BICR16 1 against 7 leaves of 5 lines | 7 | 8.352, 0.0210 | 2 | KYSE520 6, BICR16 1 against the other 36 | 13.693, 0.0005 | 3 of 43 |

- The sprs: in N96475, prune N89472 (TE14 4) onto the edge above N88659 (TE14 6), then prune the joined node
  onto the edge above N96475. In N96842 and N96808, prune N88548 (KYSE520 2) onto the edge above N90950
  (KYSE520 4, BICR16 1), then prune the joined node onto the edge above the node.
- All three new splits pass at 0.01 with the lowest possible p. Their `rho` equals the joint best flat cut's
  `rho` of the earlier run (7.651, 10.305, 13.693): the side is the node's best flat cut, so the reference,
  which is the null of the best flat cut, covers the choice of the split. It does not cover the choice of k or
  of the node.
- The line cut share goes from 0.053, 0.013 and 0.000 to 0.000 in all three.
- Only 3 or 4 leaf sets change per node; every other node keeps its identity. The node ids do not follow: after
  the sprs the leaf set of N96475 is held by the node id N94554, and the id N96475 holds the 104 other leaves.
- The split of the new rest child is refused in N96475 (95 against 9, p 1.0) and in N96842 (36 against 33, SNU46
  against TE10, k 7, `rho` 4.664, p 1.0); in N96808 it has p 0.0135 (29 against 7).

Not measured: the edge gate, the BH step (17 of 48 splits at p 0.0005 were not accepted by BH in the earlier
real-data run) and pass-through, because the fit was not rerun; branch lengths after the sprs are the stated
shortcut of `subtree_prune_regraft`.

## Limits

One dataset and one tree; dependent splits; blocks are chosen from the same cells that are then cut, and the
reference does not cover that choice; Bonferroni covers only the choice of block within one node; 781
neighbouring pairs are near-equal, so their single directions are arbitrary; one direction per block, so a
separation along a mixture of blocks is not seen by a single block; splits with more than 10 invariant
directions are not covered; nothing was applied to the fit.
