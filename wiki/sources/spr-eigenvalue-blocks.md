---
title: Sprs per Eigenvalue Block of the Invariant Directions
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-20261111.md
  - raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-scripts-20261111.py.txt
  - tree_break_selection/eigenvalue_groups.py
  - tree_break_selection/tree/spr.py
tags:
  - real-data
  - subtree
  - spectral
  - topology
---

# Sprs per Eigenvalue Block of the Invariant Directions

## Summary

The between/total ratio of a split and the transfer cost of an spr are sums
over the eigenvalue blocks of the node's invariant directions. On the 1,285
measured splits of the Markov/Ward CHOIR tree every invariant direction is its
own numerical block. A passing split lives in the first block (share 0.89);
failing tree splits spread over the blocks, and the rescuing sprs move the
split into the first block. Two passing blocks in one node are rare (9
nodes): in about half the second cut is a second cell-line boundary, in the
others it runs through the lines.

## Key Points

- 9,250 blocks, median 7 per node, no numerical tie; 781 of 7,965
  neighbouring pairs have an eigenvalue ratio below 1.1.
- Additivity: r2 = sum of r2_block and dW = sum of dW_block, exact to
  rounding error.
- First block's share of r2: 0.89 where tree split and best flat cut pass;
  in the rescued class 0.42 for the tree split and 0.73 for the greedy sprs
  and the best flat cut; 0.23 to 0.25 for tree splits elsewhere.
- Own cut of a block passing at 0.01 / k: 139 in block 1, 25 in block 2, 6
  in block 3, 2 later; 172 cuts in 162 nodes, median ARI 0.854 with the
  joint best flat cut and 0.242 with the tree split.
- 28 of the 1,007 nodes whose joint best flat cut is refused have one block
  whose own cut passes at 0.01 / k; about 10 are allowed under one Gaussian
  population per node.
- 11 pairs of passing blocks in 9 nodes, all with ARI below 0.1. Inspected
  against the cell lines: in 4 nodes the two cuts isolate different lines
  (both real, held by a tree as two successive splits); in 4 nodes the
  second cut divides every line in similar proportions (a second partition
  the lines do not explain); one node is mixed and one has two cuts inside
  one line. A low ARI alone does not show that two cuts exclude each other.
- N97405 written as leaf sets (`find_node_by_split`): the tree splits it
  into N96707 (HCC1419 219 of 220 leaves) and N97129, and N97129 into
  N96510 (BT474 122 of 136) and N96138 (MDAMB361 122 of 122); both splits
  are accepted. Block 1's cut is the first tree split with 31 leaves on the
  other side, block 2's cut is the second tree split seen from the parent,
  and block 3's cut is a union of 96 and 66 nodes, held nowhere in the tree.
- All 9 nodes as leaf sets, 19 block cuts: 3 are exactly the node's tree
  split, 6 are compatible with it (unions of 2 to 16 tree nodes), 10 are not
  compatible. The 7 cuts through the lines are unions of 12 to 196 tree
  nodes per side. Where the fit refuses a split although a line boundary
  exists (TE14 in N96475, KYSE520 in N96842 and N96808), that boundary is a
  union of 2 tree nodes.
- Sprs applied with `subtree_prune_regraft` in N96475, N96842 and N96808
  (2 sprs each, 3 or 4 leaf sets changed): the sibling reference of the
  node's split goes from p 1.0, 0.0445 and 0.021 to 0.0005 in all three,
  with `rho` equal to the joint best flat cut's (7.651, 10.305, 13.693) and
  a line cut share of 0. The fit was not rerun, so edge gate, BH step and
  pass-through are not measured.
- An spr that raises `rho` is helped by 55.4% of the blocks; all blocks
  agree in 0.6%; the first block agrees in 61.3% (75.3% for clades of 10
  cells or more).
- Cells on the wrong side per block: 27.4% in block 1, 41.7% in later
  blocks; tied cells 2.8% to 1.4%.
- Limits: blocks are chosen from the same cells; near-equal eigenvalues make
  single directions arbitrary; one direction per block; one dataset and tree.

## Evidence

`raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-20261111.md` holds the
trace, tables, reading and limits; protocol, scripts, self-check output and
summary output are in
`raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-scripts-20261111.py.txt`.
The block rule is `tied_eigenvalue_neighbors` in
`tree_break_selection/eigenvalue_groups.py`.

## Links

- [[spr-subspaces-equidistances]]
- [[spr-alternatives-best-flat-cut]]
- [[spr-refit]]
- [[subtree-descriptions-and-statistical-tests]]
- [[subtree-distance-topology-families]]
- [[subspace-distance-ties]]
