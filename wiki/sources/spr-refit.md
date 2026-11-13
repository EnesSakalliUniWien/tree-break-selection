---
title: Fit Rerun after sprs at the Refused Splits with a Passing Best Flat Cut
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/spr-refit-20261111.md
  - raw/inbox/choir-cell-lines/spr-refit-scripts-20261111.py.txt
  - raw/inbox/choir-cell-lines/spr-third-version-20261112.md
  - raw/inbox/choir-cell-lines/spr-third-version-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/spr-third-version-targeted-20261112.md
  - raw/inbox/choir-cell-lines/spr-third-version-targeted-scripts-20261112.py.txt
  - tree_break_selection/tree/spr.py
tags:
  - real-data
  - sibling-gate
  - topology
  - subtree
---

# Fit Rerun after sprs at the Refused Splits with a Passing Best Flat Cut

## Summary

A top-down pass applies sprs at every node of the Markov/Ward CHOIR tree whose
own split is refused by the sibling reference while its best flat cut passes
(166 nodes, 4,371 sprs), and the fit is rerun on the changed tree. Accepted
splits rise from 31 to 104 and found clusters from 138 to 257; the ARI against
the cell lines rises from 0.2975 to 0.3425. It stays low because the nodes
with more than 10 invariant directions, which include every node with 1,000
cells or more, are still refused.

## Key Points

- Harness: the fit runs on a tree given as a parent array; the unchanged tree
  reproduces the baseline exactly (138 clusters, 1,838 tested splits,
  assignment ARI 1.0).
- Spr pass: 5,654 nodes evaluated, 194 with more than 10 invariant directions
  left alone, 166 candidates (108 with several lines, 58 with one line);
  2,831 of the 4,371 sprs move a single leaf; 8,760 of 97,757 leaf sets
  change and none with 1,000 leaves or more.
- Refit: lines merged only 91 to 51, split only 9 to 47, recovered 76 to 82,
  largest cluster 7,678 to 6,346, clusters at least 90% pure 98 to 202, NMI
  0.8512 to 0.8846.
- Of the 166 changed splits the refit tested 145, all at p <= 0.01; BH
  accepted 76 (58 of 97 with several lines, 18 of 48 with one line); 17 at
  the lowest possible p were not accepted.
- Lines that are the majority of no cluster fall from 75 to 41.
- Checks: array sprs equal `subtree_prune_regraft` on three nodes; all 166
  changed splits have the best flat cut's `rho`; a separate script confirms
  the output tree and re-tests 110 candidates.
- Limits: one dataset, tree and run; candidates chosen by the tested data; no
  null calibration of the complete procedure.
- Under the third version of the draft definition of k (12 November,
  [[k-split-direction-draft]]) the same tree after the sprs gives 258
  clusters and ARI 0.7691 against 206 and 0.8359 on the unchanged tree:
  merged-only lines 30 to 13, split-only lines 16 to 42, recovered 127 in
  both.
- At the 110 spr nodes that are nodes of the unchanged tree the third
  version has p <= 0.01 at 2 before and 60 after the sprs; inside one cell
  line 0 of 38 before and 10 after.
- Room for sprs under the third version itself: the best cut along the
  node's direction passes at 32 of 1,302 refused nodes, 22 of them at nodes
  with several lines.
- Sprs at those nodes only (third version, refit): at 17 nodes with
  several lines, chosen with the labels, ARI 0.8266; at 25 nodes without
  labels 0.8257; unchanged tree 0.8359 (0.8349 to 0.8411 over five row
  orders). Merged-only lines 30 to 28 and 27, split-only 16 to 17 and 18.
- So on this tree the sprs do not improve the result under either rule.

## Evidence

`raw/inbox/choir-cell-lines/spr-refit-20261111.md` holds the setup, the spr
pass, the comparison table, the reading and the limits; protocol, scripts,
self-check and verification output and metrics are in
`raw/inbox/choir-cell-lines/spr-refit-scripts-20261111.py.txt`.
The run under the third version of the draft is in
`raw/inbox/choir-cell-lines/spr-third-version-20261112.md` with protocol,
scripts and outputs in
`raw/inbox/choir-cell-lines/spr-third-version-scripts-20261112.py.txt`.
The sprs at selected refused nodes are in
`raw/inbox/choir-cell-lines/spr-third-version-targeted-20261112.md` with
protocol, scripts and outputs in
`raw/inbox/choir-cell-lines/spr-third-version-targeted-scripts-20261112.py.txt`.

## Links

- [[k-split-direction-draft]]
- [[spr-alternatives-best-flat-cut]]
- [[spr-eigenvalue-blocks]]
- [[best-flat-cut-reference-real]]
- [[best-flat-cut-divisive]]
- [[poset-tree]]
- [[best-flat-cut-reference-k15]]
