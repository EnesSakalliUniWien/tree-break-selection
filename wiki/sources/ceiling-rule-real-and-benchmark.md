---
title: The Ceiling-3 Rule on Real Cell Lines and the Continuous Benchmark
type: source
status: draft
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/ceiling-rule-real-and-benchmark-20261109.md
  - raw/inbox/choir-cell-lines/ceiling-rule-real-and-benchmark-scripts-20261109.py.txt
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
tags:
  - calibration
  - real-data
  - choir
  - benchmark
---

# The Ceiling-3 Rule on Real Cell Lines and the Continuous Benchmark

## Summary

Swapping the rule-out for the ceiling-3 rule ([[single-population-split-bound]]) on all
48,879 CHOIR cells moves Markov/Ward much closer to the 190 cell lines: 664 to 217
clusters, ARI 0.765 to 0.833, lines recovered 80 to 119. More lines are merged (11 to
33 merged only). The default tree barely changes. On the 9 continuous benchmark cases
the ceiling-3 rule gives exactly the rule-out result, because the losses there come
from subtrees without invariant directions.

## Key Points

- Independent rerun: [[three-tests-construction-audit]] reproduces the A/B cell
  assignments and all 54 benchmark rows exactly. Its new full-cell C application
  gives Markov/Ward 173 clusters, ARI 0.815201 and 111 recovered lines, versus
  B's 217, 0.832948 and 119; the audit rejects a general conservative-tail claim.
- Replays: production rule-out gives the 2026-11-08 harness result (69 / 664 clusters);
  benchmark today and rule-out rows match 36/36.
- Markov/Ward: NMI 0.901 to 0.940, ARI against published CHOIR_1 0.514 to 0.538,
  clusters under 20 cells 253 to 25; lines split only 72 to 16.
- The ceiling-3 rule keeps far fewer small splits (32-159 cells: 17.4% to 1.6%) and
  more large ones (1,000+: 34.0% to 40.8%); of 1,838 tested splits, 258 vs 121 kept.
- Default tree: 64 clusters, ARI 0.028, lines 29; its giant cluster (28,611 cells) comes
  from the topology.
- Benchmark mean ARI 0.695 for both rules (today 0.85-0.88); 60 cells in 40 features,
  45 in 30 and 40 in 20,000 have no invariant directions under either rule.
- Merges (55 of 190 lines): mostly from the tree, not the ceiling. Of 39 merged line
  pairs with invariant directions, 27 pass ceiling 3 on their true labels (rho 5.5-27.7);
  12 fail it (median smaller share 0.13), 9 of which pass the Gaussian ceiling. The very
  smallest lines merge more often (5 of the 10 smallest).
- For those 27 pairs traversal stops at the cluster's top split: 8 top splits follow the
  line boundary but reach rho 1.6-4.5 where the flat-cut count needs 4.3-7.4; 2 stop at
  the edge gate; 17 cut elsewhere (median rho 1.22) and hide deeper nodes.
- Not measured: error control on real data, long-tailed populations. No production
  change.

## Evidence

`raw/inbox/choir-cell-lines/ceiling-rule-real-and-benchmark-20261109.md` holds the setup,
replay checks, tables and reading; protocol, deviations, code and outputs are in
`raw/inbox/choir-cell-lines/ceiling-rule-real-and-benchmark-scripts-20261109.py.txt`.

## Links

- [[single-population-split-bound]]
- [[split-count-exclusion-real]]
- [[rule-out-continuous-benchmark]]
- [[three-tests-construction-audit]]
