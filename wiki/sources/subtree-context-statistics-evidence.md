---
title: Subtree Context and Eigenvalue Statistics Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-context-statistics-20261109.md
  - raw/inbox/choir-cell-lines/subtree-context-statistics-20261109-manifest.json
  - raw/inbox/choir-cell-lines/subtree-context-statistics-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/subtree-context-statistics-20261109-wolfram.json
  - reports/subtree_context_statistics_20261109/protocol.md
  - reports/subtree_context_statistics_20261109/native_summary.json
  - reports/subtree_context_statistics_20261109/native_summary.md
  - reports/subtree_context_statistics_20261109/relationships/README.md
  - reports/subtree_context_statistics_20261109/math/README.md
  - tests/statistics/79_test_subtree_topology_context.py
  - tests/statistics/80_test_subtree_context_native_replay.py
tags:
  - subtree
  - eigenspaces
  - calibration
  - validation
---

# Subtree Context and Eigenvalue Statistics Evidence

## Summary

The confirmed follow-up traces four overlapping node-presence cases across
subspace-generated trees and relates them to native edge/sibling statistics,
branch factors, eigenvalue magnitude and traversal. Diagnostics and reports
were added; no production algorithm or default changed.

## Key Points

- Match clades by exact labeled leaf sets. Record the original parent's and
  children's presence separately from direct adjacency. Missing clades have no
  corresponding native comparison; root incoming comparisons are not applicable.
- The fixed panel reuses 16 matrices, with two native references and eight
  covariance alternatives each. All 160 fits complete, yielding 16,128 dependent
  node comparisons. Counts for cases 1/2/3/4 are 15,316/5,271/401/102. Case 2
  has 399 present focal nodes, including 200 roots. Counts overlap.
- All 19 nonroot same-parent/same-children comparisons change raw sibling T
  and branch factor g, with unchanged df. Across 296 same-child comparisons,
  median absolute T*g difference is 2.53e-11. This is approximate branch scaling
  under the native ridge, not a universal identity. Changed children compare
  different contrasts; T*g changes in 257/260 such nonroot comparisons.
- Among 812 present alternatives, 479 have finite collected raw p-values but
  unavailable active BH decisions. Traversal visits 523 and leaves 289 unvisited.
  All 32 native references reproduce prior root splits, decisions and partitions.
- The 5,120-row saved-data reanalysis finds that PC2 with median variance share
  .2885 recovers the omitted crossing factor in 20/20 datasets, while the leading
  span with median share .9300 omits it. Lower common-target residual and larger
  retained variance do not identify that factor.
- Exact Wolfram and live-kernel controls distinguish covariance whitening,
  diffusion eigenvalue distance weights and native PCA rank/tie decisions.
  The edge-half basis differs from the all-leaf sibling context in the inspected
  fixture. Neither the basis split nor nominal chi-square tails proves selected
  calibration.
- Review found three local rebuilt bases cutting numerical tie blocks. A saved
  data sidecar recomputes all 8,064 comparisons with whole positive blocks,
  correcting those three ranks and leaving aggregate span classifications intact.
  Original traces and the exact executed runner are preserved; no native fit
  was rerun for this correction.

## Evidence

The report retains the fixed protocol, all attempts, data/source hashes, exact
structural witness, node masks, actual spectral contexts, branch fits, collected
raw records, annotations and traversal. The final focused suite passes 69 tests;
style, archive and wiki checks are recorded in the report verification artifact.
The source manifest distinguishes the executed runner from the corrected span
helper. Common-mean null roots still reject on all four reused null matrices
under each method; the panel supplies no new production calibration guarantee.

## Links

- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-spectral-alignment-evidence]]
- [[subtree-clustering-research-evidence]]
- [[gaussian-root-calibration-and-cut-geometry]]
