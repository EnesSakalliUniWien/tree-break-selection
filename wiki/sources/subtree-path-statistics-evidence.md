---
title: Surviving Subtree Paths and Missing Neighborhood Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-path-statistics-20261109.md
  - raw/inbox/choir-cell-lines/subtree-path-statistics-20261109-manifest.json
  - raw/inbox/choir-cell-lines/subtree-path-statistics-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/subtree-path-statistics-20261109-wolfram.json
  - reports/subtree_path_statistics_20261109/protocol.md
  - reports/subtree_path_statistics_20261109/panel/summary.json
  - reports/subtree_path_statistics_20261109/panel/endpoint_summary.json
  - reports/subtree_path_statistics_20261109/panel/same_path_example.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_path_context.py
  - tests/statistics/81_test_subtree_path_context.py
tags:
  - subtree
  - paths
  - neighborhood
  - statistics
---

# Surviving Subtree Paths and Missing Neighborhood Evidence

## Summary

The user confirmed both ancestor–descendant and cross-branch paths and defined
the missing environment as disappearing neighboring parents, children and
siblings. The diagnostic reuses the prior native traces and relates path changes
to mean contrasts, actual edge-half contexts, branch lengths, gates and traversal.

## Key Points

- All 256 reference/alternative tree pairs were examined, using 160 saved fits
  on 16 matrices of 64 rows and six features. There are 569,954 paths between
  surviving clades, including leaves, and 1,890 with two internal endpoints.
  These are dependent comparisons on reused data, not new independent tests.
- Of 816 internal ancestor–descendant paths, 635 lose an intermediate clade;
  of 1,074 internal cross-branch paths, 993 do. On cross-branch paths the original
  LCA disappears in 416 comparisons and survives off the new path in 147.
- Missing neighborhood clades occur on 1,844 internal-endpoint paths. Parent,
  child and sibling roles overlap: a missing neighboring parent occurs in 341,
  child in 1,821 and sibling in 1,241. Surviving neighbors displaced elsewhere
  are recorded separately from absent clades.
- Full-data endpoint contrasts are identical for matched endpoints. Exact
  Wolfram checks show why standardized edge scores cannot be added: a four-leaf
  path has summed edge scores 42 and endpoint score 4.5. Native half-row context
  discrepancies are measured separately from this full-data identity.
- Endpoint traversal records change on 150 internal-endpoint paths, including
  16 reachability changes. Of 201 paths with unchanged labeled vertex sequences,
  six change a boundary/split decision; all six have changed neighborhoods.
  This is an observed association, not isolated causal attribution.
- One preserved 32-to-16-leaf edge changes its lower endpoint's children from
  14|2 to 8|8. Raw sibling T changes from 4.2517 to 78.7587 at df 6; an unavailable
  gate/boundary becomes active rejection/split. The 16 leaves are one known
  simulated group. Whole-partition AMI falls from 1 to .9366.
- Across 608 unique matched directed edges, raw T changes in 91, with unchanged
  df and Tree-BH-tested flags. All compared raw scores are finite. Path counts
  reuse these edge records and are not additional independent evidence.

## Evidence

The path CSV retains original and alternative native IDs as pointers into
unchanged source traces; matching itself uses exact leaf sets. Input/source
hashes, all counts, examples, exact algebra and the original 160 fits are
preserved. The focused context/path suite passes 31 tests and four Python files
pass Ruff/style checks. Final wiki and provenance checks are recorded with the
report. No native analysis fit, production preset or statistical reference changed.

## Links

- [[subtree-paths-and-missing-neighborhoods]]
- [[subtree-context-statistics-evidence]]
- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-clustering-production-gates]]
