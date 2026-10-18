---
title: Legacy c2ef9a69 Method Package 2026-06-16
type: source
status: reviewed
updated: 2026-08-10
sources:
  - benchmarks/shared/runners/method_registry.py
  - benchmarks/shared/util/method_sets.py
tags:
  - source
  - diagnostics
  - legacy
  - benchmark
---

# Legacy c2ef9a69 Method Package 2026-06-16

## Summary

The full old `tree_break_selection` method package from commit
`c2ef9a69e0888168950bdee4a41ae8ab9996e32f` was available as an isolated nested
package until it was retired on 2026-06-25:

`tree_break_selection.legacy_methods.commit_c2ef9a69.tree_break_selection`

This was a whole-package method snapshot, not only the internal-node spectral
patch. Its absolute imports were mechanically rewritten into the nested
namespace so the old implementation can be imported in the same Python process
as the current code.
The active package, runner, registry entries, and comparison panels have now
been removed; this page is retained as historical context for the raw legacy
comparison artifacts.

## Key Points

- The package provenance is recorded in
  `tree_break_selection/legacy_methods/commit_c2ef9a69/METADATA.md`.
- The standard benchmark method id is `tbs_legacy_c2ef9a69`.
- The runner in `legacy_commit_runner.py` calls the old `PosetTree` and
  old `TreeDecomposition` directly.
- The old method intentionally ignores modern-only options such as typed
  `FeatureSpace`, phylogenetic builders, root stability guards, and spectral
  transport. Using those would turn the snapshot into a hybrid rather than the
  full old method.
- The dispatcher smoke test confirms that `tbs_legacy_c2ef9a69` runs through the
  current benchmark contract and records the legacy commit hash in result
  metadata.

## Evidence

- Git history preserves the retired `legacy_commit_runner.py`, method-registry
  entry, runner-method set, and dispatch smoke contract.
- The raw legacy comparison artifacts remain the primary behavioral evidence
  for the retired package.

## Links

- [[legacy-internal-spectral-comparison-panel-20260616]]
- [[old-vs-current-method-stack-comparison-20260615]]
- [[local-marchenko-pastur-rule]]
