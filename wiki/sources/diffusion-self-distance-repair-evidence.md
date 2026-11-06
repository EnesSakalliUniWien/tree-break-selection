---
title: Diffusion Self Distance Repair Evidence
type: source
status: reviewed
updated: 2026-11-03
sources:
  - raw/inbox/package-build-verification-20261103.md
  - raw/inbox/diffusion-self-distance-and-region-tightening-20261102.md
  - raw/inbox/diffusion-factor-bound-tightening-20261102.md
  - raw/assets/diffusion-self-distance-20261102/baseline_manifest.json
  - raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/comparison_manifest.json
  - raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/comparison_summary.csv
  - raw/assets/diffusion-self-distance-20261102/kernel_bounds_before.json
  - raw/assets/diffusion-self-distance-20261102/kernel_bounds_after.json
  - raw/assets/diffusion-region-tightening-20261102/spectral_postkernel_comparison.json
  - raw/assets/diffusion-region-tightening-20261102/final-small-controls/manifest.json
  - raw/assets/diffusion-region-tightening-20261102/final-attribution-first-root/manifest.json
  - tree_break_selection/space_separation/diffusion.py
  - tree_break_selection/space_separation/README.md
  - CHANGELOG.md
  - benchmarks/shared/runners/dispatch.py
  - tests/pipeline/dispatch/test_parameter_forwarding.py
  - benchmarks/diagnostics/calibration/diffusion_self_distance_comparison.py
  - tests/core/test_adaptive_diffusion_self_distances.py
  - tests/core/test_diffusion_self_distance_comparison.py
  - tests/statistics/66_test_diffusion_root_kernel_tightening.py
tags:
  - diffusion
  - selection
  - numerical-analysis
---

# Diffusion Self Distance Repair Evidence

## Summary

The subsequent user-confirmed work reconciled identity self distances before
epsilon, compared 55 matched datasets and tightened exact-real selection bounds.
The local numerical discrepancy is corrected; selected-null calibration remains
unresolved because substantial selection-region mass is unclassified.

## Key Points

- A minimal production change zeros only the identity diagonal before epsilon
  and kernel computation, preserving every off-diagonal entry and stored zero.
- Genuine pre-fix captures match all 50 previous Gaussian factorial histories.
  Eight histories change ordering after correction; no clade or root partition
  changes. The other five datasets also preserve their histories.
- Root statistics are identical; unselected chi tails still reject 50/50
  selected Gaussian roots. This separates the numerical correction from bias repair.
- Exact radial quadratic bounds preserve shared-motion cancellation; neighbor
  certificates require support membership instead of irrelevant internal ordering.
- Restricted diffusion-factor bounds certify tested scalar half-width 1e-5
  versus 1e-12 previously, and a positive selected interval for the corrected
  80-row observation. Interval expansion and subdivision retain unknown mass.
- All four final root decisions remain unresolved. Exact-real certificates
  continue to differ from a proof of every floating-point execution.

## Evidence

The primary notes contain the derivations, Wolfram algebra check, numerical
comparisons, reproduction settings and limits. The final versioned topology
manifest discloses an initial post-fix output replacement while preserving the
genuine baseline. Tests cover self-identity masking, stored duplicate affinity,
history versus clade comparisons, raw-output preservation, quadratic extrema,
neighbor-support ambiguity and complete-domain probability bookkeeping.

Production integration was verified on 2026-11-03 with two Sol agents. The
shared constructor already supplied the correction to direct, feature-family
and block adaptive diffusion. Expanded regression coverage exercises all three
entry points with forced positive self distances and a stored off-diagonal zero.
Feature-family dispatch now forwards geometry metadata under the existing
`adaptive_diffusion` key; its new assertion failed with `KeyError` before the
fix and passes afterward. The fallback path carries no diffusion metadata when
no diffusion geometry is built. The existing runner preserves this metadata for
successful results; its pre-existing unsupported-result metadata omission was
not changed.

Fresh integration verification passed 77 tests spanning the constructor,
wrappers, dispatch, runner contracts and affected applications. Targeted Ruff
lint and formatting pass. The full suite and distribution build were not rerun;
`flit_core` is absent from the existing environment. Selected-root calibration
remains diagnostic, with no threshold or statistical-default changes.

Subsequent package verification installed the declared build backend, built a
source distribution and its wheel, and confirmed all 109 production Python files
were present. All 114 package files match the checkout and both artifacts.
The installed wheel passed the three self-distance regressions outside the
checkout with isolated imports. All 26 installed packages were compatible;
the temporary environment was removed. The raw package-verification note records
artifact hashes, dependency-cache limitations and the library-only package scope.

## Links

- [[matched-full-benchmark-evidence]]
- [[diffusion-selected-root-diagnostic]]
- [[diffusion-root-diagnostic-evidence]]
- [[diffusion-topology-branch-rank-isolation]]
- [[diffusion-selection-covariance-attribution]]
