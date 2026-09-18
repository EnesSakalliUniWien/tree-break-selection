---
title: Diffusion NNLS Versions and Support Review 20260909
type: source
status: reviewed
updated: 2026-09-09
sources:
  - reports/diffusion_nnls_versions_20260909/README.md
  - reports/diffusion_nnls_versions_20260909/configuration_summary.csv
  - reports/diffusion_nnls_versions_20260909/outcomes.csv
  - reports/diffusion_nnls_versions_20260909/historical_matched_cases.csv
  - reports/diffusion_nnls_versions_20260909/result_verification.json
  - reports/diffusion_nnls_versions_20260909/august_version_code_comparison.json
  - reports/diffusion_nnls_versions_20260909/support_review_results.csv
  - reports/diffusion_nnls_versions_20260909/support_review_manifest.json
  - reports/diffusion_nnls_versions_20260909/support_result_verification.json
  - reports/diffusion_nnls_versions_20260909/continuous_geometry_check.json
tags:
  - diffusion
  - nnls
  - benchmark
  - calibration
---

# Diffusion NNLS Versions and Support Review 20260909

## Summary

Internal empirical-null support and a validated selected-tail reference law
are different requirements. Neither unsupported reason is a direct measure
of clustering quality. A matched diagnostic panel demonstrates useful
diffusion NNLS clustering with zero internal support, while also exposing
alpha sensitivity and false splits on null cases.

The counts below describe the pre-restoration `6441b2e6` source. Subsequent
restoration and verification are recorded in [[calibration-rule-restoration-20260909]].

## Key Points

- Git history distinguishes the August-3 addition of a no-support
  `unsupported` result from the August-11 removal of the earlier adjusted
  chi-square sibling p-value for same-selected-hierarchy calibration. The
  older code already closed gates without internal support; August 11 also
  closes positive-dimensional focal gates when internal support exists.
- Diffusion, NNLS solver and diffusion runner implementations match the
  preceding July-30 committed snapshot. The August-1 capture lacks an exact
  commit/working-tree manifest, so that snapshot is supporting history rather
  than a certified checkout for the historical run.
- The full study completes three registered diffusion NNLS methods and ten
  tree configurations over 122 cases: 137 `ok`, 856 `unsupported`, 21 `skip`
  and 206 execution errors. All 137 returned partitions have K=1.
- Of 856 unsupported outcomes, 813 have internal support but an unvalidated
  reference law; 43 have no internal support. The execution errors comprise
  190 centroid/median merge-height inversions and 16 MAD zero-distance failures.
- Among 119 shared August-1/current pydiffmap cases with identical recorded
  parameters and matching input metadata, 94 change from `ok` to `unsupported`;
  five remain `ok` with unchanged ARI/K and 20 remain skipped. All 122 fresh
  pydiffmap statuses agree with the August-11 capture.
- The gate review completes 168 outcomes: 141 `ok`, 11 `unsupported` and
  16 input-compatibility skips. Fourteen skips are downstream cells propagated
  from the pydiffmap duplicate-geometry contract; those gates were not run.
- With zero internal support on dense four-cluster signal, pydiffmap NNLS
  with `fixed_coordinate_bh` at alpha 0.01 recovers K=4 and ARI=1. Removing
  only the outer unsupported status returns K=1 because sibling gates remain
  closed. Both graphtools variants also return K=1 on that dense case with
  either fixed-coordinate gate at all three tested alpha values.
- The same diagnostic gates split the small Gaussian-derived null for all
  three methods and the large null for both graphtools methods at alpha 0.01.
  Good signal ARI does not validate their selected-tail p-values.
- The historical 2026-08-01 pydiffmap NNLS preset has mean ARI 0.773363 on
  100 successful cases. Against TBS on those same cases, it wins 43, ties 47
  and loses 10, with mean ARI gain 0.121122. This is a combined diffusion/NNLS
  preset comparison, not an isolated NNLS effect or a historical code rerun.
- The review recommends retaining diagnostic cluster assignments alongside
  the calibration reason. It does not change the production rule or promote
  either fixed-coordinate gate to calibrated status.

## Evidence

The diagnostic study freezes diffusion distances per case/method and verifies
identical tree and branch-length hashes across the eight gate/alpha variants.
All 137 full-suite and 141 diagnostic saved label files were independently
rescored and match their ARI and cluster-count records. All 1,220 requested
outcomes are unique and all 398 captured source hashes remain unchanged.
Four null seeds are diagnostic observations, not an error-rate study.
Original input seeds and parameters are preserved.

The graphtools Hamming presets accept native continuous input: on one
60-sample continuous signal case, all 1,770 pair distances equal 1. Such rows
require separate interpretation when comparing coverage and quality.

## Links

- [[empirical-null-calibration-reference-law-contract]]
- [[adaptive-diffusion-nnls-method-library-audit]]
- [[latest-changes-benchmark-review-20260909]]
- [[benchmark-pipeline-contract]]
