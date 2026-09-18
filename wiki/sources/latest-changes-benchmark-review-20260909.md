---
title: Latest Changes and Benchmark Review 20260909
type: source
status: reviewed
updated: 2026-09-09
sources:
  - reports/latest_changes_benchmark_review_20260909.md
  - reports/latest_changes_benchmark_review_20260909_methods.csv
  - reports/latest_changes_benchmark_review_20260909_evidence.json
tags:
  - review
  - benchmark
  - calibration
---

# Latest Changes and Benchmark Review 20260909

## Summary

Review of `6441b2e6`, the uncommitted selected-Gaussian prototype and the latest
August benchmark captures found four reproducible implementation issues and a
benchmark summary covering only a subset of its cited CSV.

## Key Points

- Polynomial cancellation admits an incorrect hierarchy interval: selected
  p-value 0.8561354 instead of approximately 0.8. Independent reconstructed
  distance replay disagrees with the polynomial replay inside that interval.
- Chi-tail underflow rejects a positive selected event with conditional
  p-value 0.6701192.
- Production annotation bypasses support-threshold decisions and discards
  fitted diagnostic evidence; tuple-valued dependency groups crash decisions.
- Canonical TBS returns 11 successful one-cluster outcomes and 111 unsupported
  outcomes over 122 cases. Four successful cases recover exact K.
- The variants CSV has 2,790 rows and 25 methods; its saved performance summary
  describes only 732 rows and six NNLS variants.
- Fixed-coordinate BY leads fully covered methods by mean ARI, 0.7123, but
  over-splits 74 cases and recovers exact K in 30/122. These diagnostic results
  do not establish selected-null calibration or a uniform NNLS improvement.

## Evidence

The report preserves code locations, reproductions, matched comparisons and
verification limits. The method CSV and JSON capture re-aggregated evidence
and source hashes without rewriting original benchmark artifacts. All 134
targeted tests passed; separate counterexamples exposed the four issues.
The complete project gate and full benchmarks were not rerun.

## Links

- [[empirical-null-calibration-reference-law-contract]]
- [[root-selected-region-model]]
- [[benchmark-pipeline-contract]]
- [[repository-execution-and-benchmark-map]]
