---
title: Matched Full Benchmark Evidence
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/full-benchmark-lbfgsb-20261106.md
  - raw/inbox/full-benchmark-lbfgsb-summary-20261106.csv
  - raw/inbox/full-benchmark-lbfgsb-changed-20261106.csv
  - reports/matched_full_benchmark_20261103/README.md
  - reports/matched_full_benchmark_20261103/matched_method_summary.csv
  - reports/matched_full_benchmark_20261103/matched_case_results.csv
  - reports/matched_full_benchmark_20261103/changed_case_results.csv
  - reports/matched_full_benchmark_20261103/comparison_checks.json
  - reports/matched_full_benchmark_20261103/matched-run-manifest.json
  - reports/matched_full_benchmark_20261103/input-fingerprints.csv
  - reports/matched_full_benchmark_20261103/source-changes.patch
  - reports/matched_full_benchmark_20261103/completion.json
  - reports/matched_full_benchmark_20261103/compare_matched.py
  - benchmarks/full/run.py
tags:
  - benchmark
  - diffusion
  - validation
---

# Matched Full Benchmark Evidence

## Summary

The current integrated implementation completed the full 122-case, 11-method
benchmark against the 2026-10-19 baseline. All 1,342 recorded case/run settings
match. Eight rows improve in ARI and none decline, but the two adaptive-diffusion
variants lose one exact-cluster-count success by oversplitting an overdispersed
categorical case. This is not an isolated self-distance-effect estimate.

## Key Points

- Full execution completed with exit code 0 in 654.13 seconds, including plots
  and report generation. No case failed; skip counts and statuses are unchanged.
- Family-diffusion TBS: mean ARI 0.824874 to 0.825843; exact K 99/122 to 98/122.
  Its sole changed case, `cat_dirichlet_overdispersed_20cat`, improves ARI
  0.214239 to 0.332460 but returns seven clusters instead of the true four.
- Canonical TBS: mean ARI 0.794135 to 0.805953; exact K 91/122 to 93/122.
  Improvements occur on high-cardinality categorical and two phylogenetic cases.
- Older diffusion: mean ARI 0.530765 to 0.537087 on 109 supported cases.
  Adaptive diffusion with NNLS: 0.821914 to 0.823085 on 101 supported cases,
  with the same four-to-seven-cluster regression as family diffusion.
- Non-TBS methods are unchanged. Leiden remains at mean ARI 0.807762.
  K-means and spectral clustering receive true K, unlike the TBS methods.
- Canonical and family-diffusion TBS retain one cluster on all six suite nulls.
  The family-diffusion outlier-shell case remains at ARI 0. These limited null
  cases do not resolve the separate selected-root calibration failure.
- A 2026-11-06 rerun after the singleton-peel merge and the L-BFGS-B NNLS switch
  gives `tbs_family_diffusion` mean ARI 0.828107 (exact K 98) and `tbs` 0.810390
  (94) on the same 122 cases. Only the three singleton-peel cases and three
  `tbs_diffusion` phylogeny cases change; `tbs_diffusion_adaptive_nnls` and all
  statuses are unchanged, so every NNLS fit was certified
  (`raw/inbox/full-benchmark-lbfgsb-20261106.md`).

## Evidence

The comparison verifies one-to-one pairing, identical recorded parameters,
seeds, sample/feature counts, true K, case categories and benchmark metadata.
All 210 source hashes and the original CSV hash remain unchanged during the run.
Case/generator definitions and the lockfile match baseline commit `8e0093d3`.
Historical matrix hashes were not recorded, so byte equality to historical
arrays is not directly certified; reconstructed input hashes are now retained.

The source snapshot includes both the self-distance correction and the prior
empirical-inflation change. Consequently, gains cannot be attributed to
self-distance handling alone. Recorded stage timings omit parts of diffusion
and NNLS and are not complete per-method runtime measurements.

Compact results, the reproduction script and provenance are retained under
`reports/matched_full_benchmark_20261103/`. Full CSV snapshots are catalogued as
local-only generated outputs. The original full run and PDFs remain under
`benchmarks/results/run_20261103_matched_full_ixx1jf0m/`.

## Links

- [[diffusion-self-distance-repair-evidence]]
- [[diffusion-selected-root-diagnostic]]
- [[failing-case-subspace-separation]]
