# Calibration restoration — 2026-09-09

The historical context-weighted empirical-null calibration rule is restored on
`dev` as an uncommitted change over `6441b2e6`. The verified pydiffmap NNLS
benchmark reproduces **both ARI and cluster count on all 97 successful cases
shared with the August 1 capture**. Their mean ARI is **0.7972814937** in both
runs. Restoration recovers 92 previously unsupported results.

## Benchmark results

Method: `tbs_diffusion_adaptive_nnls`, adaptive pydiffmap, Hamming, average
linkage, fixed-topology NNLS, sibling alpha 0.01 and edge alpha 0.001. All
parameters and seeds are recorded in `run_manifest.json` and the result CSV.

| Capture | Cases | Results | Unsupported | Input skips | Mean ARI on results |
|---|---:|---:|---:|---:|---:|
| August 1 historical | 121 | 100 | 0 | 21 | 0.773363 |
| Before restoration, `6441b2e6` | 122 | 5 | 96 | 21 | 0.200000 |
| Restored, final source | 122 | 97 | 4 | 21 | 0.797281 |

These aggregate means have different denominators. **The larger restored
aggregate is not an improvement over August 1**: on the 97 jointly successful
cases, every ARI and cluster count matches. The August 1 mean includes three
additional zero-ARI results outside that common successful set.

There are 119 shared historical/current case IDs, all with identical recorded
parameters and matching sample counts, feature counts and true K. Among them:
97 return the same ARI/K, 20 retain their skip, and two now have explicit
no-support outcomes (`cat_highcard_20cat_4c`, `bar_binary_unbalanced_95_5`), where
the old capture returned K=1 and ARI=0. The August 3 no-support outcome remains.
The August 1 run lacks an exact commit/working-tree or input-hash manifest;
`260b11f9` is the nearest preceding committed implementation, not a certified
checkout for that capture.

The final 122-case run has zero execution errors, 70 exact-K results, 12
over-splits and 15 under-splits among 97 successful outcomes. Its four
no-support cases are:

- `gauss_dense_signal_highd`
- `gauss_sparse_signal_highd_noise`
- `cat_highcard_20cat_4c`
- `bar_binary_unbalanced_95_5`

The two scored one-cluster Gaussian-derived null cases return K=1 and ARI=1.
The other named null inputs are skipped; these few seeds do not measure an
error rate. `phylo_brownian_null_16taxa` has benchmark true K=16 and must not be
mistaken for a one-cluster null just because its name includes “null”.

## Restored behavior

- Supported sibling records estimate the same feature-family/context-weighted
  inflation scale, clipped below at one. The rule returns
  `chi2.sf(T / (reference_scale * c_hat), df)` with `internal_admissible`, then
  runs the existing traversal-aligned sibling BH.
- The August 11 production bypass that withheld those p-values is removed.
  Requested support thresholds now reach the existing decision function.
- The newer immutable calibration sample, dependency-group diagnostics,
  support policy and dedicated exact independent common-scale F API remain.
  Tuple dependency-group IDs work in leave-one-group scale sensitivity.
  Exact-F models cannot enter the empirical p-value API.
- Missing internal support still closes sibling gates and yields an unsupported
  benchmark result. No neutral inflation value or replacement gate is inserted.
- Mixed calibration sweep artifacts use v3; NNLS null calibration artifacts
  use v5, preventing cached pre-restoration decisions from being reused.
- README, statistical documentation, method manuscript and current wiki state
  now distinguish the active empirical rule from an unproved selected-tail law.

## Necessary calibration validation

Run independent held-out null and signal simulations of the **whole method**,
refitting diffusion, hierarchy, NNLS, projections, edge selection and empirical
calibration on every dataset. Measure local tail error, false splits, final K
and mixed-null false-discovery proportion with uncertainty based on independent
datasets. Check power simultaneously. Do not count overlapping nodes as
independent replicates or tune and evaluate on the same seeds.

Conservative nonuniform p-values can be valid: a uniformity-test rejection alone
does not justify disabling the rule. High ARI and stable mean inflation likewise
do not validate its false-split tail. The exact targets, simulation precision,
selected-null contamination checks and repair decision are in
[validation_protocol.md](validation_protocol.md), with primary references and
rechecked local evidence. That new validation study has **not** been run here.

## Verification and provenance

- Final-source benchmark: 122 unique outcomes, zero errors, 276.14 seconds.
- Independently rescored all 97 saved label files; ARI and K match the CSV.
- All 399 captured source/runner hashes match at final verification.
- All 122 new data/truth signatures are saved. The earlier September 9 run has
  signatures for 88 shared cases; all 88 match. The other 34 original signatures
  were not captured. Historical comparisons use recorded seeds, parameters and
  input metadata, not an unsupported claim of verified historical input bytes.
- Baseline before edits: 62 focused tests passed. Restoration checks: 103
  targeted tests passed. After the exact-API guard, 60 focused tests passed;
  all seven ordered stages passed, totaling 624 tests.
- Final `make check` exited zero: all 624 tests, Ruff, dependency audit,
  Vulture and wiki lint (193 pages) passed. `git diff --check` also passed.
  The manuscript built successfully through `make -C manuscript pdf`.
- No commit or push was made. Earlier selected-Gaussian prototype work remains
  outside this restoration.

The first restoration benchmark was a pilot before an additional exact-API
guard. It is retained separately under
`benchmarks/results/run_20260909_calibration_restoration/`. All reported final
results use the complete rerun under
`benchmarks/results/run_20260909_calibration_restoration_verified/`.

## Standards review

The independent review found one reference-law boundary bug: an exact-F model
could receive a chi-square p-value through the restored empirical API. An early
guard and regressions now reject that misuse. A separate projection-dimension
assertion was also retained alongside restored cluster/split tests. No remaining
material standards issue was reported in the bounded restoration diff.

## Spec review

The independent review confirmed that the empirical formula, weighted fitting,
context kernel and inflation floor match `260b11f9` and `6441b2e6^` for equivalent
records. It independently identified the same exact-F misuse, then verified
the fix and reported no residual restoration mismatch. The unsupported rule is
reviewed explicitly above; statistical validation is specified, not claimed.

## Reproduce

From the repository root, the existing runner accepts a fresh directory with
the saved manifest as `run_manifest.json`; its source hashes must match:

```sh
uv run --no-sync python -m reports.diffusion_nnls_versions_20260909.run_panel --run-dir PATH_TO_FRESH_RUN
```

Recheck and regenerate the compact comparison from the saved final run:

```sh
uv run --no-sync python -m reports.calibration_restoration_20260909.summarize
```

`summary.json`, `matched_cases.csv`, `full_benchmark_comparison.csv`,
`null_cases.csv`, `verification.json` and `run_manifest.json` contain the
machine-readable evidence. The raw benchmark directory retains labels,
annotations, timing output and source-pinned run logs.
