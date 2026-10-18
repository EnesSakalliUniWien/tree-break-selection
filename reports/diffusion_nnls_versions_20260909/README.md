# Diffusion NNLS versions and unsupported-rule review

Study date: 2026-09-09. All 1,220 requested full-suite outcomes and 168 focused
diagnostic outcomes are recorded. All three final benchmark workers exited 0;
this means the study completed, including recorded method execution errors.

The historical diffusion NNLS quality gain is supported by the saved results.
Current production behavior primarily withholds assignments under the newer
calibration contract: every returned full-suite partition contains one
cluster. Zero internal support can still coexist with useful diagnostic
clustering, but the same diagnostic gates also produce null splits.

## Full benchmark results

Each configuration below covers all 122 current cases. ARI is conditional on
returned labels and is not an unconditional quality ranking.

| ID | Diffusion NNLS configuration | OK | Unsupported | Skip | Execution error | Mean ARI on OK |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| C01 | Pydiffmap, average | 5 | 96 | 21 | 0 | 0.2000 |
| C02 | Graphtools fixed K, average | 15 | 107 | 0 | 0 | 0.2667 |
| C03 | Graphtools adaptive K, average | 15 | 107 | 0 | 0 | 0.2667 |
| C04 | Graphtools adaptive K, complete | 14 | 108 | 0 | 0 | 0.2857 |
| C05 | Graphtools adaptive K, weighted | 16 | 106 | 0 | 0 | 0.2500 |
| C06 | Graphtools adaptive K, single | 28 | 94 | 0 | 0 | 0.2143 |
| C07 | Graphtools adaptive K, centroid | 14 | 9 | 0 | 99 | 0.2857 |
| C08 | Graphtools adaptive K, median | 14 | 17 | 0 | 91 | 0.2857 |
| C09 | Graphtools adaptive K, ward | 14 | 108 | 0 | 0 | 0.2857 |
| C10 | Graphtools adaptive K, neighbor joining / MAD | 2 | 104 | 0 | 16 | 0.5000 |

Totals are **137 OK, 856 unsupported, 21 skipped and 206 execution errors**.
Every one of the 137 OK rows returns K=1. Of the unsupported outcomes, 813
have internal support but an unvalidated reference law; 43 have no internal
support. The constructor errors are 99 centroid and 91 median merge-height
inversions, plus 16 MAD zero-distance failures.

The apparently higher coverage of graphtools average is partly native
continuous input: 11 of its 15 OK cases use the Hamming-on-continuous preset.
Pydiffmap skips all 13 continuous cases. `representation_summary.csv` preserves
this distinction; mean ARI on a different returned subset does not identify
the best diffusion geometry. The per-configuration and per-case tables also
retain all execution errors in the denominator.

## Scope and design

This study tests the three registered diffusion NNLS methods on all 122 cases
of the current `full` suite at `dev` commit `6441b2e6`, with the existing
uncommitted working tree captured in `run_manifest.json`. It uses the existing
benchmark generators, dispatch, metrics, NNLS implementation and method
presets. No production source, method default, dependency or historical
benchmark output is changed.

- `tbs_diffusion_adaptive_nnls`: pydiffmap variable-bandwidth diffusion,
  average linkage.
- `tbs_diffusion_graphtools_nnls`: fixed-K graphtools diffusion, average linkage.
- `tbs_diffusion_graphtools_adaptive_nnls`: graphtools `fragmentation_guard`
  adaptive K, with average, complete, weighted, single, centroid, median and
  ward linkage, plus neighbor joining with MAD rooting.

There are ten configurations and 1,220 requested case/configuration outcomes.
The original case seeds and one repeat are retained. Sibling alpha is 0.01;
edge alpha is 0.001. Diffusion time is 3 and the component count is 30. The
NNLS target is `squared_standardized_euclidean` on the original distributional
features, with pair sample size 50,000, random state 0, tolerance 1e-5 and
maximum iterations 1,000. The variance policy is `normalized_branch_length`.
Complete parameters and source hashes are in the manifest.

### Execution recovery

The established `python -m benchmarks.full.run` invocation stopped on the
second case: neighbor joining with MAD rooting raised
`ValueError: MAD rooting requires positive distances between leaf pairs.`
It exited 1 after 173.69 seconds, leaving only the first case's ten rows.
`initial_execution_status.json` preserves that failed attempt.

`run_panel.py` resumes those ten recorded cells using the existing
`run_single_method_once` API. It appends each returned row immediately and
records constructor/solver exceptions separately, with their tracebacks.
An exception is an execution error, not an unsupported calibration outcome,
an ARI of zero, or an implicitly successful fallback. The standard returned-row
summary excludes these exceptions; `configuration_summary.csv` and
`outcomes.csv` include all requested outcomes and are authoritative here.

The serial continuation was intentionally interrupted with exit 130 after
preserving its completed outputs in `serial_stage/`. `parallel_finish.py`
partitions the remaining work by original case index across three isolated
workers, retaining all configurations and original indices within each case.
It appends only newly completed rows/errors to the original files and runs the
same completeness checks after merging. `parallel_manifest.json` records the
commands and hashes; this changes execution scheduling only.

Successful cells have labels and decomposition audits in separate directories
to avoid overwriting the eight adaptive-K configurations. Input and source
hashes, exact commands, logs and per-cell wall times remain under
`benchmarks/results/run_20260909_diffusion_nnls_versions_full/`. The commands
set BLAS/OpenMP thread limits to one; pydiffmap's existing nearest-neighbor
`n_jobs=-1` setting is retained. Concurrent execution and the
interrupted first attempt mean the recorded times are execution diagnostics,
not a controlled speed ranking. The original JSON manifest serializes the
adaptive-K tuple as a list; dispatch converts both to the same integer tuple.

## What the unsupported rule does

The rule concerns internal empirical-null calibration and the validity of the
selected-tail reference law. It does not measure clustering accuracy.

There are two distinct reasons:

1. `empirical_null_no_internal_support`: focal sibling tests exist, but no
   admissible internal empirical-null calibration records remain. Strong
   signal can produce this condition when the child-edge tests reject
   throughout the tree; lack of null calibration is not proof of absent
   cluster structure.
2. `empirical_null_unvalidated_reference_law`: internal calibration records
   exist, but their same-selected-hierarchy conditional reference law is
   unvalidated. Increasing the support count alone does not resolve this.

The implementation has two effects. First, `pipeline.py` marks unresolved
positive-dimensional sibling tests as skipped, invalid and closed. Second,
`tbs_support.py` detects any such focal record in the tree, and
`tbs_runner.py` returns `labels=None` before traversal. This applies to the
whole result, even if the unresolved record would not be reached in traversal.
Removing only this outer status check leaves the closed sibling gates intact.

The latest commit changes this calibration contract and its reporting; its
diff does not change the diffusion construction or NNLS solver. The earlier
implementation review also found that the unresolved-law branch bypasses
support-threshold diagnostics and discards fitted diagnostic information. See
`reports/latest_changes_benchmark_review_20260909.md` for the independent
reproductions and other outstanding prototype issues.

## Changes since the August 1 capture

The recorded diffusion, linkage, NNLS and alpha parameters match on all 119
shared cases. The diffusion implementation, NNLS solver and diffusion runner
files are also byte-identical between the preceding committed snapshot
`260b11f9` (July 30) and the current working tree. This is checked in
`august_version_code_comparison.json`. The historical run has no code-commit
or working-tree manifest, so that preceding commit must not be presented as
a certified exact snapshot of the August 1 benchmark execution.

The relevant behavior changes are:

- **August 3, `20362435`:** the benchmark runner gains a whole-result
  `unsupported` outcome for absent internal support and returns no labels.
  The earlier sibling annotation already closed gates without support;
  previously traversal could still return the one-cluster boundary and an
  ordinary benchmark score. This status change did not itself introduce the
  no-support closed-gate behavior.
- **August 11, `6441b2e6`:** previously the empirical-null model estimated a
  context-weighted inflation factor, computed
  `chi2.sf(stat / (reference_scale * inflation_factor), df)` and supplied that
  p-value to sibling BH. It could therefore open sibling gates when internal
  support existed. The current same-selected-hierarchy model instead reports
  `undefined_unvalidated_reference_law` with `p_value=None`; positive-dimensional
  focal gates are closed, and the runner also returns `unsupported` for this
  condition even with internal support present.
- That August 11 commit also adds dependency-group accounting and a separate
  exact-F decision for independent, disjoint, fixed-unit-weight, common-scale
  calibration. This restricted API is not a replacement tail for the current
  selected-hierarchy production path. Other intervening changes include
  support-role metadata, spectral row-count diagnostics, retired optional
  profiles and benchmark case/schema changes.

The reason for the statistical restriction is that the estimated scale alone
does not account for finite calibration uncertainty, shared observations and
same-data hierarchy/edge selection. Its practical effect is broader than
withholding a significance claim: it suppresses assignments which previously
had high empirical clustering accuracy. For example, `gauss_clear_small`
previously returned three clusters with ARI 1; the current pydiffmap method
reports unsupported despite 27 internal support records.

## Matched diagnostic gate comparison

`review_support.py` evaluates seven predeclared cases: clear signal, dense
high-dimensional signal, extreme overlap, and four global-null cases. It
uses each of the three diffusion NNLS methods with average linkage and eight
variants: production at alpha 0.01; only the outer status rule bypassed at
alpha 0.01; and the existing `fixed_coordinate_bh` and `fixed_coordinate_by`
gates at alpha 0.005, 0.01 and 0.05. This is 168 diagnostic outcomes.

Each method/case pair constructs diffusion geometry once. All variants reuse
those exact distances and refit the same NNLS tree. Assertions check identical
topology and branch-length hashes across the eight variants. Truth labels
are used only for scoring, never gate selection. The status bypass is a
process-local patch in the study; fixed-coordinate gates are explicit changes
of sibling test, not merely disabling the unsupported rule.

These fixed-coordinate results evaluate clustering utility on the selected
hierarchy. They do not establish calibrated selected-tail p-values, FDR or
Type I error control. Four null seeds and three alpha values are a focused
diagnostic panel, not a repeated null-calibration study or a full-suite ranking
of replacement sibling gates.

The diagnostic panel completed with 141 `ok`, 11 `unsupported` and 16 `skip`
outcomes. The pydiffmap duplicate-geometry contract skips both binary-null
cases before constructing a tree. Fourteen downstream gate cells consequently
carry `execution_kind=input_contract_skip`; they were not executed and have
no fabricated clustering scores. The initial study script stopped on its
missing skip handler after 120 completed cells; that attempt is preserved in
`support_review/initial_attempt/`, and the continuation reused those cells.

At alpha 0.01:

| Case | Method | Production | Status rule bypassed K | BH K / ARI | BY K / ARI |
| --- | --- | --- | ---: | --- | --- |
| Clear signal, true K=3 | All three | Unsupported; 27 internal support records | 1 | 3 / 1.0 | 3 / 1.0 |
| Dense signal, true K=4 | Pydiffmap | Unsupported; zero internal support | 1 | 4 / 1.0 | 3 / 0.697674 |
| Dense signal, true K=4 | Both graphtools variants | Unsupported; zero internal support | 1 | 1 / 0.0 | 1 / 0.0 |
| Extreme overlap, true K=4 | All three | OK, K=1 | 1 | 1 / 0.0 | 1 / 0.0 |
| Small Gaussian-derived null, true K=1 | All three | Unsupported | 1 | 2 / 0.0 | 2 / 0.0 |
| Large Gaussian-derived null, true K=1 | Pydiffmap | OK, K=1 | 1 | 1 / 1.0 | 1 / 1.0 |
| Large Gaussian-derived null, true K=1 | Both graphtools variants | Unsupported | 1 | 2 / 0.0 | 2 / 0.0 |
| Two binary nulls, true K=1 | Both graphtools variants | OK, K=1 | 1 | 1 / 1.0 | 1 / 1.0 |

The Gaussian-derived cases above use the benchmark's binary representation;
they are not the native continuous cases discussed below. Pydiffmap BH on
dense signal changes from K=3 / ARI 0.697674 at alpha 0.005 to K=4 / ARI 1.0
at 0.01 and 0.05. BY changes from K=1 / ARI 0 at 0.005 to K=3 / ARI 0.697674
at 0.01 and 0.05. Both graphtools methods stay at K=1 for every tested dense
signal gate/alpha. All nineteen compatible case/method pairs return K=1 when
only the status rule is bypassed; this is the unchanged closed-gate behavior.

At alpha 0.01, both diagnostic gates split one of two scored pydiffmap null
cases and two of four scored null cases for each graphtools variant. These
are individual false-split outcomes, not error-rate estimates. BY at alpha
0.005 avoids the large-null split for graphtools but still splits the small
null. This prevents recommending either gate as a generally validated
replacement solely from signal-case ARI.

## Historical evidence and comparison limits

The 2026-08-01 `full_post_crossfit` capture has 100 successful and 21 skipped
`tbs_diffusion_adaptive_nnls` rows out of 121, with mean ARI 0.773363 on
successful rows. Against TBS on those same 100 cases, diffusion NNLS wins
43, ties 47 and loses 10; mean ARI gain is 0.121122. Recorded case seeds,
input dimensions, truth K, source/representation and both alpha values match.
This compares the combined preset; it does not isolate diffusion from NNLS.
The historical pydiffmap preset returns one cluster on both scored global-null
cases; four other one-cluster cases are skipped. Thus the historical empirical
gate is not interchangeable with the newly tested fixed-coordinate diagnostics:
the latter split the small null. Neither the old two null outcomes nor its
good aggregate ARI establishes selected-tail calibration.

The old graphtools captures are also preserved in the configuration summary:
fixed K from 2026-06-30 and eight adaptive-K tree configurations from
2026-07-06. Their NNLS iteration limit was 300 rather than the current 1,000;
older schema/default differences and absent input hashes prevent treating
them as controlled reruns. The case suite also changed from 121 to 122 cases.
Exact case IDs are matched without silently renaming old extreme-noise cases
to current dense-signal cases. `historical_matched_cases.csv` records parameter
differences and status transitions individually.

For pydiffmap, all 119 shared August-1/current case IDs have identical recorded
parameters, seeds, sample/feature counts, truth K and source/representation.
The transitions are 94 `ok` to `unsupported`, five `ok` to `ok`, and 20 `skip`
to `skip`. All five jointly scored cases have unchanged ARI and K. The fresh
run agrees with the August-11 capture on all 122 statuses and recorded
parameters, with unchanged ARI on its five scored cases. This supports a
loss of returned assignments, rather than evidence that the diffusion
geometry itself became weaker.

For completeness, the older graphtools mean ARI / successful-case counts are:
fixed K average 0.741645 / 117; adaptive K average 0.743267 / 117, complete
0.728065 / 115, weighted 0.749952 / 118, single 0.716481 / 119, centroid
0.748246 / 118, median 0.729795 / 119, ward 0.751030 / 113, and neighbor joining
0.799965 / 106. Each historical configuration has 121 total rows. These are
descriptive historical results under the different solver/schema conditions
above; averaging the eight tree configurations together would not represent
a single selected method.

An unsupported row has no returned partition to score. Mean ARI conditional
on `ok` must always be read with coverage; it cannot be compared as an
unqualified full-suite quality score when most assignments are withheld.

## Additional benchmark limitations

- The graphtools presets use Hamming distance on native continuous inputs,
  whereas pydiffmap explicitly skips incompatible continuous geometry. On
  `gauss_clear_medium_continuous` (seed 101, 60 samples, 40 features), all 1,770
  Hamming pair distances equal 1. The fresh input hash and reproduction are
  in `continuous_geometry_check.json`. Such rows do not establish useful
  continuous diffusion geometry; the existing presets were tested unchanged.
- Centroid and median linkage can have decreasing merge heights. The current
  tree constructor rejects these before NNLS. MAD rooting rejects zero
  leaf-pair distances. These execution errors are separate from calibration
  support and must not disappear from the denominator.
- The uncommitted selected-Gaussian prototype is not connected to these
  production method presets and is not evaluated by their benchmark ARI.
- This review does not infer valid p-values from good ARI, use known K to cut
  a tree, tune on benchmark truth, or substitute raw scores for calibrated
  significance claims.

## Recommendation

The unsupported rule is appropriate for withholding an unvalidated calibrated
p-value. Applying that status to the whole clustering result makes it
unsuitable as the sole verdict on clustering utility. Keep the calibration
reason visible while exposing cluster assignments from an explicitly
diagnostic method for quality evaluation. Benchmark its coverage, ARI, cluster
counts and null splits independently of the calibration claim.

Do not simply delete the outer status check: the matched experiment shows that
this yields the same one-cluster closed-gate result. Recovering the older
useful behavior requires an explicit sibling decision rule and a matched
benchmark of that rule. The fixed-coordinate examples prove that useful
clustering can coexist with zero internal support, but also demonstrate null
splits and alpha sensitivity. They do not justify silently promoting BH or BY
to a calibrated default or asserting that every diffusion variant improved.

The earlier empirical method should remain a distinct historical comparator;
its lower false-split count on the two scored null cases must not be attributed
to the fixed-coordinate gate. A repeated null study with regeneration of the
selected hierarchy would be needed to assess inference claims. No production
rule has been changed by this review.

## Reproduction and verification

From the repository root, using the existing environment:

```bash
uv run --no-sync python -u -m reports.diffusion_nnls_versions_20260909.run_panel --run-dir benchmarks/results/run_20260909_diffusion_nnls_versions_full
uv run --no-sync python -u -m reports.diffusion_nnls_versions_20260909.review_support --run-dir benchmarks/results/run_20260909_diffusion_nnls_versions_full
uv run --no-sync python -m reports.diffusion_nnls_versions_20260909.summarize --run-dir benchmarks/results/run_20260909_diffusion_nnls_versions_full
```

The first command resumes the recorded run. The second rewrites only this
study's diagnostic outputs. Use the exact initial invocation in the manifest
to reproduce the established full runner's failure from an empty run directory.
The focused benchmark readiness suite passed 24 tests, and the unsupported
outcome/sibling annotation suites passed 11 tests. The latter output is saved
in `support_tests.log`. Final checks verified 1,220 unique requested outcomes,
all 122 case IDs in every configuration, and unchanged hashes for 398 source
files. All 137 full-suite label files and all 141 diagnostic label files were
independently rescored and match the stored ARI and K. Repeated input hashes
match; 121 cases have input signatures, while the first case's ten unsupported
rows predate that continuation instrumentation and retain source/seed records.
Ruff, `git diff --check` and wiki lint passed. No commit or push was made; the
complete `make check` gate was not run for this review.
