# Naming consistency audit and fixes

What this repo does: Tree-Break Selection builds a hierarchy, tests child-parent
and sibling contrasts, and selects a cluster partition. This audit assumed
local research scripts and batch benchmark runs, including local checks of the
AWS shard merger.

The user requested a naming audit and implementation of verified
corrections. The audit covered active package APIs, benchmark and
application callers, CLI options, output labels, tests, and relevant manuscript
and wiki terminology. The initial working tree was clean. Findings below are
resolved in the working tree; no commit or push was made.

## Should fix — applied

1. **A method named divergences only populated distributions.**
   (`tree_break_selection/tree/poset_tree.py`, `PosetTree`)
   - **What this is:** This method stores node distributions, leaf counts,
     leaf flags, and the initial annotation table.
   - **Problem:** `populate_node_divergences` implied that divergence tests
     had already run, although it performs none.
   - **Fix:** Renamed it `populate_node_distributions` throughout active
     callers, tests, quick start, ASV timing names, and package documentation.
   - **If we skip it:** A caller can mistake distribution preparation for
     completed statistical annotation.

2. **An FDR API also implemented family-wise error adjustments.**
   (`tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py`)
   - **What this is:** Four coordinate p-value adjustments share an API.
     FDR means false discovery rate; FWER means family-wise error rate.
   - **Problem:** Its `fdr` names included Holm and Bonferroni, while the
     generic interface also supports the BH and BY FDR procedures.
   - **Fix:** Renamed the type and helpers to `CoordinateAdjustmentMethod`,
     `coordinate_adjusted_minimum`, and `fixed_coordinate_adjusted_p_value`.
     Removed the redundant `fixed_coordinate_bh_p_value` forwarding wrapper;
     callers now pass `method="bh"` explicitly.
   - **If we skip it:** The interface misidentifies the error criterion
     associated with two supported methods.

3. **Binary weighting accepted arbitrary raw values.**
   (`tree_break_selection/space_separation/adaptive_cosine.py`, `weight_feature_matrix`)
   - **What this is:** Spectral separation accepts raw or TF-IDF-weighted
     feature values.
   - **Problem:** `weighting="binary"` returned the input unchanged, even
     for signed continuous features; it did not binarize them.
   - **Fix:** Renamed this option to `"raw"` in functions, CLI choices,
     callers, and tests. Relevant saved-input readers reject retired spectral
     weighting labels explicitly, while preserving the whole-matrix diffusion
     producer's established `weighting="whole"` input.
   - **If we skip it:** Users can infer a binary transformation that never
     happened, and renamed reader filters could silently discard old rows.

4. **The NNLS distance name omitted its average.**
   (`tree_break_selection/tree/optimized_branch_lengths.py`)
   - **What this is:** Non-negative least squares (NNLS) fits branch lengths
     to pairwise distances on a fixed tree topology.
   - **Problem:** `squared_standardized_euclidean` divided its coordinate
     sum by the total input feature count, unlike the unaveraged distance.
   - **Fix:** Renamed the target to `mean_squared_standardized_euclidean`,
     its constant to `BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN`,
     and its private helper to `_mean_squared_standardized_targets`.
   - **If we skip it:** The target name hides a scale factor that matters
     when interpreting fitted branch lengths.

5. **Tree captions described fitted length as calibrated time.**
   (`applications/scrna/plots/selected_nnls_report.py`;
   `tree_break_selection/tree/branch_lengths.py`)
   - **What this is:** Rectangular coordinates sum fitted branch lengths;
     radial coordinates further divide by their maximum.
   - **Problem:** Captions said cumulative NNLS branch time and omitted
     radial normalization. Linkage documentation also implied stochastic time.
   - **Fix:** Used cumulative fitted branch length and normalized cumulative
     fitted branch length, preserving the topological-depth fallback.
     Linkage documentation now describes normalized merge-height increments.
   - **If we skip it:** Readers can assign time units or a calibration
     guarantee that the plotted quantities do not establish.

6. **A simulation summary obscured its event and overstated control.**
   (`benchmarks/validation/statistics/traversal_sibling_fdr_null.py`;
   `benchmarks/cloud/aws_traversal_sibling_fdr_null.py`)
   - **What this is:** The diagnostic reports mean false-discovery proportion
     and the frequency of at least one false rejection across successful runs.
   - **Problem:** `false_rejection_rate` did not name that event, and
     `algorithmic_fdr_control` overstated a finite simulation result.
   - **Fix:** Used `any_false_rejection_rate` and `empirical_control_observed`;
     documented the successful-run denominator. Advanced the output schema
     to `traversal_sibling_fdr_null/v2`; the merger rejects old or mixed schemas.
   - **If we skip it:** Readers can confuse different error summaries or
     treat a simulation observation as a proof of FDR control.

## Nice to have — applied

7. **A metric helper's name described only three of its outputs.**
   (`benchmarks/shared/metrics.py`)
   - **What this is:** The helper returns clustering metrics, including
     external agreement, fragmentation, internal geometry, and outlier metrics.
   - **Problem:** `calculate_ari_nmi_purity_metrics` described only three
     measures despite returning the wider `ClusteringMetrics` record.
   - **Fix:** Renamed it `calculate_clustering_metrics` and updated all
     active callers and tests.
   - **If we skip it:** Callers looking for the other supported metrics
     cannot identify the existing helper from its name.

8. **Reference-result descriptions did not identify the returned values.**
   (`hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py`;
   `hierarchy_analysis/statistics/multiple_testing/tree_bh/` under `tree_break_selection/`)
   - **What this is:** One API returns `ProjectedQuadraticReference`; the
     Tree-BH API returns local adjusted p-values and hierarchical decisions.
   - **Problem:** The first documented a three-item tuple. The second did
     not explicitly distinguish within-family adjustments from ancestor-scaled
     rejection levels.
   - **Fix:** Corrected the return type and explained the Tree-BH values,
     including placeholder values for untested edges.
   - **If we skip it:** A caller can unpack the wrong type or compare local
     adjusted values with base alpha and fail to reproduce Tree-BH decisions.

## Migration

The old Python names and the old `binary` weighting and
`squared_standardized_euclidean` option values have no compatibility aliases.
External callers must adopt the names above. The coordinate test file is now
`tests/statistics/44_test_fixed_coordinate_adjustment_candidates.py`.

Existing raw captures, saved reports, benchmark outputs, and wiki chronology
retain their original names and schemas. To process historical weighting
outputs with current strict readers, regenerate the input under current
option names or use the matching historical code. Two report helpers,
`reports/eigenvalue_block_comparison_20260926/summarize_results.py` and
`reports/eigenspace_reassessment_20260928/spectral_blocks/compare_before_after.py`,
keep `binary`: each first requires source hashes from a recorded working tree
that matches no commit, so neither can reach a weighting call on later source.
Their code comments and READMEs state this. Old traversal diagnostic
shards cannot be merged with v2 shards. Existing ASV timing history retains
the former method identifier.

## Evidence and verification

- Readiness: macOS arm64, zsh 5.9, existing Python 3.11.12 environment,
  uv 0.9.5, SciPy 1.16.3, pytest 8.4.2, Ruff 0.15.14, deptry 0.25.1,
  and Vulture 2.16. Manifests, lockfile, configured checks, test fixtures,
  and affected caller paths were inspected before edits.
- Baselines and focused suites passed. Exact comparisons preserved all four
  coordinate-adjustment results for checked vector sizes, all 171 NNLS
  pair targets in a deterministic 19-by-7 case including a constant feature,
  and the unweighted feature matrix.
- A local two-shard diagnostic run exercised current production and merging,
  schema rejection, the new output column, and both CLI parsers without AWS.
- `UV_OFFLINE=true make check` passed: Ruff lint and formatting, deptry,
  Vulture, wiki lint (213 pages), generated-file checks (305 classified files),
  and all seven test stages (895 passed). One existing scikit-learn MDS
  `FutureWarning` concerns a future `n_init` default.
- `UV_OFFLINE=true uv run --no-sync python quick_start.py` completed:
  80 samples, six returned clusters, ARI 0.9530. This is a smoke run,
  not new evidence of calibration or improved clustering.
- Integration review preserved the producer's `whole` weighting label for
  whole-matrix diffusion. Three reader tests cover that valid input and
  explicit rejection of retired spectral labels. Final reference searches
  found no retired Python identifiers in active code; `git diff --check`
  passed, and the Git index remains unchanged.

The multiple-testing distinction follows the
[official R documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
The distance definition follows the
[official SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.seuclidean.html).
Captured paraphrases and access limits are in
`raw/inbox/naming-audit-reference-notes-20261017.md`.

Verdict: All eight confirmed naming findings are corrected, and the full
project gate passes. External callers and historical inputs need the documented
name/schema migration.

Lean: one forwarding wrapper and its export removed; no dependencies removed.

Not checked: every historical report or notebook, external consumers, full
research benchmark campaigns, AWS execution, and wheel packaging. The
configured wheel backend `flit_core` is absent from the installed environment;
no dependency installation was attempted. This audit does not validate
selected-tail calibration or prove whole-pipeline error control.
