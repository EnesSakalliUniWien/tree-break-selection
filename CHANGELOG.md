# Changelog

All notable changes to this project are documented in this file.

## [Unreleased] - 2026-10-14

### Changed
- Pass-through traversal counts a deeper split only when every node on the path would split or pass through. Runs with pass-through enabled (the default) can return fewer clusters than before.
- `build_tree` builds linkage trees only and no longer takes `builder` or `rooting`; `TreeBuildResult` drops `builder`, `rooting` and `phylogenetic_rooting`. Benchmark methods no longer accept `tree_builder` or `tree_rooting`.
- `run_clustering_result` rejects unknown parameter keys for every registered method; `benchmark_*` bookkeeping keys remain allowed.
- Decomposition results report `sibling_alpha_effective` (sibling alpha divided by `sibling_gate_alpha_penalty`) next to `sibling_alpha`.
- Core Benjamini-Hochberg and Tree-BH use `scipy.stats.false_discovery_control`. The scipy floor is 1.11; `statsmodels` is no longer a core dependency and remains only in the `all` extra.
- Full-benchmark resume tracks completed cases by `case_id` instead of suite position.
- Diffusion coordinates accept sparse kernels; connected graphs above 2000 samples use `eigsh` for the leading modes.
- The graphtools adaptive-K tree grid and tree consensus cover the seven linkage methods only, named by `SUPPORTED_LINKAGE_METHODS`; the `TREE_CONSENSUS_STRATEGIES` and `TREE_PRIORITY` aliases were removed.
- `decompose_invariant_equivariant_space` / `InvariantEquivariantSpace` (`invariant_equivariant.py`) were renamed `decompose_leading_axes_space(following_dim=...)` / `LeadingAxesSpace` (`leading_axes.py`); the scRNA space-decomposition script now uses `leading_*` / `following_*` column names and `--following-dim`.
- `run_single_method_once` returns `(row, computed)`; `TreeBuildResult.linkage_matrix` and the distance diagnostics are always present.
- The wiki lint requires exactly one H1 outside code fences.
- `make check` runs the new `generated-files-check`. The 36 large report outputs are archived outside the repository and git-ignored by exact path; the four run logs stay local-only.
- `plotly` moved from the core dependencies to the `viz` extra.
- Benchmark diagnostics read boolean cells through `bool_series`/`bool_value` only (`root_split`, `false_split`, guard flags, `required_for_production`, edge flags); `1.0`, `t` and `y` now count as true everywhere.
- The selected-hierarchy null audit picks the root target by `parent_depth == 0` and non-root targets by `parent_depth > 0`; a case whose root record is not selected raises instead of substituting the largest open parent.
- The manifold comparison plot raises `ImportError` when `umap` is missing instead of plotting PCA under a UMAP label.
- The artifact catalog no longer stores archive locations; `check --verify-archive` and `restore` use `~/kl-te-cluster-archive` unless `--archive-root` is given.
- `make lint` also runs `ruff format --check .`; the `[tool.ruff] src` entry that hid first-party import ordering is gone, and imports are sorted project-wide.


### Removed
- `benchmarks/shared/util/pdf/image_concat.py` and `figure_split.py` (`concat_*_plots`, `split_collected_figs_to_pdfs`): no caller; the pipeline writes to `PdfPages` directly.
- The `diffusion`, `trajectory` and `brancharchitect` extras with deeptime, PyGSP, PHATE, kaleido, biopython and orjson; none was imported.
- BranchArchitect: the `vendor/BranchArchitect` submodule, `.gitmodules`, and `benchmarks/validation/tree/brancharchitect_tree_comparison.py`. Existing clones should delete the leftover `vendor/BranchArchitect` folder and `.git/modules/vendor`.
- The unreachable 3D-UMAP plot chain (`plot_umap_3d`, `create_umap_3d_plots_from_results`, `create_clustering_comparison_plot_3d`, `concat_umap_3d_plots`, the `umap3d` category) and the `TBS_EMBEDDING_BACKEND_3D` variable only it read.
- `benchmarks/diagnostics/analysis/distributional_action.py` and the scRNA edge-export loop that copied the `Distributional_*` columns only it wrote.
- The second `is_leaf` (`projection.spectral.tree_helpers.is_leaf`, re-exported as `projection.is_leaf`); the spectral estimator uses `tree.topology.is_leaf`.
- Three tests that could not fail on their own (the near-threshold decomposer test and two `REPO_ROOT` checks) and the orphan `skills-lock.json`.
- IQ-TREE 3 tree building (`iqtree3_tree_from_alignment`, `tbs_iqtree3`, `iqtree_*` parameters).
- Neighbor joining and MAD rooting (`tree/construction/phylogenetic.py`, `tbs_neighbor_joining`, the adaptive-K neighbor-joining cell), and `SUPPORTED_TREE_BUILDERS`, `TreeBuilderName`, `LINKAGE_TREE_BUILDER` and `DEFAULT_LINKAGE_TREE_ROOTING`.
- Unused modules: `hierarchy_analysis/bootstrap_consensus.py`, `multiple_testing/stopping_edge_recovery/`, and `sibling_divergence/neighborhood_bandwidth/`.
- The all-pairs tree-distance cache built for `neighborhood_bandwidth_profile`, which only recorded a pair count.
- The `matrix_audit` option and its TensorFlow matrix export, which no caller enabled and no dependency group could run, plus the unread `TBS_MATRIX_AUDIT_ROOT` variable.
- Duplicate helpers: the guard's copy of the pass-through rule (now `GateEvaluator`), `TreeDecomposition._annotation_bool_value`, diagnostic copies of the core inflation formulas, `_node_pvalue`, and repeated `bool_value`, `parse_csv_list`, JSON-default, `finite_float`, git-provenance, timestamp and AWS shard helpers; the copied traversal skip summary and projection-input helpers; and the `read_dataframe`, `_as_bool`, `_compute_js_divergence` and annotation-PDF `load_binary_matrix` pass-through wrappers.
- More copied helpers (`_select_cases`, `_selected_cases`, `_df_bin`, `_rank_bounds`, `_majority_label`, `_branch_length_to_parent`, `_bool_sum`, `_root`, `safe_float`, `sha256_file`, the pancreas `descendant_leaves` tree helper, the renamed copies `parse_alpha_grid`, `_median` and `_load_cases`, the AWS `iter_expected_shard_dirs`, S3 sync and JSON writers (now in `benchmarks/cloud/sharding.py`), the root diagnostics' `_require_columns` and `_finite_median` (now `root_values.finite_median`), `_coerce_bool_series` (now `values.bool_series`), and the spectral profiler's copy of the Bernoulli whitening rule); each module imports the one kept copy.
- The `tbs_bandwidth_context_v1` benchmark method and its `neighborhood_bandwidth_profile` parameter: since the bandwidth module was removed it ran plain `tbs` and only added metadata. The parameter is now rejected as an unknown key.
- Test-only distributional-action API: `edge_distributional_action`, `binary_split_distributional_action`, `binary_split_distributional_action_summary` and `BinarySplitDistributionalActionSummary`; use `edge_distributional_action_summary` and `split_distributional_action_summary`.
- Unread data: the GO dataset inventory's `ResultKind` folder templates and descriptions (`RESULT_KINDS` is now a tuple of codes) and the unused `DEFAULT_RESULT_ROOT` of the root hard-negative replay panel.
- The endotypes KAK page column `common_invariant_axis_score` is now `common_axis_score` (and `is_common_invariant_axis` is `is_common_axis`).

### Fixed
- The full benchmark's relationship analysis no longer fails with `KeyError: 'gaussian_sparse_signal_highd_noise'`; the category is registered in the cover-page group map, and a test checks every default case category.
- Reject duplicate `leaf_data` sample labels instead of reusing the last matching row.
- Permutation guards keep fractional Bernoulli and soft categorical values instead of truncating them to zero, and no longer test splits that traversal cannot reach.
- A non-root permutation-guard scope with zero replicates now raises instead of silently skipping the guard.
- Isolated benchmark cases read their results before waiting for the subprocess, so large results no longer hit the timeout, and a child that lingers after sending its result is stopped at the deadline.
- Hamming diffusion with `k_neighbors=1` keeps a non-self neighbor instead of returning all-zero distances.
- Endotype gene lookups no longer store Ensembl IDs as Entrez IDs.
- The GO annotation level audit lists unreadable config and summary files in a `read_errors` column instead of skipping them silently.
- Diagnostic panels read boolean cells the same way: missing values are false (two copies treated them as true) and text such as `"false"` is false.
- The benchmark relationship analysis and failure debug trace read a 0/1 column that has a blank cell correctly; before, every `1.0` was read as false.

## [Unreleased] - 2026-09-09

### Changed
- Restore the historical context-weighted empirical-null chi-square p-value rule and traversal-aligned sibling BH when internal support exists. `internal_admissible` denotes availability under the support policy, not a validated selected-tail law.
- Retain dependency-group support diagnostics, optional support thresholds, the exact independent common-scale F API, and explicit no-internal-support outcomes.
- Version calibration sweep artifacts to prevent reuse of cached withheld-p-value outcomes after restoration.

### Fixed
- Allow hashable tuple dependency-group IDs in leave-one-group scale-sensitivity diagnostics.
- Pass production sibling records through the requested support-threshold decision instead of bypassing it.

## [Unreleased] - 2026-08-11

### Added
- A separately named exact-F sibling calibration path for independent focal and calibration statistics with pairwise-disjoint observation ownership, independent calibration records, unit fixed weights, and a common chi-square scale.
- Immutable calibration-sample and versioned support-policy types, including parent provenance and stopping-event dependency groups.

### Changed
- Internal support diagnostics now count dependency groups rather than nested records for supported-group counts, effective support, maximum group weight share, and leave-one-group stability.
- Positive-dimensional same-selected-hierarchy sibling calibration now preserves raw evidence but returns `undefined_unvalidated_reference_law`; it no longer reports a calibrated chi-square p-value.

### Fixed
- Reject invalid support thresholds, including booleans as counts, nonfinite values, invalid weight shares, and negative stability limits.
- Reject the invalid sibling-record role `is_edge_blocked=True, is_null_like=False` and remove the redundant stopped-or-null calibration count and threshold.

## [Unreleased] - 2026-02-17

### Added
- **Spectral dimension estimation** (`spectral_dimension.py`): Per-node eigendecomposition of the local null-whitened tangent correlation matrix provides the projection dimension and parent PCA basis for both edge and sibling projected-Wald tests.
- **Dual-form eigendecomposition** (`spectral_dimension.py`): When `n_desc < d_active`, computes the `n×n` Gram matrix instead of the `d×d` correlation matrix — O(n²d + n³) vs O(d³). For subtrees with n=10 leaves and d=2000 features, this is 10×10 eigh instead of 2000×2000, eliminating the performance bottleneck on high-dimensional cases.
- **Internal node distributions in spectral decomposition** (`spectral_dimension.py`): The data matrix for eigendecomposition now includes both leaf rows AND internal descendant node distribution vectors. This enriches the covariance estimate, especially for nodes high in the tree where internal descendants capture intermediate subtree structure.

### Refactored
- **Step 3.4 — `gates.py` extraction**: Gate logic (`should_split`, `should_split_v2`, `_check_edge_significance`) extracted from `TreeDecomposition` into `GateEvaluator` class in `gates.py` (343 lines). `tree_decomposition.py` reduced from 992 to 757 lines (−234 lines). `decompose_tree()` and `decompose_tree_v2()` now delegate to free functions `iterate_worklist`, `process_node`, `process_node_v2` in `gates.py`. `GateEvaluator` constructor accepts injected `children_map`, `descendant_leaf_sets`, `root` to decouple from `PosetTree` internals. Five dead inline methods removed.
- **Power guard in `should_split_v2`** (`gates.py`): When signal localization finds zero significant difference pairs after BH correction, returns `(True, None)` — trusting the aggregate sibling-divergence gate SPLIT but discarding the powerless localization result. Prevents false cross-boundary merges via similarity-only edges.

### Fixed
- **Benchmark hang on high-dimensional cases**: `gaussian_extreme_noise_2` (n=300, d=2000, K=30) previously hung due to ~598 × O(2000³) eigendecompositions at every internal node. Dual-form optimization and information cap eliminate the hang.
- **Test import error** (`59_test_pipeline_pdf_naming.py`): Updated import from removed `benchmarks.shared.pipeline._resolve_pdf_output_path` to `benchmarks.shared.util.pdf.session.resolve_pdf_output_path`.
- **`merge_similarity_graphs` p-value direction** (`signal_localization.py`): When duplicate similarity edges existed across localization levels, kept the **lower** p-value. For similarity edges, higher p-value = stronger evidence of similarity (fail to reject H₀). Fixed to keep the **higher** p-value. Previously caused false cross-boundary merges by under-reporting similarity strength.

## [Unreleased] - 2026-02-14

### Changed
- Remove obsolete `n_permutations` parameter from the decomposition API; callers passing it will no longer be accepted. Tests updated accordingly.
- Pipeline plotting behavior: default no longer writes intermediate PNGs; when `concat_plots_pdf=True` the pipeline collects Figures and writes categorized PDFs (`k_distance_plots.pdf`, `tree_plots.pdf`, `umap_plots.pdf`) instead of emitting PNGs.
- PDF utilities: improved diagnosis, figure classification (manifold plots grouped with UMAP), and robust headless handling (Agg backend).

### Removed
- `kl_clustering_analysis.threshold` package removed — helper functionality that was required by debug scripts now lives in those scripts or should be implemented directly where needed.

### Fixed
- Various import-time and optional-dependency issues (runners now import optional dependencies lazily and return skip results when missing) to improve test/CI stability.
- Resolved multiple test failures and tightened integration behavior for PDF generation.


(For full history, add older entries here in reverse chronological order.)
