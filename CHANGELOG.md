# Changelog

All notable changes to this project are documented in this file.

## [Unreleased] - 2026-10-14

### Changed
- Fixed-coordinate BH and BY minima use `scipy.stats.false_discovery_control`; values are unchanged.
- Trees built from sample data go through `build_tree`, which fixes the sample order before linkage, in the quick start, the selected-root permutation and root-stability guards, the MNIST and UMAP experiments and five diagnostics, including the oracle cut of the data-independent gate panel. Tied distances now resolve the same way on these paths, so they can return a different tree on tie-heavy data. `run_feature_matrix_with_umap.py` keeps row-order linkage because its dendrogram and `cut_tree` labels reuse it. `tree_from_linkage` stays the low-level constructor.
- `benchmark_cluster_algorithm` returns the results DataFrame only instead of `(df, None)`. The plot export helpers write PDF pages only and no longer take `output_dir`, `save`, `collect`, `collected` or `timestamp`; report PDFs no longer include an empty "Real Data" section page, and the evaluation-metrics text now ends the `method_proof` section page.
- `build_benchmark_result_row` takes one `metrics: ClusteringMetrics` instead of 27 metric keywords, and the result columns are derived from `BenchmarkResultRow` with the same names and order.
- Adaptive diffusion metadata records the eigensolver actually used (`scipy.linalg.eigh` or `scipy.sparse.linalg.eigsh`).
- Renames: the decomposition result key `independence_analysis` is `gate_alphas`; MNIST `assignment_key` is `sweep_assignment_key`; same-name helpers with different behavior are `format_counts`, `cluster_size_ranks`, `load_relabeled_assignment` and `plot_polar_linkage_tree`. `PosetTree.get_leaves` and `topology.get_leaf_values` lose their never-used id/label flags.
- `independent_gaussian` writes edge BH values to `Child_Parent_Divergence_P_Value_BH`, the column traversal traces read.
- The NNLS null-calibration sweep schema is `nnls_null_calibration_sweep/v6` without 15 unread parent-spectrum columns. Selected-hierarchy covariates report `eigenvalue_effective_rank` as NaN when no eigenvalue is positive (it was 1.0); one `effective_rank`/`spectral_entropy` helper uses `scipy.stats.entropy`.
- The GO tree-matrix pipeline drops the single-value `--method-versions` flag; saved `method_version` values and `current__` run ids are unchanged.
- Shared code replaces copies: the KAK lens sweep uses the production block diffusion geometries, the AWS drivers share one shard layout in `benchmarks/cloud/sharding.py`, the MNIST and UMAP experiments share `tbs_decompose`, the scRNA audits share their dataset and method tables, and the subspace roster exporter uses the PDF builder's folder-naming rule.
- `SpectralContext` is gone; the edge gate returns `SpectralDecompositionResult` with its own timing added.
- Naming audit: `PosetTree.populate_node_divergences` is now `populate_node_distributions`; `calculate_ari_nmi_purity_metrics` is now `calculate_clustering_metrics`. Active callers, tests, quick start and ASV timing names use the new names without compatibility aliases.
- Coordinate adjustment APIs use `CoordinateAdjustmentMethod`, `coordinate_adjusted_minimum` and `fixed_coordinate_adjusted_p_value`, because the accepted methods include the Bonferroni FWER procedure as well as BH/BY FDR procedures. The redundant `fixed_coordinate_bh_p_value` wrapper is removed; callers pass `method="bh"`.
- Unweighted spectral input is named `raw` instead of `binary`; the averaged NNLS target is named `mean_squared_standardized_euclidean` instead of `squared_standardized_euclidean`. CLI choices, callers and current output metadata use these names; strict spectral readers reject retired labels. Historical outputs retain their original names.
- The traversal sibling diagnostic uses `any_false_rejection_rate` and `empirical_control_observed`, with schema `traversal_sibling_fdr_null/v2`. Its local/AWS merger rejects old or mixed-schema shards. Probability formulas and decision thresholds are unchanged.
- NNLS tree labels distinguish fitted branch length from calibrated time and disclose radial normalization. Projected-Wald and Tree-BH docstrings identify the actual return type and within-family adjusted values.
- Pass-through traversal counts a deeper split only when every node on the path would split or pass through. Runs with pass-through enabled (the default) can return fewer clusters than before.
- `build_tree` builds linkage trees only and no longer takes `builder` or `rooting`; `TreeBuildResult` drops `builder`, `rooting` and `phylogenetic_rooting`. Benchmark methods no longer accept `tree_builder` or `tree_rooting`.
- `run_clustering_result` rejects unknown parameter keys for every registered method; `benchmark_*` bookkeeping keys remain allowed.
- Decomposition results report `sibling_alpha_effective` (sibling alpha divided by `sibling_gate_alpha_penalty`) next to `sibling_alpha`.
- Core Benjamini-Hochberg and Tree-BH use `scipy.stats.false_discovery_control`. The scipy floor is 1.11; `statsmodels` is no longer a core dependency and remains only in the `all` extra.
- Full-benchmark resume tracks completed cases by `case_id` instead of suite position.
- Diffusion coordinates accept sparse kernels; connected graphs above 2000 samples use `eigsh` for the leading modes.
- The graphtools adaptive-K tree grid and tree consensus cover the supported linkage methods only, named by `SUPPORTED_LINKAGE_METHODS`; the `TREE_CONSENSUS_STRATEGIES` and `TREE_PRIORITY` aliases were removed.
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
- `tree_break_selection/core_utils/` is gone: `tree_root`, `bottom_up_nodes` and `compute_node_depths` live in `tree/topology.py`, and `data_utils.py` moved to `hierarchy_analysis/statistics/data_utils.py`.
- `hierarchy_analysis`, `hierarchy_analysis.statistics`, `statistics.projection` and `statistics.sibling_divergence` no longer re-export names; import from the defining module. Importing one statistics module no longer loads the whole statistics stack.
- The rectangular tree layout, the internal-node halo, the linkage tree builder and three diagnostics find the root with `tree.topology.tree_root`. The layout no longer lays out forests, and the NNLS null-calibration sweep raises instead of reporting an empty `tree_root`. Tree-BH `find_root_nodes` and sibling-FDR `_root_nodes` keep returning every parentless node: the gate fixtures use a detached calibration component.
- `benjamini_hochberg_correction` and `apply_tree_bh_correction` require `alpha`; the unused `0.05` defaults are gone.
- Gate and guard root-stability and permutation-guard defaults, `fit_independent_gaussian_tree(linkage_method=...)` and the TBS runner mirrors use `DEFAULT_BINARY_TREE_DISTANCE_METRIC` and `DEFAULT_TREE_LINKAGE_METHOD` instead of retyped `"hamming"`/`"average"`.
- Seven scalar coercion copies use `values.finite_float` or `root_values.finite_median`: the alpha-structure sweep and gate-path trace map `inf` to NaN, and the selected-edge geometry panel maps unparseable cells to NaN instead of raising.
- `benchmarks/README.md` documents every `TBS_*` variable; the dated 2026-04-27 and 2026-07-27 repository audits moved from `docs/` to `raw/inbox/`.
- `tree_from_linkage` requires `leaf_names` and raises when their count differs from the leaf count; it no longer invents `leaf_{i}` names.
- `benjamini_hochberg_correction` returns `(rejected, adjusted)`; the third value only echoed `alpha`.
- `extract_bool_column_dict` always keeps the index type; the `coerce_index_to_str` option is gone.
- Cluster metadata must carry `leaf_signature`; `build_sample_cluster_assignments` no longer derives a missing one.
- `BranchLengthOptimizationResult` requires `target_definition`, `target_normalization` and `target_normalization_divisor`.
- Full-benchmark resume raises on a results CSV without `case_id`/`run_id` instead of rerunning every case into it.
- Spectral branch-length state reads `branch_length` only; the `length` and `weight` edge aliases are gone.
- `go_ic_tree_summary_plots.py` requires `--method-matrix-summary`; it reads method-matrix rows only.
- The benchmark Docker image installs `scipy>=1.11.0`, matching `pyproject.toml`.
- Tree-BH runs one top-down pass: each family's level is its parent family's level times that family's rejection share. Results are unchanged; a 4,000-leaf single-linkage tree takes 0.08 s instead of about 2 s. `ChildParentEdgeTreeBHResult` drops `tree_bh_base_alpha_by_depth` and `sibling_group_outcomes`, and `TreeBHSiblingGroupOutcome` is gone; nothing read them.
- The fixed-topology NNLS design is built from a sparse leaf-by-edge incidence matrix in chunks; the matrix is identical, about 4x faster to build and half the peak memory on chain-shaped trees.
- Fixed-subspace evidence columns compute all block p-values with one vectorized `chi2.sf` call and write each column once; on a 1,500 x 600 binary case they take 1.1 s instead of 13 s, with identical values.
- The full benchmark refuses to append rows whose columns differ from the existing results file.
- The AWS selected-tail merge raises when a shard's `case_summary.csv` reports a case that did not run, and the KAK lens merge raises for a lens/linkage group with no `ok` row and records status counts in its manifest.
- `tabulate` is declared in the `benchmark` and `all` extras; `DataFrame.to_markdown` reports needed it.
- Config cleanup: the `all` extra is `tree-break-selection[dev,benchmark,benchmark-infra,viz,scrna]` plus `statsmodels` (same resolved packages); `make wiki-lint` runs through `uv run --no-sync`; `asv.conf.json` points `project_url` at the real repository; the AWS image no longer installs the unused `statsmodels`; the Batch job role no longer grants `s3:DeleteObject`; six `.gitignore` entries for paths nothing creates are gone.
- `decompose_leading_axes_space` always standardizes columns; the `standardize` option is gone. `build_cluster_color_spec` loses `base_cmap` and `plot_tree_with_clusters` loses `colormap` and `figsize`.
- Diagnostics share `finite_int`, `is_finite`, `finite_min` (finite values only), `numeric_column`, `split_path`, `input_data_role` and `output_data_role` from `benchmarks/shared/util/values.py`, and `load_shard_table`, `load_shard_manifest` and `resolve_synced_input` from `cloud/sharding.py`.
- `SUPPORTED_LINKAGE_METHODS` drops `centroid` and `median`, and `build_tree` rejects any linkage method outside it by name. Both produce nonmonotone merge heights on almost every input (all 100 random binary and Euclidean inputs tried), which `build_tree` rejects, and a failing cell of the graphtools adaptive-K tree-strategy grid takes its whole benchmark case down. The grid, tree consensus and the family-metric sweep now cover average, complete, weighted, single and Ward.
- The selected-root permutation guards build sibling contrasts with `sibling_fixed_subspace_contrast`, the helper the fixed-subspace gate uses; the guard-only discrete fast path is gone.
- `sibling_divergence_test` requires keyword-only `branch_length_sum` and `mean_branch_length` (pass `0.0` and `None` for a tree without branch lengths). The projected-Wald sibling test always scales its contrast variance by `1 + (L_left + L_right) / (2 * mean positive branch length)`; the fixed-subspace sibling gates and guards use no branch time, and `edge_branch_length_variance_policy` governs only the edge test. Both rules are documented and pinned by a test.

- `TreeDecomposition` no longer takes `selected_family_passthrough_guard`; only a gate bundle whose `root_selective_permutation_guard_scope` is a pass-through scope turns the selected-family guard on.
- Tests now check production outcomes: the root-stability and selective-permutation guard tests assert that the root is blocked, the pipeline parity test compares every non-evidence column and adds a benchmark case where decisions depend on alpha, three tests build trees with `build_tree` (average linkage) or `run_tbs_on_distance` instead of hand-built complete-linkage trees, and the supplied-distance dispatch tests assert the runner receives the supplied matrix.

- Ambiguous names: `SiblingPairRecord.has_empirical_null_support` (role only) is `has_null_role`; `CalibrationDecision.exact_context` is `decision_context`.
- One implementation per helper name: `parse_go_term`/`GO_ID_RE` and `default_output_dir` live in `applications/endotypes/_shared.py`, one `load_binary_matrix` in `applications/endotypes/io/`, one `markdown_table` in `benchmarks/shared/performance_grid.py`, and `aws_alpha_grid_search` uses `values.parse_csv_list`. A GO label without an id now comes back stripped, with an empty id.
- Research-tool tests moved from `tests/core` and `tests/pipeline` into `tests/validation`; timestamp tests restore `sys.modules`; tests of the null-whitened tangent and spectral jobs go through `build_null_whitened_tangent_matrix` and `compute_spectral_decomposition`.

- `benchmarks/diagnostics/calibration/values.py` moved to `benchmarks/shared/util/values.py`, so applications and shared benchmark code read these helpers outside the diagnostics tree. `parse_csv_list` is the one comma-list parser (it accepts `None`); `alpha_grid_search.select_cases` is the one case selector (the selected-hierarchy diagnostics reject an empty `--case-names`); `family_metric_nnls_grid._bool_sum` is `values.bool_sum`.
- About 30 helpers imported across `applications/` and `benchmarks/` modules lost their leading underscore (for example `normalize_labels`, `ok_result_from_labels`, `case_contract`, `run_benchmarks`); module-private helpers and the `_run_*_method` registry runners keep theirs.


### Removed
- The Holm coordinate adjustment and its two benchmark methods (`tbs_fixed_coordinate_holm` and `tbs_fixed_coordinate_holm_nnls`): Holm's minimum adjusted p-value always equals Bonferroni's.
- Write-only or test-only code: `TreeBuildDiagnostics` and the post-build tree validator, `benchmarks/shared/tbs_tree_context.py`, `predict_empirical_inflation_factor`, the `math_trace` failure classifier and trace validators, the `preloaded` case generator, seven `values.py` helpers, `SUMMARY_METRICS`, four `write_diagnostic_bundle` options, `ClusterAnnotation.representative_gene_proteins`, 15 parent spectral-summary fields on sibling records, the always-false edge invalid flag with its `Child_Parent_Divergence_Invalid` column, and `full_component_eigenvalues_by_node`.
- `scripts/rest_request.py`; the GO annotation PDF builder calls the QuickGO and UniProt REST APIs directly with `requests`.
- Verification code only research or tests reached: `feature_covariance_calibration.py`, `selected_pca_projected_wald_calibration.py` and `alpha_structure_sweep.py` in `benchmarks/validation/`; `traversal/cli.py`, `mixed_internal_calibration_sweep.py`, `internal_support_threshold_validation.py`, `sibling_null_weight_rule_validation.py`, `selected_edge_geometry_analysis.py`, `case_geometry_audit.py` and `independent_projection_control.py` in `benchmarks/diagnostics/`, with the tests that only covered them; `calibration_support.summarize_p_value_calibration`.
- 71 past-study files in `benchmarks/diagnostics/` (32 entry scripts, 29 helpers only they imported, 10 empty packages, plus `runner_support.py`): nothing imported or tested them, and their results are recorded in the wiki. They remain in git history at 7cc4e1c0; `spectral/apply_vandermonde_diagnostic.py` and `path_b/recursive_pvalue_geometry.py` stay because retained reports rerun or import them.
- Six duplicate tests, the test-only `EmpiricalNullInflationModel.sample_is_strict_null`/`sample_is_edge_blocked` properties, `LeadingAxesSpace.leading_coordinates`/`following_coordinates`, and four unused `project_root()` copies.
- Code only tests reached: the `selected_gaussian_{radial,hierarchy,diffusion}` prototypes and tests 51–54; the exact-F calibration mode (`fit_independent_unweighted_common_scale_inflation_model`, `decide_independent_unweighted_common_scale_calibration`, `IndependentUnweightedCommonScaleContract`, the model's `reference_law`/`independent_contract` fields and the `independent_exact_admissible` status); `split_distributional_action_summary`.
- Code only tests and `reports/` scripts reached: the independent Gaussian discovery/test API (`independent_gaussian.py` with `fit_independent_gaussian_tree`, `test_independent_gaussian_tree` and `IndependentGaussian{Tree,Result,Diffusion}`), its conditional energy reference `energy_selected_reference.py` and tests 56–57. The four `reports/*_20260928/` scripts that imported it no longer run against current code; their recorded outputs stay.
- The never-produced `undefined_unvalidated_reference_law` calibration status and `UnsupportedReasonCode.EMPIRICAL_NULL_UNVALIDATED_REFERENCE_LAW`.
- `benchmarks/external/` (clustbench and OpenML loaders) and the `clustering-benchmarks` and `openml` dependencies of the `benchmark-infra` extra.
- Duplicate guard outputs: the `Root_Selective_Permutation_*` columns (read `Selective_Permutation_*` on the root row) and the gate-config keys echoed into TBS `result.extra` (read `extra["gate_bundle"].metadata.config`; `sibling_gate_method` stays).
- Package re-exports nobody imports, in `plot`, `space_separation`, `tree.construction`, `projected_wald`, `projection_dimension_estimation`, `multiple_testing`, `tree_bh` and `benchmarks.shared.util.pdf`.
- `SMALL_TEST_CASES` (moved into its one test); the empty defaults of six `SpectralContext`/`SpectralDecompositionResult` fields, which are now required; the `old_display_rank` ranking column; relationship-analysis support for audit CSVs without gate or branch-length columns.
- Compatibility leftovers: `report_helpers` re-exports of `artifact_record`/`sha256_file`; the ignored `save_png`/`collect_figs` parameters of `generate_benchmark_plots`; `_resolve_requested_cluster_count` and its ignored `n_samples` (k-means and spectral read `params["n_clusters"]`); the private `_labels_from_decomposition` wrapper; the `specificity_aware_rank` and `_rank` fallbacks in the subspace annotation PDF; the `run_id` back-fill for artifact indexes written before 2026-06-18.
- Tests that only guarded removed parameters or a re-export: the two `build_contrast_covariance` branch-length `TypeError` tests, the `tree_builder` absence asserts, the scRNA distributional-action absence test, and `test_goncalves_benchmark_reuses_scrna_method_configs` (it failed when run after the manifest timestamp tests).
- `sibling_z_scores.py` (`_compute_sibling_z_scores`): a pass-through to `compute_whitened_wald_contrast`, now called directly by `sibling_divergence_test`.
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
- Test-only distributional-action API: `edge_distributional_action`, `binary_split_distributional_action`, `binary_split_distributional_action_summary` and `BinarySplitDistributionalActionSummary`; use `edge_distributional_action_summary`.
- Unread data: the GO dataset inventory's `ResultKind` folder templates and descriptions (`RESULT_KINDS` is now a tuple of codes) and the unused `DEFAULT_RESULT_ROOT` of the root hard-negative replay panel.
- The endotypes KAK page column `common_invariant_axis_score` is now `common_axis_score` (and `is_common_invariant_axis` is `is_common_axis`).
- Dead symbols: `write_cover_pages_to_pdf`, `list_clustbench_batteries`, `list_clustbench_datasets`, `projected_quadratic_directional_derivative` and `EmpiricalNullInflationModel.sample_parent_ids`.
- Unread parameters: `plot_multi_split_heatmap(branch_lengths, n_groups_list)`, `_run_tbs_distance_tree_method(method_id)`, `build_lens_coordinates(data)` and `geometry_model_specs(table)`.
- Write-only fields: `TreeBuildDiagnostics.zero_branch_length_count`, `ClusterColorSpec.unassigned_color`, `SiblingInflationInputs.non_binary_nodes` and `SiblingInflationTrace.raw_degrees_of_freedom`.
- The `TBS_FORCE_UMAP_FOR_LARGE` variable, which the full runner set and nothing read.
- The `.mex/` compatibility bridge for the `mex-agent` CLI (router, setup and sync notes, an empty pattern index, `config.json`) and its two contract tests. It only pointed back into `wiki/`, which is now the single project memory system; `AGENTS.md` says so and a wiki contract test fails if `.mex/` returns.
- Compatibility paths nothing still needs:
  - In the package: `PosetTree.build_sample_cluster_assignments`, the decomposition result's constant `decision_mode` key, the `kernel_backend` metadata copy of `backend`, the orchestrator re-exports of `SIBLING_GATE_PROFILES`, `resolve_sibling_gate_profile` and the two permutation p-value helpers (import them from `gates.profiles` / `gates.guards`), and the missing-column guards on `Sibling_Divergence_P_Value`.
  - In the benchmarks: the `suite="default"` alias for `full`, the `resolve_shard_index` alias, the AWS selected-edge `--resume/--no-resume` flag that changed nothing, the `row_status`/`lens_tree_id` and `getattr` defaults in the KAK lens panels, the `debug_trace` fallbacks for missing result columns, and `hasattr(tree, "root")` lookups that returned `None` (six diagnostics now call `tree_root`).
  - In the applications: the `--raw-kak-dir`, `--adaptive-diffusion-kak-dir` and `--whole-adaptive-dir` inputs of `go_ic_tree_summary_plots.py`, its `*_kak` family renames and assignment-column aliases, the artifact-index `display_rank` and `method_run_id` fallbacks, and the roster reader for unnamed linkage columns.
- Dead code: `applications/endotypes/analysis/assess_method_correctness.py`, `benchmarks/diagnostics/calibration/statistics/same_z_functional_decomposition.py`, `benchmarks/diagnostics/calibration/cli.py`, the unused `strict` flag of `benchmark_cluster_algorithm` and the case runners, the always-true `compute_eigendecomposition_outputs` spectral flag, `FeatureBlock.signature`/`FeatureSpace.signature`, the `TREE_CONSENSUS_PRIORITY_ORDER` export, the `--progress` pytest plugin, and local copies of `parse_csv_values`, `parse_names`, `dataset_slug`, `select_cases` and the oracle runner's runtime and output-directory setup.

### Fixed
- `get_env_bool` error messages show the rejected value instead of the literal text `{raw!r}`.
- `Sibling_Role_Supported` and the empirical-null gate share one rule, `SiblingPairRecord.has_internal_empirical_null_support` (positive degrees of freedom, positive null weight, null role). The stamp used to count df-0 merges of identical rows, so `binary_perfect_2c`, `_4c` and `_8c` came back `ok` with one cluster and ARI 0; they now report `unsupported` (`EMPIRICAL_NULL_NO_INTERNAL_SUPPORT`). No other default-benchmark row changes.
- A sibling-gate profile rejects an explicit setting that differs from it even when the setting equals the unprofiled default; before, `sibling_gate_profile="fixed_coordinate_guarded_v1"` with `sibling_gate_method="projected_wald_inflation"` silently ran the profile's method. `None` now means "not given": `resolve_sibling_gate_profile_config`, `run_gate_annotation_pipeline`, `run_tbs_on_distance` and benchmark dispatch default the ten profile-controlled settings to `None`, and the fixed-profile validation CLI leaves unset guard options to the profile.
- `--resume` in the alpha-grid, alpha-structure and mixed internal-calibration sweeps refuses saved results written for other cases, methods or settings instead of mixing them into the new run. Alpha-grid summaries record a `case_set_sha256`; summaries without it are not reused.
- Fixed-topology NNLS adds one crossing leaf pair for each edge the sampled pairs miss; those edges used to keep the solver's starting value of 0.1 under status `ok`. `n_zero_design_columns` now counts the edges that needed such a pair.
- Fixed-topology NNLS splits the identified root-pair sum in the ratio of the two root edges' linkage lengths instead of the solver's incidental even split. On 122 default cases with three NNLS methods, this changed one result: `gauss_outlier_cluster_4c_continuous` under `tbs_fixed_coordinate_bh_nnls` goes from 1 to 2 clusters (ARI 0 to 0.46).
- The MP projection-dimension sweep passes each sibling pair's branch-length sum and the tree's mean branch length, so its sibling statistics match production instead of omitting the branch-time factor.
- `build_cluster_color_spec(0)`, `overlap_selected_pass_through_fixture_miner` and `mixed_null_signal_geometry_validation` no longer crash; the previous cleanup left a removed field and two removed helpers referenced. A test now imports every `benchmarks/` module.
- The tree plot maps sample labels to leaves only, so a sample named like an internal node no longer breaks colouring; `lowest_common_ancestor_for_set` works on trees without a cached root.
- The GO-IC tree pages seed their subspace UMAP from a CRC32 of the run id instead of the per-process salted `hash()`, so layouts repeat across runs.
- `profile_spectral_backends.py` calls `_process_node` with its current signature.
- `make check` runs every test stage again: the stage-2 non-finite-handling fixture lacked `is_leaf`, so `scripts/run_tests_ordered.py` stopped before stages 3-7.
- The full benchmark keeps going when a case times out or crashes: the case is skipped and not recorded, the remaining cases and reports run, and the run exits non-zero naming the failed cases. Rerunning with the same `TBS_RUN_DIR` retries them.
- The endotype feature-matrix pipeline stores the gate-annotated frame on the tree, so its default tree plot no longer raises `KeyError: 'Sibling_Divergence_Skipped'`.
- The adaptive-diffusion GO-IC ranking keeps all-failed subspaces as unranked rows instead of raising `KeyError: 'display_rank'`, and an empty block selection raises a clear error.
- `discrete_covariance_guard.py` requires `--output-dir` and `local_mp_identity_law_diagnostic.py` defaults to `benchmarks/results/`; both previously wrote into tracked `raw/assets/` evidence. The subspace roster export, the per-method GO-IC split and the GO-IC tree pages delete earlier outputs only after their inputs load, and the GO annotation pipeline's `--dry-run` prints the plan without writing logs, manifests or audits.
- The scRNA and MNIST scripts no longer write into retained `raw/assets/benchmark-results/` evidence by default. Producers (both benchmarks, `analyze_space_decomposition.py`, the selected-NNLS reports, the distributional-action and branch-length audits, and both MNIST renderers) default to `results/analyses/<report>_<YYYYMMDD_HHMMSS>/` via `applications/output_dirs.py`; scripts that read and write one benchmark run in place (`plot_pipeline.py`, the UMAP and R tree renderers, `analyze_goncalves_progenitors.py`, `analyze_pancreas_inner_nodes.py`, `compare_pancreas_progenitor_signatures.py`, `goncalves_progenitor_coherence.py`) require `--output-dir`.
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
