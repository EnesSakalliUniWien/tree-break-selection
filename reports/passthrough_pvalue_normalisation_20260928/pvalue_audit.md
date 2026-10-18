# Passthrough p-value audit

## Finding

Passthrough does not recompute p-values or change the tested sibling family in the application pipeline. The endotypes app completes `run_gate_annotation_pipeline` before constructing `TreeDecomposition` ([application path](../../applications/endotypes/pipelines/run_feature_matrix_with_umap.py#L271-L281)); the annotation function has no passthrough argument ([orchestrator signature](../../tree_break_selection/hierarchy_analysis/decomposition/gates/orchestrator.py#L300-L336)). `TreeDecomposition` passes the flag to `GateEvaluator` ([constructor](../../tree_break_selection/hierarchy_analysis/tree_decomposition.py#L66-L75), [wiring](../../tree_break_selection/hierarchy_analysis/tree_decomposition.py#L146-L156)), whose decision branch uses it to choose `BOUNDARY`, `SPLIT`, or `PASS_THROUGH` from already annotated gates ([decision](../../tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py#L213-L237)). Traversal then visits children only for `SPLIT` and `PASS_THROUGH` ([loop](../../tree_break_selection/hierarchy_analysis/tree_decomposition.py#L635-L645)).

For the projected-Wald inflation path, `Sibling_Divergence_P_Value` is the empirical-inflation plug-in p-value supplied to BH, and `Sibling_Divergence_P_Value_Corrected` is its traversal-aligned BH adjustment ([calibration and FDR call](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/pipeline.py#L184-L220), [BH fields](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py#L134-L177)). This adjustment does not establish selected-tail validity. The nominal record p-value is not retained in a separate production annotation field. The sibling BH frontier is formed from edge-reachable candidates, before clustering traversal, so it does not depend on passthrough ([frontier construction](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py#L153-L182)).

## Real-data production probe

The reproducible probe uses the endotypes app's default binary matrix (626 samples, 456 features), Hamming/average linkage, edge alpha 0.001, sibling alpha 0.01, `projected_wald_inflation`, and support contract v2. Its one gate run produced 1,251 annotation rows and 37 fields. There are 11 fields with `P_Value` in their names, including two categorical calibration/role fields; the remaining nine are numeric p-value fields. The production sibling p-value and BH-adjusted p-value were present for the same 345 nodes, all marked `empirical_null_inflation`; the edge raw and BH fields were present for 1,250 and 718 nodes.

The exact same completed `GateAnnotationBundle` was passed to both traversals. Assertions compare each decomposer's `annotations_df` and the shared bundle against a deep snapshot: all 37 annotation fields and all nine numeric p-value fields were identical. Passthrough disabled reached 1 internal node and produced 1 cluster. Enabled reached 182 internal nodes and 35 leaves (217 total), with 32 splits and 76 pass-throughs, producing 109 clusters. There were 216 nodes reached only when enabled, and the leaf partitions differed.

The p-values among reached nodes have a different observed distribution because the reached set changes: the disabled run has one reached node with a sibling p-value (median 0.0269), versus 182 in the enabled run (median 0.4461). This is conditional selection over fixed per-node values, not a change to those values. The real matrix has no generator truth, so this probe cannot label its additional splits as false splits or count signal recovery.

Artifacts retain exact per-node before/after fields, per-mode reachability and decisions, source/input SHA-256 hashes, and configuration:

- [`pvalue_probe.py`](pvalue_probe.py)
- [`production_probe.json`](production_probe.json)
- [`production_node_pvalues.csv`](production_node_pvalues.csv)

## Relation to controlled passthrough results

The matched known-branch Gaussian benchmark also reuses the same exact known-covariance local sibling tests under both flag values ([benchmark protocol](../passthrough_benchmark_20260928/README.md#L3-L5)). Its synthetic truth shows the clustering consequence: passthrough recovers deep signals blocked by a null root, while global-null datasets have eight additional reached-false-split cases per 1,000 under each covariance profile ([outcomes](../passthrough_benchmark_20260928/README.md#L7-L18)). That is added reachability and opportunity for downstream rejection, not altered local p-values. The independent projection study separately shows that same-data PCA can distort null p-values with rank fixed, so projection/calibration concerns must not be attributed to passthrough ([projection control](../independent_projection_vandermonde_20260928/README.md#L30-L52)).
