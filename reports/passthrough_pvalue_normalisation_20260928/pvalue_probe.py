"""Run one real-data production gate annotation and compare traversal modes."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from applications.endotypes.pipelines.run_feature_matrix_with_umap import (
    _build_linkage_tree,
    _load_binary_matrix,
)
from pandas.testing import assert_frame_equal
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_EDGE_ALPHA,
    DEFAULT_SIBLING_ALPHA,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.construction import tree_from_linkage

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(__file__).resolve().parent
INPUT_PATH = ROOT / "data/feature_matrices/feature_matrix.tsv"
SOURCE_PATHS = (
    "applications/endotypes/pipelines/run_feature_matrix_with_umap.py",
    "tree_break_selection/hierarchy_analysis/decomposition/gates/orchestrator.py",
    "tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py",
    "tree_break_selection/hierarchy_analysis/tree_decomposition.py",
    "tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/pipeline.py",
    "tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py",
    "tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/inflation_adjusted_sibling_tests.py",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _partition(result: dict[str, object]) -> frozenset[frozenset[str]]:
    clusters = result["cluster_assignments"]
    assert isinstance(clusters, dict)
    return frozenset(frozenset(info["leaves"]) for info in clusters.values())


def _same_value(left: object, right: object) -> bool:
    if pd.isna(left) and pd.isna(right):
        return True
    value = left == right
    return bool(value) if not pd.isna(value) else False


def main() -> None:
    data = _load_binary_matrix(INPUT_PATH)
    linkage_matrix, distance_metric, linkage_method, _ = _build_linkage_tree(
        data_df=data,
        tree_method="tbs",
        tree_distance_metric="hamming",
        tree_linkage_method="average",
        diffusion_k_neighbors=15,
        diffusion_time=3,
        diffusion_components=30,
        adaptive_neighbor_k=15,
        adaptive_bandwidth_type="-1/(d+2)",
        adaptive_epsilon="median",
        adaptive_metric="hamming",
    )
    tree = tree_from_linkage(linkage_matrix, leaf_names=data.index.tolist())
    tree.populate_node_divergences(data)
    bundle = run_gate_annotation_pipeline(
        tree,
        tree.annotations_df.copy(),
        edge_alpha=DEFAULT_EDGE_ALPHA,
        sibling_alpha=DEFAULT_SIBLING_ALPHA,
        leaf_data=data,
    )

    annotations_before = bundle.annotated_df.copy(deep=True)
    pvalue_columns = [
        column
        for column in annotations_before.columns
        if "p_value" in column.lower() or "pvalue" in column.lower()
    ]
    results: dict[bool, dict[str, object]] = {}
    decomposers: dict[bool, TreeDecomposition] = {}
    snapshots: dict[bool, pd.DataFrame] = {}
    for passthrough in (False, True):
        decomposers[passthrough] = TreeDecomposition(
            tree=tree,
            gate_annotation_bundle=bundle,
            passthrough=passthrough,
            trace_level="full",
        )
        results[passthrough] = decomposers[passthrough].decompose_tree()
        snapshots[passthrough] = bundle.annotated_df.copy(deep=True)
        assert_frame_equal(annotations_before, snapshots[passthrough], check_exact=True)
        assert_frame_equal(
            annotations_before,
            decomposers[passthrough].annotations_df,
            check_exact=True,
        )

    traces = {
        flag: result["traversal_trace"] for flag, result in results.items()
    }
    visited = {
        flag: {row["node_id"] for row in trace}
        for flag, trace in traces.items()
    }
    internal_nodes = {node for node in tree.nodes if tree.out_degree(node) > 0}
    leaf_nodes = set(tree.nodes) - internal_nodes
    fields_equal = {
        column: bool(annotations_before[column].equals(snapshots[True][column]))
        for column in pvalue_columns
    }
    numeric_pvalue_columns = [
        column
        for column in pvalue_columns
        if pd.api.types.is_numeric_dtype(annotations_before[column].dtype)
    ]

    node_rows: list[dict[str, object]] = []
    for node in tree.nodes:
        row: dict[str, object] = {
            "node_id": node,
            "is_internal": node in internal_nodes,
            "reached_passthrough_false": node in visited[False],
            "reached_passthrough_true": node in visited[True],
        }
        for flag in (False, True):
            label = "true" if flag else "false"
            decision = next(
                (item["decision"] for item in traces[flag] if item["node_id"] == node),
                "not_reached",
            )
            row[f"decision_passthrough_{label}"] = decision
        for column in pvalue_columns:
            left = annotations_before.at[node, column]
            right = snapshots[True].at[node, column]
            row[f"before__{column}"] = left
            row[f"after__{column}"] = right
            row[f"equal__{column}"] = _same_value(left, right)
        node_rows.append(row)

    node_table = pd.DataFrame(node_rows)
    node_table.to_csv(OUTPUT_DIR / "production_node_pvalues.csv", index=False)

    sibling_tested = annotations_before["Sibling_Divergence_P_Value"].notna()
    corrected_tested = annotations_before["Sibling_Divergence_P_Value_Corrected"].notna()
    calibration_counts = {
        str(key): int(value)
        for key, value in annotations_before[
            "Sibling_Gate_P_Value_Calibration"
        ].value_counts(dropna=False).items()
    }
    mode_summary: dict[str, object] = {}
    for flag, result in results.items():
        trace = traces[flag]
        mode_summary[str(flag).lower()] = {
            "clusters": int(result["num_clusters"]),
            "visited_nodes_all": len(trace),
            "visited_internal_nodes": sum(row["node_id"] in internal_nodes for row in trace),
            "visited_leaf_nodes": sum(row["node_id"] in leaf_nodes for row in trace),
            "split_decisions": sum(row["decision"] == "split" for row in trace),
            "pass_through_decisions": sum(row["decision"] == "pass_through" for row in trace),
            "reached_nodes_with_sibling_p": int(
                annotations_before.loc[list(visited[flag]), "Sibling_Divergence_P_Value"]
                .notna()
                .sum()
            ),
        }

    payload = {
        "purpose": "Compare production p-value annotations and partitions with passthrough disabled/enabled on one real feature matrix.",
        "run_command": ".venv/bin/python reports/passthrough_pvalue_normalisation_20260928/pvalue_probe.py",
        "input": {
            "path": str(INPUT_PATH.relative_to(ROOT)),
            "sha256": _sha256(INPUT_PATH),
            "shape": [int(value) for value in data.shape],
            "all_binary": bool(np.isin(data.to_numpy(), [0, 1]).all()),
        },
        "source_sha256": {
            source: _sha256(ROOT / source) for source in SOURCE_PATHS
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "configuration": {
            "tree_method": "tbs",
            "tree_distance_metric": distance_metric,
            "tree_linkage_method": linkage_method,
            "node_count": len(tree.nodes),
            "edge_alpha": bundle.metadata.edge.alpha,
            "sibling_alpha": bundle.metadata.sibling.alpha,
            "sibling_gate_method": bundle.metadata.config.sibling_gate_method,
            "calibration_support_contract_version": (
                bundle.metadata.config.internal_support_policy.support_contract_version
            ),
        },
        "annotation_comparison": {
            "annotation_rows": len(annotations_before),
            "annotation_fields": len(annotations_before.columns),
            "pvalue_named_fields": pvalue_columns,
            "pvalue_named_field_count": len(pvalue_columns),
            "numeric_pvalue_field_count": len(numeric_pvalue_columns),
            "numeric_pvalue_fields": numeric_pvalue_columns,
            "all_annotation_fields_equal_after_both_traversals": bool(
                annotations_before.equals(snapshots[False])
                and annotations_before.equals(snapshots[True])
            ),
            "all_pvalue_fields_equal": fields_equal,
            "changed_pvalue_fields": [
                field for field, equal in fields_equal.items() if not equal
            ],
            "sibling_calibrated_p_nonmissing": int(sibling_tested.sum()),
            "sibling_BH_adjusted_p_nonmissing": int(corrected_tested.sum()),
            "sibling_calibration_status_counts": calibration_counts,
            "edge_raw_p_nonmissing": int(
                annotations_before["Child_Parent_Divergence_P_Value"].notna().sum()
            ),
            "edge_BH_p_nonmissing": int(
                annotations_before["Child_Parent_Divergence_P_Value_BH"].notna().sum()
            ),
        },
        "traversal": mode_summary,
        "reachability": {
            "nodes_reached_only_with_passthrough": len(visited[True] - visited[False]),
            "nodes_reached_only_without_passthrough": len(visited[False] - visited[True]),
            "partitions_equal": _partition(results[False]) == _partition(results[True]),
        },
        "node_level_csv": "production_node_pvalues.csv",
    }
    (OUTPUT_DIR / "production_probe.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
