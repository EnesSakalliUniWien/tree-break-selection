"""Decomposition helpers for benchmark runners and pipeline."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from benchmarks.shared.types.method_run_result import MethodRunResult
from tree_break_selection.hierarchy_analysis.cluster_assignments import (
    build_sample_cluster_assignments,
)
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.construction import build_tree


def _empty_report_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["cluster_id", "cluster_size"],
        index=pd.Index([], name="sample_id"),
    )


def _report_dataframe_from_assignments(assignments_table: pd.DataFrame) -> pd.DataFrame:
    if assignments_table.empty:
        return _empty_report_dataframe()
    return assignments_table[["cluster_id", "cluster_size"]].copy()


def _labels_from_assignment_table(
    assignments_table: pd.DataFrame,
    sample_index: Sequence[object] | pd.Index,
) -> np.ndarray:
    ordered_index = _sample_index(sample_index).rename("sample_id")
    if assignments_table.empty:
        return np.full(len(ordered_index), -1, dtype=int)
    return assignments_table["cluster_id"].reindex(ordered_index).fillna(-1).astype(int).to_numpy()


def _sample_index(sample_names: Sequence[object] | pd.Index) -> pd.Index:
    return sample_names if isinstance(sample_names, pd.Index) else pd.Index(sample_names)


def labels_and_report_from_decomposition(
    decomposition: dict,
    sample_index: Sequence[object] | pd.Index,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Convert decomposition output once into ordered labels and a report table."""
    assignments_table = build_sample_cluster_assignments(decomposition)
    labels = _labels_from_assignment_table(assignments_table, sample_index)
    report_df = _report_dataframe_from_assignments(assignments_table)
    return labels, report_df


def tbs_decompose(
    data: pd.DataFrame,
    distance_condensed: np.ndarray,
    *,
    linkage_method: str,
    alpha: float,
) -> tuple[np.ndarray, dict]:
    """Build the TBS tree, run both gates at ``alpha`` and decompose it.

    Returns labels in ``data`` row order (``-1`` for unassigned samples) and the
    decomposition result.
    """
    tree = build_tree(data, distance_condensed, linkage_method=linkage_method).tree
    tree.populate_node_distributions(data)
    gate_bundle = run_gate_annotation_pipeline(
        tree,
        tree.annotations_df.copy(),
        leaf_data=data,
        edge_alpha=alpha,
        sibling_alpha=alpha,
    )
    decomposition = TreeDecomposition(
        tree=tree,
        gate_annotation_bundle=gate_bundle,
    ).decompose_tree()
    labels, _report = labels_and_report_from_decomposition(decomposition, data.index)
    return labels, decomposition


def create_report_dataframe_from_labels(
    labels: list[int] | np.ndarray, sample_names: Sequence[object] | pd.Index
) -> pd.DataFrame:
    """Create a report dataframe from flat cluster labels."""
    labels_arr = np.asarray(labels, dtype=int)
    if labels_arr.size == 0:
        return _empty_report_dataframe()
    series = pd.Series(labels_arr, index=_sample_index(sample_names), name="cluster_id")
    sizes = series.value_counts()
    report = pd.DataFrame(
        {
            "cluster_id": series,
            "cluster_size": series.map(sizes),
        }
    )
    report.index.name = "sample_id"
    return report


def _count_non_noise_clusters(labels: np.ndarray) -> int:
    return int(len({int(label) for label in labels if int(label) >= 0}))


def ok_result_from_labels(
    labels: list[int] | np.ndarray,
    sample_names: Sequence[object] | pd.Index,
    *,
    extra: dict | None = None,
) -> MethodRunResult:
    """Build a successful runner result from flat labels."""
    labels_array = np.asarray(labels, dtype=int)
    return MethodRunResult(
        labels=labels_array,
        found_clusters=_count_non_noise_clusters(labels_array),
        report_df=create_report_dataframe_from_labels(labels_array, sample_names),
        status="ok",
        skip_reason=None,
        extra=extra,
    )


__all__ = [
    "tbs_decompose",
    "labels_and_report_from_decomposition",
    "create_report_dataframe_from_labels",
    "ok_result_from_labels",
]
