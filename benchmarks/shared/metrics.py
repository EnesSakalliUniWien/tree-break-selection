"""Metric helpers for benchmarking pipeline.

Contains helpers for computing clustering metrics and outlier-isolation metrics.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import entropy
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    f1_score,
    fowlkes_mallows_score,
    homogeneity_completeness_v_measure,
    normalized_mutual_info_score,
    precision_score,
    recall_score,
    silhouette_score,
)
from sklearn.metrics.cluster import contingency_matrix
from sklearn.preprocessing import LabelEncoder


@dataclass(frozen=True)
class ClusteringMetrics:
    """Clustering metrics computed for a single method run."""

    ari: float
    nmi: float
    ami: float
    purity: float
    homogeneity: float
    completeness: float
    v_measure: float
    fowlkes_mallows: float
    macro_recall: float
    macro_f1: float
    worst_cluster_recall: float
    n_singleton_clusters: float
    singleton_fraction: float
    median_cluster_size: float
    largest_cluster_fraction: float
    effective_cluster_count: float
    cluster_size_entropy: float
    cluster_size_gini: float
    noise_label_fraction: float
    silhouette_score: float
    davies_bouldin_index: float
    calinski_harabasz_index: float
    outlier_precision: float
    outlier_recall: float
    outlier_f1: float
    singleton_outlier_isolated: float
    grouped_outlier_cluster_recovered: float
    matched_macro_f1: float = np.nan
    matched_rare_f1: float = np.nan
    pair_split_share: float = np.nan
    pair_merge_share: float = np.nan


def nan_clustering_metrics() -> ClusteringMetrics:
    """Return the all-NaN metrics recorded for runs without usable labels."""
    return ClusteringMetrics(**{field.name: np.nan for field in fields(ClusteringMetrics)})


def _unit_interval_score(value: float) -> float:
    """Normalize bounded scores against tiny floating-point overshoots."""
    score = float(value)
    if not np.isfinite(score):
        return score
    return float(np.clip(score, 0.0, 1.0))


def _is_noise_label(value: object) -> bool:
    if pd.isna(value):
        return False
    try:
        return bool(float(value) == -1.0)
    except (TypeError, ValueError):
        return str(value).strip() == "-1"


def _effective_cluster_count_from_sizes(sizes: pd.Series) -> tuple[float, float]:
    if sizes.empty or float(sizes.sum()) <= 0.0:
        return np.nan, np.nan
    size_entropy = float(entropy(sizes.to_numpy(dtype=float)))
    return float(np.exp(size_entropy)), size_entropy


def _cluster_size_gini(sizes: pd.Series) -> float:
    values = np.sort(sizes.to_numpy(dtype=float))
    if values.size == 0 or np.isclose(values.sum(), 0.0):
        return np.nan
    index = np.arange(1, values.size + 1)
    numerator = ((2.0 * index - values.size - 1.0) * values).sum()
    denominator = values.size * values.sum()
    return float(numerator / denominator)


def _cluster_purity(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Return majority-label purity across predicted clusters."""
    matrix = contingency_matrix(y_true, y_pred)
    if matrix.size == 0 or matrix.sum() <= 0:
        return np.nan
    return float(np.max(matrix, axis=0).sum() / matrix.sum())


def _calculate_fragmentation_metrics(
    pred_cluster: pd.Series,
) -> tuple[float, float, float, float, float, float, float, float]:
    """Return cluster-size diagnostics excluding DBSCAN-style noise label -1."""
    if pred_cluster.empty or pred_cluster.isna().any():
        return (np.nan,) * 8

    noise_mask = pred_cluster.map(_is_noise_label).to_numpy(dtype=bool)
    noise_label_fraction = float(noise_mask.mean()) if noise_mask.size else np.nan
    non_noise = pred_cluster.loc[~noise_mask]
    if non_noise.empty:
        return (
            0.0,
            0.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            noise_label_fraction,
        )

    sizes = non_noise.value_counts()
    n_samples = float(len(pred_cluster))
    effective_count, size_entropy = _effective_cluster_count_from_sizes(sizes)
    return (
        float((sizes == 1).sum()),
        float(sizes[sizes == 1].sum() / n_samples),
        float(sizes.median()),
        float(sizes.max() / n_samples),
        effective_count,
        size_entropy,
        _cluster_size_gini(sizes),
        noise_label_fraction,
    )


def _calculate_internal_geometry_metrics(
    feature_matrix: np.ndarray | pd.DataFrame | None,
    y_pred: np.ndarray,
) -> tuple[float, float, float]:
    """Return silhouette, Davies-Bouldin, and Calinski-Harabasz when defined."""
    if feature_matrix is None:
        return np.nan, np.nan, np.nan
    values = np.asarray(feature_matrix, dtype=float)
    labels = np.asarray(y_pred)
    if values.ndim != 2 or values.shape[0] != labels.shape[0] or values.shape[0] == 0:
        return np.nan, np.nan, np.nan
    if not np.isfinite(values).all():
        return np.nan, np.nan, np.nan

    non_noise = np.array([not _is_noise_label(label) for label in labels], dtype=bool)
    values = values[non_noise]
    labels = labels[non_noise]
    n_samples = labels.shape[0]
    n_labels = np.unique(labels).size
    if not 1 < n_labels < n_samples:
        return np.nan, np.nan, np.nan

    silhouette = np.nan
    davies_bouldin = np.nan
    calinski_harabasz = np.nan
    try:
        sample_size = min(2000, n_samples)
        silhouette = float(
            silhouette_score(
                values,
                labels,
                sample_size=sample_size if sample_size < n_samples else None,
                random_state=0,
            )
        )
    except ValueError:
        silhouette = np.nan
    try:
        davies_bouldin = float(davies_bouldin_score(values, labels))
    except ValueError:
        davies_bouldin = np.nan
    try:
        calinski_harabasz = float(calinski_harabasz_score(values, labels))
    except ValueError:
        calinski_harabasz = np.nan
    return silhouette, davies_bouldin, calinski_harabasz


def _calculate_outlier_isolation_metrics(
    report_df: pd.DataFrame,
    sample_names: pd.Index,
    outlier_indices: list[int] | tuple[int, ...] | np.ndarray,
) -> tuple[float, float, float, float, float]:
    """Compute outlier-isolation metrics from predicted clusters.

    A predicted cluster counts as an outlier cluster when a strict majority of its
    members are true outliers. This makes the metric meaningful even when methods
    do not emit a dedicated noise label.
    """

    if len(outlier_indices) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    sample_index = pd.Index(sample_names)
    if sample_index.empty:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    outlier_mask_by_position = np.zeros(len(sample_index), dtype=bool)
    outlier_positions = np.asarray(outlier_indices, dtype=int)
    if np.any(outlier_positions < 0) or np.any(outlier_positions >= len(sample_index)):
        raise ValueError("outlier_indices contains positions outside sample_names")
    outlier_mask_by_position[outlier_positions] = True

    outlier_lookup = pd.Series(outlier_mask_by_position, index=sample_index)
    true_outlier = report_df.index.to_series().map(outlier_lookup)
    pred_cluster = report_df["cluster_id"]

    if true_outlier.isna().any() or pred_cluster.isna().any():
        return np.nan, np.nan, np.nan, np.nan, np.nan

    cluster_outlier_fraction = true_outlier.groupby(pred_cluster).mean()
    predicted_outlier_clusters = cluster_outlier_fraction[cluster_outlier_fraction > 0.5].index
    predicted_outlier = pred_cluster.isin(predicted_outlier_clusters).to_numpy(dtype=bool)
    true_outlier_arr = true_outlier.to_numpy(dtype=bool)

    outlier_precision = float(precision_score(true_outlier_arr, predicted_outlier, zero_division=0))
    outlier_recall = float(recall_score(true_outlier_arr, predicted_outlier, zero_division=0))
    outlier_f1 = float(f1_score(true_outlier_arr, predicted_outlier, zero_division=0))

    singleton_outlier_isolated = np.nan
    if len(outlier_positions) == 1:
        outlier_sample_name = sample_index[int(outlier_positions[0])]
        if outlier_sample_name in pred_cluster.index:
            outlier_cluster = pred_cluster.loc[outlier_sample_name]
            singleton_outlier_isolated = float((pred_cluster == outlier_cluster).sum() == 1)

    grouped_outlier_cluster_recovered = np.nan
    if len(outlier_positions) > 1:
        outlier_sample_names = sample_index[outlier_positions]
        outlier_cluster_ids = pred_cluster.loc[outlier_sample_names]
        if outlier_cluster_ids.nunique(dropna=False) == 1:
            grouped_cluster_id = outlier_cluster_ids.iloc[0]
            grouped_cluster_mask = pred_cluster == grouped_cluster_id
            grouped_outlier_cluster_recovered = float(
                true_outlier[grouped_cluster_mask].all()
                and int(grouped_cluster_mask.sum()) == len(outlier_positions)
            )
        else:
            grouped_outlier_cluster_recovered = 0.0

    return (
        outlier_precision,
        outlier_recall,
        outlier_f1,
        singleton_outlier_isolated,
        grouped_outlier_cluster_recovered,
    )


def _remap_predicted_labels_to_truth(
    y_true: np.ndarray,
    y_pred: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map predicted cluster ids to true-label ids for per-class metrics."""
    le_true = LabelEncoder()
    le_pred = LabelEncoder()
    y_true_enc = le_true.fit_transform(y_true)
    y_pred_enc = le_pred.fit_transform(y_pred)

    cm = contingency_matrix(y_true_enc, y_pred_enc)
    if cm.size == 0:
        return y_true_enc, y_pred_enc

    row_ind, col_ind = linear_sum_assignment(-cm)  # maximize agreement
    mapping_enc = {pred: true for true, pred in zip(row_ind, col_ind)}

    # Map any extra predicted clusters to their majority true class.
    for pred in range(cm.shape[1]):
        if pred not in mapping_enc:
            mapping_enc[pred] = int(np.argmax(cm[:, pred]))

    y_pred_remapped_enc = np.asarray([mapping_enc[p] for p in y_pred_enc], dtype=int)
    return y_true_enc, y_pred_remapped_enc


def _matched_recovery_from_table(
    table: pd.DataFrame, rare_labels
) -> tuple[pd.DataFrame, float, float]:
    counts = table.to_numpy(dtype=np.int64)
    truth_n, pred_n = counts.sum(axis=1), counts.sum(axis=0)
    f1 = 2 * counts / (truth_n[:, None] + pred_n)
    rows, cols = linear_sum_assignment(f1, maximize=True)
    shared = counts[rows, cols]
    keep = shared > 0
    rows, cols, shared = rows[keep], cols[keep], shared[keep]
    matched = pd.DataFrame({"group": table.index, "truth_n": truth_n})
    matched["predicted"] = pd.Series(index=matched.index, dtype=object)
    matched[["precision", "recall", "f1"]] = 0.0
    matched.loc[rows, "predicted"] = table.columns[cols].to_numpy(dtype=object)
    matched.loc[rows, "precision"] = shared / pred_n[cols]
    matched.loc[rows, "recall"] = shared / truth_n[rows]
    matched.loc[rows, "f1"] = f1[rows, cols]
    if rare_labels is None:
        rare_labels = matched.sort_values(["truth_n", "group"]).head((len(matched) + 3) // 4).group
    matched["rare"] = matched.group.isin(rare_labels)
    matched["matched_90_90"] = (matched.precision >= 0.9) & (matched.recall >= 0.9)
    biggest = counts.argmax(axis=1)
    matched["legacy_90_90"] = (counts.max(axis=1) / truth_n >= 0.9) & (
        counts.max(axis=0)[biggest] / pred_n[biggest] >= 0.9
    )
    shared_pairs = (counts * (counts - 1)).sum()
    truth_pairs = (truth_n * (truth_n - 1)).sum()
    predicted_pairs = (pred_n * (pred_n - 1)).sum()
    split = 1 - shared_pairs / truth_pairs if truth_pairs else np.nan
    merge = 1 - shared_pairs / predicted_pairs if predicted_pairs else np.nan
    return matched, float(split), float(merge)


def matched_group_recovery(true_labels, predicted_labels, *, rare_labels=None) -> pd.DataFrame:
    """One-to-one maximum-total-F1 recovery with unmatched truth groups scored zero.

    Two Series are aligned by their unique indexes. Other inputs are aligned
    by position. Rare defaults to the smallest ceil(K/4) truth groups, with
    label order breaking size ties. Pass an empty collection for no rare groups.
    """
    if isinstance(true_labels, pd.Series) and isinstance(predicted_labels, pd.Series):
        if (
            not true_labels.index.is_unique
            or not predicted_labels.index.is_unique
            or len(true_labels) != len(predicted_labels)
            or not true_labels.index.isin(predicted_labels.index).all()
        ):
            raise ValueError("Truth and prediction Series must share a unique sample index.")
        predicted_labels = predicted_labels.reindex(true_labels.index)
    truth, predicted = np.asarray(true_labels), np.asarray(predicted_labels)
    if (
        truth.ndim != 1
        or predicted.shape != truth.shape
        or not len(truth)
        or pd.isna(truth).any()
        or pd.isna(predicted).any()
    ):
        raise ValueError(
            "Expected aligned nonempty truth and predicted labels without missing values."
        )
    return _matched_recovery_from_table(pd.crosstab(truth, predicted), rare_labels)[0]


def calculate_clustering_metrics(
    report_df: pd.DataFrame | None,
    sample_names: pd.Index,
    true_labels: np.ndarray,
    metadata: dict[str, object] | None = None,
    feature_matrix: np.ndarray | pd.DataFrame | None = None,
) -> ClusteringMetrics:
    """Calculate external, fragmentation, internal, and outlier metrics."""

    if report_df is None or report_df.empty:
        return nan_clustering_metrics()
    if "cluster_id" not in report_df.columns:
        return nan_clustering_metrics()

    sample_index = pd.Index(sample_names)
    true_arr = np.asarray(true_labels)

    # Metrics are only valid when predicted labels and truth align by sample.
    if len(sample_index) != len(true_arr):
        raise ValueError(
            f"true_labels has {len(true_arr)} entries for {len(sample_index)} samples."
        )
    if len(report_df.index) != len(sample_index):
        raise ValueError(
            f"report_df has {len(report_df.index)} rows for {len(sample_index)} samples."
        )
    if sample_index.has_duplicates or report_df.index.has_duplicates:
        raise ValueError("Duplicate sample ids in sample_names or report_df index.")

    true_lookup = pd.Series(true_arr, index=sample_index)
    true_cluster = report_df.index.to_series().map(true_lookup)
    pred_cluster = report_df["cluster_id"]

    if true_cluster.isna().any() or pred_cluster.isna().any():
        raise ValueError(
            "Unmapped or missing labels: report_df ids not in sample_names, or NaN cluster_id."
        )

    y_true = true_cluster.to_numpy()
    y_pred = pred_cluster.to_numpy()
    (
        n_singleton_clusters,
        singleton_fraction,
        median_cluster_size,
        largest_cluster_fraction,
        effective_cluster_count,
        cluster_size_entropy,
        cluster_size_gini,
        noise_label_fraction,
    ) = _calculate_fragmentation_metrics(pred_cluster)
    (
        internal_silhouette,
        davies_bouldin,
        calinski_harabasz,
    ) = _calculate_internal_geometry_metrics(feature_matrix, y_pred)
    outlier_precision = np.nan
    outlier_recall = np.nan
    outlier_f1 = np.nan
    singleton_outlier_isolated = np.nan
    grouped_outlier_cluster_recovered = np.nan
    try:
        ari = adjusted_rand_score(y_true, y_pred)
        nmi = normalized_mutual_info_score(y_true, y_pred)
        ami = adjusted_mutual_info_score(y_true, y_pred)
        purity = _cluster_purity(y_true, y_pred)
        homogeneity, completeness, v_measure = homogeneity_completeness_v_measure(
            y_true,
            y_pred,
        )
        fowlkes_mallows = fowlkes_mallows_score(y_true, y_pred)
        matched, pair_split_share, pair_merge_share = _matched_recovery_from_table(
            pd.crosstab(y_true, y_pred), None if metadata is None else metadata.get("rare_labels")
        )
        y_true_enc, y_pred_remapped_enc = _remap_predicted_labels_to_truth(y_true, y_pred)
        labels = np.arange(len(np.unique(y_true_enc)), dtype=int)
        if labels.size == 0:
            macro_recall = np.nan
            macro_f1 = np.nan
            worst_cluster_recall = np.nan
        else:
            per_class_recall = recall_score(
                y_true_enc,
                y_pred_remapped_enc,
                labels=labels,
                average=None,
                zero_division=0,
            )
            macro_recall = float(np.mean(per_class_recall))
            macro_f1 = float(
                f1_score(
                    y_true_enc,
                    y_pred_remapped_enc,
                    labels=labels,
                    average="macro",
                    zero_division=0,
                )
            )
            worst_cluster_recall = float(np.min(per_class_recall))
        outlier_indices_raw = None if metadata is None else metadata.get("outlier_indices")
        if (
            isinstance(outlier_indices_raw, (list, tuple, np.ndarray))
            and len(outlier_indices_raw) > 0
        ):
            (
                outlier_precision,
                outlier_recall,
                outlier_f1,
                singleton_outlier_isolated,
                grouped_outlier_cluster_recovered,
            ) = _calculate_outlier_isolation_metrics(
                report_df,
                sample_names,
                outlier_indices_raw,
            )
    except ValueError:
        return nan_clustering_metrics()

    return ClusteringMetrics(
        ari=float(ari),
        nmi=_unit_interval_score(nmi),
        ami=float(ami),
        purity=_unit_interval_score(purity),
        homogeneity=_unit_interval_score(homogeneity),
        completeness=_unit_interval_score(completeness),
        v_measure=_unit_interval_score(v_measure),
        fowlkes_mallows=_unit_interval_score(fowlkes_mallows),
        macro_recall=_unit_interval_score(macro_recall),
        macro_f1=_unit_interval_score(macro_f1),
        worst_cluster_recall=_unit_interval_score(worst_cluster_recall),
        n_singleton_clusters=float(n_singleton_clusters),
        singleton_fraction=_unit_interval_score(singleton_fraction),
        median_cluster_size=float(median_cluster_size),
        largest_cluster_fraction=_unit_interval_score(largest_cluster_fraction),
        effective_cluster_count=float(effective_cluster_count),
        cluster_size_entropy=float(cluster_size_entropy),
        cluster_size_gini=_unit_interval_score(cluster_size_gini),
        noise_label_fraction=_unit_interval_score(noise_label_fraction),
        silhouette_score=float(internal_silhouette),
        davies_bouldin_index=float(davies_bouldin),
        calinski_harabasz_index=float(calinski_harabasz),
        outlier_precision=float(outlier_precision),
        outlier_recall=float(outlier_recall),
        outlier_f1=float(outlier_f1),
        singleton_outlier_isolated=float(singleton_outlier_isolated),
        grouped_outlier_cluster_recovered=float(grouped_outlier_cluster_recovered),
        matched_macro_f1=float(matched.f1.mean()),
        matched_rare_f1=float(matched.loc[matched.rare, "f1"].mean()),
        pair_split_share=pair_split_share,
        pair_merge_share=pair_merge_share,
    )
