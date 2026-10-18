"""Per-method execution helper for benchmark pipeline runs."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from benchmarks.shared.benchmark_grid import (
    benchmark_param_metadata,
    strip_benchmark_metadata,
)
from benchmarks.shared.metrics import calculate_clustering_metrics, nan_clustering_metrics
from benchmarks.shared.result_records import (
    BenchmarkResultRow,
    ComputedResultRecord,
    build_benchmark_result_row,
    build_computed_result_record,
)
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.types import MethodSpec
from benchmarks.shared.util.decomposition import create_report_dataframe_from_labels
from benchmarks.shared.util.execution_mode import coerce_bool_param
from benchmarks.shared.util.method_sets import TBS_DISTANCE_TREE_METHODS
from benchmarks.shared.util.time import BENCHMARK_STAGE_TIMING_KEYS
from tree_break_selection.space_separation import FEATURE_FAMILY_DIFFUSION_TREE_DISTANCE_METRIC
from tree_break_selection.tree.construction import DEFAULT_BINARY_TREE_DISTANCE_METRIC
from tree_break_selection.tree.continuous_distance import (
    CONTINUOUS_STANDARDIZED_EUCLIDEAN_TREE_DISTANCE_METRIC,
    CONTINUOUS_TREE_DISTANCE_METRIC,
    continuous_time_distance_condensed,
    standardized_euclidean_distance_condensed,
)
from tree_break_selection.tree.optimized_branch_lengths import (
    BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
    BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC,
    BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
)

TBS_TREE_DISTANCE_SOURCE_KEY = "tree_distance_source"
TBS_TREE_DISTANCE_SOURCE_FEATURE_METRIC = "feature_metric"
TBS_TREE_DISTANCE_SOURCE_PRECOMPUTED = "precomputed"


def require_precomputed_tbs_distance_metric(
    *,
    meta: dict[str, object],
    case_name: str,
) -> str:
    """Return the named metric/source for a required precomputed TBS distance."""
    if "distance_metric" not in meta:
        raise ValueError(
            f"Case '{case_name}' requires precomputed TBS distance metadata but "
            "does not define 'distance_metric'."
        )
    distance_metric = str(meta["distance_metric"])
    if not distance_metric:
        raise ValueError(
            f"Case '{case_name}' requires non-empty precomputed TBS distance metadata."
        )
    return distance_metric


def _report_for_metric_evaluation(
    report_df: pd.DataFrame | None,
    labels: np.ndarray,
    sample_index: pd.Index,
) -> pd.DataFrame:
    """Return a report table aligned to the sample index used for metrics."""
    if report_df is None:
        return create_report_dataframe_from_labels(labels, sample_index)
    if "cluster_id" not in report_df.columns:
        raise ValueError("Runner report_df must include a 'cluster_id' column.")
    if len(report_df.index) != len(sample_index):
        raise ValueError(
            "Runner report_df must have one row per sample. "
            f"Got {len(report_df.index)} rows for {len(sample_index)} samples."
        )
    if report_df.index.has_duplicates:
        raise ValueError("Runner report_df index contains duplicate sample ids.")
    if sample_index.has_duplicates:
        raise ValueError("Sample index contains duplicate sample ids.")
    if report_df.index.equals(sample_index):
        return report_df

    missing = sample_index.difference(report_df.index)
    extras = report_df.index.difference(sample_index)
    if missing.empty and extras.empty:
        aligned_report = report_df.loc[sample_index].copy()
        aligned_report.index.name = "sample_id"
        return aligned_report
    raise ValueError(
        "Runner report_df index must match sample ids. "
        f"Missing={list(missing[:5])}, extra={list(extras[:5])}."
    )


def _extract_stage_timings(
    *,
    method_id: str,
    result_status: str,
    result_extra: dict | None,
) -> dict[str, object] | None:
    """Return normalized stage timings for TBS-family benchmark rows."""
    stage_timings = (
        result_extra.get("stage_timings")
        if result_extra and isinstance(result_extra.get("stage_timings"), dict)
        else None
    )
    if method_id.startswith("tbs") and result_status == "ok":
        if stage_timings is None:
            raise ValueError("Successful TBS-family method results must include stage_timings.")
        missing = sorted(set(BENCHMARK_STAGE_TIMING_KEYS) - set(stage_timings))
        if missing:
            raise ValueError(
                "Successful TBS-family method results must include all stage_timings; "
                f"missing={missing}."
            )
    return stage_timings


def _should_fail_closed_hard_overlap_internal_filter(
    *,
    method_id: str,
    case_name: str,
    true_clusters: int,
    result: object,
    recorded_run_params: dict[str, object],
) -> bool:
    """Return whether a guarded internal-filter hard-overlap OK row should skip."""
    if not str(case_name).startswith("overlap_extreme_4c"):
        return False
    if not str(method_id).startswith("tbs"):
        return False
    if not coerce_bool_param(
        recorded_run_params.get("enforce_internal_support_thresholds", False),
        name="enforce_internal_support_thresholds",
    ):
        return False
    if not coerce_bool_param(
        recorded_run_params.get("spectral_include_internal_barycenters", False),
        name="spectral_include_internal_barycenters",
    ):
        return False
    if true_clusters <= 1:
        return False
    return (
        getattr(result, "status", None) == "ok"
        and getattr(result, "labels", None) is not None
        and int(getattr(result, "found_clusters", 0)) <= 1
    )


def run_single_method_once(
    *,
    method_id: str,
    spec: MethodSpec,
    params: dict[str, object],
    case_idx: int,
    case_name: str,
    tc_seed: object,
    significance_level: float,
    edge_alpha: float,
    data_t: object,
    y_t: object,
    x_original: object,
    meta: dict[str, object],
    distance_matrix: np.ndarray | None,
    distance_condensed: np.ndarray | None,
) -> tuple[BenchmarkResultRow, ComputedResultRecord | None]:
    """Execute one method+params run and return typed outputs."""
    benchmark_metadata = benchmark_param_metadata(
        method_id,
        params,
        default_class=spec.benchmark_class,
        default_grid=spec.benchmark_grid,
    )
    run_params = strip_benchmark_metadata(params)
    if method_id in {"kmeans", "spectral"}:
        raw_k = run_params["n_clusters"]
        if str(raw_k).strip().lower() in {"true", "expected", "auto"}:
            run_params["n_clusters"] = int(meta["n_clusters"])

    run_seed = tc_seed
    repeat = int(benchmark_metadata["benchmark_repeat"])
    if repeat and tc_seed is not None:
        run_seed = int(tc_seed) + repeat

    meta_run = meta.copy()
    feature_space = meta.get("feature_space")
    distance_condensed_for_run = None
    recorded_run_params = dict(run_params)
    recorded_run_params["benchmark_seed"] = run_seed
    if method_id.startswith("tbs"):
        recorded_run_params["edge_alpha"] = float(edge_alpha)
        recorded_run_params["sibling_alpha"] = float(significance_level)
        branch_length_method = str(
            recorded_run_params.get(
                "branch_length_optimization_method",
                BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC,
            )
        )
        recorded_run_params["branch_length_optimization_method"] = branch_length_method
        if branch_length_method == BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS:
            recorded_run_params.setdefault(
                "branch_length_optimization_target_metric",
                BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
            )
            recorded_run_params.setdefault("branch_length_optimization_pair_sample_size", 100_000)
            recorded_run_params.setdefault("branch_length_optimization_random_state", 0)
            recorded_run_params.setdefault(
                "branch_length_optimization_apply_nonconverged",
                False,
            )
    if method_id in TBS_DISTANCE_TREE_METHODS:
        metric = str(run_params["tree_distance_metric"])
        requires_precomputed_tbs_distance = bool(meta["requires_precomputed_tbs_distance"])
        if metric == FEATURE_FAMILY_DIFFUSION_TREE_DISTANCE_METRIC:
            # The method builds its own geometry; the case distance is only its fallback.
            distance_condensed_for_run = distance_condensed
            recorded_run_params["tree_distance_metric"] = metric
            recorded_run_params[TBS_TREE_DISTANCE_SOURCE_KEY] = (
                TBS_TREE_DISTANCE_SOURCE_FEATURE_METRIC
            )
        elif requires_precomputed_tbs_distance:
            if distance_condensed is None:
                raise ValueError(
                    f"Case '{meta['name']}' requires "
                    "'precomputed_distance_condensed' for TBS but it is missing."
                )
            distance_condensed_for_run = distance_condensed
            recorded_run_params["tree_distance_metric"] = require_precomputed_tbs_distance_metric(
                meta=meta,
                case_name=str(meta["name"]),
            )
            recorded_run_params[TBS_TREE_DISTANCE_SOURCE_KEY] = TBS_TREE_DISTANCE_SOURCE_PRECOMPUTED
        else:
            if metric == CONTINUOUS_TREE_DISTANCE_METRIC:
                if feature_space is None:
                    raise ValueError(
                        "mahalanobis_time TBS tree distances require feature_space metadata."
                    )
                distance_condensed_for_run = continuous_time_distance_condensed(
                    data_t.values,
                    feature_space,
                )
            elif metric == CONTINUOUS_STANDARDIZED_EUCLIDEAN_TREE_DISTANCE_METRIC:
                if feature_space is None:
                    raise ValueError(
                        "standardized_euclidean TBS tree distances require feature_space metadata."
                    )
                distance_condensed_for_run = standardized_euclidean_distance_condensed(
                    data_t.values,
                    feature_space,
                )
            elif metric == DEFAULT_BINARY_TREE_DISTANCE_METRIC and distance_condensed is not None:
                distance_condensed_for_run = distance_condensed
            else:
                distance_condensed_for_run = pdist(data_t.values, metric=metric)
            recorded_run_params["tree_distance_metric"] = metric
            recorded_run_params[TBS_TREE_DISTANCE_SOURCE_KEY] = (
                TBS_TREE_DISTANCE_SOURCE_FEATURE_METRIC
            )

    result = run_clustering_result(
        data_df=data_t,
        method_id=method_id,
        params=run_params,
        seed=run_seed,
        significance_level=significance_level,
        edge_alpha=edge_alpha,
        distance_matrix=distance_matrix,
        distance_condensed=distance_condensed_for_run,
        feature_space=feature_space,
    )

    true_clusters_raw = meta["n_clusters"]
    true_clusters = int(true_clusters_raw)
    stage_timings = _extract_stage_timings(
        method_id=method_id,
        result_status=result.status,
        result_extra=result.extra,
    )

    if _should_fail_closed_hard_overlap_internal_filter(
        method_id=method_id,
        case_name=case_name,
        true_clusters=true_clusters,
        result=result,
        recorded_run_params=recorded_run_params,
    ):
        raise ValueError(
            "Fail-closed hard-overlap internal support guard: "
            "guarded internal-barycenter TBS returned one cluster without "
            "admissible split support."
        )

    if result.status == "ok" and result.labels is not None:
        labels = result.labels
        report_df = _report_for_metric_evaluation(result.report_df, labels, data_t.index)
        found_clusters = result.found_clusters
        labels_len = len(labels)
        metrics = calculate_clustering_metrics(
            report_df,
            data_t.index,
            y_t,
            meta,
            feature_matrix=data_t,
        )
    else:
        labels_len = 0
        found_clusters = 0
        metrics = nan_clustering_metrics()

    if result.status == "ok":
        cluster_count_abs_error = float(abs(int(found_clusters) - true_clusters))
        over_split = float(int(found_clusters > true_clusters))
        under_split = float(int(found_clusters < true_clusters))
    else:
        cluster_count_abs_error = np.nan
        over_split = np.nan
        under_split = np.nan

    result_row = build_benchmark_result_row(
        test_case=case_idx,
        case_id=case_name,
        case_category=meta["category"],
        source_family=meta["source_family"],
        feature_representation=meta["feature_representation"],
        simulation_model=meta["simulation_model"],
        observation_model=meta["observation_model"],
        benchmark_intent=meta["benchmark_intent"],
        scientific_caution=meta["scientific_caution"],
        recommended_simulation_family=meta["recommended_simulation_family"],
        method=method_id,
        run_params=recorded_run_params,
        run_id=str(benchmark_metadata["run_id"]),
        benchmark_class=str(benchmark_metadata["benchmark_class"]),
        benchmark_grid=str(benchmark_metadata["benchmark_grid"]),
        benchmark_repeat=int(benchmark_metadata["benchmark_repeat"]),
        true_clusters=true_clusters,
        found_clusters=found_clusters,
        samples=meta["n_samples"],
        features=meta["n_features"],
        noise=meta["noise"],
        metrics=metrics,
        cluster_count_abs_error=cluster_count_abs_error,
        over_split=over_split,
        under_split=under_split,
        status=result.status,
        skip_reason=result.skip_reason,
        unsupported_reason=result.unsupported_reason,
        labels_length=labels_len,
        stage_timings=stage_timings,
    )

    computed_result = None
    if result.status == "ok" and result.labels is not None:
        meta_run["found_clusters"] = found_clusters
        if stage_timings is not None:
            meta_run["stage_timings"] = dict(stage_timings)
        computed_result = build_computed_result_record(
            test_case_num=case_idx,
            method=method_id,
            method_name=spec.name,
            run_id=str(benchmark_metadata["run_id"]),
            benchmark_class=str(benchmark_metadata["benchmark_class"]),
            benchmark_grid=str(benchmark_metadata["benchmark_grid"]),
            benchmark_repeat=int(benchmark_metadata["benchmark_repeat"]),
            params=recorded_run_params,
            ari=metrics.ari,
            nmi=metrics.nmi,
            purity=metrics.purity,
            outlier_precision=metrics.outlier_precision,
            outlier_recall=metrics.outlier_recall,
            outlier_f1=metrics.outlier_f1,
            singleton_outlier_isolated=metrics.singleton_outlier_isolated,
            grouped_outlier_cluster_recovered=metrics.grouped_outlier_cluster_recovered,
            labels=result.labels,
            data=data_t,
            meta=meta_run,
            x_original=x_original,
            y_true=y_t,
            tree=result.extra.get("tree") if result.extra else None,
            decomposition=result.extra.get("decomposition") if result.extra else None,
            annotations=result.extra.get("annotations") if result.extra else None,
        )

    return result_row, computed_result


__all__ = ["require_precomputed_tbs_distance_metric", "run_single_method_once"]
