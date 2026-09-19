"""Adaptive diffusion subspace clustering and GO-IC orchestration."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict

import numpy as np
import pandas as pd
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from scipy.cluster.hierarchy import linkage
from tree_break_selection.space_separation import (
    adaptive_spectral_blocks,
    block_adaptive_diffusion_geometry,
    coordinates_for_block,
    cosine_eigendecomposition,
    weight_feature_matrix,
)

from applications.endotypes._shared import matrix_slug, safe_name
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_ranking import (
    add_specificity_aware_rank,
    rank_go_ic_results,
)
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_terms import (
    axis_term_summary_lines,
    build_axis_loadings,
    build_term_metadata,
    component_feature_loadings,
)
from applications.endotypes.io.adaptive_diffusion_dataframe_io import (
    load_binary_matrix,
    write_dataframe,
)
from applications.endotypes.pipelines.spectral_records import spectral_block_record
from applications.endotypes.plots.go_ic_tree_summary_plots import (
    cluster_coherence,
    cluster_size_metrics,
    go_information_criterion,
    tfidf_within_cosine_quality,
)
from applications.endotypes.reports.go_ic.adaptive_diffusion_plot_exports import (
    plot_axis_terms,
    plot_combined_axis_terms,
    plot_dendrogram,
    plot_diffusion_embedding,
    plot_subspace_embedding,
    plot_tree_subtree_clusters,
)
from applications.endotypes.reports.go_ic.adaptive_diffusion_subspace_artifacts import (
    RESULT_PREFIX,
    attach_artifact_paths,
    default_output_dir,
    write_artifact_index,
    write_experiment_readme,
)
from applications.endotypes.reports.go_ic.adaptive_diffusion_subspace_reports import (
    write_ranking_plots,
    write_workflow_pdfs,
)


def assignments_from_labels(labels: np.ndarray, index: pd.Index) -> pd.DataFrame:
    labels = labels.astype(int)
    sizes = pd.Series(labels).value_counts().sort_index()
    out = pd.DataFrame({"gene": index.astype(str), "cluster_id": labels})
    out["cluster_size"] = out["cluster_id"].map(sizes).astype(int)
    return out


def run_current_kl(
    data: pd.DataFrame,
    distances: np.ndarray,
    *,
    edge_alpha: float,
    sibling_alpha: float,
) -> np.ndarray:
    result = run_tbs_on_distance(
        data,
        distances,
        sibling_alpha,
        tree_linkage_method="average",
        edge_alpha=edge_alpha,
    )
    if result.labels is None:
        raise RuntimeError(f"Current TBS returned no labels: {result.skip_reason or result.status}")
    return np.asarray(result.labels, dtype=int)


def run_analysis(args: argparse.Namespace) -> None:
    artifact_prefix = matrix_slug(args.input, args.dataset_label)
    output_dir = args.output_dir or default_output_dir(args.input, args.dataset_label)
    rankings_dir = output_dir / "rankings"
    plots_dir = output_dir / "plots"
    subspaces_dir = output_dir / "subspaces"
    rankings_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    subspaces_dir.mkdir(parents=True, exist_ok=True)

    data = load_binary_matrix(args.input)
    term_metadata = build_term_metadata(data)
    config = {
        "input": str(args.input),
        "output_dir": str(output_dir),
        "artifact_prefix": artifact_prefix,
        "dataset_label": args.dataset_label or artifact_prefix,
        "method_version": "current",
        "tree_geometry": "adaptive_diffusion_cosine_subspace",
        "edge_alpha": float(args.edge_alpha),
        "sibling_alpha": float(args.sibling_alpha),
        "max_rank": int(args.max_rank),
        "min_segment_length": int(args.min_segment_length),
        "max_segments": int(args.max_segments),
        "diffusion_k_neighbors": int(args.diffusion_k_neighbors),
        "diffusion_time": int(args.diffusion_time),
        "diffusion_components": int(args.diffusion_components),
        "adaptive_bandwidth_type": args.adaptive_bandwidth_type,
        "adaptive_epsilon": args.adaptive_epsilon,
        "adaptive_metric": args.adaptive_metric,
        "weightings": list(args.weightings),
        "block_names": list(args.block_names or []),
        "top_terms_per_axis": int(args.top_terms_per_axis),
    }
    (output_dir / "experiment_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )

    summary_rows: list[dict[str, object]] = []
    cluster_coherence_frames: list[pd.DataFrame] = []
    tfidf_quality_frames: list[pd.DataFrame] = []
    block_rows: list[dict[str, object]] = []
    spectrum_rows: list[dict[str, object]] = []
    allow_blocks = set(args.block_names or [])

    for weighting in args.weightings:
        print(f"[{weighting}] cosine eigendecomposition", flush=True)
        values = weight_feature_matrix(data, weighting)
        eigvals, eigvecs = cosine_eigendecomposition(values, args.max_rank)
        feature_loadings = component_feature_loadings(values, eigvals, eigvecs)
        total_energy = float(np.sum(eigvals))
        blocks, diagnostics = adaptive_spectral_blocks(
            eigvals,
            min_segment_length=args.min_segment_length,
            max_segments=args.max_segments,
        )
        for component, eigval in enumerate(eigvals, start=1):
            spectrum_rows.append(
                {
                    "weighting": weighting,
                    "component": int(component),
                    "eigenvalue": float(eigval),
                    "fraction_of_kept_operator_energy": float(eigval / total_energy)
                    if total_energy > 0
                    else math.nan,
                }
            )

        for block in blocks:
            if allow_blocks and block.block_name not in allow_blocks:
                continue
            subspace_id = f"{weighting}__{block.block_name}"
            subspace_safe = safe_name(subspace_id)
            print(f"[{subspace_id}] adaptive diffusion current TBS", flush=True)
            subspace_dir = subspaces_dir / safe_name(weighting) / safe_name(block.block_name)
            subspace_dir.mkdir(parents=True, exist_ok=True)
            coords = coordinates_for_block(eigvals, eigvecs, block)
            block_energy = (
                float(np.sum(eigvals[block.block_start - 1 : block.block_end]) / total_energy)
                if total_energy > 0
                else math.nan
            )
            block_record = {
                "run_id": f"current__adaptive_diffusion_cosine_subspace__{subspace_id}",
                "method_version": "current",
                "tree_geometry": "adaptive_diffusion_cosine_subspace",
                "weighting": weighting,
                **spectral_block_record(
                    block,
                    block_energy=block_energy,
                    diagnostics=diagnostics,
                ),
                "subspace_dir": str(subspace_dir),
            }
            block_rows.append(block_record)
            (subspace_dir / f"{subspace_safe}__metadata.json").write_text(
                json.dumps({**block_record, "spectral_block": asdict(block)}, indent=2),
                encoding="utf-8",
            )
            write_dataframe(
                pd.DataFrame(
                    coords,
                    index=data.index.astype(str),
                    columns=[
                        f"mode_{idx:02d}" for idx in range(block.block_start, block.block_end + 1)
                    ],
                )
                .rename_axis("gene")
                .reset_index(),
                subspace_dir / f"{subspace_safe}__subspace_coordinates.csv",
            )
            all_loadings, top_loadings = build_axis_loadings(
                term_metadata=term_metadata,
                loadings=feature_loadings,
                block_start=block.block_start,
                block_end=block.block_end,
                top_n=args.top_terms_per_axis,
            )
            write_dataframe(
                all_loadings, subspace_dir / f"{subspace_safe}__axis_term_loadings_all.csv"
            )
            write_dataframe(top_loadings, subspace_dir / f"{subspace_safe}__axis_top_terms.csv")
            plot_axis_terms(
                all_loadings,
                subspace_dir / f"{subspace_safe}__axis_terms.png",
                top_n=min(args.top_terms_per_axis, 20),
            )
            plot_combined_axis_terms(
                all_loadings,
                subspace_dir / f"{subspace_safe}__axis_terms_combined.png",
                terms_per_axis=min(8, args.top_terms_per_axis),
                max_terms=60,
            )
            axis_term_lines = axis_term_summary_lines(top_loadings)
            plot_subspace_embedding(
                coords,
                None,
                subspace_dir / f"{subspace_safe}__subspace_embedding_terms_only.png",
                subspace_id,
                axis_term_lines=axis_term_lines,
            )

            distances: np.ndarray | None = None
            linkage_path = ""
            try:
                geometry = block_adaptive_diffusion_geometry(
                    coords,
                    k_neighbors=args.diffusion_k_neighbors,
                    diffusion_time=args.diffusion_time,
                    n_components=args.diffusion_components,
                    metric=args.adaptive_metric,
                    bandwidth_type=args.adaptive_bandwidth_type,
                    epsilon=args.adaptive_epsilon,
                )
                distances = geometry.distance_condensed
                diffusion_metadata = geometry.metadata
                z = linkage(distances, method="average")
                linkage_path = str(subspace_dir / f"{subspace_safe}__linkage_matrix.csv")
                write_dataframe(
                    pd.DataFrame(z, columns=["left", "right", "distance", "count"]), linkage_path
                )
                plot_dendrogram(
                    z,
                    subspace_dir / f"{subspace_safe}__tree_dendrogram.png",
                    f"{subspace_id}: adaptive diffusion average-linkage tree",
                )
            except Exception as exc:
                diffusion_metadata = {}
                status = "failed_diffusion"
                error = repr(exc)
                labels = None
            else:
                try:
                    labels = run_current_kl(
                        data,
                        distances,
                        edge_alpha=args.edge_alpha,
                        sibling_alpha=args.sibling_alpha,
                    )
                    status = "ok"
                    error = ""
                except Exception as exc:
                    status = "failed_gate"
                    error = repr(exc)
                    labels = None

            (subspace_dir / f"{subspace_safe}__diffusion_metadata.json").write_text(
                json.dumps(diffusion_metadata, indent=2, default=str),
                encoding="utf-8",
            )

            if status == "ok":
                assert labels is not None
                assert distances is not None
                assignments = assignments_from_labels(labels, data.index)
                assignments_path = subspace_dir / f"{subspace_safe}__cluster_assignments.csv"
                write_dataframe(assignments, assignments_path)
                write_dataframe(
                    assignments["cluster_id"]
                    .value_counts()
                    .sort_index()
                    .rename_axis("cluster_id")
                    .reset_index(name="cluster_size"),
                    subspace_dir / f"{subspace_safe}__cluster_sizes.csv",
                )
                plot_tree_subtree_clusters(
                    z,
                    assignments,
                    subspace_dir / f"{subspace_safe}__tree_subtree_clusters.png",
                    f"{subspace_id}: adaptive diffusion tree with current TBS cluster assignments",
                )
                plot_subspace_embedding(
                    coords,
                    labels,
                    subspace_dir / f"{subspace_safe}__subspace_embedding_clusters.png",
                    subspace_id,
                )
                plot_subspace_embedding(
                    coords,
                    labels,
                    subspace_dir
                    / f"{subspace_safe}__subspace_embedding_clusters_annotated_terms.png",
                    subspace_id,
                    axis_term_lines=axis_term_lines,
                )
                plot_diffusion_embedding(
                    distances,
                    labels,
                    subspace_dir / f"{subspace_safe}__adaptive_diffusion_embedding_clusters.png",
                    subspace_id,
                )
                plot_diffusion_embedding(
                    distances,
                    labels,
                    subspace_dir
                    / f"{subspace_safe}__adaptive_diffusion_embedding_clusters_annotated_terms.png",
                    subspace_id,
                    axis_term_lines=axis_term_lines,
                )
                info = go_information_criterion(data, labels)
                coherence = cluster_coherence(data, labels)
                tfidf_quality = tfidf_within_cosine_quality(data, labels)
                sizes = cluster_size_metrics(labels)
                coherence_frame = coherence["coherence_table"].copy()
                coherence_frame.insert(0, "run_id", block_record["run_id"])
                write_dataframe(
                    coherence_frame, subspace_dir / f"{subspace_safe}__cluster_coherence.csv"
                )
                cluster_coherence_frames.append(coherence_frame)
                tfidf_frame = tfidf_quality["tfidf_quality_table"].copy()
                tfidf_frame.insert(0, "run_id", block_record["run_id"])
                write_dataframe(
                    tfidf_frame, subspace_dir / f"{subspace_safe}__tfidf_cluster_quality.csv"
                )
                tfidf_quality_frames.append(tfidf_frame)
                rank_fields = {
                    **info,
                    **sizes,
                    "coherent_cluster_count": coherence["coherent_cluster_count"],
                    "coherent_cluster_fraction": coherence["coherent_cluster_fraction"],
                    "median_significant_terms_q05": coherence["median_significant_terms_q05"],
                    "min_cluster_q_value": coherence["min_cluster_q_value"],
                    "median_within_tfidf_cosine": tfidf_quality["median_within_tfidf_cosine"],
                    "weighted_mean_within_tfidf_cosine": tfidf_quality[
                        "weighted_mean_within_tfidf_cosine"
                    ],
                }
            else:
                assignments_path = None
                write_dataframe(
                    pd.DataFrame(
                        [
                            {
                                "status": status,
                                "error": error,
                                "note": "No cluster assignments were produced because the current TBS gate did not complete for this subspace.",
                            }
                        ]
                    ),
                    subspace_dir / f"{subspace_safe}__failure_status.csv",
                )
                rank_fields = {
                    "go_log_likelihood": math.nan,
                    "go_bic_active": math.nan,
                    "go_aic_active": math.nan,
                    "go_bic_full": math.nan,
                    "go_active_parameters": pd.NA,
                    "go_full_parameters": pd.NA,
                    "go_bic_active_per_gene": math.nan,
                    "n_clusters": pd.NA,
                    "largest_cluster_size": pd.NA,
                    "largest_cluster_fraction": math.nan,
                    "singleton_clusters": pd.NA,
                    "singleton_fraction": math.nan,
                    "singleton_gene_fraction": math.nan,
                    "median_cluster_size": math.nan,
                    "coherent_cluster_count": pd.NA,
                    "coherent_cluster_fraction": math.nan,
                    "median_significant_terms_q05": math.nan,
                    "min_cluster_q_value": math.nan,
                    "median_within_tfidf_cosine": math.nan,
                    "weighted_mean_within_tfidf_cosine": math.nan,
                }

            write_dataframe(
                pd.DataFrame([rank_fields]),
                subspace_dir / f"{subspace_safe}__go_ic_quality_summary.csv",
            )
            summary_rows.append(
                {
                    **block_record,
                    "status": status,
                    "error": error,
                    "assignments_path": str(assignments_path)
                    if assignments_path is not None
                    else "",
                    "linkage_path": linkage_path,
                    "axis_top_terms_path": str(
                        subspace_dir / f"{subspace_safe}__axis_top_terms.csv"
                    ),
                    **rank_fields,
                }
            )

    ranking = rank_go_ic_results(summary_rows)
    coherence_long = (
        pd.concat(cluster_coherence_frames, ignore_index=True) if cluster_coherence_frames else None
    )
    ranking = attach_artifact_paths(ranking)
    ranking = add_specificity_aware_rank(ranking, coherence_long)
    ranking_path = rankings_dir / f"{RESULT_PREFIX}_ranking.csv"
    write_dataframe(ranking, ranking_path)
    write_dataframe(ranking, rankings_dir / f"{RESULT_PREFIX}_specificity_aware_ranking.csv")
    write_dataframe(pd.DataFrame(block_rows), rankings_dir / f"{RESULT_PREFIX}_subspace_blocks.csv")
    write_dataframe(pd.DataFrame(spectrum_rows), rankings_dir / f"{RESULT_PREFIX}_spectrum.csv")
    if coherence_long is not None:
        write_dataframe(
            coherence_long, rankings_dir / f"{RESULT_PREFIX}_cluster_coherence_long.csv"
        )
    if tfidf_quality_frames:
        write_dataframe(
            pd.concat(tfidf_quality_frames, ignore_index=True),
            rankings_dir / f"{RESULT_PREFIX}_tfidf_quality_long.csv",
        )

    write_ranking_plots(ranking, output_dir)

    workflow_outputs = write_workflow_pdfs(output_dir, ranking, artifact_prefix)
    write_artifact_index(output_dir, ranking, artifact_prefix, workflow_outputs)

    write_experiment_readme(args.input, output_dir, artifact_prefix, ranking)
    print(f"Wrote current adaptive diffusion subspace tree experiment: {output_dir}", flush=True)
