"""Adaptive diffusion subspace artifact export."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

from applications.endotypes._shared import matrix_slug, safe_name
from applications.endotypes.io.adaptive_diffusion_dataframe_io import write_dataframe

# Retain the artifact naming contract used by report readers and saved runs.
RESULT_PREFIX = "current_adaptive_diffusion_subspace_tree"


def attach_artifact_paths(ranking: pd.DataFrame) -> pd.DataFrame:
    out = ranking.copy()
    rows: list[dict[str, object]] = []
    for row in out.itertuples(index=False):
        subspace_id = f"{row.weighting}__{row.block_name}"
        subspace_safe = safe_name(subspace_id)
        subspace_dir = Path(row.subspace_dir)
        rows.append(
            {
                "cluster_assignments": getattr(row, "assignments_path", ""),
                "failure_status": str(subspace_dir / f"{subspace_safe}__failure_status.csv")
                if getattr(row, "status", "") != "ok"
                else "",
                "linkage_matrix": getattr(row, "linkage_path", ""),
                "tree_dendrogram_png": str(subspace_dir / f"{subspace_safe}__tree_dendrogram.png"),
                "tree_subtree_clusters_png": str(
                    subspace_dir / f"{subspace_safe}__tree_subtree_clusters.png"
                ),
                "subspace_embedding_png": str(
                    subspace_dir / f"{subspace_safe}__subspace_embedding_clusters.png"
                ),
                "adaptive_diffusion_embedding_png": str(
                    subspace_dir / f"{subspace_safe}__adaptive_diffusion_embedding_clusters.png"
                ),
                "axis_terms_combined_png": str(
                    subspace_dir / f"{subspace_safe}__axis_terms_combined.png"
                ),
                "axis_top_terms": str(subspace_dir / f"{subspace_safe}__axis_top_terms.csv"),
                "axis_term_loadings_all": str(
                    subspace_dir / f"{subspace_safe}__axis_term_loadings_all.csv"
                ),
                "go_ic_quality_summary": str(
                    subspace_dir / f"{subspace_safe}__go_ic_quality_summary.csv"
                ),
                "cluster_coherence": str(subspace_dir / f"{subspace_safe}__cluster_coherence.csv"),
                "tfidf_cluster_quality": str(
                    subspace_dir / f"{subspace_safe}__tfidf_cluster_quality.csv"
                ),
                "subspace_embedding_terms_only_png": str(
                    subspace_dir / f"{subspace_safe}__subspace_embedding_terms_only.png"
                ),
                "subspace_embedding_annotated_terms_png": str(
                    subspace_dir
                    / f"{subspace_safe}__subspace_embedding_clusters_annotated_terms.png"
                ),
                "adaptive_diffusion_embedding_annotated_terms_png": str(
                    subspace_dir
                    / f"{subspace_safe}__adaptive_diffusion_embedding_clusters_annotated_terms.png"
                ),
            }
        )
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def write_artifact_index(
    output_dir: Path,
    ranking: pd.DataFrame,
    artifact_prefix: str,
    workflow_outputs: dict[str, str],
) -> None:
    write_dataframe(ranking, output_dir / "artifact_index.csv")
    plot_rows = []
    for row in ranking.itertuples(index=False):
        plot_rows.append(
            {
                "display_rank": getattr(row, "display_rank", ""),
                "weighting": row.weighting,
                "block_name": row.block_name,
                "status": row.status,
                "tree_dendrogram_png": getattr(row, "tree_dendrogram_png", ""),
                "subspace_embedding_png": getattr(row, "subspace_embedding_png", ""),
                "adaptive_diffusion_embedding_png": getattr(
                    row, "adaptive_diffusion_embedding_png", ""
                ),
                "axis_terms_combined_png": getattr(row, "axis_terms_combined_png", ""),
            }
        )
    write_dataframe(pd.DataFrame(plot_rows), output_dir / "subspace_plot_index.csv")
    manifest = {
        "artifact_prefix": artifact_prefix,
        "result_prefix": RESULT_PREFIX,
        "ranking_note": (
            "specificity_aware_rank is the preferred reading order: quality tier, "
            "specificity score, specificity fraction, weighted specificity delta, "
            "then GO-BIC/gene. old_display_rank preserves the earlier quality-tiered GO-IC order."
        ),
        "ranking_csv": str(output_dir / "rankings" / f"{RESULT_PREFIX}_ranking.csv"),
        "specificity_ranking_csv": str(
            output_dir / "rankings" / f"{RESULT_PREFIX}_specificity_aware_ranking.csv"
        ),
        "workflow_outputs": workflow_outputs,
        "subspaces": ranking.to_dict(orient="records"),
    }
    (output_dir / "connected_results_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str),
        encoding="utf-8",
    )
    lines = [
        "# Artifact Index",
        "",
        "This experiment uses the full current adaptive-diffusion cosine-subspace tree pipeline.",
        "Every input matrix should produce this same directory structure.",
        "",
        "## Main Files",
        "",
        f"- [ranking CSV](rankings/{RESULT_PREFIX}_ranking.csv)",
        f"- [specificity-aware ranking CSV](rankings/{RESULT_PREFIX}_specificity_aware_ranking.csv)",
        "- [connected manifest](connected_results_manifest.json)",
        "- [artifact table](artifact_index.csv)",
        "- [subspace plot table](subspace_plot_index.csv)",
        "",
        "## PDFs",
        "",
        f"- [method tree-pages PDF]({Path(workflow_outputs['method_tree_pages_pdf']).relative_to(output_dir)})",
        f"- [all tree-pages PDF]({Path(workflow_outputs['all_tree_pages_pdf']).relative_to(output_dir)})",
        "",
        "## Subspaces",
        "",
        "- `subspaces/<weighting>/<block_name>/` contains coordinates, linkage tree, cluster assignments or failure status, GO-IC quality summary, coherence tables, TF-IDF quality, and axis term-loading files.",
    ]
    (output_dir / "ARTIFACT_INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_experiment_readme(
    input_path: Path, output_dir: Path, artifact_prefix: str, ranking: pd.DataFrame
) -> None:
    ranking_path = output_dir / "rankings" / f"{RESULT_PREFIX}_ranking.csv"
    ok_mask = ranking["status"].eq("ok")
    readme_lines = [
        "# Current Adaptive Diffusion Subspace Tree Experiment",
        "",
        f"Input: `{input_path}`",
        "Method version: `current`",
        "Tree geometry: `adaptive_diffusion_cosine_subspace`",
        f"Experiment directory: `{output_dir}`",
        "",
        "Directory layout:",
        "- `rankings/`: experiment-level ranking, spectrum, block metadata, and long quality tables.",
        "- `plots/`: experiment-level summary plots.",
        "- `subspaces/<weighting>/<block_name>/`: one folder per subspace with assignments, tree, quality CSVs, coordinates, and axis term-loading plots.",
        f"- `{artifact_prefix}_quality_aware_go_ic_by_method/current__adaptive_diffusion_cosine_subspace/`: method-separated PDF, PNG pages, and copied ranking CSV.",
        f"- `{artifact_prefix}_quality_aware_go_ic_plots/`: all-tree/ordered PDFs for the experiment.",
        "- `ARTIFACT_INDEX.md`, `artifact_index.csv`, `subspace_plot_index.csv`, and `connected_results_manifest.json`: connected artifact tables.",
        "",
        "Ranking:",
        "- `display_rank` equals `specificity_aware_rank` for completed rows.",
        "- `specificity_aware_rank` sorts by quality tier, cluster specificity score, specific-cluster fraction, weighted specificity delta, then lower GO-BIC active/gene.",
        "- `old_display_rank` preserves the older quality-tier then GO-BIC order.",
        "- `raw_go_ic_rank` preserves the raw GO-IC order for audit.",
        "- `go_bic_active_per_gene` is lower-is-better only within comparable quality tiers.",
        "",
        "Axis term loadings:",
        "- `axis_term_loadings_all.csv` stores every GO term loading for every cosine mode in the subspace.",
        "- `axis_top_terms.csv` stores the top positive, negative, and absolute GO-term loadings per axis.",
        "- Positive and negative signs are orientation-dependent; the absolute loading is the stable strength score.",
        "",
        f"Ranking CSV: `{ranking_path}`",
        f"Artifact index: `{output_dir / 'ARTIFACT_INDEX.md'}`",
    ]
    if ok_mask.any():
        best = ranking[ranking["status"].eq("ok")].sort_values("display_rank").iloc[0]
        readme_lines.extend(
            [
                "",
                "Top ranked subspace:",
                f"- `{best['weighting']} / {best['block_name']}`",
                f"- clusters: `{int(best['n_clusters'])}`",
                f"- quality tier: `{best['quality_tier_label']}`",
                f"- GO-BIC active/gene: `{float(best['go_bic_active_per_gene']):.6f}`",
                f"- coherent clusters: `{int(best['coherent_cluster_count'])}/{int(best['n_clusters'])}`",
            ]
        )
    (output_dir / "README.md").write_text("\n".join(readme_lines) + "\n", encoding="utf-8")


def default_output_dir(input_path: Path, dataset_label: str | None = None) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return (
        Path("results/analyses")
        / f"{matrix_slug(input_path, dataset_label)}_{RESULT_PREFIX}_{stamp}"
    )
