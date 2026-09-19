"""Experiment-wide ranking plots and adaptive diffusion subspace PDFs."""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import matplotlib

from applications.endotypes.io.adaptive_diffusion_dataframe_io import (
    read_dataframe,
    write_dataframe,
)

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from tree_break_selection.plot.image_panel import draw_image_panel

from applications.endotypes._shared import safe_name
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_terms import (
    axis_term_summary_lines,
)
from applications.endotypes.plots.embeddings.embedding_comparison_plot import (
    create_embedding_comparison,
)
from applications.endotypes.plots.ranking.subspace_ranking_plot import create_subspace_ranking
from applications.endotypes.plots.terms.axis_term_bar_plot import create_axis_term_bars
from applications.endotypes.plots.terms.axis_term_heatmap import create_axis_term_heatmap
from applications.endotypes.plots.trees.tree_cluster_plot import create_tree_clusters
from applications.endotypes.reports.go_ic.adaptive_diffusion_plot_exports import write_figure
from applications.endotypes.reports.go_ic.adaptive_diffusion_subspace_artifacts import RESULT_PREFIX


def write_combined_axis_terms_pdf(ranking: pd.DataFrame, output_dir: Path) -> Path | None:
    ok = ranking[ranking["status"].eq("ok")].sort_values("display_rank")
    if ok.empty:
        return None
    pdf_path = output_dir / f"{RESULT_PREFIX}_axis_terms_combined_by_subspace.pdf"
    with PdfPages(pdf_path) as pdf:
        for row in ok.itertuples(index=False):
            fig, ax = plt.subplots(figsize=(16, 10))
            draw_image_panel(
                ax,
                getattr(row, "axis_terms_combined_png", ""),
                f"Rank {int(row.display_rank)}: {row.weighting} / {row.block_name}\n"
                "Combined signed GO-term feature loadings across subspace eigenmodes",
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return pdf_path


def write_workflow_pdfs(
    output_dir: Path, ranking: pd.DataFrame, artifact_prefix: str
) -> dict[str, str]:
    ok = ranking[ranking["status"].eq("ok")].sort_values("display_rank")
    plots_dir = output_dir / f"{artifact_prefix}_quality_aware_go_ic_plots"
    method_dir = (
        output_dir
        / f"{artifact_prefix}_quality_aware_go_ic_by_method"
        / "current__adaptive_diffusion_cosine_subspace"
    )
    pages_dir = method_dir / "tree_pages"
    plots_dir.mkdir(parents=True, exist_ok=True)
    method_dir.mkdir(parents=True, exist_ok=True)
    pages_dir.mkdir(parents=True, exist_ok=True)

    method_prefix = (
        f"{artifact_prefix}_quality_aware_go_ic_current__adaptive_diffusion_cosine_subspace"
    )
    write_dataframe(ranking, method_dir / f"{method_prefix}_tree_ranking.csv")

    pdf_path = method_dir / f"{method_prefix}_tree_pages.pdf"
    with PdfPages(pdf_path) as pdf:
        for row in ok.itertuples(index=False):
            rank = int(row.display_rank)
            title = (
                f"Rank {rank}: current adaptive diffusion cosine subspace, "
                f"{row.weighting} / {row.block_name}"
            )
            metrics = (
                "Method: current TBS gate on an average-linkage tree built from adaptive diffusion "
                "distances inside this cosine eigenspace subspace.\n"
                "Values: GO-BIC/gene is lower-is-better within quality tier; specificity score "
                "combines the fraction and strength of clusters with specific enriched GO terms.\n"
                f"Clusters: {int(row.n_clusters)}; specific clusters: {int(row.specific_cluster_count)}; "
                f"specificity score: {float(row.specificity_score):.3f}; "
                f"GO-BIC/gene: {float(row.go_bic_active_per_gene):.3f}."
            )

            linkage_path = Path(str(getattr(row, "linkage_matrix", "")))
            assignments_path = Path(str(getattr(row, "cluster_assignments", "")))
            if linkage_path.exists() and assignments_path.exists():
                fig = create_tree_clusters(
                    read_dataframe(linkage_path).to_numpy(),
                    read_dataframe(assignments_path),
                    f"{title}: tree with cluster assignments",
                    tight_layout=False,
                )
            else:
                fig = plt.figure(figsize=(18, 10.5))
                ax = fig.add_subplot(111)
                draw_image_panel(
                    ax,
                    getattr(row, "tree_subtree_clusters_png", "")
                    or getattr(row, "tree_dendrogram_png", ""),
                    f"{title}: tree with cluster assignments",
                )
            fig.text(
                0.02,
                0.02,
                textwrap.fill(metrics.replace("\n", " "), width=180),
                ha="left",
                va="bottom",
                fontsize=8.5,
            )
            fig.tight_layout(rect=(0.0, 0.055, 1.0, 1.0))
            page_path = pages_dir / (
                f"{rank:02d}_current__adaptive_diffusion_cosine_subspace__"
                f"{safe_name(row.weighting)}__{safe_name(row.block_name)}__tree_page.png"
            )
            fig.savefig(page_path, dpi=180)
            pdf.savefig(fig)
            plt.close(fig)

            images = []
            for column in ("subspace_embedding_png", "adaptive_diffusion_embedding_png"):
                path = Path(str(getattr(row, column, "")))
                images.append(plt.imread(path) if path.is_file() else None)
            terms_path = Path(str(getattr(row, "axis_top_terms", "")))
            term_lines = (
                axis_term_summary_lines(read_dataframe(terms_path)) if terms_path.is_file() else []
            )
            fig = create_embedding_comparison(images[0], images[1], title, term_lines)
            page_path = pages_dir / (
                f"{rank:02d}_current__adaptive_diffusion_cosine_subspace__"
                f"{safe_name(row.weighting)}__{safe_name(row.block_name)}__embedding_page.png"
            )
            fig.savefig(page_path, dpi=180)
            pdf.savefig(fig)
            plt.close(fig)

            loadings_path = Path(str(getattr(row, "axis_term_loadings_all", "")))
            if loadings_path.exists():
                fig = create_axis_term_heatmap(
                    read_dataframe(loadings_path),
                    f"{title}: strongest GO-term feature loadings",
                    terms_per_axis=6,
                    max_terms=36,
                    report=True,
                )
            else:
                fig = plt.figure(figsize=(18, 10.5))
                ax = fig.add_subplot(111)
                draw_image_panel(
                    ax, getattr(row, "axis_terms_combined_png", ""), f"{title}: GO-term loadings"
                )
                fig.tight_layout()
            page_path = pages_dir / (
                f"{rank:02d}_current__adaptive_diffusion_cosine_subspace__"
                f"{safe_name(row.weighting)}__{safe_name(row.block_name)}__terms_page.png"
            )
            fig.savefig(page_path, dpi=180)
            pdf.savefig(fig)
            plt.close(fig)

    shutil.copyfile(
        pdf_path, plots_dir / f"{artifact_prefix}_quality_aware_go_ic_all_tree_pages.pdf"
    )
    write_dataframe(ranking, plots_dir / f"{artifact_prefix}_quality_aware_go_ic_tree_ranking.csv")

    (method_dir / "README.md").write_text(
        "\n".join(
            [
                "# Current Adaptive Diffusion Cosine-Subspace Results",
                "",
                f"- Main tree/subspace PDF: `{pdf_path.name}`",
                "- `tree_pages/`: three PNG pages per ranked subspace.",
                "- Ranking is specificity-aware within quality tier.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (plots_dir / "README.md").write_text(
        "\n".join(
            [
                "# Quality-Aware GO-IC Plots",
                "",
                f"- All tree/subspace pages PDF: `{artifact_prefix}_quality_aware_go_ic_all_tree_pages.pdf`",
                "- Method-separated results are under `../"
                f"{artifact_prefix}_quality_aware_go_ic_by_method/`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "method_dir": str(method_dir),
        "method_tree_pages_pdf": str(pdf_path),
        "all_tree_pages_pdf": str(
            plots_dir / f"{artifact_prefix}_quality_aware_go_ic_all_tree_pages.pdf"
        ),
    }


def write_ranking_plots(ranking: pd.DataFrame, output_dir: Path) -> None:
    plots_dir = output_dir / "plots"
    ok_mask = ranking["status"].eq("ok")
    if ok_mask.any():
        write_figure(
            create_subspace_ranking(ranking),
            plots_dir / f"{RESULT_PREFIX}_top_ranked_subspaces.png",
            dpi=180,
        )

        pdf_path = output_dir / f"{RESULT_PREFIX}_axis_terms_by_subspace.pdf"
        with PdfPages(pdf_path) as pdf:
            for row in (
                ranking[ranking["status"].eq("ok")]
                .sort_values("display_rank")
                .itertuples(index=False)
            ):
                subspace_id = f"{row.weighting}__{row.block_name}"
                subspace_safe = safe_name(subspace_id)
                axis_path = Path(row.subspace_dir) / f"{subspace_safe}__axis_top_terms.csv"
                top_terms = read_dataframe(axis_path)
                for axis, frame in top_terms[top_terms["selection"].eq("top_absolute")].groupby(
                    "axis", sort=True
                ):
                    fig = create_axis_term_bars(
                        frame,
                        f"Rank {int(row.display_rank)}: {subspace_id}, mode {int(axis):02d}\n"
                        "Top absolute GO-term loadings for this subspace axis",
                        top_n=15,
                        figsize=(11, 7),
                    )
                    pdf.savefig(fig)
                    plt.close(fig)
        write_combined_axis_terms_pdf(ranking, output_dir)
