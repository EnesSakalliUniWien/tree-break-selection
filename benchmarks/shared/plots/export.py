"Batch helpers to render and save plots from benchmark results."

from __future__ import annotations

import logging
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from benchmarks.shared.result_records import ComputedResultRecord
from benchmarks.shared.util.params import format_params_for_display
from benchmarks.shared.util.pdf.layout import (
    PDF_PAGE_SIZE_INCHES,
    PDF_WIDE_PAGE_SIZE_INCHES,
    prepare_pdf_figure,
    set_pdf_page_size,
)
from tree_break_selection.plot.cluster_tree_visualization import plot_tree_with_clusters

from .embedding import create_clustering_comparison_plots

logger = logging.getLogger(__name__)

DEFAULT_BENCHMARK_TREE_LAYOUT = "radial"


def _benchmark_tree_layout() -> str:
    """Return the tree layout used by benchmark report pages."""
    layout = os.getenv("TBS_BENCHMARK_TREE_LAYOUT", DEFAULT_BENCHMARK_TREE_LAYOUT)
    layout = layout.strip().lower()
    if layout not in {"radial", "rectangular"}:
        raise ValueError(
            f"TBS_BENCHMARK_TREE_LAYOUT must be 'radial' or 'rectangular', got {layout!r}."
        )
    return layout


def _format_ari_nmi(result: ComputedResultRecord, *, compact: bool = False) -> str:
    """Format primary benchmark metrics when available on a computed result record."""
    ari = result.ari
    nmi = result.nmi
    parts: list[str] = []
    found_clusters = result.meta.get("found_clusters")
    if isinstance(found_clusters, (int, float, np.integer, np.floating)) and np.isfinite(
        float(found_clusters)
    ):
        parts.append(f"K={int(found_clusters)}" if compact else f"Found_K={int(found_clusters)}")
    if isinstance(ari, (int, float, np.floating)) and np.isfinite(float(ari)):
        parts.append(f"A={float(ari):.3f}" if compact else f"ARI={float(ari):.3f}")
    if isinstance(nmi, (int, float, np.floating)) and np.isfinite(float(nmi)):
        parts.append(f"N={float(nmi):.3f}" if compact else f"NMI={float(nmi):.3f}")
    outlier_f1 = result.outlier_f1
    singleton_hit = result.singleton_outlier_isolated
    grouped_recovery = result.grouped_outlier_cluster_recovered
    if isinstance(outlier_f1, (int, float, np.floating)) and np.isfinite(float(outlier_f1)):
        parts.append(
            f"OF1={float(outlier_f1):.3f}" if compact else f"Outlier_F1={float(outlier_f1):.3f}"
        )
    if isinstance(singleton_hit, (int, float, np.floating)) and np.isfinite(float(singleton_hit)):
        parts.append(
            f"SOI={float(singleton_hit):.0f}"
            if compact
            else f"Singleton_Isolated={float(singleton_hit):.0f}"
        )
    if isinstance(grouped_recovery, (int, float, np.floating)) and np.isfinite(
        float(grouped_recovery)
    ):
        parts.append(
            f"GOR={float(grouped_recovery):.0f}"
            if compact
            else f"Grouped_Outlier_Recovered={float(grouped_recovery):.0f}"
        )
    return ", ".join(parts)


def _shorten_title_line(text: str, width: int) -> str:
    text = " ".join(str(text).split())
    if len(text) <= width:
        return text
    return textwrap.shorten(text, width=width, placeholder="...")


def _wrap_title_line(text: str, width: int) -> list[str]:
    return textwrap.wrap(
        _shorten_title_line(text, width * 2),
        width=width,
        break_long_words=True,
        break_on_hyphens=False,
    )


def _tree_plot_style(tree, decomposition: dict) -> dict[str, object]:
    """Return size and legend defaults scaled for tree density."""
    node_count = len(list(tree.nodes()))
    num_clusters = int(decomposition.get("num_clusters", 0))
    dense = node_count > 80 or num_clusters > 20
    return {
        "figsize": PDF_WIDE_PAGE_SIZE_INCHES if dense else PDF_PAGE_SIZE_INCHES,
        "node_size": 5 if dense else 12,
        "font_size": 7 if dense else 9,
        "show_legend": not dense,
        "max_cluster_legend_entries": 0 if num_clusters > 20 else 20,
        "subplots_right": 0.82 if num_clusters <= 20 else 0.88,
    }


def _group_results_by_case(
    test_results: list[ComputedResultRecord],
) -> dict[int, list[ComputedResultRecord]]:
    results_by_case: dict[int, list[ComputedResultRecord]] = {}
    for result in test_results:
        case_num = int(result.test_case_num)
        results_by_case.setdefault(case_num, []).append(result)
    return results_by_case


def _save_pdf_page(pdf: PdfPages, fig: plt.Figure) -> None:
    prepare_pdf_figure(fig)
    pdf.savefig(fig)
    plt.close(fig)


def _create_tree_figures_for_case(
    *,
    case_num: int,
    case_results: list[ComputedResultRecord],
) -> list[plt.Figure]:
    # Keep only tree-capable runs.
    tree_results = [r for r in case_results if r.tree is not None and r.decomposition is not None]
    if not tree_results:
        return []

    meta = tree_results[0].meta
    meta_text = (
        f"samples={meta['n_samples']}, "
        f"features={meta['n_features']}, "
        f"generator={meta['generator']}, "
        f"noise={meta['noise']}"
    )

    figs: list[plt.Figure] = []
    n_items = len(tree_results)

    for idx, result in enumerate(tree_results, start=1):
        plot_style = _tree_plot_style(result.tree, result.decomposition)
        layout = _benchmark_tree_layout()
        fig, ax = plt.subplots(1, 1, figsize=plot_style["figsize"])
        set_pdf_page_size(fig, plot_style["figsize"])
        fig.suptitle(
            f"Tree Comparisons – Test Case {case_num}\n{meta_text}",
            fontsize=16,
            weight="bold",
            y=0.97,
        )

        tree_t = result.tree
        decomp_t = result.decomposition
        annotations_df = result.annotations
        method_name = result.method_name
        params = result.params
        param_str_display = format_params_for_display(params)
        metrics_text = _format_ari_nmi(result)

        page_tag = f" ({idx}/{n_items})" if n_items > 1 else ""
        title_lines = _wrap_title_line(f"{method_name}{page_tag}", 72)
        if param_str_display:
            title_lines.extend(_wrap_title_line(f"({param_str_display})", 72))
        if metrics_text:
            title_lines.extend(_wrap_title_line(metrics_text, 72))
        title = "\n".join(title_lines)

        plot_tree_with_clusters(
            tree=tree_t,
            decomposition_results=decomp_t,
            annotations_df=annotations_df,
            node_size=int(plot_style["node_size"]),
            font_size=int(plot_style["font_size"]),
            title=title,
            layout=layout,
            ax=ax,
            show_legend=bool(plot_style["show_legend"]),
            max_cluster_legend_entries=int(plot_style["max_cluster_legend_entries"]),
            legend_outside=True,
        )
        fig.subplots_adjust(
            top=0.82,
            bottom=0.06,
            left=0.03,
            right=float(plot_style["subplots_right"]),
        )
        figs.append(fig)

    return figs


def create_case_report_pages_from_results(
    test_results: list[ComputedResultRecord],
    *,
    pdf: PdfPages,
    verbose: bool = True,
) -> None:
    """Write each case's tree pages followed by its UMAP comparison pages."""
    for case_num, case_results in sorted(_group_results_by_case(test_results).items()):
        first_result = case_results[0]
        meta = first_result.meta
        if verbose:
            print(f"  Creating Tree→UMAP plots for test case {case_num}...")

        for tree_fig in _create_tree_figures_for_case(case_num=case_num, case_results=case_results):
            _save_pdf_page(pdf, tree_fig)

        labels_to_plot = {"Ground Truth": first_result.y_true}
        for res in case_results:
            param_str = format_params_for_display(res.params)
            unique_key = f"{res.method_name} ({param_str})" if param_str else res.method_name
            metrics_text = _format_ari_nmi(res, compact=True)
            if metrics_text:
                unique_key = f"{unique_key} [{metrics_text}]"
            if res.labels is not None:
                labels_to_plot[unique_key] = res.labels

        for umap_fig in create_clustering_comparison_plots(
            X_original=first_result.x_original,
            labels_dict=labels_to_plot,
            test_case_num=case_num,
            meta=meta,
            cache_key=f"{meta['name']}_n{meta['n_samples']}",
        ):
            _save_pdf_page(pdf, umap_fig)
        plt.close("all")


def create_manifold_plots_from_results(
    test_results: list[ComputedResultRecord],
    *,
    pdf: PdfPages,
    verbose: bool = True,
) -> None:
    """Write one manifold-alignment page per result that has labels."""
    from .manifold import create_manifold_alignment_plot

    for result in test_results:
        if result.labels is None:
            continue

        i = result.test_case_num
        method_name = result.method_name
        param_str_display = format_params_for_display(result.params)
        metrics_text = _format_ari_nmi(result)

        title = f"Manifold Alignment - {method_name} ({param_str_display})\nTest Case {i}"
        if metrics_text:
            title = f"{title}\n{metrics_text}"

        if verbose:
            print(f"  Creating manifold plot for test case {i}, {method_name}...")

        fig, _mantel_r, _mantel_p = create_manifold_alignment_plot(
            result.x_original,
            np.asarray(result.labels),
            test_case_num=i,
            meta=result.meta,
            y_true=result.y_true,
            title=title,
        )
        _save_pdf_page(pdf, fig)
        plt.close("all")


def create_tree_plots_from_results(
    test_results: list[ComputedResultRecord],
    *,
    pdf: PdfPages,
) -> None:
    """Write one tree page per tree-capable result, grouped by test case."""
    for case_num, case_results in sorted(_group_results_by_case(test_results).items()):
        for fig in _create_tree_figures_for_case(case_num=case_num, case_results=case_results):
            _save_pdf_page(pdf, fig)


__all__ = [
    "create_case_report_pages_from_results",
    "create_manifold_plots_from_results",
    "create_tree_plots_from_results",
    "DEFAULT_BENCHMARK_TREE_LAYOUT",
]
