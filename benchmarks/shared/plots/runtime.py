"""
Runtime helpers to render and save benchmark plots."""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from benchmarks.shared.result_records import (
    ComputedResultRecord,
)
from benchmarks.shared.util.pdf.layout import prepare_pdf_figure

from .cover_page import GROUP_ORDER, generate_overview_page, generate_section_page
from .export import (
    create_case_report_pages_from_results,
    create_manifold_plots_from_results,
    create_tree_plots_from_results,
)
from .summary import create_validation_plot

logger = logging.getLogger(__name__)


def generate_benchmark_plots(
    df_results: pd.DataFrame,
    computed_results: list[ComputedResultRecord],
    *,
    pdf: PdfPages,
    verbose: bool,
    plot_umap: bool,
    plot_manifold: bool,
    include_cover_pages: bool = True,
    include_validation_page: bool = True,
) -> None:
    """Write benchmark report pages into an open PDF."""
    if df_results.empty:
        return

    # Full benchmark reports add their own global overview/section pages during
    # top-level PDF assembly. Allow per-case runs to skip these to avoid
    # duplicated/misaligned explanation pages in concatenated reports.
    if include_cover_pages:
        n_cases = df_results["test_case"].nunique()
        cover_figs = [generate_overview_page(n_cases=n_cases)]
        for group in GROUP_ORDER:
            section_fig = generate_section_page(group)
            if section_fig is not None:
                cover_figs.append(section_fig)
        for cfig in cover_figs:
            pdf.savefig(cfig)
            plt.close(cfig)
        logger.info("Wrote %d cover/setup pages to PDF.", len(cover_figs))

    if include_validation_page:
        fig = create_validation_plot(df_results)
        prepare_pdf_figure(fig)
        pdf.savefig(fig)
        plt.close(fig)

    if plot_umap:
        if verbose:
            logger.info("Generating Tree→UMAP comparison pages...")
        create_case_report_pages_from_results(computed_results, pdf=pdf, verbose=verbose)
    else:
        if verbose:
            logger.info("Generating tree plots...")
        create_tree_plots_from_results(computed_results, pdf=pdf)

    if plot_manifold:
        if verbose:
            logger.info("Generating manifold diagnostics...")
        create_manifold_plots_from_results(computed_results, pdf=pdf, verbose=verbose)
