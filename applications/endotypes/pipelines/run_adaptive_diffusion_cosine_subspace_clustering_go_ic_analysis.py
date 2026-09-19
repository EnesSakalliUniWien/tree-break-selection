#!/usr/bin/env python3
"""Run adaptive diffusion cosine subspace clustering and GO-IC analysis.

Separate the feature matrix into adaptive cosine subspaces, construct
adaptive diffusion average-linkage trees, apply Tree-Break Selection clustering,
and evaluate GO information criteria, coherence, TF-IDF quality, and
specificity-aware ranking.

The output layout is intentionally subspace-first:

    <output-dir>/
      rankings/
      plots/
      subspaces/<weighting>/<block_name>/

Each subspace folder contains the cluster assignments, linkage tree, GO-IC
scores, cluster coherence, TF-IDF quality, subspace coordinates, diffusion
metadata, and per-axis GO-term loading tables/plots.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_EDGE_ALPHA,
    DEFAULT_SIBLING_ALPHA,
)

from applications.endotypes.pipelines.adaptive_diffusion_go_ic_workflow import run_analysis
from applications.endotypes.pipelines.tree_analysis_args import add_tree_analysis_arguments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/feature_matrices/feature_matrix_julia_allGO_new.tsv"),
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    add_tree_analysis_arguments(
        parser,
        edge_alpha_default=DEFAULT_EDGE_ALPHA,
        sibling_alpha_default=DEFAULT_SIBLING_ALPHA,
    )
    parser.add_argument("--top-terms-per-axis", type=int, default=20)
    parser.add_argument(
        "--dataset-label",
        default=None,
        help="Optional label used in experiment directory and reader-facing artifact names.",
    )
    return parser.parse_args()


def main() -> None:
    run_analysis(parse_args())


if __name__ == "__main__":
    main()
