"""File export adapters for the independent GO-IC plotting modules."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from applications.endotypes.plots.embeddings.diffusion_embedding_plot import (
    create_diffusion_embedding,
)
from applications.endotypes.plots.embeddings.subspace_embedding_plot import (
    create_subspace_embedding,
)
from applications.endotypes.plots.terms.axis_term_bar_plot import create_axis_term_bars
from applications.endotypes.plots.terms.axis_term_heatmap import create_axis_term_heatmap
from applications.endotypes.plots.trees.dendrogram_plot import create_dendrogram
from applications.endotypes.plots.trees.tree_cluster_plot import create_tree_clusters


def write_figure(fig: plt.Figure | None, path: Path, *, dpi: int) -> None:
    """Save and release a figure, including when the destination is unwritable."""
    if fig is None:
        return
    try:
        fig.savefig(path, dpi=dpi)
    finally:
        plt.close(fig)


def plot_axis_terms(axis_frame: pd.DataFrame, output_path: Path, *, top_n: int) -> None:
    for axis, frame in axis_frame.groupby("axis", sort=True):
        fig = create_axis_term_bars(
            frame, f"Mode {int(axis):02d}: strongest GO-term loadings", top_n=top_n
        )
        write_figure(
            fig, output_path.parent / f"{output_path.stem}__mode_{int(axis):02d}.png", dpi=180
        )


def plot_combined_axis_terms(
    axis_frame: pd.DataFrame, output_path: Path, *, terms_per_axis: int = 8, max_terms: int = 60
) -> None:
    write_figure(
        create_axis_term_heatmap(axis_frame, terms_per_axis=terms_per_axis, max_terms=max_terms),
        output_path,
        dpi=220,
    )


def plot_subspace_embedding(
    coords: np.ndarray,
    labels: np.ndarray | None,
    output_path: Path,
    title: str,
    *,
    axis_term_lines: list[str] | None = None,
) -> None:
    write_figure(
        create_subspace_embedding(coords, labels, title, axis_term_lines=axis_term_lines),
        output_path,
        dpi=180,
    )


def plot_diffusion_embedding(
    distances: np.ndarray,
    labels: np.ndarray,
    output_path: Path,
    title: str,
    *,
    axis_term_lines: list[str] | None = None,
) -> None:
    write_figure(
        create_diffusion_embedding(distances, labels, title, axis_term_lines=axis_term_lines),
        output_path,
        dpi=180,
    )


def plot_dendrogram(linkage_matrix: np.ndarray, output_path: Path, title: str) -> None:
    write_figure(create_dendrogram(linkage_matrix, title), output_path, dpi=180)


def plot_tree_subtree_clusters(
    linkage_matrix: np.ndarray, assignments: pd.DataFrame, output_path: Path, title: str
) -> None:
    write_figure(create_tree_clusters(linkage_matrix, assignments, title), output_path, dpi=240)
