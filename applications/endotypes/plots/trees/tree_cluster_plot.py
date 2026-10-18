"""Render tree cluster strips in dendrogram leaf order."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from applications.endotypes.plots.shared._cluster_legend import (
    cluster_color_spec,
    cluster_legend_height,
    draw_cluster_legend,
)
from applications.endotypes.plots.trees.dendrogram_plot import draw_dendrogram


def draw_tree_clusters(
    fig: plt.Figure,
    linkage_matrix: np.ndarray,
    assignments: pd.DataFrame,
    title: str,
) -> None:
    labels = assignments["cluster_id"].to_numpy(dtype=int)
    cluster_sizes = assignments["cluster_id"].value_counts().sort_values(ascending=False)
    n_leaves = len(labels)
    n_clusters = int(assignments["cluster_id"].nunique())

    legend_height = cluster_legend_height(labels, columns=12)
    grid = fig.add_gridspec(4, 1, height_ratios=[7.4, 0.5, legend_height, 0.8])
    ax_tree = fig.add_subplot(grid[0])
    ax_strip = fig.add_subplot(grid[1], sharex=ax_tree)
    ax_legend = fig.add_subplot(grid[2])
    ax_text = fig.add_subplot(grid[3])

    leaf_order = draw_dendrogram(ax_tree, linkage_matrix)
    ordered_labels = labels[leaf_order]
    spec = cluster_color_spec(labels)
    ax_strip.imshow(
        ordered_labels[np.newaxis, :],
        aspect="auto",
        interpolation="nearest",
        cmap=spec.cmap,
        norm=spec.norm,
        extent=(0, 10 * n_leaves, 0, 1),
    )
    ax_tree.set_title(title, fontsize=15, pad=10)
    ax_tree.set_ylabel("adaptive diffusion distance", fontsize=12)
    ax_tree.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_tree.tick_params(axis="y", labelsize=10)
    ax_strip.set_yticks([])
    ax_strip.set_ylabel("cluster", rotation=0, ha="right", va="center", labelpad=30, fontsize=12)
    ax_strip.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax_strip.spines.values():
        spine.set_linewidth(0.6)

    draw_cluster_legend(ax_legend, labels, columns=12)
    ax_text.axis("off")
    top_sizes = ", ".join(
        f"C{int(cid)}={int(size)}" for cid, size in cluster_sizes.head(16).items()
    )
    ax_text.text(
        0.0,
        0.98,
        "Tree cluster assignments: the colored strip is the final current TBS cluster id in dendrogram leaf order.\n"
        f"Clusters: {n_clusters}; genes: {n_leaves}; largest cluster sizes: {top_sizes}",
        va="top",
        fontsize=10,
        linespacing=1.25,
        transform=ax_text.transAxes,
    )


def create_tree_clusters(
    linkage_matrix: np.ndarray, assignments: pd.DataFrame, title: str, *, tight_layout: bool = True
) -> plt.Figure:
    legend_height = cluster_legend_height(assignments["cluster_id"].to_numpy(), columns=12)
    fig = plt.figure(figsize=(18, 10.5 + max(0.0, legend_height - 0.65)))
    draw_tree_clusters(fig, linkage_matrix, assignments, title)
    if tight_layout:
        fig.tight_layout()
    return fig
