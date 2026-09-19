"""Shared scatter rendering for subspace and diffusion embeddings."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from applications.endotypes.plots.shared._cluster_legend import (
    cluster_color_spec,
    cluster_legend_height,
    draw_cluster_legend,
)
from applications.endotypes.plots.shared._plot_text import wrap_annotation_lines


def embedding_figure(
    axis_term_lines: list[str] | None, labels: np.ndarray | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    annotation = wrap_annotation_lines(axis_term_lines or [])
    plot_height = max(6.0, 0.2 * len(annotation.splitlines()) + 0.8)
    legend_height = cluster_legend_height(labels, columns=8) if labels is not None else 0.0
    fig = plt.figure(figsize=(16 if annotation else 8, plot_height + legend_height + 0.5))
    grid = fig.add_gridspec(
        2 if labels is not None else 1,
        2 if annotation else 1,
        height_ratios=[plot_height, legend_height] if labels is not None else [plot_height],
    )
    ax = fig.add_subplot(grid[0, 0])
    if annotation:
        text_ax = fig.add_subplot(grid[0, 1])
        text_ax.axis("off")
        text_ax.set_title("Highest GO-term loadings by eigenmode", loc="left", fontsize=12)
        text_ax.text(
            0.0, 0.98, annotation, va="top", ha="left", fontsize=11,
            linespacing=1.2, transform=text_ax.transAxes,
        )
    if labels is not None:
        draw_cluster_legend(fig.add_subplot(grid[1, :]), labels, columns=8)
    return fig, ax


def draw_embedding(
    embedding: np.ndarray,
    labels: np.ndarray | None,
    title: str,
    *,
    xlabel: str,
    ylabel: str,
    axis_term_lines: list[str] | None = None,
) -> plt.Figure:
    fig, ax = embedding_figure(axis_term_lines, labels)
    if labels is None:
        ax.scatter(
            embedding[:, 0], embedding[:, 1], color="#4c78a8", s=18, alpha=0.72, linewidths=0
        )
    else:
        spec = cluster_color_spec(labels)
        ax.scatter(
            embedding[:, 0], embedding[:, 1], c=labels, cmap=spec.cmap, norm=spec.norm,
            s=18, linewidths=0,
        )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    return fig
