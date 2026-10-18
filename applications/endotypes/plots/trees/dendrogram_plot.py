"""Render SciPy dendrograms shared by plain and cluster-annotated trees."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram


def draw_dendrogram(ax: plt.Axes, linkage_matrix: np.ndarray) -> np.ndarray:
    result = dendrogram(
        linkage_matrix,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#333333",
        link_color_func=lambda _node_id: "#333333",
        ax=ax,
    )
    return np.asarray(result["leaves"], dtype=int)


def create_dendrogram(linkage_matrix: np.ndarray, title: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(12, 6))
    draw_dendrogram(ax, linkage_matrix)
    ax.set_title(title)
    ax.set_ylabel("adaptive diffusion distance")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    fig.tight_layout()
    return fig
