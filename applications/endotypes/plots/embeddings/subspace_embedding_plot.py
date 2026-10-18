"""Create a cosine-subspace embedding figure without file I/O."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

from applications.endotypes.plots.embeddings._embedding_scatter import draw_embedding


def create_subspace_embedding(
    coords: np.ndarray,
    labels: np.ndarray | None,
    title: str,
    *,
    axis_term_lines: list[str] | None = None,
) -> plt.Figure:
    if coords.shape[1] == 1:
        embedding = np.column_stack([coords[:, 0], np.zeros(coords.shape[0])])
        method = "axis"
    else:
        try:
            import umap

            embedding = umap.UMAP(
                n_components=2,
                n_neighbors=min(18, max(2, coords.shape[0] - 1)),
                min_dist=0.05,
                metric="euclidean",
                random_state=1729,
            ).fit_transform(coords)
            method = "UMAP"
        except Exception:
            embedding = PCA(n_components=2, random_state=1729).fit_transform(coords)
            method = "PCA"
    description = "no cluster assignments" if labels is None else "by cluster"
    return draw_embedding(
        embedding,
        labels,
        f"{title}\nsubspace embedding; {description} ({method})"
        if labels is None
        else f"{title}\nsubspace embedding {description} ({method})",
        xlabel="axis 1",
        ylabel="axis 2",
        axis_term_lines=axis_term_lines,
    )
