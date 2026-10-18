"""Create an adaptive diffusion distance embedding without file I/O."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS

from applications.endotypes.plots.embeddings._embedding_scatter import draw_embedding


def create_diffusion_embedding(
    distances: np.ndarray,
    labels: np.ndarray,
    title: str,
    *,
    axis_term_lines: list[str] | None = None,
) -> plt.Figure:
    matrix = squareform(distances)
    embedding = MDS(
        n_components=2,
        dissimilarity="precomputed",
        random_state=1729,
        normalized_stress="auto",
    ).fit_transform(matrix)
    return draw_embedding(
        embedding,
        labels,
        f"{title}\nadaptive diffusion distance embedding by cluster",
        xlabel="MDS axis 1",
        ylabel="MDS axis 2",
        axis_term_lines=axis_term_lines,
    )
