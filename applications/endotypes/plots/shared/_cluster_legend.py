"""Shared categorical colors and explicit integer keys for GO-IC plots."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from tree_break_selection.plot.cluster_color_mapping import (
    ClusterColorSpec,
    build_cluster_color_spec,
)


def cluster_color_spec(labels: np.ndarray) -> ClusterColorSpec:
    return build_cluster_color_spec(max(int(np.max(labels)) + 1, 0))


def cluster_legend_height(labels: np.ndarray, *, columns: int) -> float:
    rows = int(np.ceil(len(np.unique(labels)) / columns))
    return 0.4 + 0.25 * rows


def draw_cluster_legend(ax: plt.Axes, labels: np.ndarray, *, columns: int) -> None:
    spec = cluster_color_spec(labels)
    ids = np.unique(labels)
    ax.axis("off")
    ax.legend(
        handles=[Patch(facecolor=spec.id_to_color[int(cid)], label=str(int(cid))) for cid in ids],
        title="Cluster ID",
        ncol=min(columns, len(ids)),
        loc="center",
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        handlelength=1.2,
        columnspacing=1.2,
    )
