from __future__ import annotations

import pandas as pd
import pytest

from applications.endotypes.pipelines.run_feature_matrix_with_umap import build_linkage_tree


def _small_feature_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        index=["A", "B", "C", "D"],
    )


def test_build_linkage_tree_rejects_ward_with_non_euclidean_metric() -> None:
    with pytest.raises(ValueError, match="Ward linkage requires euclidean distance"):
        build_linkage_tree(
            data_df=_small_feature_matrix(),
            tree_method="tbs",
            tree_distance_metric="cosine",
            tree_linkage_method="ward",
            diffusion_k_neighbors=10,
            diffusion_time=1,
            diffusion_components=10,
            adaptive_neighbor_k=None,
            adaptive_bandwidth_type="local",
            adaptive_epsilon="auto",
            adaptive_metric="hamming",
        )


def test_build_linkage_tree_allows_ward_with_euclidean_metric() -> None:
    linkage_matrix, distance_metric, linkage_method, adaptive_metadata = build_linkage_tree(
        data_df=_small_feature_matrix(),
        tree_method="tbs",
        tree_distance_metric="euclidean",
        tree_linkage_method="ward",
        diffusion_k_neighbors=10,
        diffusion_time=1,
        diffusion_components=10,
        adaptive_neighbor_k=None,
        adaptive_bandwidth_type="local",
        adaptive_epsilon="auto",
        adaptive_metric="hamming",
    )

    assert linkage_matrix.shape == (3, 4)
    assert distance_metric == "euclidean"
    assert linkage_method == "ward"
    assert adaptive_metadata is None


def test_decomposition_tree_carries_gate_columns_for_the_tree_plot() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from applications.endotypes.pipelines.run_feature_matrix_with_umap import _run_decomposition
    from tree_break_selection.plot.cluster_tree_visualization import plot_tree_with_clusters

    rng = np.random.default_rng(0)
    centers = np.array([[0.9] * 15 + [0.1] * 15, [0.1] * 15 + [0.9] * 15])
    values = (rng.random((40, 30)) < np.repeat(centers, 20, axis=0)).astype(float)
    data = pd.DataFrame(values, index=[f"g{i}" for i in range(40)])

    tree, decomposition, _assignments = _run_decomposition(
        data_df=data,
        edge_alpha=0.05,
        sibling_alpha=0.05,
        tree_method="tbs",
        tree_distance_metric="hamming",
        tree_linkage_method="average",
        diffusion_k_neighbors=10,
        diffusion_time=1,
        diffusion_components=10,
        adaptive_neighbor_k=None,
        adaptive_bandwidth_type="local",
        adaptive_epsilon="auto",
        adaptive_metric="hamming",
    )

    assert "Sibling_Divergence_Skipped" in tree.annotations_df.columns
    fig, _ax = plot_tree_with_clusters(
        tree, decomposition, annotations_df=tree.annotations_df, layout="rectangular"
    )
    plt.close(fig)
