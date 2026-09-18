"""Diffusion-hierarchy boundaries for the selected Gaussian radial path."""

from __future__ import annotations

import numpy as np
import pytest
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_diffusion import (
    replay_adaptive_pydiffmap_average_linkage,
    replay_hamming_knn_average_linkage,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    build_gaussian_selected_projection_path,
)


def _continuous_path():
    feature_matrix = np.random.default_rng(0).normal(size=(12, 3))
    contrast = np.array(
        [1.0, -1.0, 1.0, -1.0, 0.5, -0.5, 0.25, -0.25, 0.0, 0.0, 0.0, 0.0],
        dtype=np.float64,
    )
    return build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        np.eye(3, dtype=np.float64),
        projection_dimension=2,
    )


def _binary_path():
    feature_matrix = np.array(
        [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [1, 0, 0, 1],
            [0, 1, 1, 0],
            [0, 1, 0, 1],
        ],
        dtype=np.float64,
    )
    contrast = np.array(
        [1.0, -1.0, 1.0, -1.0, 0.5, -0.5, 0.25, -0.25, 0.0, 0.0],
        dtype=np.float64,
    )
    return build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        np.eye(4, dtype=np.float64),
        projection_dimension=2,
    )


def _adaptive_replay(path, radius: float, *, metric: str = "euclidean"):
    return replay_adaptive_pydiffmap_average_linkage(
        path,
        radius=radius,
        k_neighbors=10 if metric == "euclidean" else 9,
        diffusion_time=3,
        n_components=30 if metric == "euclidean" else 5,
        metric=metric,
        bandwidth_type="-1/(d+2)",
        epsilon="median",
    )


def test_adaptive_pydiffmap_replay_tracks_data_dependent_distance_and_topology() -> None:
    path = _continuous_path()

    low = _adaptive_replay(path, 0.2 * path.observed_radius)
    observed = _adaptive_replay(path, path.observed_radius)
    high = _adaptive_replay(path, 2.0 * path.observed_radius)

    assert low.geometry_metadata["backend"] == "pydiffmap"
    assert observed.geometry_metadata["metric"] == "euclidean"
    assert low.merge_signature[0] == ((9,), (10,))
    assert observed.merge_signature[0] == ((6,), (11,))
    assert high.merge_signature[0] == ((0,), (11,))
    assert len({low.merge_signature, observed.merge_signature, high.merge_signature}) == 3
    assert np.linalg.norm(low.distance_condensed - observed.distance_condensed) > 1e-4
    assert np.linalg.norm(observed.distance_condensed - high.distance_condensed) > 1e-4
    np.testing.assert_allclose(
        observed.distance_condensed[:3],
        np.array(
            [
                0.00025969006048411433,
                0.0027715216759704817,
                0.003746763194109221,
            ]
        ),
        rtol=1e-10,
        atol=1e-14,
    )


def test_hamming_neighbor_replay_is_defined_at_the_observed_binary_matrix() -> None:
    path = _binary_path()

    replay = replay_hamming_knn_average_linkage(
        path,
        radius=path.observed_radius,
        k_neighbors=3,
        diffusion_time=3,
        n_components=5,
    )

    assert replay.geometry_metadata["backend"] == "sklearn_hamming_knn"
    assert replay.geometry_metadata["metric"] == "hamming"
    assert replay.merge_signature[0] == ((0,), (5,))
    assert len(replay.merge_signature) == 9


def test_adaptive_hamming_replay_preserves_the_benchmark_binary_contract() -> None:
    path = _binary_path()

    replay = _adaptive_replay(path, path.observed_radius, metric="hamming")

    assert replay.geometry_metadata["backend"] == "pydiffmap"
    assert replay.geometry_metadata["metric"] == "hamming"
    assert replay.merge_signature[0] == ((0,), (5,))


def test_hamming_neighbor_replay_tolerates_observed_reconstruction_roundoff() -> None:
    feature_matrix = (np.random.default_rng(0).random((10, 4)) < 0.5).astype(np.float64)
    contrast = np.array(
        [1.0, -1.0, 1.0, -1.0, 0.5, -0.5, 0.25, -0.25, 0.0, 0.0],
        dtype=np.float64,
    )
    path = build_gaussian_selected_projection_path(
        feature_matrix,
        contrast,
        np.eye(4, dtype=np.float64),
        projection_dimension=2,
    )
    reconstructed = path.reconstruct(path.observed_radius)
    assert not np.isin(reconstructed, (0.0, 1.0)).all()
    np.testing.assert_allclose(reconstructed, feature_matrix, atol=1e-12)

    replay = replay_hamming_knn_average_linkage(
        path,
        radius=path.observed_radius,
        k_neighbors=3,
        diffusion_time=3,
        n_components=5,
    )

    assert replay.geometry_metadata["metric"] == "hamming"
    assert len(replay.merge_signature) == 9


@pytest.mark.parametrize(
    "replay",
    [
        pytest.param(
            lambda path, radius: replay_hamming_knn_average_linkage(
                path,
                radius=radius,
                k_neighbors=3,
                diffusion_time=3,
                n_components=5,
            ),
            id="hamming-neighbor",
        ),
        pytest.param(
            lambda path, radius: _adaptive_replay(path, radius, metric="hamming"),
            id="adaptive-pydiffmap-hamming",
        ),
    ],
)
def test_hamming_diffusion_rejects_the_continuous_radial_path_away_from_binary_values(
    replay,
) -> None:
    path = _binary_path()

    with pytest.raises(
        ValueError,
        match="continuous Gaussian radial candidates are unsupported",
    ):
        replay(path, 0.5 * path.observed_radius)
