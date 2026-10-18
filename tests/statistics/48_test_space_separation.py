"""Tests for reusable space-separation method interfaces."""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

from tree_break_selection.space_separation import (
    SpectralBlock,
    adaptive_spectral_blocks,
    block_diffusion_geometry,
    coordinates_for_block,
    cosine_eigendecomposition,
    decompose_leading_axes_space,
    hamming_knn_diffusion_geometry,
    separate_adaptive_cosine_space,
)


@pytest.mark.parametrize("perturbation", [0.0, 1e-12])
def test_spectral_blocks_keep_repeated_eigenspace_whole(perturbation: float) -> None:
    eigenvalues = np.array([10, 9, 8, 7, 6, 6, 6, 6, 6, 6, 2, 1], dtype=float)
    eigenvalues[4:10] += np.linspace(perturbation, 0, 6)
    blocks, _ = adaptive_spectral_blocks(eigenvalues, min_segment_length=4, max_segments=8)
    assert [i for b in blocks for i in range(b.block_start, b.block_end + 1)] == list(range(1, 13))
    assert all(b.block_end - b.block_start + 1 >= 4 for b in blocks)
    assert any(b.block_start <= 5 and b.block_end >= 10 for b in blocks)
    eigenvectors = np.eye(12)
    rotated = eigenvectors.copy()
    rotation, _ = np.linalg.qr(np.random.default_rng(42).normal(size=(6, 6)))
    rotated[:, 4:10] = eigenvectors[:, 4:10] @ rotation
    for block in blocks:
        np.testing.assert_allclose(
            pdist(coordinates_for_block(eigenvalues, eigenvectors, block)),
            pdist(coordinates_for_block(eigenvalues, rotated, block)),
            atol=1e-11,
        )


@pytest.mark.parametrize("max_rank,expected_rank", [(1, 1), (2, 4), (4, 4), (8, 4)])
@pytest.mark.parametrize("perturbation", [0.0, 1e-12])
def test_rank_cutoff_preserves_repeated_cosine_eigenspace(
    max_rank: int, expected_rank: int, perturbation: float
) -> None:
    # A valid cosine operator with eigenvalues 2.5, 0.5, 0.5, 0.5.
    operator = 0.5 * np.eye(4) + 0.5 * np.ones((4, 4))
    operator[0, 1] += perturbation
    operator[1, 0] += perturbation
    values = np.linalg.cholesky(operator)
    eigenvalues, eigenvectors = cosine_eigendecomposition(values, max_rank)
    assert len(eigenvalues) == expected_rank
    assert eigenvectors.shape == (4, expected_rank)
    if expected_rank == 4:
        np.testing.assert_allclose((eigenvectors * eigenvalues) @ eigenvectors.T, operator)


def test_tied_leading_modes_are_not_isolated_as_common_mode() -> None:
    eigenvalues = np.ones(12)
    eigenvalues[0] += 1e-12
    blocks, diagnostics = adaptive_spectral_blocks(
        eigenvalues, min_segment_length=4, max_segments=8
    )
    assert [(b.block_start, b.block_end) for b in blocks] == [(1, 12)]
    assert not diagnostics["common_mode_isolated"]
    values, vectors = cosine_eigendecomposition(np.eye(12), max_rank=1)
    assert values.shape == (12,)
    np.testing.assert_allclose((vectors * values) @ vectors.T, np.eye(12))


def test_distinct_spectrum_retains_common_mode_and_decay_regimes() -> None:
    blocks, diagnostics = adaptive_spectral_blocks(
        np.array([100, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5]),
        min_segment_length=4,
        max_segments=8,
    )
    assert [(b.block_start, b.block_end) for b in blocks] == [(1, 1), (2, 8), (9, 12)]
    assert diagnostics["common_mode_isolated"]


def test_leading_axes_decomposition_is_finite_and_energy_ordered() -> None:
    values = np.array(
        [
            [1.0, 2.0, 3.0],
            [2.0, 2.0, 4.0],
            [3.0, 2.0, 5.0],
            [4.0, 2.0, 7.0],
        ]
    )

    result = decompose_leading_axes_space(values, following_dim=2)

    assert result.coordinates.shape == (4, 3)
    assert result.axes.shape == (3, 3)
    assert np.isfinite(result.coordinates).all()
    assert np.all(np.diff(result.singular_values) <= 0.0)


def test_adaptive_cosine_space_exposes_selected_block_coordinates() -> None:
    data = pd.DataFrame(
        [
            [1.0, 0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
        ]
    )

    result = separate_adaptive_cosine_space(
        data,
        weighting="raw",
        max_rank=4,
        min_segment_length=2,
        max_segments=2,
    )

    assert result.blocks
    coordinates = result.coordinates(result.blocks[0])
    assert coordinates.shape[0] == len(data)
    assert np.isfinite(coordinates).all()


def test_coordinates_for_block_uses_one_based_inclusive_modes() -> None:
    eigenvalues = np.array([9.0, 4.0, 1.0])
    eigenvectors = np.eye(3)
    block = SpectralBlock(0, "modes_02_03", 2, 3, "test")

    coordinates = coordinates_for_block(eigenvalues, eigenvectors, block)

    np.testing.assert_allclose(
        coordinates,
        np.array(
            [
                [0.0, 0.0],
                [2.0, 0.0],
                [0.0, 1.0],
            ]
        ),
    )


def test_block_diffusion_geometry_is_a_public_separation_interface() -> None:
    coordinates = np.array(
        [[0.0, 0.0], [0.1, 0.0], [1.0, 1.0], [1.1, 1.0]],
        dtype=float,
    )

    geometry = block_diffusion_geometry(
        coordinates,
        k_neighbors=2,
        diffusion_time=2,
        n_components=3,
    )

    assert geometry.distance_condensed.shape == (6,)
    assert geometry.coordinates.shape == (4, 3)
    assert np.isfinite(geometry.distance_condensed).all()
    assert geometry.metadata["kernel"] == "knn_gaussian"


def test_hamming_diffusion_geometry_is_a_public_separation_interface() -> None:
    binary = pd.DataFrame(
        [
            [1, 0, 0, 1],
            [1, 0, 1, 1],
            [0, 1, 1, 0],
            [0, 1, 0, 0],
        ]
    )

    geometry = hamming_knn_diffusion_geometry(
        binary,
        k_neighbors=2,
        diffusion_time=2,
        n_components=3,
    )

    assert geometry.distance_condensed.shape == (6,)
    assert np.isfinite(geometry.distance_condensed).all()


def test_hamming_diffusion_geometry_rejects_continuous_values() -> None:
    with np.testing.assert_raises_regex(ValueError, "binary or one-hot"):
        hamming_knn_diffusion_geometry(
            np.array([[0.0, 0.5], [1.0, 0.0]]),
            k_neighbors=1,
            diffusion_time=1,
            n_components=1,
        )


def test_hamming_diffusion_with_one_neighbor_keeps_a_non_self_edge():
    from tree_break_selection.space_separation import hamming_knn_diffusion_geometry

    rng = np.random.default_rng(0)
    values = (rng.random((40, 20)) < 0.5).astype(float)

    geometry = hamming_knn_diffusion_geometry(
        values, k_neighbors=1, diffusion_time=3, n_components=10
    )

    assert geometry.metadata["neighbor_search_k"] == 2
    assert float(np.max(geometry.distance_condensed)) > 0.0
