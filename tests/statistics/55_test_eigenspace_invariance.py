"""Regression contracts for complete repeated eigenspaces at rank boundaries."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.spatial.distance import pdist
from scipy.stats import chi2

from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_projection_basis import (
    build_pca_projection_basis,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projection_dimension_estimation.projection_dimension_estimators import (
    estimate_marchenko_pastur_dimension,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.marchenko_pastur import (
    _process_node,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.node_spectral_task import (
    NodeSpectralTask,
)
from tree_break_selection.space_separation import diffusion
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def _six_cycle() -> np.ndarray:
    weights = np.zeros((6, 6))
    for index in range(6):
        weights[index, (index + 1) % 6] = weights[(index + 1) % 6, index] = 1.0
    return weights


@pytest.mark.parametrize("diffusion_time", [0, 1, 2])
def test_diffusion_tied_cutoff_is_sample_permutation_invariant(diffusion_time: int) -> None:
    weights = _six_cycle()
    permutation = np.array([2, 5, 0, 4, 1, 3])
    coordinates = diffusion.compute_diffusion_coordinates(
        weights,
        diffusion_time=diffusion_time,
        n_components=1,
    )
    permuted = diffusion.compute_diffusion_coordinates(
        weights[np.ix_(permutation, permutation)],
        diffusion_time=diffusion_time,
        n_components=1,
    )[np.argsort(permutation)]

    assert coordinates.shape == (6, 2)
    assert np.max(pdist(coordinates)) > 0.1
    np.testing.assert_allclose(pdist(coordinates), pdist(permuted), atol=1e-12)


def test_diffusion_is_invariant_to_rotation_of_retained_tied_basis(monkeypatch) -> None:
    weights = _six_cycle()
    original_eigh = diffusion.eigh
    reference = diffusion.compute_diffusion_coordinates(weights, diffusion_time=2, n_components=1)

    def rotated_eigh(matrix):
        eigenvalues, eigenvectors = original_eigh(matrix)
        tied = np.flatnonzero(np.isclose(eigenvalues, 0.5))
        rotation = np.array([[1.0, -1.0], [1.0, 1.0]]) / np.sqrt(2.0)
        eigenvectors[:, tied] = eigenvectors[:, tied] @ rotation
        np.testing.assert_allclose(matrix @ eigenvectors, eigenvectors * eigenvalues, atol=1e-12)
        return eigenvalues, eigenvectors

    monkeypatch.setattr(diffusion, "eigh", rotated_eigh)
    rotated = diffusion.compute_diffusion_coordinates(weights, diffusion_time=2, n_components=1)
    np.testing.assert_allclose(pdist(reference), pdist(rotated), atol=1e-12)


def test_diffusion_keeps_a_repeated_leading_group_whole() -> None:
    weights = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    coordinates = diffusion.compute_diffusion_coordinates(weights, diffusion_time=2, n_components=1)
    assert coordinates.shape == (4, 2)
    expected_projector = 0.5 * np.array(
        [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]],
    )
    np.testing.assert_allclose(coordinates @ coordinates.T, expected_projector, atol=1e-12)


def test_diffusion_distinct_cutoff_preserves_existing_coordinates() -> None:
    weights = np.array([[0, 1, 0, 0], [1, 0, 2, 0], [0, 2, 0, 3], [0, 0, 3, 0]])
    degrees = weights.sum(axis=1)
    operator = weights / np.sqrt(degrees[:, None] * degrees[None, :])
    eigenvalues, eigenvectors = np.linalg.eigh(operator)
    expected = eigenvectors[:, -2:-1] * eigenvalues[-2] ** 2
    coordinates = diffusion.compute_diffusion_coordinates(weights, diffusion_time=2, n_components=1)
    assert coordinates.shape == (4, 1)
    np.testing.assert_allclose(pdist(coordinates), pdist(expected), atol=1e-12)


@pytest.mark.parametrize("use_dual", [False, True])
@pytest.mark.parametrize("perturbation", [0.0, 1e-12])
def test_pca_tied_cutoff_is_feature_permutation_invariant(
    use_dual: bool, perturbation: float
) -> None:
    if use_dual:
        centered_identity = np.eye(4) - np.ones((4, 4)) / 4.0
        values = np.column_stack([centered_identity, centered_identity])
        expected_rank = 3
    else:
        values = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]])
        expected_rank = 2
    values[:, 0] *= 1.0 + perturbation
    permutation = np.arange(values.shape[1])[::-1]
    projectors = []
    for order in (np.arange(values.shape[1]), permutation):
        eig = eigendecompose_covariance(values[:, order], compute_eigenvectors=True)
        assert eig.use_dual is use_dual
        basis, eigenvalues = build_pca_projection(
            eig,
            projection_dimension=1,
            n_features_total=values.shape[1],
        )
        assert basis.shape == (expected_rank, values.shape[1])
        assert len(eigenvalues) == expected_rank
        np.testing.assert_allclose(basis @ basis.T, np.eye(expected_rank), atol=1e-12)
        basis = basis[:, np.argsort(order)]
        projectors.append(basis.T @ basis)
    np.testing.assert_allclose(projectors[0], projectors[1], atol=1e-12)


def test_pca_projection_does_not_depend_on_tied_eigenvector_rotation() -> None:
    axes = np.diag(np.sqrt(np.array([16.0, 4.0, 4.0, 1.0])))
    eig = eigendecompose_covariance(np.vstack([axes, -axes]), compute_eigenvectors=True)
    rotated_vectors = eig.eigenvectors_active.copy()
    rotation = np.array([[1.0, -1.0], [1.0, 1.0]]) / np.sqrt(2.0)
    rotated_vectors[:, 1:3] = rotated_vectors[:, 1:3] @ rotation
    original, _ = build_pca_projection(eig, projection_dimension=2, n_features_total=4)
    rotated, _ = build_pca_projection(
        replace(eig, eigenvectors_active=rotated_vectors),
        projection_dimension=2,
        n_features_total=4,
    )
    np.testing.assert_allclose(original.T @ original, rotated.T @ rotated, atol=1e-12)
    np.testing.assert_allclose(original.T @ original, np.diag([1.0, 1.0, 1.0, 0.0]), atol=1e-12)


def test_small_positive_eigenvalues_do_not_become_artificial_wald_ties() -> None:
    expected_eigenvalues = np.array([1e-11, 1e-13, 5e-14])
    axes = np.diag(np.sqrt(3.0 * expected_eigenvalues))
    eig = eigendecompose_covariance(np.vstack([axes, -axes]), compute_eigenvectors=True)
    basis, eigenvalues = build_pca_projection(eig, projection_dimension=3, n_features_total=3)
    np.testing.assert_allclose(eigenvalues, expected_eigenvalues, rtol=1e-12, atol=0)
    selected_basis, selected_eigenvalues = build_pca_projection_basis(
        k=2,
        pca_projection=basis,
        pca_eigenvalues=eigenvalues,
    )
    assert selected_basis.shape == (2, 3)
    np.testing.assert_allclose(selected_eigenvalues, expected_eigenvalues[:2], rtol=1e-12, atol=0)


@pytest.mark.parametrize(
    "eigenvalues,minimum_dimension,expected_dimension",
    [([0.5, 0.5], 1, 2), ([1.0, 1.0, 0.0], 1, 2), ([0.5, 0.5], 0, 0), ([0.5, 0.1], 1, 1)],
)
def test_mp_dimension_floor_keeps_ties_without_adding_null_directions(
    eigenvalues,
    minimum_dimension,
    expected_dimension,
) -> None:
    result = estimate_marchenko_pastur_dimension(
        np.array(eigenvalues),
        n_samples=10,
        n_features=len(eigenvalues),
        minimum_projection_dimension=minimum_dimension,
    )
    assert result.raw_mp_signal_count == 0
    assert result.test_projection_dimension == expected_dimension


def test_node_metadata_matches_the_complete_whitened_eigenspace() -> None:
    values = np.random.default_rng(104).normal(size=(20, 3))
    space = continuous_feature_space_from_columns(["x", "y", "z"])
    task = NodeSpectralTask(
        node_id="root",
        row_indices=tuple(range(20)),
        internal_distributions=(),
        null_distribution=values.mean(axis=0),
        feature_space=space,
        continuous_covariance_by_block={
            space.continuous_blocks[0].name: np.cov(values, rowvar=False)
        },
    )
    result = _process_node(task, values, 1, None, 3, True)
    assert result.raw_mp_signal_count == 0
    assert result.test_projection_dimension == result.projection_matrix.shape[0] == 3
    assert len(result.eigenvalues) == 3
    np.testing.assert_allclose(
        result.projection_matrix.T @ result.projection_matrix, np.eye(3), atol=1e-12
    )


def test_mp_boundary_tie_extends_test_dimension_but_preserves_raw_count() -> None:
    threshold = (1.0 + np.sqrt(3.0 / 100.0)) ** 2
    result = estimate_marchenko_pastur_dimension(
        np.array([threshold + 1e-12, threshold - 1e-12, 0.1]),
        n_samples=100,
        n_features=3,
        minimum_projection_dimension=0,
    )
    assert result.raw_mp_signal_count == 1
    assert result.test_projection_dimension == 2


def test_wald_basis_extends_through_transitively_tied_values() -> None:
    eigenvalues = np.array([2.0, 1.0 + 1.5e-10, 1.0 + 0.8e-10, 1.0])
    basis, selected = build_pca_projection_basis(
        k=2, pca_projection=np.eye(4), pca_eigenvalues=eigenvalues
    )
    assert basis.shape == (4, 4)
    np.testing.assert_array_equal(selected, eigenvalues)


def test_wald_basis_rejects_known_tied_directions_missing_from_parent_basis() -> None:
    with pytest.raises(ValueError, match="PCA projection has"):
        build_pca_projection_basis(k=1, pca_projection=np.eye(2)[:1], pca_eigenvalues=np.ones(2))


@pytest.mark.parametrize("adaptive_fraction", [None, 0.9])
@pytest.mark.parametrize("angle", [np.pi / 4, np.pi / 2])
def test_wald_statistic_rank_and_rejection_are_tied_basis_invariant(
    adaptive_fraction, angle
) -> None:
    z = np.array([3.0, 0.0])
    rotation = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
    outputs = [
        run_projected_wald_kernel(
            z,
            spectral_k=1,
            pca_projection=basis,
            pca_eigenvalues=np.ones(2),
            adaptive_dimension_energy_fraction=adaptive_fraction,
        )
        for basis in (np.eye(2), rotation)
    ]
    for output in outputs:
        assert output.projection_dimension == 2
        assert output.degrees_of_freedom == 2
        np.testing.assert_allclose(output.statistic, 9.0, atol=1e-12)
        np.testing.assert_allclose(output.p_value, chi2.sf(9.0, 2), rtol=1e-12)
        assert (output.p_value <= 0.01) is False


def test_adaptive_wald_selects_complete_blocks_with_distinct_surrounding_values() -> None:
    z = np.array([0.0, 3.0, 0.0, 0.0])
    basis = np.eye(4)
    rotation = basis.copy()
    rotation[1:3, 1:3] = np.array([[1.0, -1.0], [1.0, 1.0]]) / np.sqrt(2.0)
    for candidate in (basis, rotation):
        output = run_projected_wald_kernel(
            z,
            spectral_k=1,
            pca_projection=candidate,
            pca_eigenvalues=np.array([4.0, 1.0, 1.0, 0.25]),
            adaptive_dimension_energy_fraction=0.9,
        )
        assert output.projection_dimension == 3
        np.testing.assert_allclose(output.statistic, 9.0, atol=1e-12)
        np.testing.assert_allclose(output.p_value, chi2.sf(9.0, 3), rtol=1e-12)


def _knn_weights(points: np.ndarray, k: int = 8):
    from scipy import sparse
    from sklearn.neighbors import NearestNeighbors

    distances, indices = NearestNeighbors(n_neighbors=k).fit(points).kneighbors(points)
    rows = np.repeat(np.arange(len(points)), k)
    return sparse.csr_matrix(
        (np.exp(-(distances.ravel() ** 2)), (rows, indices.ravel())),
        shape=(len(points), len(points)),
    )


@pytest.mark.parametrize("offset", [0.0, 1000.0])
def test_sparse_leading_modes_match_dense_distances(monkeypatch, offset: float) -> None:
    rng = np.random.default_rng(3)
    points = rng.normal(size=(300, 4))
    points[150:] += offset  # offset 1000 makes two graph components (dense fallback)
    weights = _knn_weights(points)
    dense = diffusion.compute_diffusion_coordinates(
        weights.toarray(), diffusion_time=2, n_components=12
    )
    monkeypatch.setattr(diffusion, "DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES", 50)

    sparse_path = diffusion.compute_diffusion_coordinates(
        weights, diffusion_time=2, n_components=12
    )

    assert sparse_path.shape == dense.shape
    np.testing.assert_allclose(pdist(sparse_path), pdist(dense), atol=1e-9)
