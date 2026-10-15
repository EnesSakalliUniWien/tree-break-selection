import numpy as np
import pytest

from benchmarks.diagnostics.spectral.vandermonde import (
    decompose_psd_toeplitz,
    sinusoidal_atoms,
)


@pytest.mark.parametrize("frequencies", [[0.12, 0.37], [0.125, 0.875], [0.0]])
def test_recovers_identifiable_atoms_and_weights(frequencies):
    frequencies = np.array(frequencies)
    weights = np.arange(1, len(frequencies) + 1) / 3
    atoms = sinusoidal_atoms(8, frequencies)
    covariance = (atoms * weights) @ atoms.conj().T
    result = decompose_psd_toeplitz(covariance)
    np.testing.assert_allclose(
        np.exp(2j * np.pi * result.frequencies), np.exp(2j * np.pi * frequencies), atol=1e-10
    )
    np.testing.assert_allclose(result.weights, weights, atol=1e-10)
    assert result.unique_at_numerical_rank
    assert result.relative_reconstruction_error < 1e-12


@pytest.mark.parametrize("rho", [0.0, 0.6])
def test_full_rank_anchor_exposes_nonuniqueness_without_changing_covariance(rho):
    covariance = rho ** np.abs(np.arange(8)[:, None] - np.arange(8))
    first = decompose_psd_toeplitz(covariance, anchor_frequency=0)
    second = decompose_psd_toeplitz(covariance, anchor_frequency=0.137)
    for result in (first, second):
        assert result.numerical_rank == 8
        assert not result.unique_at_numerical_rank
        np.testing.assert_allclose(
            (result.atoms * result.weights) @ result.atoms.conj().T, covariance, atol=1e-10
        )
    assert not np.allclose(first.frequencies, second.frequencies)


@pytest.mark.parametrize(
    "matrix,reason",
    [
        (np.diag([1.0, 2.0, 3.0]), "not_toeplitz"),
        (np.array([[1.0, 2.0], [2.0, 1.0]]), "not_positive_semidefinite"),
        (np.array([[1.0, 1.0], [0.0, 1.0]]), "not_hermitian"),
        (np.array([[np.nan]]), "nonfinite_matrix"),
    ],
)
def test_rejects_unsupported_covariances_without_fitting(matrix, reason):
    with pytest.raises(ValueError, match=reason):
        decompose_psd_toeplitz(matrix)


def test_zero_and_scalar_and_scale_invariance():
    assert decompose_psd_toeplitz(np.zeros((3, 3))).numerical_rank == 0
    scalar = decompose_psd_toeplitz(np.array([[3.0]]))
    np.testing.assert_allclose(scalar.weights, [3.0])
    assert not scalar.unique_at_numerical_rank
    atoms = sinusoidal_atoms(8, np.array([0.17, 0.63]))
    covariance = (atoms * [1.0, 2.0]) @ atoms.conj().T
    for scale in (1e-15, 1e15):
        result = decompose_psd_toeplitz(scale * covariance)
        np.testing.assert_allclose(result.weights / scale, [1.0, 2.0], atol=1e-10)
