"""Guarded PSD Toeplitz Vandermonde factorization for diagnostics only.

Uses the rank-deficient shift relation and the full-rank one-atom reduction
in Yang & Xie, arXiv:1605.02431, Theorem 1. Frequencies are not eigenvalues;
atoms are not generally orthogonal eigenvectors. No Toeplitz fitting occurs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import toeplitz


@dataclass(frozen=True)
class VandermondeResult:
    frequencies: np.ndarray
    weights: np.ndarray
    atoms: np.ndarray
    numerical_rank: int
    unique_at_numerical_rank: bool
    relative_reconstruction_error: float
    atom_condition: float


def sinusoidal_atoms(dimension: int, frequencies: np.ndarray) -> np.ndarray:
    """Columns a(f)[j] = exp(2 pi i j f), with unnormalized first entry one."""
    return np.exp(2j * np.pi * np.arange(dimension)[:, None] * frequencies[None, :])


def toeplitz_residual(matrix: np.ndarray) -> float:
    """Relative distance from the matrix fixed by its first column and row."""
    matrix = np.asarray(matrix)
    scale = float(np.linalg.norm(matrix))
    return (
        float(np.linalg.norm(matrix - toeplitz(matrix[:, 0], matrix[0, :])) / scale)
        if scale
        else 0.0
    )


def _rank_deficient_frequencies(matrix: np.ndarray, rank: int) -> np.ndarray:
    if rank == 0:
        return np.empty(0)
    values, vectors = np.linalg.eigh(matrix)
    factor = vectors[:, -rank:] * np.sqrt(values[-rank:])
    if np.linalg.matrix_rank(factor[:-1]) != rank:
        raise ValueError("numerically_unresolved_shift_rank")
    shift = np.linalg.lstsq(factor[:-1], factor[1:], rcond=None)[0]
    roots = np.linalg.eigvals(shift)
    if np.max(np.abs(np.abs(roots) - 1.0)) > 1e-6:
        raise ValueError("numerically_unresolved_unit_circle_roots")
    return np.mod(np.angle(roots) / (2 * np.pi), 1.0)


def decompose_psd_toeplitz(
    matrix: np.ndarray,
    *,
    anchor_frequency: float = 0.0,
    rank_rtol: float = 1e-10,
    structure_rtol: float = 1e-10,
    reconstruction_rtol: float = 1e-7,
) -> VandermondeResult:
    """Return T = A diag(weights) A^H, or fail explicitly on unsupported input.

    Full rank needs a chosen anchor and is nonunique. Rank is numerical under
    the declared tolerance; this function does not certify exact rank or
    distinguish extremely close frequencies below numerical resolution.
    """
    matrix = np.asarray(matrix, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError("expected_nonempty_square_matrix")
    if not np.isfinite(matrix).all():
        raise ValueError("nonfinite_matrix")
    if not np.isfinite(anchor_frequency):
        raise ValueError("nonfinite_anchor")
    if any(not np.isfinite(v) or v <= 0 for v in (rank_rtol, structure_rtol, reconstruction_rtol)):
        raise ValueError("invalid_tolerance")
    n = len(matrix)
    scale = float(np.linalg.norm(matrix))
    if scale == 0:
        return VandermondeResult(
            np.empty(0), np.empty(0), np.empty((n, 0), complex), 0, True, 0.0, 0.0
        )
    normalized = matrix / scale
    if np.linalg.norm(normalized - normalized.conj().T) > structure_rtol:
        raise ValueError("not_hermitian")
    if toeplitz_residual(normalized) > structure_rtol:
        raise ValueError("not_toeplitz")
    values = np.linalg.eigvalsh(normalized)
    tolerance = rank_rtol * float(np.max(np.abs(values)))
    if values[0] < -tolerance:
        raise ValueError("not_positive_semidefinite")
    rank = int(np.count_nonzero(values > tolerance))
    if rank == n:
        anchor = float(anchor_frequency % 1.0)
        atom = sinusoidal_atoms(n, np.array([anchor]))[:, 0]
        weight = 1.0 / float(np.vdot(atom, np.linalg.solve(normalized, atom)).real)
        reduced = normalized - weight * np.outer(atom, atom.conj())
        reduced_values = np.linalg.eigvalsh(reduced)
        if (
            np.count_nonzero(reduced_values > tolerance) != n - 1
            or reduced_values[0] < -10 * tolerance
        ):
            raise ValueError("numerically_unresolved_rank_reduction")
        frequencies = np.append(_rank_deficient_frequencies(reduced, n - 1), anchor)
    else:
        frequencies = _rank_deficient_frequencies(normalized, rank)
    frequencies = np.sort(frequencies)
    atoms = sinusoidal_atoms(n, frequencies)
    # Real nonnegative weights are fitted against every complex matrix entry.
    outer_columns = np.einsum("ik,jk->ijk", atoms, atoms.conj()).reshape(n * n, rank)
    design = np.vstack((outer_columns.real, outer_columns.imag))
    target = np.concatenate((normalized.real.ravel(), normalized.imag.ravel()))
    weights = np.linalg.lstsq(design, target, rcond=None)[0]
    if np.any(weights <= 0):
        raise ValueError("nonpositive_atom_weight")
    reconstruction = (atoms * weights) @ atoms.conj().T
    error = float(np.linalg.norm(normalized - reconstruction))
    if error > reconstruction_rtol:
        raise ValueError(f"reconstruction_failed:{error:.6g}")
    return VandermondeResult(
        frequencies=frequencies,
        weights=weights * scale,
        atoms=atoms,
        numerical_rank=rank,
        unique_at_numerical_rank=rank < n,
        relative_reconstruction_error=error,
        atom_condition=float(np.linalg.cond(atoms)),
    )
