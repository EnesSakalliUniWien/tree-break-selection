"""Deterministic diagnostics for aligned covariance and diffusion subspaces."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, norm, subspace_angles

from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors
from tree_break_selection.space_separation.diffusion import DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES

type Block = tuple[int, int]
SPAN_ATOL = 1e-10


@dataclass(frozen=True, eq=False)
class SubspaceComparison:
    rank_a: int
    rank_b: int
    singular_overlaps: np.ndarray
    principal_angles: np.ndarray
    a_outside_b: float
    b_outside_a: float
    a_contained_in_b: bool
    b_contained_in_a: bool
    projector_frobenius: float
    same_span: bool


@dataclass(frozen=True, eq=False)
class EigenSystem:
    matrix: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    blocks: tuple[Block, ...]


@dataclass(frozen=True, eq=False)
class CovarianceDecomposition:
    total: np.ndarray
    within: np.ndarray
    between: np.ndarray


@dataclass(frozen=True, eq=False)
class ReflectionComparison:
    subspaces: SubspaceComparison
    covariance_frobenius: float | None
    covariance_relative_frobenius: float | None
    covariance_matches: bool | None


@dataclass(frozen=True, eq=False)
class RestrictedBasis:
    basis: np.ndarray
    parent_rank: int
    restricted_rank: int
    rank_loss: int
    singular_values: np.ndarray


@dataclass(frozen=True, eq=False)
class DiffusionSpectrum:
    weights: np.ndarray
    degrees: np.ndarray
    symmetric_operator: np.ndarray
    markov_operator: np.ndarray
    eigenvalues: np.ndarray
    symmetric_vectors: np.ndarray
    right_vectors: np.ndarray
    left_vectors: np.ndarray
    blocks: tuple[Block, ...]
    residuals: dict[str, float]

    def block_alignment(self, block: Block) -> dict[str, SubspaceComparison]:
        """Compare one whole numerical block before and after degree transforms."""
        if block not in self.blocks:
            raise ValueError("block must be a whole numerical eigenvalue block")
        start, stop = block
        symmetric, right, left = (
            basis[:, start:stop]
            for basis in (self.symmetric_vectors, self.right_vectors, self.left_vectors)
        )
        root = np.sqrt(self.degrees)[:, None]
        return {
            "symmetric_vs_right": compare_subspaces(symmetric, right),
            "right_vs_left": compare_subspaces(right, left),
            "symmetric_vs_degree_right": compare_subspaces(symmetric, root * right),
            "symmetric_vs_inverse_degree_left": compare_subspaces(symmetric, left / root),
        }


def _matrix(values) -> np.ndarray:
    x = np.asarray(values)
    if np.iscomplexobj(x):
        raise ValueError("Expected a real matrix")
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[0] == 0 or not np.isfinite(x).all():
        raise ValueError("Expected a finite real matrix with at least one row")
    return x


def _rank_tolerance(shape, rtol: float | None) -> float:
    value = np.finfo(float).eps * max(shape) if rtol is None else float(rtol)
    if not 0 <= value < 1:
        raise ValueError("rtol must be finite and lie in [0, 1)")
    return value


def orthonormal_basis(columns, *, rtol: float | None = None) -> np.ndarray:
    """Orthonormal column span, keeping singular values above rtol times the largest."""
    x = _matrix(columns)
    tolerance = _rank_tolerance(x.shape, rtol)
    scale = np.max(np.abs(x), initial=0.0)
    if scale == 0:
        return np.empty((len(x), 0))
    u, singular, _ = np.linalg.svd(x / scale, full_matrices=False)
    return u[:, singular > singular[0] * tolerance]


def compare_subspaces(a, b, *, rtol: float | None = None) -> SubspaceComparison:
    """Compare column spans in common coordinates. Angles are radians, smallest first."""
    qa, qb = (orthonormal_basis(x, rtol=rtol) for x in (a, b))
    if len(qa) != len(qb):
        raise ValueError("Subspaces must share the same ambient dimension")
    overlap = qa.T @ qb
    angles = subspace_angles(qa, qb)[::-1] if overlap.size else np.empty(0)
    singular = np.clip(np.linalg.svd(overlap, compute_uv=False), 0, 1)
    outside_a = float(np.linalg.norm(qa - qb @ overlap.T))
    outside_b = float(np.linalg.norm(qb - qa @ overlap))
    rank_a, rank_b = qa.shape[1], qb.shape[1]
    a_in_b = rank_a <= rank_b and outside_a <= SPAN_ATOL
    b_in_a = rank_b <= rank_a and outside_b <= SPAN_ATOL
    return SubspaceComparison(
        rank_a,
        rank_b,
        singular,
        angles,
        outside_a,
        outside_b,
        a_in_b,
        b_in_a,
        float(np.hypot(outside_a, outside_b)),
        a_in_b and b_in_a,
    )


def _centered(values, ddof: int) -> tuple[np.ndarray, int]:
    x = _matrix(values)
    if isinstance(ddof, (bool, np.bool_)) or ddof not in (0, 1) or len(x) <= ddof or not x.shape[1]:
        raise ValueError("Covariance needs features, ddof 0 or 1, and more rows than ddof")
    with np.errstate(over="raise", invalid="raise"):
        return x - x.mean(axis=0), len(x) - ddof


def _symmetric_spectrum(matrix: np.ndarray) -> EigenSystem:
    eigenvalues, eigenvectors = eigh(matrix)
    eigenvalues, eigenvectors = eigenvalues[::-1], eigenvectors[:, ::-1]
    ends = [
        *(np.flatnonzero(~tied_eigenvalue_neighbors(eigenvalues)) + 1).tolist(),
        len(eigenvalues),
    ]
    return EigenSystem(
        matrix, eigenvalues, eigenvectors, tuple(zip([0, *ends[:-1]], ends, strict=True))
    )


def covariance_spectrum(values, *, ddof: int = 1) -> EigenSystem:
    """Whole numerical eigenspaces of raw covariance in the supplied feature coordinates."""
    centered, divisor = _centered(values, ddof)
    with np.errstate(over="raise", invalid="raise"):
        scaled = centered / np.sqrt(divisor)
        return _symmetric_spectrum(scaled.T @ scaled)


def covariance_decomposition(values, groups, *, ddof: int = 1) -> CovarianceDecomposition:
    """Parent covariance = within-group scatter + between-group scatter, using one divisor."""
    centered, divisor = _centered(values, ddof)
    labels = np.asarray(groups)
    if labels.shape != (len(centered),):
        raise ValueError("groups must supply one label per row")
    _, indices = np.unique(labels, return_inverse=True)
    with np.errstate(over="raise", invalid="raise"):
        within = np.zeros((centered.shape[1], centered.shape[1]))
        between = within.copy()
        for index in np.unique(indices):
            child = centered[indices == index]
            shift = child.mean(axis=0)
            residual = (child - shift) / np.sqrt(divisor)
            within += residual.T @ residual
            weighted_shift = np.sqrt(len(child) / divisor) * shift
            between += np.outer(weighted_shift, weighted_shift)
        scaled = centered / np.sqrt(divisor)
        return CovarianceDecomposition(scaled.T @ scaled, within, between)


def compare_reflection(
    a, b, reflection, *, covariance_a=None, covariance_b=None, rtol: float | None = None
) -> ReflectionComparison:
    """Compare a prescribed orthogonal involution applied to A against B."""
    a, q = _matrix(a), _matrix(reflection)
    identity = np.eye(len(a))
    if (
        q.shape != identity.shape
        or not np.allclose(q, q.T, rtol=0, atol=SPAN_ATOL)
        or not np.allclose(q.T @ q, identity, rtol=0, atol=SPAN_ATOL)
    ):
        raise ValueError("reflection must be an orthogonal involution in the common coordinates")
    subspaces = compare_subspaces(q @ orthonormal_basis(a, rtol=rtol), b, rtol=rtol)
    if covariance_a is None and covariance_b is None:
        return ReflectionComparison(subspaces, None, None, None)
    if covariance_a is None or covariance_b is None:
        raise ValueError("Supply both covariance matrices or neither")
    ca, cb = _matrix(covariance_a), _matrix(covariance_b)
    if any(c.shape != q.shape or not np.allclose(c, c.T, rtol=SPAN_ATOL, atol=0) for c in (ca, cb)):
        raise ValueError("Covariance matrices must be symmetric in the common coordinates")
    scale = max(np.max(np.abs(ca)), np.max(np.abs(cb))) or 1.0
    ca_scaled, cb_scaled = ca / scale, cb / scale
    error = float(np.linalg.norm(q @ ca_scaled @ q.T - cb_scaled))
    reference = float(max(np.linalg.norm(ca_scaled), np.linalg.norm(cb_scaled))) or 1.0
    relative = error / reference
    return ReflectionComparison(subspaces, error * float(scale), relative, relative <= SPAN_ATOL)


def reversible_diffusion_spectrum(weights) -> DiffusionSpectrum:
    """Dense, unclipped reversible spectrum using production symmetrization conventions."""
    if sparse.issparse(weights):
        if (
            weights.shape[0] != weights.shape[1]
            or weights.shape[0] > DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES
        ):
            raise ValueError("Diffusion weights exceed the dense diagnostic size limit")
        weights = weights.toarray()
    w = _matrix(weights)
    n = len(w)
    if w.shape != (n, n) or not 2 <= n <= DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES:
        raise ValueError("Expected a square diffusion weight matrix within the dense size limit")
    if np.any(w < 0):
        raise ValueError("Diffusion weights must be nonnegative")
    w = np.maximum(w, w.T)
    np.fill_diagonal(w, 0)
    with np.errstate(over="raise", invalid="raise"):
        degrees = w.sum(axis=1)
        if np.any(degrees <= 0):
            raise ValueError(
                "Diffusion diagnostic requires positive degrees after removing the diagonal"
            )
        root = np.sqrt(degrees)
        symmetric = (w / root[:, None]) / root
        # Equal diagonal scaling is symmetric mathematically. Remove rounding asymmetry.
        symmetric = (symmetric + symmetric.T) / 2
        markov = w / degrees[:, None]
        eigen = _symmetric_spectrum(symmetric)
        u, values = eigen.eigenvectors, eigen.eigenvalues
        right, left = u / root[:, None], u * root[:, None]
        residuals = {
            "symmetric_eigen": float(np.linalg.norm(symmetric @ u - u * values)),
            "markov_right_eigen": float(
                norm((markov @ right - right * values).ravel()) / norm(right.ravel())
            ),
            "markov_left_eigen": float(
                norm((markov.T @ left - left * values).ravel()) / norm(left.ravel())
            ),
            "biorthogonality": float(np.linalg.norm(left.T @ right - np.eye(n))),
            "degree_similarity": float(np.linalg.norm(root[:, None] * markov / root - symmetric)),
        }
    return DiffusionSpectrum(
        w, degrees, symmetric, markov, values, u, right, left, eigen.blocks, residuals
    )


def restrict_parent_basis(
    parent_basis, child_rows, *, rtol: float | None = None
) -> RestrictedBasis:
    """Restrict parent rows in the given child order and reorthogonalize the surviving span."""
    parent = orthonormal_basis(parent_basis, rtol=rtol)
    indices = tuple(child_rows)
    if (
        not indices
        or any(
            isinstance(i, (bool, np.bool_))
            or not isinstance(i, Integral)
            or not 0 <= i < len(parent)
            for i in indices
        )
        or len(set(indices)) != len(indices)
    ):
        raise ValueError("child_rows must contain distinct valid row indices")
    restricted = parent[np.asarray(indices, dtype=np.intp)]
    u, singular, _ = np.linalg.svd(restricted, full_matrices=False)
    # Restriction is a contraction of a unit-norm parent basis. Reference its
    # scale, so roundoff leakage from a vanishing restriction stays rank zero.
    basis = u[:, singular > _rank_tolerance(parent.shape, rtol)]
    parent_rank, restricted_rank = parent.shape[1], basis.shape[1]
    return RestrictedBasis(
        basis, parent_rank, restricted_rank, parent_rank - restricted_rank, singular
    )


def compare_parent_child(
    parent_basis, child_basis, child_rows, *, rtol: float | None = None
) -> SubspaceComparison:
    """Compare restricted parent and supplied child spans, not their graph operators."""
    restricted = restrict_parent_basis(parent_basis, child_rows, rtol=rtol)
    return compare_subspaces(restricted.basis, child_basis, rtol=rtol)
