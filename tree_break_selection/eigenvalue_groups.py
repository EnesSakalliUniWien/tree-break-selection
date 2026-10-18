"""Numerical eigenvalue groups shared by spectral rank selection paths."""

from __future__ import annotations

import numpy as np


def tied_eigenvalue_neighbors(eigenvalues: np.ndarray) -> np.ndarray:
    """Identify adjacent numerical ties in an ordered spectrum.

    Relative tolerance 1e-10 and absolute tolerance 1e-12 times the full
    spectral scale match the adaptive-cosine rule. Adjacent ties connect
    transitively; they describe numerical groups, not exact multiplicities.
    """
    values = np.asarray(eigenvalues, dtype=np.float64)
    scale = float(np.max(np.abs(values))) if len(values) else 0.0
    return np.isclose(values[:-1], values[1:], rtol=1e-10, atol=1e-12 * scale)


def extend_rank_through_ties(eigenvalues: np.ndarray, target_rank: int) -> int:
    """Extend a prefix rank through its last numerical eigenvalue group.

    Zero selects no directions. Callers supply only directions admissible
    under their existing positive-rank/availability constraints.
    """
    rank = max(0, min(int(target_rank), len(eigenvalues)))
    tied = tied_eigenvalue_neighbors(eigenvalues)
    while 0 < rank < len(eigenvalues) and tied[rank - 1]:
        rank += 1
    return rank
