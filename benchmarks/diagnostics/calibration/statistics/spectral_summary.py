"""Spectral summaries used only by calibration and validation reports."""

from __future__ import annotations

import numpy as np
from scipy.stats import entropy


def spectral_entropy(eigenvalues: np.ndarray) -> float:
    """Return the Shannon entropy of the normalized positive finite eigenvalues.

    The value is ``nan`` when no eigenvalue is positive, because the spectrum
    then has no distribution to summarize.
    """
    values = np.asarray(eigenvalues, dtype=np.float64)
    positive = values[np.isfinite(values) & (values > 0.0)]
    return float(entropy(positive)) if positive.size else float("nan")


def effective_rank(eigenvalues: np.ndarray) -> float:
    """Return the exponentiated spectral entropy; ``nan`` without positive eigenvalues."""
    return float(np.exp(spectral_entropy(eigenvalues)))


__all__ = ["effective_rank", "spectral_entropy"]
