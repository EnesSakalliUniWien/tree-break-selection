"""Conditional energy-selected Gaussian tails for a fixed orthonormal basis.

This law requires an isotropic Gaussian contrast independent of the candidate
basis. It corrects contrast-energy rank selection, not same-data tree/PCA fitting
or estimated covariance. Numerical eigenvalue groups must be retained whole.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import betainc, gammaincc, logsumexp, roots_jacobi


@lru_cache(maxsize=1024)
def _conditional_shares(m: int, k: int, previous: int, q: float, order: int):
    # C_k ~ Beta(k/2, (m-k)/2); given C_k=c, C_previous/c is
    # Beta(previous/2, (k-previous)/2). K=k iff C_previous<q<=C_k.
    a, b = k / 2, (m - k) / 2
    # c=q+(1-q)u^2 also resolves the half-integer endpoint of the
    # conditional Beta CDF at c=q; the Jacobi weight is u(1-u)^(b-1).
    nodes, weights = roots_jacobi(order, b - 1, 1)
    u = (nodes + 1) / 2
    c = q + (1 - q) * u**2
    selection = np.ones_like(c) if previous == 0 else betainc(
        previous / 2, (k - previous) / 2, q / c,
    )
    with np.errstate(divide="ignore"):
        log_weights = np.log(weights) + (b - 1) * np.log1p(u) + (a - 1) * np.log(c) + np.log(selection)
    normalizer = logsumexp(log_weights)
    if not np.isfinite(normalizer):
        raise ValueError("Energy-selection event has numerically unresolved probability.")
    normalized = np.exp(log_weights - normalizer)
    c.setflags(write=False)
    normalized.setflags(write=False)
    return c, normalized


def energy_selected_gaussian_sf(
    statistic: float,
    *,
    candidate_dimension: int,
    projection_dimension: int,
    previous_group_dimension: int,
    energy_fraction: float,
) -> float:
    """P(T>=t | selected complete-group rank), under the declared Gaussian law.

    The radius squared is chi-square ``candidate_dimension`` independently of
    the directional energy shares. Integrating its survival over the selected
    share gives the exact mathematical reference, evaluated by converged
    Gauss-Jacobi quadrature. Retained rank is not an effective chi-square df.
    """
    m, k, previous = candidate_dimension, projection_dimension, previous_group_dimension
    if any(isinstance(v, bool) or int(v) != v for v in (m, k, previous)):
        raise ValueError("Energy-selection dimensions must be integers.")
    m, k, previous = int(m), int(k), int(previous)
    if not 0 <= previous < k <= m:
        raise ValueError("Require 0 <= previous group dimension < retained rank <= candidate rank.")
    q, t = float(energy_fraction), float(statistic)
    if not np.isfinite(q) or not 0 < q <= 1:
        raise ValueError("Energy fraction must be finite and in (0, 1].")
    if not np.isfinite(t) or t < 0:
        raise ValueError("Statistic must be finite and non-negative.")
    if q == 1 and k != m:
        raise ValueError("A partial-rank event at energy fraction one has probability zero.")
    if k == m:
        return float(gammaincc(m / 2, t / 2))
    if t == 0:
        return 1.0
    previous_value = None
    for order in (64, 128, 256, 512):
        shares, weights = _conditional_shares(m, k, previous, q, order)
        value = float(weights @ gammaincc(m / 2, t / (2 * shares)))
        if previous_value is not None and abs(value - previous_value) <= 1e-12 + 1e-9 * value:
            return float(np.clip(value, 0, 1))
        previous_value = value
    raise ValueError("Energy-selected tail quadrature did not converge; reference unavailable.")


__all__ = ["energy_selected_gaussian_sf"]
