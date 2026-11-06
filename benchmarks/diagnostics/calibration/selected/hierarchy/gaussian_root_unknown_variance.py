"""Scalar root diagnostic under iid N(mu, sigma²), both parameters unknown.

This is the global common-mean null, stronger than equality of two cluster
average means with arbitrary within-cluster means. Condition on sample mean,
centered norm, residual direction, contrast sign, and the full observed merge
sequence of average squared-Euclidean linkage. No production gates or PCA.
See Yun and Barber (2023), https://arxiv.org/abs/2301.12999, Section 3 (K=2).
"""

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.integrate import quad
from scipy.spatial.distance import pdist
from scipy.special import betaln, logsumexp
from scipy.stats import t

from .gaussian_root_region import merge_sequence_region


@dataclass(frozen=True)
class ScalarGaussianRootResult:
    ratio: float
    statistic: float
    residual_degrees_of_freedom: int
    intervals: tuple[tuple[float, float], ...]
    p_value: float
    log_selection_probability: float
    mean: float
    centered_norm: float
    offset: np.ndarray
    slope: np.ndarray
    linkage_matrix: np.ndarray

    def reconstruct(self, ratio):
        """Return the constant-centered-norm path at h=|contrast|/residual norm."""
        return self.mean + self.centered_norm * (self.offset + self.slope * ratio) / np.hypot(
            1, ratio
        )


def _log_half_t_mass(lo, hi, degrees):
    """Log mass for h=sqrt(F(1, degrees)/degrees), including narrow intervals."""
    if lo >= hi:
        return -np.inf
    root_df = np.sqrt(degrees)
    large, small = t.logsf(root_df * lo, degrees), t.logsf(root_df * hi, degrees)
    if np.isfinite(large) and small < large - 1e-7:
        return float(np.log(2) + large + np.log(-np.expm1(small - large)))
    if not np.isfinite(hi):
        raise ValueError("Selected F interval mass cannot be resolved numerically.")

    def logdensity(h):
        log_h = -np.inf if h == 0 else np.log(h)
        return np.log(2) - betaln(0.5, degrees / 2) - (degrees + 1) / 2 * np.logaddexp(0, 2 * log_h)

    scale = logdensity(lo)
    mass, error = quad(
        lambda u: np.exp(logdensity(lo + (hi - lo) * u) - scale), 0, 1, epsabs=1e-12, epsrel=1e-10
    )
    if not mass > 0 or error > max(1e-12, mass * 1e-7):
        raise ValueError("Selected F interval mass cannot be resolved numerically.")
    return float(scale + np.log(hi - lo) + np.log(mass))


def selected_scalar_gaussian_root_test(values):
    """Return the truncated F(1,n-2) upper tail for a selected scalar root.

    For root contrast e=v/||v||, a=e'x and w=(I-11'/n-ee')x, condition
    on sign(a), w/||w|| and R=||x-mean(x)||. Then h=|a|/||w|| has
    density 2(1+h²)^(-(n-1)/2)/Beta(1/2,(n-2)/2). The path is
    mean(x)1 + R(sign(a)e*h+w/||w||)/sqrt(1+h²). Every squared distance
    has the same R²/(1+h²) factor, so merge comparisons are quadratic in h.
    Zero residual variance and observed minimum-linkage ties are unsupported.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("values must contain at least three finite scalar observations.")
    n = len(x)
    z = linkage(pdist(x[:, None], metric="sqeuclidean"), method="average")
    clusters = {i: [i] for i in range(n)}
    for step, row in enumerate(z):
        first, second = map(int, row[:2])
        clusters[n + step] = clusters[first] + clusters[second]
    left, right = (clusters[int(i)] for i in z[-1, :2])
    weights = [Fraction(1, len(left)) if i in left else -Fraction(1, len(right)) for i in range(n)]
    v = np.asarray(weights, dtype=float)
    norm = np.linalg.norm(v)
    a = (x[left].mean() - x[right].mean()) / norm
    residual = x.copy()
    for group in (left, right):
        residual[group] -= x[group].mean()
    residual_norm = np.linalg.norm(residual)
    if not np.isfinite(residual_norm) or residual_norm <= 0:
        raise ValueError("Positive finite within-child residual norm is required.")
    if not np.isfinite(a) or a == 0:
        raise ValueError("Nonzero finite root contrast is required.")
    h = abs(a) / residual_norm
    offset = residual / residual_norm
    direction = np.array([np.sign(a) / norm])
    slope = v * direction[0]
    intervals = merge_sequence_region(offset[:, None], weights, direction, x[:, None], z)
    if not np.isfinite(h) or not any(lo < h < hi for lo, hi in intervals):
        raise ValueError("Observed ratio lies outside the resolved selection region.")
    degrees = n - 2
    denominator = logsumexp([_log_half_t_mass(lo, hi, degrees) for lo, hi in intervals])
    numerator = logsumexp(
        [_log_half_t_mass(max(lo, h), hi, degrees) for lo, hi in intervals if hi > h]
    )
    if not np.isfinite(denominator) or not np.isfinite(numerator):
        raise ValueError("Selected F probability cannot be resolved numerically.")
    return ScalarGaussianRootResult(
        float(h),
        float(degrees * h * h),
        degrees,
        intervals,
        float(np.exp(min(0, numerator - denominator))),
        float(denominator),
        float(x.mean()),
        float(np.linalg.norm(x - x.mean())),
        offset,
        slope,
        z,
    )
