"""Analytic selected root for independent Gaussian rows with covariance sigma² I.

Diagnostic only: full-dimensional contrast, average squared-Euclidean linkage,
full observed merge-sequence conditioning. No PCA, estimated covariance, edge
opening, multiplicity or traversal. Observed ties are unsupported. The conditional
null is eta.T @ E[X] = 0; iid rows with a common mean satisfy it.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.integrate import quad
from scipy.spatial.distance import pdist
from scipy.special import logsumexp
from scipy.stats import chi


@dataclass(frozen=True)
class GaussianRootResult:
    radius: float
    intervals: tuple[tuple[float, float], ...]
    p_value: float
    log_selection_probability: float
    contrast: np.ndarray
    offset: np.ndarray
    slope: np.ndarray
    linkage_matrix: np.ndarray


def _quadratic_region(coefficients):
    """Open intervals where an exact-rational quadratic is nonpositive."""
    c, b, a = coefficients
    with localcontext() as context:
        context.prec = 80

        def decimal(value):
            return Decimal(value.numerator) / Decimal(value.denominator)

        roots = []
        if a:
            discriminant = b * b - 4 * a * c
            if discriminant >= 0:
                root = decimal(discriminant).sqrt()
                roots = [
                    (-decimal(b) - root) / (2 * decimal(a)),
                    (-decimal(b) + root) / (2 * decimal(a)),
                ]
        elif b:
            roots = [-decimal(c) / decimal(b)]
        boundaries = [Decimal(0)] + sorted(set(r for r in roots if r > 0))
        boundaries.append(Decimal("Infinity"))
        intervals = []
        for lo, hi in zip(boundaries[:-1], boundaries[1:]):
            point = (lo + hi) / 2 if hi.is_finite() else lo + max(Decimal(1), abs(lo))
            if decimal(a) * point * point + decimal(b) * point + decimal(c) <= 0:
                left, right = float(lo), float(hi)
                if not left < right:
                    raise ValueError("Selected interval boundaries collapse in float64.")
                intervals.append((left, right))
        return intervals


def _intersect(left, right):
    return [(max(a, c), min(b, d)) for a, b in left for c, d in right if max(a, c) < min(b, d)]


def _log_chi_mass(lo, hi, dimension):
    if lo >= hi:
        return -np.inf
    # Use the tail away from cancellation; fall back to scaled density quadrature.
    for log_large, log_small in [
        (chi.logsf(lo, dimension), chi.logsf(hi, dimension)),
        (chi.logcdf(hi, dimension), chi.logcdf(lo, dimension)),
    ]:
        if np.isfinite(log_large) and log_small < log_large - 1e-7:
            return float(log_large + np.log(-np.expm1(log_small - log_large)))
    mode = max(lo, min(hi, np.sqrt(max(0, dimension - 1))))
    scale = chi.logpdf(mode, dimension)
    mass, error = quad(
        lambda r: np.exp(chi.logpdf(r, dimension) - scale), lo, hi, epsabs=1e-12, epsrel=1e-10
    )
    if not np.isfinite(scale) or not mass > 0 or error > max(1e-12, mass * 1e-7):
        raise ValueError("Selected chi interval mass cannot be resolved numerically.")
    return float(scale + np.log(mass))


def merge_sequence_region(offset, weights, direction, observed_values, linkage_matrix):
    """Affine-path quadratic interval union for a fixed average-linkage sequence.

    Pairwise differences are offset_i-offset_j + r*(weights_i-weights_j)*direction.
    observed_values supplies the observed tie check; a common positive scaling of
    all affine-path distances leaves the merge comparisons unchanged.
    """
    n = len(offset)
    z = linkage_matrix
    clusters = {i: [i] for i in range(n)}
    for step, row in enumerate(z):
        first, second = map(int, row[:2])
        clusters[n + step] = clusters[first] + clusters[second]
    # Shared scalar weights retain exact cancellation of equal slope energies.
    exact_offset = [[Fraction(float(value)) for value in row] for row in offset]
    exact_direction = [Fraction(float(value)) for value in direction]
    weights = [Fraction(value) for value in weights]
    distances = {}
    observed_distances = {}
    for i, j in combinations(range(n), 2):
        delta = [a - b for a, b in zip(exact_offset[i], exact_offset[j])]
        weight = weights[i] - weights[j]
        distances[i, j] = (
            sum(v * v for v in delta),
            2 * weight * sum(a * b for a, b in zip(delta, exact_direction)),
            weight * weight * sum(v * v for v in exact_direction),
        )
        observed_distances[i, j] = float(np.sum((observed_values[i] - observed_values[j]) ** 2))
    region = [(0.0, np.inf)]
    active = set(range(n))
    for step, row in enumerate(z):
        selected = tuple(sorted(map(int, row[:2])))
        chosen = distances[selected]
        for pair in combinations(sorted(active), 2):
            if pair == selected:
                continue
            if observed_distances[pair] == observed_distances[selected]:
                raise ValueError("Observed minimum linkage tie is unsupported.")
            region = _intersect(
                region, _quadratic_region(tuple(a - b for a, b in zip(chosen, distances[pair])))
            )
        first, second = selected
        merged = n + step
        for other in active - {first, second}:
            first_key, second_key = tuple(sorted((first, other))), tuple(sorted((second, other)))
            total = len(clusters[first]) + len(clusters[second])
            distances[tuple(sorted((merged, other)))] = tuple(
                (len(clusters[first]) * a + len(clusters[second]) * b) / total
                for a, b in zip(distances[first_key], distances[second_key])
            )
            observed_distances[tuple(sorted((merged, other)))] = (
                len(clusters[first]) * observed_distances[first_key]
                + len(clusters[second]) * observed_distances[second_key]
            ) / total
        active.difference_update(selected)
        active.add(merged)
    return tuple(region)


def selected_gaussian_root_test(values, *, sigma):
    """Condition on nuisance, contrast direction and the complete merge sequence.

    Under the declared Gaussian null, the conditional radius has chi_d density
    restricted to the interval union returned here. Linkage is recomputed here
    to fix squared-Euclidean average linkage rather than trusting an external tree.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or min(x.shape) < 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("values must be a finite n-by-d matrix with n >= 2 and d >= 1.")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be a known positive finite scalar.")
    n, dimension = x.shape
    z = linkage(pdist(x, metric="sqeuclidean"), method="average")
    clusters = {i: [i] for i in range(n)}
    for step, row in enumerate(z):
        left, right = map(int, row[:2])
        clusters[n + step] = clusters[left] + clusters[right]
    left, right = map(int, z[-1, :2])
    eta = np.zeros(n)
    eta[clusters[left]] = 1 / len(clusters[left])
    eta[clusters[right]] = -1 / len(clusters[right])
    norm = np.linalg.norm(eta)
    contrast = x.T @ eta / norm
    radius = np.linalg.norm(contrast) / sigma
    if not radius > 0:
        raise ValueError("Observed contrast direction is undefined at zero radius.")
    offset = x - np.outer(eta / norm, contrast)
    direction = sigma * contrast / np.linalg.norm(contrast) / norm
    slope = np.outer(eta, direction)

    weights = [
        Fraction(1, len(clusters[left]))
        if i in clusters[left]
        else -Fraction(1, len(clusters[right]))
        for i in range(n)
    ]
    region = merge_sequence_region(offset, weights, direction, x, z)
    if not any(lo < radius < hi for lo, hi in region):
        raise ValueError("Observed radius lies outside the resolved selection region.")
    denominator = logsumexp([_log_chi_mass(lo, hi, dimension) for lo, hi in region])
    numerator = logsumexp(
        [_log_chi_mass(max(lo, radius), hi, dimension) for lo, hi in region if hi > radius]
    )
    if not np.isfinite(denominator) or not np.isfinite(numerator):
        raise ValueError("Selected chi probability cannot be resolved numerically.")
    return GaussianRootResult(
        float(radius),
        tuple(region),
        float(np.exp(min(0.0, numerator - denominator))),
        float(denominator),
        eta,
        offset,
        slope,
        z,
    )
