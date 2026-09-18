"""Nuisance-derived radial path for continuous Gaussian selected inference.

This internal prototype isolates one standardized Gaussian contrast from its
contrast-orthogonal nuisance coordinates.  Its PCA projection is learned only
from those nuisance coordinates, so the projection remains fixed while the
projected contrast radius is varied.  It is not connected to the production
gate or traversal defaults.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import linalg
from scipy.integrate import quad
from scipy.special import gammaln, log_ndtr, logsumexp
from scipy.stats import chi

from ...decomposition.backends.eigen.decomposition import eigendecompose_covariance
from ...decomposition.backends.eigen.projection import build_pca_projection


@dataclass(frozen=True)
class GaussianSelectedProjectionPath:
    """One conditional radial path through a continuous Gaussian data matrix."""

    nuisance_matrix: NDArray[np.float64]
    contrast: NDArray[np.float64]
    contrast_norm: float
    covariance_factor: NDArray[np.float64]
    projection: NDArray[np.float64]
    projected_direction: NDArray[np.float64]
    orthogonal_whitened_contrast: NDArray[np.float64]
    observed_radius: float

    @property
    def projection_dimension(self) -> int:
        """Return the dimension of the nuisance-derived projected contrast."""
        return int(self.projection.shape[0])

    @property
    def observed_statistic(self) -> float:
        """Return the Gaussian score/likelihood-ratio statistic ``R^2``."""
        return float(self.observed_radius * self.observed_radius)

    def reconstruct(self, radius: float) -> NDArray[np.float64]:
        """Reconstruct data at ``radius`` with all conditioned coordinates fixed."""
        radius_value = float(radius)
        if not np.isfinite(radius_value) or radius_value < 0.0:
            raise ValueError(
                f"Gaussian selected radial path requires a finite non-negative radius; "
                f"got {radius!r}."
            )
        whitened_contrast = (
            self.projection.T @ (radius_value * self.projected_direction)
            + self.orthogonal_whitened_contrast
        )
        raw_contrast = self.covariance_factor @ whitened_contrast
        return self.nuisance_matrix + np.outer(
            self.contrast / self.contrast_norm,
            raw_contrast,
        )


@dataclass(frozen=True)
class RadialSelectionInterval:
    """One positive-length interval in an exact selected radial region."""

    lower: float
    upper: float

    def __post_init__(self) -> None:
        lower = float(self.lower)
        upper = float(self.upper)
        if not np.isfinite(lower) or lower < 0.0:
            raise ValueError(
                "Radial selection interval lower bound must be finite and non-negative; "
                f"got {self.lower!r}."
            )
        if np.isnan(upper) or upper == -np.inf or upper <= lower:
            raise ValueError(
                "Radial selection interval upper bound must exceed its lower bound; "
                f"got ({self.lower!r}, {self.upper!r})."
            )
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True)
class SelectedGaussianTailResult:
    """Selected radial tail; absolute masses may underflow, but the ratio need not."""

    statistic: float
    degrees_of_freedom: int
    selection_probability: float
    selected_tail_probability: float
    p_value: float
    selection_intervals: tuple[RadialSelectionInterval, ...]


def _canonicalize_selection_intervals(
    intervals: Sequence[RadialSelectionInterval],
) -> tuple[RadialSelectionInterval, ...]:
    if not intervals:
        raise ValueError("Selected Gaussian tail requires at least one selection interval.")
    validated: list[RadialSelectionInterval] = []
    for interval in intervals:
        if not isinstance(interval, RadialSelectionInterval):
            raise TypeError(
                "Selected Gaussian tail intervals must be RadialSelectionInterval objects."
            )
        validated.append(interval)
    ordered = sorted(validated, key=lambda interval: (interval.lower, interval.upper))
    merged = [ordered[0]]
    for interval in ordered[1:]:
        previous = merged[-1]
        if interval.lower <= previous.upper:
            merged[-1] = RadialSelectionInterval(
                previous.lower,
                max(previous.upper, interval.upper),
            )
        else:
            merged.append(interval)
    return tuple(merged)


def _chi_log_survival(radius: float, degrees_of_freedom: int) -> float:
    value = float(chi.logsf(radius, df=degrees_of_freedom))
    if np.isfinite(value):
        return value
    # chi.logsf also underflows in SciPy. For integer dimensions, the upper
    # incomplete-gamma recurrence starts from an exponential or normal tail.
    x = 0.5 * radius * radius
    base_dimension = 1 if degrees_of_freedom % 2 else 2
    value = float(np.log(2.0) + log_ndtr(-radius)) if base_dimension == 1 else -x
    for dimension in range(base_dimension, degrees_of_freedom, 2):
        shape = 0.5 * dimension
        value = float(np.logaddexp(value, shape * np.log(x) - x - gammaln(shape + 1)))
    return value


def _chi_radial_interval_log_probability(
    interval: RadialSelectionInterval,
    *,
    degrees_of_freedom: int,
) -> float:
    if interval.upper == np.inf:
        return _chi_log_survival(interval.lower, degrees_of_freedom)
    if chi.cdf(interval.upper, df=degrees_of_freedom) < 0.5:
        larger = float(chi.logcdf(interval.upper, df=degrees_of_freedom))
        smaller = float(chi.logcdf(interval.lower, df=degrees_of_freedom))
    else:
        larger = _chi_log_survival(interval.lower, degrees_of_freedom)
        smaller = _chi_log_survival(interval.upper, degrees_of_freedom)
    if np.isfinite(larger) and smaller - larger < -1e-5:
        return float(larger + np.log(-np.expm1(smaller - larger)))

    # Integrate narrow intervals directly to avoid subtracting nearly equal
    # tails. Scaling by the maximum density also handles underflowed CDFs.
    width = interval.upper - interval.lower
    mode = float(np.clip(np.sqrt(degrees_of_freedom - 1), interval.lower, interval.upper))
    log_scale = float(chi.logpdf(mode, df=degrees_of_freedom))
    mass, _error = quad(
        lambda fraction: np.exp(
            chi.logpdf(interval.lower + fraction * width, df=degrees_of_freedom) - log_scale
        ),
        0.0, 1.0, epsabs=1e-12, epsrel=1e-12,
    )
    if not np.isfinite(log_scale) or mass <= 0.0:
        raise ValueError("Selected Gaussian interval probability could not be resolved numerically.")
    return float(log_scale + np.log(width) + np.log(mass))


def compute_selected_gaussian_tail(
    *,
    observed_radius: float,
    projection_dimension: int,
    selection_intervals: Sequence[RadialSelectionInterval],
) -> SelectedGaussianTailResult:
    r"""Integrate the chi-radial tail over an exact selected-region interval union.

    Conditional on nuisance coordinates, projected direction, and an exact
    selection event ``S`` along the radial path, ``R`` has a chi distribution
    with ``projection_dimension`` degrees of freedom truncated to ``S``.  This
    function evaluates

    ``P(R >= r_observed, R in S) / P(R in S)``

    deterministically through log-domain chi interval integrals. Absolute
    probability fields may underflow to zero; the p-value is normalized before
    converting back from logarithms.
    Constructing the exact interval union from a hierarchy replay is a separate
    responsibility and is intentionally outside this prototype.
    """
    radius = float(observed_radius)
    if not np.isfinite(radius) or radius < 0.0:
        raise ValueError(
            "Selected Gaussian tail requires a finite non-negative observed_radius; "
            f"got {observed_radius!r}."
        )
    if isinstance(projection_dimension, bool) or not isinstance(
        projection_dimension,
        (int, np.integer),
    ):
        raise ValueError("Selected Gaussian tail projection_dimension must be an integer.")
    degrees_of_freedom = int(projection_dimension)
    if degrees_of_freedom <= 0:
        raise ValueError("Selected Gaussian tail projection_dimension must be positive.")

    intervals = _canonicalize_selection_intervals(selection_intervals)
    if not any(interval.lower <= radius <= interval.upper for interval in intervals):
        raise ValueError(
            "Selected Gaussian tail observed radius lies outside the declared selection event."
        )

    log_selection_probability = float(
        logsumexp([
            _chi_radial_interval_log_probability(
                interval,
                degrees_of_freedom=degrees_of_freedom,
            )
            for interval in intervals
        ])
    )
    if not np.isfinite(log_selection_probability):
        raise ValueError(
            "Selected Gaussian tail requires a declared selection event with positive "
            "chi-radial probability."
        )

    log_tail_masses = []
    for interval in intervals:
        tail_lower = max(interval.lower, radius)
        if tail_lower >= interval.upper:
            continue
        log_tail_masses.append(_chi_radial_interval_log_probability(
            RadialSelectionInterval(tail_lower, interval.upper),
            degrees_of_freedom=degrees_of_freedom,
        ))
    log_tail_probability = min(float(logsumexp(log_tail_masses)), log_selection_probability)
    selection_probability = float(np.exp(log_selection_probability))
    selected_tail_probability = float(np.exp(log_tail_probability))
    p_value = float(np.exp(log_tail_probability - log_selection_probability))
    return SelectedGaussianTailResult(
        statistic=float(radius * radius),
        degrees_of_freedom=degrees_of_freedom,
        selection_probability=selection_probability,
        selected_tail_probability=selected_tail_probability,
        p_value=p_value,
        selection_intervals=intervals,
    )


def _validate_gaussian_path_inputs(
    feature_matrix: np.ndarray,
    contrast: np.ndarray,
    covariance: np.ndarray,
    projection_dimension: int,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    int,
]:
    matrix = np.asarray(feature_matrix, dtype=np.float64)
    if matrix.ndim != 2 or min(matrix.shape, default=0) <= 0:
        raise ValueError(
            "Gaussian selected radial path requires a non-empty two-dimensional "
            f"feature matrix; got shape {matrix.shape}."
        )
    if not np.isfinite(matrix).all():
        raise ValueError("Gaussian selected radial feature matrix must be finite.")

    contrast_vector = np.asarray(contrast, dtype=np.float64)
    if contrast_vector.shape != (matrix.shape[0],):
        raise ValueError(
            "Gaussian selected radial contrast must have one value per matrix row; "
            f"got shape {contrast_vector.shape} for {matrix.shape[0]} rows."
        )
    if not np.isfinite(contrast_vector).all():
        raise ValueError("Gaussian selected radial contrast must be finite.")
    if not np.any(contrast_vector != 0.0):
        raise ValueError("Gaussian selected radial path requires a non-zero contrast.")

    covariance_matrix = np.asarray(covariance, dtype=np.float64)
    expected_covariance_shape = (matrix.shape[1], matrix.shape[1])
    if covariance_matrix.shape != expected_covariance_shape:
        raise ValueError(
            "Gaussian selected radial covariance must be square with the feature width; "
            f"got shape {covariance_matrix.shape}, expected {expected_covariance_shape}."
        )
    if not np.isfinite(covariance_matrix).all():
        raise ValueError("Gaussian selected radial covariance must be finite.")
    if not np.allclose(covariance_matrix, covariance_matrix.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Gaussian selected radial covariance must be symmetric.")

    if isinstance(projection_dimension, bool) or not isinstance(
        projection_dimension,
        (int, np.integer),
    ):
        raise ValueError("Gaussian selected radial projection_dimension must be an integer.")
    dimension = int(projection_dimension)
    if dimension <= 0 or dimension > matrix.shape[1]:
        raise ValueError(
            "Gaussian selected radial projection_dimension must lie between 1 and "
            f"the feature width {matrix.shape[1]}; got {projection_dimension!r}."
        )

    return matrix, contrast_vector, covariance_matrix, dimension


def build_gaussian_selected_projection_path(
    feature_matrix: np.ndarray,
    contrast: np.ndarray,
    covariance: np.ndarray,
    *,
    projection_dimension: int,
) -> GaussianSelectedProjectionPath:
    r"""Build a nuisance-conditioned Gaussian projected-score radial path.

    For row contrast ``eta`` and ``h = ||eta||``, the data are decomposed as
    ``X = X_perp + (eta / h) t^T``.  The known feature covariance is factored
    as ``Sigma = L L^T`` and ``z = L^-1 t``.  PCA rows are learned from the
    whitened ``X_perp`` only.  Holding ``X_perp``, the projected direction, and
    the projection-orthogonal part of ``z`` fixed leaves only the projected
    radius ``R`` free.
    """
    matrix, contrast_vector, covariance_matrix, dimension = _validate_gaussian_path_inputs(
        feature_matrix,
        contrast,
        covariance,
        projection_dimension,
    )
    covariance_matrix = 0.5 * (covariance_matrix + covariance_matrix.T)
    try:
        covariance_factor = np.linalg.cholesky(covariance_matrix)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Gaussian selected radial covariance must be positive-definite."
        ) from exc

    contrast_norm = float(linalg.norm(contrast_vector))
    if not np.isfinite(contrast_norm) or contrast_norm <= 0.0:
        raise ValueError(
            "Gaussian selected radial path requires a contrast with finite non-zero norm."
        )
    unit_contrast = contrast_vector / contrast_norm
    nuisance_matrix = matrix - np.outer(
        unit_contrast,
        unit_contrast @ matrix,
    )
    whitened_nuisance = linalg.solve_triangular(
        covariance_factor,
        nuisance_matrix.T,
        lower=True,
        check_finite=False,
    ).T
    nuisance_eigendecomposition = eigendecompose_covariance(
        whitened_nuisance,
        compute_eigenvectors=True,
    )
    if nuisance_eigendecomposition is None:
        raise ValueError(
            "Gaussian selected radial nuisance PCA requires at least one varying "
            "contrast-orthogonal feature direction."
        )
    try:
        projection, _eigenvalues = build_pca_projection(
            nuisance_eigendecomposition,
            projection_dimension=dimension,
            n_features_total=matrix.shape[1],
        )
    except ValueError as exc:
        raise ValueError(
            "Gaussian selected radial nuisance PCA cannot supply the requested "
            f"projection dimension {dimension}."
        ) from exc
    if projection.shape != (dimension, matrix.shape[1]):
        raise ValueError(
            "Gaussian selected radial nuisance PCA cannot supply the requested "
            f"projection dimension {dimension}; available={projection.shape[0]}."
        )
    if not np.allclose(projection @ projection.T, np.eye(dimension), atol=1e-10):
        raise ValueError("Gaussian selected radial nuisance PCA basis is not orthonormal.")

    raw_contrast = unit_contrast @ matrix
    whitened_contrast = linalg.solve_triangular(
        covariance_factor,
        raw_contrast,
        lower=True,
        check_finite=False,
    )
    projected_contrast = projection @ whitened_contrast
    observed_radius = float(np.linalg.norm(projected_contrast))
    if not np.isfinite(observed_radius) or observed_radius <= 0.0:
        raise ValueError(
            "Gaussian selected radial conditioning requires a non-zero observed "
            "projected radius."
        )
    projected_direction = projected_contrast / observed_radius
    orthogonal_whitened_contrast = (
        whitened_contrast - projection.T @ projected_contrast
    )

    return GaussianSelectedProjectionPath(
        nuisance_matrix=np.asarray(nuisance_matrix, dtype=np.float64),
        contrast=np.asarray(contrast_vector, dtype=np.float64),
        contrast_norm=contrast_norm,
        covariance_factor=np.asarray(covariance_factor, dtype=np.float64),
        projection=np.asarray(projection, dtype=np.float64),
        projected_direction=np.asarray(projected_direction, dtype=np.float64),
        orthogonal_whitened_contrast=np.asarray(
            orthogonal_whitened_contrast,
            dtype=np.float64,
        ),
        observed_radius=observed_radius,
    )


__all__ = [
    "GaussianSelectedProjectionPath",
    "RadialSelectionInterval",
    "SelectedGaussianTailResult",
    "build_gaussian_selected_projection_path",
    "compute_selected_gaussian_tail",
]
