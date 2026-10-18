"""Sibling pair record dataclass for the pair_testing package."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SiblingPairRecord:
    """Raw per-parent sibling-test record used by calibration pipelines."""

    parent: object
    left: object
    right: object
    stat: float
    reference_scale: float
    """Scale in the projected quadratic reference law; unit for orthonormal PCA."""
    degrees_of_freedom: float
    p_value: float
    branch_length_sum: float
    n_parent: int
    is_null_like: bool
    is_edge_blocked: bool = False
    sibling_null_weight: float = 0.0
    """Weight that this sibling pair represents empirical-null structure."""
    calibration_dependency_group: object | None = None
    """Stopping event that owns this record's calibration support."""
    calibration_dependency_weight: float | None = None
    """Support weight of the owning stopping event."""
    sibling_projection_dimension: float = 0.0
    """Projection dimension used by the sibling test and inflation context."""
    parent_positive_eigenvalue_count: float = 0.0
    """Number of positive parent eigenvalues after finite/non-positive filtering."""
    feature_family: str = "bernoulli"
    """Feature-space family label for the sibling contrast covariance model."""

    def __post_init__(self) -> None:
        if self.is_edge_blocked and not self.is_null_like:
            raise ValueError(
                "SiblingPairRecord requires is_edge_blocked implies is_null_like; "
                f"parent={self.parent!r}."
            )
        if self.calibration_dependency_weight is not None and (
            not np.isfinite(self.calibration_dependency_weight)
            or self.calibration_dependency_weight <= 0.0
        ):
            raise ValueError(
                "SiblingPairRecord calibration_dependency_weight must be finite and positive "
                f"when provided; parent={self.parent!r}."
            )

    @property
    def has_null_role(self) -> bool:
        """Return whether this record has a null-like or edge-blocked role."""
        return bool(self.is_null_like or self.is_edge_blocked)

    @property
    def has_internal_empirical_null_support(self) -> bool:
        """Return whether this record can calibrate the empirical null.

        It needs positive degrees of freedom, positive null weight and a null role.
        """
        return bool(
            self.degrees_of_freedom > 0.0 and self.sibling_null_weight > 0.0 and self.has_null_role
        )


__all__ = ["SiblingPairRecord"]
