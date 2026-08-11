"""Empirical-null inflation model result type."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Literal

import numpy as np

CalibrationDecisionStatus = Literal[
    "internal_admissible",
    "independent_exact_admissible",
    "undefined_no_internal_support",
    "undefined_no_family_support",
    "undefined_sparse_context",
    "undefined_unvalidated_reference_law",
]

CalibrationReferenceLaw = Literal[
    "unresolved_same_selected_hierarchy",
    "exact_independent_unweighted_common_scale_f",
]


@dataclass(frozen=True)
class CalibrationDecision:
    """Focal sibling calibration decision and support evidence."""

    status: CalibrationDecisionStatus
    c_hat: float | None
    p_value: float | None
    estimator: str
    support: dict[str, float | int | str | bool]
    exact_context: dict[str, object]
    descriptive_strata: dict[str, object]


@dataclass(frozen=True)
class CalibrationSupportThresholds:
    """Validation thresholds for the internal calibration support contract."""

    min_supported_groups: int = 30
    min_family_supported_groups: int = 15
    min_family_effective_sample_size: float = 10.0
    min_local_effective_sample_size: float = 8.0
    max_weight_share: float = 0.25
    max_leave_one_group_delta_log_c: float = float(np.log(1.25))

    def __post_init__(self) -> None:
        count_thresholds = (
            "min_supported_groups",
            "min_family_supported_groups",
        )
        for field_name in count_thresholds:
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError(
                    f"{field_name} must be an integer, excluding booleans, and at least 1; "
                    f"got {value!r}."
                )

        effective_sample_size_thresholds = (
            "min_family_effective_sample_size",
            "min_local_effective_sample_size",
        )
        for field_name in effective_sample_size_thresholds:
            value = getattr(self, field_name)
            if not np.isfinite(value) or value < 1.0:
                raise ValueError(f"{field_name} must be finite and at least 1; got {value!r}.")

        if not np.isfinite(self.max_weight_share) or not 0.0 < self.max_weight_share <= 1.0:
            raise ValueError(
                "max_weight_share must be finite and lie in (0, 1]; "
                f"got {self.max_weight_share!r}."
            )
        if (
            not np.isfinite(self.max_leave_one_group_delta_log_c)
            or self.max_leave_one_group_delta_log_c < 0.0
        ):
            raise ValueError(
                "max_leave_one_group_delta_log_c must be finite and non-negative; "
                f"got {self.max_leave_one_group_delta_log_c!r}."
            )


DEFAULT_INTERNAL_SUPPORT_THRESHOLDS = CalibrationSupportThresholds()
INTERNAL_SUPPORT_CONTRACT_VERSION = 2


@dataclass(frozen=True)
class CalibrationSupportPolicySnapshot:
    """Immutable support thresholds paired with their interpretation version."""

    support_contract_version: int
    thresholds: CalibrationSupportThresholds


DEFAULT_INTERNAL_SUPPORT_POLICY = CalibrationSupportPolicySnapshot(
    support_contract_version=INTERNAL_SUPPORT_CONTRACT_VERSION,
    thresholds=DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
)


def _read_only_array(values: np.ndarray, *, dtype: type) -> np.ndarray:
    array = np.array(values, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class CalibrationSample:
    """One authoritative, aligned empirical-null calibration sample."""

    contexts: np.ndarray
    feature_families: tuple[str, ...]
    weights: np.ndarray
    parent_ids: tuple[object, ...]
    dependency_group_ids: tuple[object, ...]
    dependency_group_weights: np.ndarray
    is_strict_null: np.ndarray
    is_edge_blocked: np.ndarray
    statistics: np.ndarray
    reference_scales: np.ndarray
    degrees_of_freedom: np.ndarray

    def __post_init__(self) -> None:
        contexts = _read_only_array(self.contexts, dtype=float)
        vectors = {
            "weights": _read_only_array(self.weights, dtype=float),
            "dependency_group_weights": _read_only_array(
                self.dependency_group_weights,
                dtype=float,
            ),
            "is_strict_null": _read_only_array(self.is_strict_null, dtype=bool),
            "is_edge_blocked": _read_only_array(self.is_edge_blocked, dtype=bool),
            "statistics": _read_only_array(self.statistics, dtype=float),
            "reference_scales": _read_only_array(self.reference_scales, dtype=float),
            "degrees_of_freedom": _read_only_array(self.degrees_of_freedom, dtype=float),
        }
        n_records = len(self.parent_ids)
        if contexts.ndim != 2 or contexts.shape != (n_records, 2):
            raise ValueError(
                "CalibrationSample.contexts must have shape (n_records, 2); "
                f"got {contexts.shape!r} for {n_records} records."
            )
        for field_name, values in vectors.items():
            if values.ndim != 1 or values.shape[0] != n_records:
                raise ValueError(
                    f"CalibrationSample.{field_name} must be one-dimensional and aligned "
                    f"to parent_ids; got shape={values.shape!r}, n_records={n_records}."
                )
        if len(self.feature_families) != n_records:
            raise ValueError("CalibrationSample.feature_families must be aligned to parent_ids.")
        if len(self.dependency_group_ids) != n_records:
            raise ValueError("CalibrationSample.dependency_group_ids must be aligned to parent_ids.")
        try:
            unique_parent_count = len(set(self.parent_ids))
            set(self.dependency_group_ids)
        except TypeError as exc:
            raise ValueError(
                "CalibrationSample parent_ids and dependency_group_ids must be hashable."
            ) from exc
        if unique_parent_count != n_records:
            raise ValueError("CalibrationSample.parent_ids must be unique.")
        if not np.all(np.isfinite(vectors["weights"])) or np.any(vectors["weights"] <= 0.0):
            raise ValueError("CalibrationSample.weights must be finite and positive.")
        if not np.all(np.isfinite(vectors["dependency_group_weights"])) or np.any(
            vectors["dependency_group_weights"] <= 0.0
        ):
            raise ValueError(
                "CalibrationSample.dependency_group_weights must be finite and positive."
            )
        group_weights: dict[object, float] = {}
        for group_id, weight in zip(
            self.dependency_group_ids,
            vectors["dependency_group_weights"],
            strict=True,
        ):
            numeric_weight = float(weight)
            if group_id in group_weights and not np.isclose(
                group_weights[group_id],
                numeric_weight,
                rtol=1e-12,
                atol=1e-15,
            ):
                raise ValueError(
                    "CalibrationSample dependency-group weights must agree within each group; "
                    f"group={group_id!r}."
                )
            group_weights[group_id] = numeric_weight

        object.__setattr__(self, "contexts", contexts)
        object.__setattr__(self, "feature_families", tuple(self.feature_families))
        object.__setattr__(self, "parent_ids", tuple(self.parent_ids))
        object.__setattr__(self, "dependency_group_ids", tuple(self.dependency_group_ids))
        for field_name, values in vectors.items():
            object.__setattr__(self, field_name, values)

    @property
    def n_records(self) -> int:
        return len(self.parent_ids)

    @property
    def unique_dependency_group_ids(self) -> tuple[object, ...]:
        return tuple(dict.fromkeys(self.dependency_group_ids))

    @property
    def unique_dependency_group_weights(self) -> np.ndarray:
        first_index_by_group: dict[object, int] = {}
        for index, group_id in enumerate(self.dependency_group_ids):
            first_index_by_group.setdefault(group_id, index)
        values = np.array(
            [self.dependency_group_weights[index] for index in first_index_by_group.values()],
            dtype=float,
        )
        values.setflags(write=False)
        return values


@dataclass(frozen=True)
class CalibrationFitDiagnostics:
    """Counts about candidate records excluded before fitting the sample."""

    n_positive_weight_records: int

    def __post_init__(self) -> None:
        value = self.n_positive_weight_records
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
            raise ValueError(
                "n_positive_weight_records must be a non-negative integer; "
                f"got {value!r}."
            )


@dataclass(frozen=True)
class IndependentUnweightedCommonScaleContract:
    """Runtime-checkable provenance required by the restricted exact F law."""

    calibration_observation_ids_by_parent: tuple[tuple[object, frozenset[object]], ...]
    common_reference_scale: float
    denominator_degrees_of_freedom: float

    def __post_init__(self) -> None:
        if not self.calibration_observation_ids_by_parent:
            raise ValueError("Independent calibration observation IDs must be non-empty.")
        parent_ids = tuple(parent_id for parent_id, _ in self.calibration_observation_ids_by_parent)
        if len(set(parent_ids)) != len(parent_ids):
            raise ValueError("Independent calibration parent IDs must be unique.")
        seen_observation_ids: set[object] = set()
        for parent_id, observation_ids in self.calibration_observation_ids_by_parent:
            if not observation_ids:
                raise ValueError(
                    "Every independent calibration parent requires observation IDs; "
                    f"parent={parent_id!r}."
                )
            if seen_observation_ids.intersection(observation_ids):
                raise ValueError("Independent calibration record observation IDs must be pairwise disjoint.")
            seen_observation_ids.update(observation_ids)
        if not np.isfinite(self.common_reference_scale) or self.common_reference_scale <= 0.0:
            raise ValueError("Independent calibration common_reference_scale must be positive.")
        if (
            not np.isfinite(self.denominator_degrees_of_freedom)
            or self.denominator_degrees_of_freedom <= 0.0
        ):
            raise ValueError(
                "Independent calibration denominator_degrees_of_freedom must be positive."
            )

    @property
    def calibration_observation_ids(self) -> frozenset[object]:
        return frozenset(
            observation_id
            for _, observation_ids in self.calibration_observation_ids_by_parent
            for observation_id in observation_ids
        )


def _effective_sample_size_from_weights(weights: np.ndarray) -> float:
    log_weights = np.log(weights)
    return float(
        np.exp(
            2.0 * np.logaddexp.reduce(log_weights)
            - np.logaddexp.reduce(2.0 * log_weights)
        )
    )


@dataclass(frozen=True)
class EmpiricalNullInflationModel:
    """Fitted scale state separated from its authoritative calibration sample."""

    method: str
    sample: CalibrationSample
    baseline_scale_estimate: float
    fit_diagnostics: CalibrationFitDiagnostics
    reference_law: CalibrationReferenceLaw = "unresolved_same_selected_hierarchy"
    independent_contract: IndependentUnweightedCommonScaleContract | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.baseline_scale_estimate) or self.baseline_scale_estimate < 0.0:
            raise ValueError(
                "EmpiricalNullInflationModel.baseline_scale_estimate must be finite "
                f"and non-negative; got {self.baseline_scale_estimate!r}."
            )
        if self.fit_diagnostics.n_positive_weight_records < self.sample.n_records:
            raise ValueError(
                "n_positive_weight_records cannot be smaller than the calibration sample."
            )
        has_exact_law = self.reference_law == "exact_independent_unweighted_common_scale_f"
        if has_exact_law != (self.independent_contract is not None):
            raise ValueError(
                "The exact independent F reference law and its provenance contract "
                "must be present together."
            )

    @property
    def n_calibration(self) -> int:
        return self.sample.n_records

    @property
    def n_positive_weight_records(self) -> int:
        return self.fit_diagnostics.n_positive_weight_records

    @property
    def n_selected_nonnull_positive_weight_records(self) -> int:
        return self.n_positive_weight_records - self.n_calibration

    @property
    def n_strict_null_calibration(self) -> int:
        return int(np.sum(self.sample.is_strict_null))

    @property
    def n_edge_blocked_calibration(self) -> int:
        return int(np.sum(self.sample.is_edge_blocked))

    @property
    def baseline_empirical_inflation_factor(self) -> float:
        return float(max(self.baseline_scale_estimate, 1.0))

    @property
    def effective_sample_size(self) -> float:
        return _effective_sample_size_from_weights(self.sample.unique_dependency_group_weights)

    @property
    def context_center(self) -> np.ndarray:
        center = np.sum(self.sample.contexts * self.sample.weights[:, None], axis=0) / np.sum(
            self.sample.weights
        )
        center.setflags(write=False)
        return center

    @property
    def context_bandwidth(self) -> np.ndarray:
        center = self.context_center
        variances = np.sum(
            ((self.sample.contexts - center) ** 2) * self.sample.weights[:, None],
            axis=0,
        ) / np.sum(self.sample.weights)
        bandwidth = np.sqrt(np.maximum(variances, 0.0))
        bandwidth = np.where(bandwidth <= 1e-12, 0.0, bandwidth)
        bandwidth.setflags(write=False)
        return bandwidth

    @property
    def sample_contexts(self) -> np.ndarray:
        return self.sample.contexts

    @property
    def sample_feature_families(self) -> tuple[str, ...]:
        return self.sample.feature_families

    @property
    def sample_weights(self) -> np.ndarray:
        return self.sample.weights

    @property
    def sample_parent_ids(self) -> tuple[object, ...]:
        return self.sample.parent_ids

    @property
    def sample_is_strict_null(self) -> np.ndarray:
        return self.sample.is_strict_null

    @property
    def sample_is_edge_blocked(self) -> np.ndarray:
        return self.sample.is_edge_blocked

    @property
    def sample_statistics(self) -> np.ndarray:
        return self.sample.statistics

    @property
    def sample_reference_scales(self) -> np.ndarray:
        return self.sample.reference_scales

    @property
    def sample_degrees_of_freedom(self) -> np.ndarray:
        return self.sample.degrees_of_freedom


__all__ = [
    "CalibrationDecision",
    "CalibrationDecisionStatus",
    "CalibrationFitDiagnostics",
    "CalibrationReferenceLaw",
    "CalibrationSample",
    "CalibrationSupportPolicySnapshot",
    "CalibrationSupportThresholds",
    "DEFAULT_INTERNAL_SUPPORT_POLICY",
    "DEFAULT_INTERNAL_SUPPORT_THRESHOLDS",
    "EmpiricalNullInflationModel",
    "IndependentUnweightedCommonScaleContract",
    "INTERNAL_SUPPORT_CONTRACT_VERSION",
]
