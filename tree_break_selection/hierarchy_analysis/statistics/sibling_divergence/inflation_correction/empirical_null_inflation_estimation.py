"""Context-weighted empirical-null inflation estimation."""

from __future__ import annotations

from collections.abc import Collection, Mapping

import numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2, f

from ..pair_testing.types.sibling_pair_record import SiblingPairRecord
from .types.inflation_model import (
    DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
    CalibrationDecision,
    CalibrationFitDiagnostics,
    CalibrationReferenceLaw,
    CalibrationSample,
    CalibrationSupportThresholds,
    EmpiricalNullInflationModel,
    IndependentUnweightedCommonScaleContract,
)

# =============================================================================
# Context-weighted empirical-null inflation estimation
# =============================================================================


def _validate_calibration_record(record: SiblingPairRecord) -> None:
    """Validate one candidate record for empirical-null inflation estimation."""
    if not np.isfinite(record.stat):
        raise ValueError(
            f"Sibling inflation calibration requires finite statistics; parent={record.parent!r}."
        )
    if record.degrees_of_freedom < 0:
        raise ValueError(
            "Sibling inflation calibration requires non-negative degrees of freedom; "
            f"parent={record.parent!r}."
        )
    if not np.isfinite(record.reference_scale) or record.reference_scale <= 0:
        raise ValueError(
            "Sibling inflation calibration requires finite positive reference_scale; "
            f"parent={record.parent!r}."
        )
    if record.stat < 0:
        raise ValueError(
            "Sibling inflation calibration requires non-negative statistics; "
            f"parent={record.parent!r}."
        )
    if not np.isfinite(record.sibling_null_weight) or not (
        0.0 <= record.sibling_null_weight <= 1.0
    ):
        raise ValueError(
            "Sibling inflation calibration requires finite sibling_null_weight in [0, 1]; "
            f"parent={record.parent!r}."
        )
    if record.degrees_of_freedom > 0 and (
        not np.isfinite(record.sibling_projection_dimension)
        or record.sibling_projection_dimension <= 0
    ):
        raise ValueError(
            "Positive-degree sibling calibration records require positive "
            f"sibling_projection_dimension; parent={record.parent!r}."
        )
    if record.degrees_of_freedom > 0 and record.n_parent <= 0:
        raise ValueError(
            "Positive-degree sibling calibration records require positive "
            f"parent sample size; parent={record.parent!r}."
        )
    if record.feature_family not in {"bernoulli", "categorical", "continuous", "mixed"}:
        raise ValueError(
            "Sibling inflation calibration requires feature_family to be "
            f"'bernoulli', 'categorical', 'continuous', or 'mixed'; "
            f"parent={record.parent!r}, feature_family={record.feature_family!r}."
        )


def _inflation_mle(
    statistics: np.ndarray,
    reference_expectations: np.ndarray,
    weights: np.ndarray,
) -> float:
    return float(np.sum(weights * statistics) / np.sum(weights * reference_expectations))


def _has_internal_empirical_null_support(record: SiblingPairRecord) -> bool:
    """Return whether a record is admissible empirical-null calibration."""
    return record.has_empirical_null_support


def _effective_sample_size(weights: np.ndarray) -> float:
    positive_weights = weights[weights > 0.0]
    if positive_weights.size == 0:
        raise ValueError("Effective sample size requires positive weights.")
    log_weights = np.log(positive_weights)
    return float(np.exp(2.0 * logsumexp(log_weights) - logsumexp(2.0 * log_weights)))


def _max_weight_share(weights: np.ndarray) -> float:
    positive_weights = weights[weights > 0.0]
    if positive_weights.size == 0:
        return 0.0
    return float(np.max(positive_weights) / np.sum(positive_weights))


def _leave_one_group_max_delta_log_c(
    statistics: np.ndarray,
    reference_expectations: np.ndarray,
    weights: np.ndarray,
    dependency_group_ids: tuple[object, ...],
    baseline_c_hat: float,
) -> float:
    unique_groups = tuple(dict.fromkeys(dependency_group_ids))
    if len(unique_groups) <= 1:
        return float("inf")
    deltas: list[float] = []
    baseline_log_c = float(np.log(max(baseline_c_hat, 1e-300)))
    for group_id in unique_groups:
        keep = np.array([candidate != group_id for candidate in dependency_group_ids])
        if float(np.sum(weights[keep])) <= 0.0:
            return float("inf")
        c_hat = max(
            _inflation_mle(
                statistics[keep],
                reference_expectations[keep],
                weights[keep],
            ),
            1.0,
        )
        deltas.append(abs(float(np.log(max(c_hat, 1e-300))) - baseline_log_c))
    return float(max(deltas)) if deltas else float("inf")


def _dependency_group_support_weights(
    model: EmpiricalNullInflationModel,
    *,
    family_mask: np.ndarray | None = None,
    local_record_weights: np.ndarray | None = None,
) -> np.ndarray:
    sample = model.sample
    selected_indices = (
        np.arange(sample.n_records, dtype=int)
        if family_mask is None
        else np.flatnonzero(family_mask)
    )
    if local_record_weights is not None and local_record_weights.shape != selected_indices.shape:
        raise ValueError("Local calibration weights must align to the selected sample records.")

    representative_by_group: dict[object, int] = {}
    selected_position_by_index = {
        int(sample_index): position for position, sample_index in enumerate(selected_indices)
    }
    for sample_index in selected_indices:
        index = int(sample_index)
        group_id = sample.dependency_group_ids[index]
        representative_by_group.setdefault(group_id, index)
        if sample.parent_ids[index] == group_id:
            representative_by_group[group_id] = index

    group_weights: list[float] = []
    for sample_index in representative_by_group.values():
        weight = float(sample.dependency_group_weights[sample_index])
        if local_record_weights is not None:
            selected_position = selected_position_by_index[sample_index]
            kernel_multiplier = float(
                local_record_weights[selected_position] / sample.weights[sample_index]
            )
            weight *= kernel_multiplier
        group_weights.append(weight)
    return np.asarray(group_weights, dtype=float)


def _fit_empirical_null_inflation_model(
    records: list[SiblingPairRecord],
    *,
    _require_unit_reference_scale: bool = True,
    _method: str = "context_weighted_supported_empirical_null_inflation",
    _reference_law: CalibrationReferenceLaw = "unresolved_same_selected_hierarchy",
    _independent_contract: IndependentUnweightedCommonScaleContract | None = None,
) -> EmpiricalNullInflationModel:
    """Fit the context-weighted empirical-null inflation model."""
    if not records:
        raise ValueError("Cannot fit sibling inflation model: no sibling calibration records.")
    for record in records:
        _validate_calibration_record(record)
    if _require_unit_reference_scale and any(
        not np.isclose(record.reference_scale, 1.0, rtol=1e-12, atol=1e-12)
        for record in records
    ):
        first_non_unit = next(
            record
            for record in records
            if not np.isclose(record.reference_scale, 1.0, rtol=1e-12, atol=1e-12)
        )
        raise ValueError(
            "Sibling empirical-null inflation currently requires unit reference_scale "
            "from the orthonormal projected-Wald reference law; "
            f"parent={first_non_unit.parent!r}, "
            f"reference_scale={first_non_unit.reference_scale!r}."
        )

    positive_df_records = [record for record in records if record.degrees_of_freedom > 0]
    if not positive_df_records:
        raise ValueError(
            "Cannot fit sibling inflation model: no positive-degree calibration records."
        )

    positive_weight_records = [
        record for record in positive_df_records if record.sibling_null_weight > 0.0
    ]
    if not positive_weight_records:
        raise ValueError(
            "Cannot fit sibling inflation model: no positive sibling-null calibration weight."
        )

    role_supported_records = [
        record for record in positive_weight_records if _has_internal_empirical_null_support(record)
    ]
    if not role_supported_records:
        selected_nonnull_count = sum(
            not _has_internal_empirical_null_support(record) for record in positive_weight_records
        )
        raise ValueError(
            "Cannot fit sibling inflation model: no strict-null or stopped-edge "
            "empirical-null calibration records with positive weight. "
            f"Found {selected_nonnull_count} selected non-null positive-weight "
            "record(s), which are not valid empirical-null calibration support."
        )
    supported_records = role_supported_records

    statistics = np.array([record.stat for record in supported_records], dtype=float)
    reference_scales = np.array(
        [record.reference_scale for record in supported_records],
        dtype=float,
    )
    degrees_of_freedom = np.array(
        [record.degrees_of_freedom for record in supported_records],
        dtype=float,
    )
    null_weights = np.array(
        [record.sibling_null_weight for record in supported_records],
        dtype=float,
    )
    is_strict_null = np.array(
        [
            record.is_null_like and not record.is_edge_blocked
            for record in supported_records
        ],
        dtype=bool,
    )
    is_edge_blocked = np.array(
        [record.is_edge_blocked for record in supported_records],
        dtype=bool,
    )
    calibration_scales = np.array(
        [record.sibling_projection_dimension for record in supported_records],
        dtype=float,
    )
    calibration_parent_sample_sizes = np.array(
        [record.n_parent for record in supported_records],
        dtype=float,
    )
    calibration_feature_families = tuple(record.feature_family for record in supported_records)
    calibration_parent_ids = tuple(record.parent for record in supported_records)
    calibration_dependency_group_ids = tuple(
        record.parent
        if record.calibration_dependency_group is None
        else record.calibration_dependency_group
        for record in supported_records
    )
    calibration_dependency_group_weights = np.array(
        [
            record.sibling_null_weight
            if record.calibration_dependency_weight is None
            else record.calibration_dependency_weight
            for record in supported_records
        ],
        dtype=float,
    )
    reference_expectations = reference_scales * degrees_of_freedom

    baseline_scale_estimate = _inflation_mle(
        statistics,
        reference_expectations,
        null_weights,
    )

    sample_contexts = np.column_stack(
        [
            np.log(calibration_scales),
            np.log(calibration_parent_sample_sizes),
        ]
    )
    return EmpiricalNullInflationModel(
        method=_method,
        sample=CalibrationSample(
            contexts=sample_contexts,
            feature_families=calibration_feature_families,
            weights=null_weights,
            parent_ids=calibration_parent_ids,
            dependency_group_ids=calibration_dependency_group_ids,
            dependency_group_weights=calibration_dependency_group_weights,
            is_strict_null=is_strict_null,
            is_edge_blocked=is_edge_blocked,
            statistics=statistics,
            reference_scales=reference_scales,
            degrees_of_freedom=degrees_of_freedom,
        ),
        baseline_scale_estimate=baseline_scale_estimate,
        fit_diagnostics=CalibrationFitDiagnostics(
            n_positive_weight_records=int(len(positive_weight_records)),
        ),
        reference_law=_reference_law,
        independent_contract=_independent_contract,
    )


def fit_empirical_null_inflation_model(
    records: list[SiblingPairRecord],
) -> EmpiricalNullInflationModel:
    """Fit the empirical-null scale model; its selected-tail law remains unvalidated."""
    return _fit_empirical_null_inflation_model(records)


def _validated_observation_ids(
    observation_ids: Collection[object],
    *,
    role: str,
) -> frozenset[object]:
    values = tuple(observation_ids)
    if not values:
        raise ValueError(f"{role} observation IDs must be non-empty.")
    try:
        unique_values = frozenset(values)
    except TypeError as exc:
        raise ValueError(f"{role} observation IDs must be hashable.") from exc
    if len(unique_values) != len(values):
        raise ValueError(f"{role} observation IDs must be unique.")
    return unique_values


def fit_independent_unweighted_common_scale_inflation_model(
    records: list[SiblingPairRecord],
    *,
    calibration_observation_ids_by_parent: Mapping[object, Collection[object]],
) -> EmpiricalNullInflationModel:
    """Fit the denominator for the restricted independent exact F reference law."""
    if any(not record.has_empirical_null_support for record in records):
        raise ValueError("Independent exact calibration records must all have a null role.")
    if any(record.degrees_of_freedom <= 0.0 for record in records):
        raise ValueError("Independent exact calibration records must have positive degrees of freedom.")
    if any(
        record.sibling_null_weight != 1.0
        or (
            record.calibration_dependency_weight is not None
            and record.calibration_dependency_weight != 1.0
        )
        for record in records
    ):
        raise ValueError("Independent exact calibration requires unit fixed weights.")
    dependency_group_ids = tuple(
        record.parent
        if record.calibration_dependency_group is None
        else record.calibration_dependency_group
        for record in records
    )
    if len(set(dependency_group_ids)) != len(dependency_group_ids):
        raise ValueError("Independent exact calibration records must have distinct dependency groups.")
    reference_scales = np.asarray([record.reference_scale for record in records], dtype=float)
    if reference_scales.size == 0:
        raise ValueError("Independent exact calibration requires calibration records.")
    common_reference_scale = float(reference_scales[0])
    if not np.allclose(reference_scales, common_reference_scale, rtol=1e-12, atol=1e-12):
        raise ValueError("Independent exact calibration requires one common reference_scale.")
    parent_ids = tuple(record.parent for record in records)
    if set(calibration_observation_ids_by_parent) != set(parent_ids):
        raise ValueError(
            "Independent calibration observation-ID keys must match calibration parent IDs."
        )
    observation_ids_by_parent: list[tuple[object, frozenset[object]]] = []
    seen_observation_ids: set[object] = set()
    for parent_id in parent_ids:
        observation_ids = _validated_observation_ids(
            calibration_observation_ids_by_parent[parent_id],
            role=f"Calibration parent {parent_id!r}",
        )
        if seen_observation_ids.intersection(observation_ids):
            raise ValueError(
                "Independent calibration record observation IDs must be pairwise disjoint."
            )
        seen_observation_ids.update(observation_ids)
        observation_ids_by_parent.append((parent_id, observation_ids))
    denominator_degrees_of_freedom = float(
        np.sum([record.degrees_of_freedom for record in records])
    )
    contract = IndependentUnweightedCommonScaleContract(
        calibration_observation_ids_by_parent=tuple(observation_ids_by_parent),
        common_reference_scale=common_reference_scale,
        denominator_degrees_of_freedom=denominator_degrees_of_freedom,
    )
    return _fit_empirical_null_inflation_model(
        records,
        _require_unit_reference_scale=False,
        _method="exact_independent_unweighted_common_scale_f",
        _reference_law="exact_independent_unweighted_common_scale_f",
        _independent_contract=contract,
    )


def _decision_support(
    *,
    model: EmpiricalNullInflationModel,
    record: SiblingPairRecord,
    family_mask: np.ndarray | None = None,
    local_weights: np.ndarray | None = None,
) -> dict[str, float | int | str | bool]:
    support: dict[str, float | int | str | bool] = {
        "n_positive_weight_records": int(model.n_positive_weight_records),
        "n_supported_records": int(model.n_calibration),
        "n_supported_groups": len(model.sample.unique_dependency_group_ids),
        "n_selected_nonnull_positive_weight_records": int(
            model.n_selected_nonnull_positive_weight_records
        ),
        "n_strict_null_records": int(model.n_strict_null_calibration),
        "n_edge_blocked_records": int(model.n_edge_blocked_calibration),
        "model_effective_sample_size": float(model.effective_sample_size),
        "model_max_group_weight_share": _max_weight_share(
            model.sample.unique_dependency_group_weights
        ),
        "feature_family": str(record.feature_family),
    }
    if family_mask is not None:
        support["n_family_supported_records"] = int(np.sum(family_mask))
        if np.any(family_mask):
            family_weights = model.sample_weights[family_mask]
            family_group_weights = _dependency_group_support_weights(
                model,
                family_mask=family_mask,
            )
            support["n_family_supported_groups"] = int(family_group_weights.size)
            support["family_effective_sample_size"] = _effective_sample_size(
                family_group_weights
            )
            family_dependency_group_ids = tuple(
                group_id
                for group_id, include in zip(
                    model.sample.dependency_group_ids,
                    family_mask,
                    strict=True,
                )
                if include
            )
            support["leave_one_group_max_delta_log_c"] = _leave_one_group_max_delta_log_c(
                model.sample_statistics[family_mask],
                (
                    model.sample_reference_scales[family_mask]
                    * model.sample_degrees_of_freedom[family_mask]
                ),
                family_weights,
                family_dependency_group_ids,
                max(
                    _inflation_mle(
                        model.sample_statistics[family_mask],
                        (
                            model.sample_reference_scales[family_mask]
                            * model.sample_degrees_of_freedom[family_mask]
                        ),
                        family_weights,
                    ),
                    1.0,
                ),
            )
    if local_weights is not None:
        local_group_weights = _dependency_group_support_weights(
            model,
            family_mask=family_mask,
            local_record_weights=local_weights,
        )
        support["local_effective_sample_size"] = (
            _effective_sample_size(local_group_weights)
            if float(np.sum(local_group_weights)) > 0.0
            else 0.0
        )
        support["local_max_group_weight_share"] = _max_weight_share(local_group_weights)
    return support


def _support_contract_failures(
    support: dict[str, float | int | str | bool],
    *,
    thresholds: CalibrationSupportThresholds,
) -> tuple[str, ...]:
    failures: list[str] = []
    if int(support.get("n_supported_groups", 0)) < thresholds.min_supported_groups:
        failures.append("supported_groups_below_threshold")
    if int(support.get("n_family_supported_groups", 0)) < (
        thresholds.min_family_supported_groups
    ):
        failures.append("family_supported_groups_below_threshold")
    if float(support.get("family_effective_sample_size", 0.0)) < (
        thresholds.min_family_effective_sample_size
    ):
        failures.append("family_effective_sample_size_below_threshold")
    if float(support.get("local_effective_sample_size", 0.0)) < (
        thresholds.min_local_effective_sample_size
    ):
        failures.append("local_effective_sample_size_below_threshold")
    if float(support.get("local_max_group_weight_share", 1.0)) > thresholds.max_weight_share:
        failures.append("local_max_group_weight_share_above_threshold")
    if float(support.get("leave_one_group_max_delta_log_c", float("inf"))) > (
        thresholds.max_leave_one_group_delta_log_c
    ):
        failures.append("leave_one_group_delta_log_c_above_threshold")
    return tuple(failures)


def _with_support_contract(
    support: dict[str, float | int | str | bool],
    *,
    thresholds: CalibrationSupportThresholds,
) -> tuple[dict[str, float | int | str | bool], tuple[str, ...]]:
    failures = _support_contract_failures(support, thresholds=thresholds)
    annotated_support = dict(support)
    annotated_support["support_contract_status"] = (
        "passes_internal_support_thresholds"
        if not failures
        else "below_internal_support_thresholds"
    )
    annotated_support["support_contract_failure_reasons"] = ";".join(failures)
    return annotated_support, failures


def _decision_context(record: SiblingPairRecord) -> dict[str, object]:
    return {
        "feature_family": record.feature_family,
        "sibling_projection_dimension": float(record.sibling_projection_dimension),
        "n_parent": int(record.n_parent),
    }


def decide_independent_unweighted_common_scale_calibration(
    model: EmpiricalNullInflationModel,
    record: SiblingPairRecord,
    *,
    focal_observation_ids: Collection[object],
) -> CalibrationDecision:
    """Return the restricted exact F decision for disjoint focal observations."""
    if model.reference_law != "exact_independent_unweighted_common_scale_f":
        raise ValueError("Exact F calibration requires the independently fitted exact model.")
    contract = model.independent_contract
    if contract is None:
        raise ValueError("Exact F calibration model is missing its provenance contract.")
    focal_ids = _validated_observation_ids(focal_observation_ids, role="Focal")
    if focal_ids.intersection(contract.calibration_observation_ids):
        raise ValueError("Focal and calibration observation IDs must be disjoint.")
    if record.parent in model.sample.parent_ids:
        raise ValueError("Focal and calibration parent IDs must be disjoint.")
    _validate_calibration_record(record)
    if record.degrees_of_freedom <= 0.0:
        raise ValueError("Exact F calibration requires positive focal degrees of freedom.")
    if not np.isclose(
        record.reference_scale,
        contract.common_reference_scale,
        rtol=1e-12,
        atol=1e-12,
    ):
        raise ValueError("Exact F calibration requires the focal common reference_scale.")
    expected_unadjusted_p_value = float(
        chi2.sf(
            record.stat / record.reference_scale,
            df=float(record.degrees_of_freedom),
        )
    )
    if not np.isclose(record.p_value, expected_unadjusted_p_value, rtol=1e-10, atol=1e-14):
        raise ValueError(
            "Exact F calibration requires the focal p_value to match its unadjusted "
            "common-scale chi-square law."
        )
    scale_estimate = float(model.baseline_scale_estimate)
    f_statistic = (
        float("inf")
        if scale_estimate <= 0.0
        else float(
            record.stat
            / (
                record.reference_scale
                * record.degrees_of_freedom
                * scale_estimate
            )
        )
    )
    exact_p_value = float(
        f.sf(
            f_statistic,
            dfn=float(record.degrees_of_freedom),
            dfd=contract.denominator_degrees_of_freedom,
        )
    )
    one_sided_p_value = max(exact_p_value, expected_unadjusted_p_value)
    return CalibrationDecision(
        status="independent_exact_admissible",
        c_hat=model.baseline_empirical_inflation_factor,
        p_value=one_sided_p_value,
        estimator="exact_independent_unweighted_common_scale_f",
        support=_decision_support(model=model, record=record),
        exact_context={
            **_decision_context(record),
            "reference_law": "exact_independent_unweighted_common_scale_f",
            "denominator_degrees_of_freedom": contract.denominator_degrees_of_freedom,
            "common_reference_scale": contract.common_reference_scale,
        },
        descriptive_strata={
            "f_statistic": f_statistic,
            "unadjusted_p_value": expected_unadjusted_p_value,
            "one_sided_p_value_floor_applied": bool(
                expected_unadjusted_p_value > exact_p_value
            ),
        },
    )


def _adjusted_p_value(record: SiblingPairRecord, c_hat: float) -> float:
    """Evaluate the historical plug-in chi-square tail after empirical inflation."""
    if not np.isfinite(record.reference_scale) or record.reference_scale <= 0.0:
        raise ValueError("Empirical-null adjustment requires a finite positive reference scale.")
    if record.degrees_of_freedom == 0.0:
        return 1.0
    return float(
        chi2.sf(
            record.stat / (record.reference_scale * c_hat),
            df=record.degrees_of_freedom,
        )
    )


def decide_empirical_null_calibration(
    model: EmpiricalNullInflationModel,
    record: SiblingPairRecord,
    *,
    enforce_support_thresholds: bool = False,
    support_thresholds: CalibrationSupportThresholds = DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
) -> CalibrationDecision:
    """Return the historical empirical rule, without claiming selected-tail validity.

    ``internal_admissible`` denotes availability under the internal support policy.
    The plug-in chi-square p-value does not establish selective error control.
    """
    if model.reference_law != "unresolved_same_selected_hierarchy":
        raise ValueError(
            "Empirical-null calibration requires a same-selected-hierarchy model; "
            "use the independent calibration decision with focal observation ownership "
            "for an exact F model."
        )
    exact_context = _decision_context(record)
    if record.degrees_of_freedom == 0:
        return CalibrationDecision(
            status="internal_admissible",
            c_hat=1.0,
            p_value=1.0,
            estimator="zero_dimensional_sibling_record",
            support=_decision_support(model=model, record=record),
            exact_context=exact_context,
            descriptive_strata={"degrees_of_freedom": 0.0},
        )
    if (
        not np.isfinite(record.sibling_projection_dimension)
        or record.sibling_projection_dimension <= 0
    ):
        raise ValueError(
            "Positive-degree sibling records require positive sibling_projection_dimension "
            f"for empirical-null inflation prediction; parent={record.parent!r}."
        )
    if model.n_calibration == 0:
        return CalibrationDecision(
            status="undefined_no_internal_support",
            c_hat=None,
            p_value=None,
            estimator=model.method,
            support=_decision_support(model=model, record=record),
            exact_context=exact_context,
            descriptive_strata={"reason": "empty_calibration_model"},
        )
    if record.feature_family not in {"bernoulli", "categorical", "continuous", "mixed"}:
        raise ValueError(
            "Positive-degree sibling records require feature_family to be "
            f"'bernoulli', 'categorical', 'continuous', or 'mixed'; "
            f"parent={record.parent!r}, feature_family={record.feature_family!r}."
        )

    family_mask = np.array(
        [
            feature_family == record.feature_family
            for feature_family in model.sample_feature_families
        ],
        dtype=bool,
    )
    if not np.any(family_mask):
        return CalibrationDecision(
            status="undefined_no_family_support",
            c_hat=None,
            p_value=None,
            estimator=model.method,
            support=_decision_support(
                model=model,
                record=record,
                family_mask=family_mask,
            ),
            exact_context=exact_context,
            descriptive_strata={"reason": "no_matching_feature_family"},
        )
    family_reference_expectations = (
        model.sample_reference_scales[family_mask] * model.sample_degrees_of_freedom[family_mask]
    )
    family_baseline_inflation_factor = max(
        _inflation_mle(
            model.sample_statistics[family_mask],
            family_reference_expectations,
            model.sample_weights[family_mask],
        ),
        1.0,
    )

    active_context_axes = np.array([True, record.feature_family != "bernoulli"])
    active_context_axes = active_context_axes & (model.context_bandwidth > 0.0)
    descriptive_strata: dict[str, object] = {
        "active_projection_dimension_axis": bool(active_context_axes[0]),
        "active_parent_size_axis": bool(active_context_axes[1]),
        "context_bandwidth_projection_dimension": float(model.context_bandwidth[0])
        if model.context_bandwidth.size > 0
        else 0.0,
        "context_bandwidth_parent_size": float(model.context_bandwidth[1])
        if model.context_bandwidth.size > 1
        else 0.0,
    }
    if not np.any(active_context_axes):
        c_hat = float(max(family_baseline_inflation_factor, 1.0))
        local_weights = model.sample_weights[family_mask]
        support, failures = _with_support_contract(
            _decision_support(
                model=model,
                record=record,
                family_mask=family_mask,
                local_weights=local_weights,
            ),
            thresholds=support_thresholds,
        )
        if enforce_support_thresholds and failures:
            return CalibrationDecision(
                status="undefined_sparse_context",
                c_hat=None,
                p_value=None,
                estimator=f"{model.method}:family_baseline",
                support=support,
                exact_context=exact_context,
                descriptive_strata={
                    **descriptive_strata,
                    "reason": "internal_support_thresholds_failed",
                },
            )
        return CalibrationDecision(
            status="internal_admissible",
            c_hat=c_hat,
            p_value=_adjusted_p_value(record, c_hat),
            estimator=f"{model.method}:family_baseline",
            support=support,
            exact_context=exact_context,
            descriptive_strata={
                **descriptive_strata,
                "reason": "same_selected_hierarchy_reference_law_unvalidated",
            },
        )

    if record.n_parent <= 0:
        raise ValueError(
            "Positive-degree sibling records require positive parent sample size "
            f"for empirical-null inflation prediction; parent={record.parent!r}."
        )
    target_context = np.array(
        [
            np.log(record.sibling_projection_dimension),
            np.log(float(record.n_parent)),
        ],
        dtype=float,
    )
    family_contexts = model.sample_contexts[family_mask]
    family_weights = model.sample_weights[family_mask]
    family_statistics = model.sample_statistics[family_mask]
    scaled_offsets = (
        family_contexts[:, active_context_axes] - target_context[active_context_axes]
    ) / model.context_bandwidth[active_context_axes]
    log_kernel_weights = -0.5 * np.sum(scaled_offsets**2, axis=1)
    log_kernel_weights = log_kernel_weights - float(np.max(log_kernel_weights))
    local_weights = family_weights * np.exp(log_kernel_weights)
    if float(np.sum(local_weights)) <= 0.0:
        return CalibrationDecision(
            status="undefined_sparse_context",
            c_hat=None,
            p_value=None,
            estimator=model.method,
            support=_decision_support(
                model=model,
                record=record,
                family_mask=family_mask,
                local_weights=local_weights,
            ),
            exact_context=exact_context,
            descriptive_strata={
                **descriptive_strata,
                "reason": "zero_local_calibration_weight",
            },
        )

    reference_expectations = family_reference_expectations
    local_inflation_factor = _inflation_mle(
        family_statistics,
        reference_expectations,
        local_weights,
    )
    c_hat = float(max(local_inflation_factor, 1.0))
    support, failures = _with_support_contract(
        _decision_support(
            model=model,
            record=record,
            family_mask=family_mask,
            local_weights=local_weights,
        ),
        thresholds=support_thresholds,
    )
    if enforce_support_thresholds and failures:
        return CalibrationDecision(
            status="undefined_sparse_context",
            c_hat=None,
            p_value=None,
            estimator=f"{model.method}:local_kernel",
            support=support,
            exact_context=exact_context,
            descriptive_strata={
                **descriptive_strata,
                "reason": "internal_support_thresholds_failed",
            },
        )
    return CalibrationDecision(
        status="internal_admissible",
        c_hat=c_hat,
        p_value=_adjusted_p_value(record, c_hat),
        estimator=f"{model.method}:local_kernel",
        support=support,
        exact_context=exact_context,
        descriptive_strata={
            **descriptive_strata,
            "reason": "same_selected_hierarchy_reference_law_unvalidated",
        },
    )


def predict_empirical_inflation_factor(
    model: EmpiricalNullInflationModel,
    record: SiblingPairRecord,
    *,
    enforce_support_thresholds: bool = False,
    support_thresholds: CalibrationSupportThresholds = DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
) -> float:
    """Predict the empirical post-selection inflation for one sibling record."""
    decision = decide_empirical_null_calibration(
        model,
        record,
        enforce_support_thresholds=enforce_support_thresholds,
        support_thresholds=support_thresholds,
    )
    diagnostic_statuses = {
        "internal_admissible",
        "undefined_unvalidated_reference_law",
    }
    if decision.status not in diagnostic_statuses or decision.c_hat is None:
        reason = decision.descriptive_strata.get("reason", decision.status)
        raise ValueError(
            "Empirical-null diagnostic inflation prediction is unavailable: "
            f"{decision.status}; parent={record.parent!r}; reason={reason!r}."
        )
    return decision.c_hat


__all__ = [
    "CalibrationDecision",
    "CalibrationSupportThresholds",
    "EmpiricalNullInflationModel",
    "DEFAULT_INTERNAL_SUPPORT_THRESHOLDS",
    "decide_empirical_null_calibration",
    "decide_independent_unweighted_common_scale_calibration",
    "fit_empirical_null_inflation_model",
    "fit_independent_unweighted_common_scale_inflation_model",
    "predict_empirical_inflation_factor",
]
