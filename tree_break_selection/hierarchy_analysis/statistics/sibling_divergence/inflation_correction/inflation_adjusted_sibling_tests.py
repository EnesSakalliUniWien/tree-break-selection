"""Inflation-adjusted sibling-test summaries."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from ..pair_testing.types.sibling_pair_record import SiblingPairRecord
from .empirical_null_inflation_estimation import (
    DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
    decide_empirical_null_calibration,
)
from .types.inflation_model import (
    CalibrationDecision,
    CalibrationSupportThresholds,
    EmpiricalNullInflationModel,
)

InflationAdjustedSiblingTestSummary = tuple[float, float, float]


def _compute_inflation_adjusted_sibling_test(
    sibling_test_record: SiblingPairRecord,
    *,
    model: EmpiricalNullInflationModel | None,
    enforce_support_thresholds: bool = False,
    support_thresholds: CalibrationSupportThresholds = DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
    family_support_cache: dict[str, dict[str, float | int]] | None = None,
) -> tuple[InflationAdjustedSiblingTestSummary | None, str, CalibrationDecision | None]:
    """Return one inflation-adjusted sibling-test summary and its calibration decision.

    The summary is None when the decision is not internally admissible.
    """
    if not np.isfinite(sibling_test_record.stat):
        raise ValueError(
            "Sibling record must have a finite statistic before "
            f"adjustment; parent={sibling_test_record.parent!r}."
        )
    if sibling_test_record.degrees_of_freedom == 0:
        if sibling_test_record.stat != 0.0 or sibling_test_record.p_value != 1.0:
            raise ValueError(
                "Zero-dimensional sibling records must carry statistic=0 and p_value=1; "
                f"parent={sibling_test_record.parent!r}."
            )
        method = "zero_dimensional_sibling_record" if model is None else model.method
        return (0.0, 0.0, 1.0), method, None
    if sibling_test_record.degrees_of_freedom < 0:
        raise ValueError(
            "Sibling record must have non-negative degrees of freedom before "
            f"adjustment; parent={sibling_test_record.parent!r}."
        )

    if model is None:
        raise ValueError(
            "Sibling empirical calibration has no internal model; "
            f"parent={sibling_test_record.parent!r}."
        )
    decision = decide_empirical_null_calibration(
        model,
        sibling_test_record,
        enforce_support_thresholds=enforce_support_thresholds,
        support_thresholds=support_thresholds,
        family_support_cache=family_support_cache,
    )
    if decision.status != "internal_admissible" or decision.c_hat is None:
        return None, decision.estimator, decision
    empirical_inflation_factor = decision.c_hat
    if not np.isfinite(empirical_inflation_factor) or empirical_inflation_factor < 1.0:
        raise ValueError(
            "Sibling empirical inflation factor must be finite and >= 1.0; "
            f"got {empirical_inflation_factor!r}."
        )

    if (
        not np.isfinite(sibling_test_record.reference_scale)
        or sibling_test_record.reference_scale <= 0.0
    ):
        raise ValueError(
            "Sibling reference_scale must be finite and positive before adjustment; "
            f"parent={sibling_test_record.parent!r}."
        )

    inflation_adjusted_statistic = sibling_test_record.stat / (
        sibling_test_record.reference_scale * empirical_inflation_factor
    )
    inflation_adjusted_degrees_of_freedom = float(sibling_test_record.degrees_of_freedom)
    if decision.p_value is None:
        raise ValueError(
            "Sibling empirical calibration decision did not include an adjusted "
            f"p-value; parent={sibling_test_record.parent!r}."
        )
    inflation_adjusted_p_value = decision.p_value
    return (
        (
            inflation_adjusted_statistic,
            inflation_adjusted_degrees_of_freedom,
            inflation_adjusted_p_value,
        ),
        decision.estimator,
        decision,
    )


def compute_inflation_adjusted_sibling_tests(
    sibling_test_records: Iterable[SiblingPairRecord],
    *,
    model: EmpiricalNullInflationModel | None,
    enforce_support_thresholds: bool = False,
    support_thresholds: CalibrationSupportThresholds = DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
) -> tuple[
    list[str],
    list[InflationAdjustedSiblingTestSummary],
    list[str],
    list[CalibrationDecision | None],
    list[tuple[SiblingPairRecord, CalibrationDecision]],
]:
    """Return summaries, method labels and calibration decisions for tested parents.

    A decision is None for a zero-dimensional record, which needs no calibration.
    Focal records whose decision is not internally admissible are returned last,
    with their decisions, so the caller can close their gates.
    """
    tested_parent_ids: list[str] = []
    inflation_adjusted_test_summaries: list[InflationAdjustedSiblingTestSummary] = []
    inflation_adjustment_method_labels: list[str] = []
    calibration_decisions: list[CalibrationDecision | None] = []
    inadmissible: list[tuple[SiblingPairRecord, CalibrationDecision]] = []
    family_support_cache: dict[str, dict[str, float | int]] = {}

    for sibling_test_record in sibling_test_records:
        if sibling_test_record.is_null_like:
            continue

        test_summary, method_label, decision = _compute_inflation_adjusted_sibling_test(
            sibling_test_record,
            model=model,
            enforce_support_thresholds=enforce_support_thresholds,
            support_thresholds=support_thresholds,
            family_support_cache=family_support_cache,
        )
        if test_summary is None:
            inadmissible.append((sibling_test_record, decision))
            continue
        tested_parent_ids.append(sibling_test_record.parent)
        inflation_adjusted_test_summaries.append(test_summary)
        inflation_adjustment_method_labels.append(method_label)
        calibration_decisions.append(decision)

    return (
        tested_parent_ids,
        inflation_adjusted_test_summaries,
        inflation_adjustment_method_labels,
        calibration_decisions,
        inadmissible,
    )


__all__ = [
    "InflationAdjustedSiblingTestSummary",
    "compute_inflation_adjusted_sibling_tests",
]
