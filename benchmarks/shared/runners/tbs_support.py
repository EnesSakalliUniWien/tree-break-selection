"""Recognized scientific-support outcomes for TBS benchmark runners."""

from __future__ import annotations

import pandas as pd

from benchmarks.shared.types import (
    UnsupportedEvidence,
    UnsupportedReason,
    UnsupportedReasonCode,
)

EMPIRICAL_NULL_GATE_METHOD = "projected_wald_inflation"
UNDEFINED_INTERNAL_SUPPORT = "undefined_no_internal_support"
UNVALIDATED_REFERENCE_LAW = "undefined_unvalidated_reference_law"


def _require_columns(annotations: pd.DataFrame, columns: tuple[str, ...]) -> None:
    missing = [column for column in columns if column not in annotations.columns]
    if missing:
        raise ValueError(
            "TBS unsupported-outcome detection requires production-stamped "
            f"annotation columns; missing={missing}."
        )


def unsupported_empirical_null_reason(
    annotations: pd.DataFrame,
    *,
    sibling_gate_method: str,
) -> UnsupportedReason | None:
    """Return the typed no-support reason stamped by the empirical-null gate."""
    if str(sibling_gate_method) != EMPIRICAL_NULL_GATE_METHOD:
        return None

    calibration_column = "Sibling_Gate_P_Value_Calibration"
    _require_columns(annotations, (calibration_column,))
    calibration_values = annotations[calibration_column].astype(str)
    no_support_mask = calibration_values.eq(UNDEFINED_INTERNAL_SUPPORT)
    unvalidated_law_mask = calibration_values.eq(UNVALIDATED_REFERENCE_LAW)
    focal_mask = no_support_mask | unvalidated_law_mask
    focal_record_count = int(focal_mask.sum())
    if focal_record_count == 0:
        return None

    required = (
        "Sibling_Role_Supported",
        "Sibling_Divergence_Invalid",
        "Child_Parent_Divergence_Tested",
        "Child_Parent_Divergence_Significant",
    )
    _require_columns(annotations, required)
    supported_roles = annotations["Sibling_Role_Supported"].fillna(False).astype(bool)
    focal_supported_count = int(supported_roles.loc[focal_mask].sum())
    if focal_supported_count != 0:
        raise ValueError(
            "TBS annotations are inconsistent: fail-closed empirical-null focal rows "
            f"were stamped with {focal_supported_count} supported roles."
        )
    admissible_support_count = int(supported_roles.sum())
    has_unvalidated_law = bool(unvalidated_law_mask.any())
    if not has_unvalidated_law and admissible_support_count != 0:
        return None

    invalid_record_count = int(
        annotations.loc[focal_mask, "Sibling_Divergence_Invalid"]
        .fillna(False)
        .astype(bool)
        .sum()
    )
    upstream_tested_count = int(
        annotations["Child_Parent_Divergence_Tested"].fillna(False).astype(bool).sum()
    )
    upstream_rejected_count = int(
        annotations["Child_Parent_Divergence_Significant"]
        .fillna(False)
        .astype(bool)
        .sum()
    )
    reason_code = (
        UnsupportedReasonCode.EMPIRICAL_NULL_UNVALIDATED_REFERENCE_LAW
        if has_unvalidated_law
        else UnsupportedReasonCode.EMPIRICAL_NULL_NO_INTERNAL_SUPPORT
    )
    message = (
        "The selected hierarchy has internal calibration support, but its "
        "conditional selected-tail reference law is not validated."
        if has_unvalidated_law
        else (
            "The selected hierarchy contains focal sibling tests but no "
            "admissible internal empirical-null calibration support."
        )
    )
    return UnsupportedReason(
        code=reason_code,
        stage="sibling_calibration",
        message=message,
        evidence=UnsupportedEvidence(
            focal_record_count=focal_record_count,
            admissible_support_count=admissible_support_count,
            invalid_record_count=invalid_record_count,
            upstream_tested_count=upstream_tested_count,
            upstream_rejected_count=upstream_rejected_count,
        ),
    )


__all__ = ["unsupported_empirical_null_reason"]
