"""Projection-dimension annotations for inflated projected-Wald sibling tests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..inflation_correction.types.inflation_model import CalibrationDecision
from ..pair_testing.types.sibling_pair_record import SiblingPairRecord


def write_record_projection_dimensions(
    annotations_dataframe: pd.DataFrame,
    records: list[SiblingPairRecord],
) -> None:
    """Write sibling projection dimensions into the annotations frame."""
    record_parents = [record.parent for record in records]
    annotations_dataframe.loc[record_parents, "Sibling_Projection_Dimension"] = [
        record.sibling_projection_dimension for record in records
    ]


def write_record_calibration_evidence(
    annotations_dataframe: pd.DataFrame,
    records: list[SiblingPairRecord],
) -> None:
    """Write per-record calibration-support evidence for every sibling record.

    Skipped parents are never re-annotated by the focal test path, so their
    degrees of freedom would otherwise be absent from the audit frame and the
    calibration-support population could not be reconstructed.
    """
    record_parents = [record.parent for record in records]
    annotations_dataframe.loc[record_parents, "Sibling_Degrees_of_Freedom"] = [
        float(record.degrees_of_freedom) for record in records
    ]
    annotations_dataframe.loc[record_parents, "Sibling_Role_Supported"] = [
        record.has_internal_empirical_null_support for record in records
    ]
    annotations_dataframe.loc[record_parents, "Sibling_Null_Weight"] = [
        float(record.sibling_null_weight) for record in records
    ]
    annotations_dataframe.loc[record_parents, "Sibling_Parent_Positive_Eigenvalue_Count"] = [
        float(record.parent_positive_eigenvalue_count) for record in records
    ]


def write_inflation_evidence(
    annotations_dataframe: pd.DataFrame,
    records: list[SiblingPairRecord],
    tested_parent_ids: list[str],
    decisions: list[CalibrationDecision | None],
) -> None:
    """Write the raw statistic and inflation factors behind each tested sibling p-value.

    The stored test statistic is already divided by the reference scale and applied
    factor. With these columns, chi2.sf(raw / (reference_scale * factor), df) can be
    re-evaluated under the local or family-baseline factor without rerunning the tree.
    """
    by_parent = {record.parent: record for record in records}
    for parent, decision in zip(tested_parent_ids, decisions, strict=True):
        record = by_parent[parent]
        strata = {} if decision is None else decision.descriptive_strata
        annotations_dataframe.loc[
            parent,
            [
                "Sibling_Raw_Test_Statistic",
                "Sibling_Reference_Scale",
                "Sibling_Inflation_Factor",
                "Sibling_Local_Inflation_Factor",
                "Sibling_Family_Baseline_Inflation_Factor",
            ],
        ] = [
            float(record.stat),
            float(record.reference_scale),
            np.nan if decision is None else float(decision.c_hat),
            float(strata.get("local_inflation_factor", np.nan)),
            float(strata.get("family_baseline_inflation_factor", np.nan)),
        ]


__all__ = [
    "write_inflation_evidence",
    "write_record_calibration_evidence",
    "write_record_projection_dimensions",
]
