from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from benchmarks.diagnostics.failure.debug_trace import analyze_single_case


def _write_audit(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def test_analyze_single_case_rejects_missing_sibling_decision_column(tmp_path: Path) -> None:
    audit_path = tmp_path / "case_1_tbs_stats.csv"
    _write_audit(
        audit_path,
        [
            {
                "node_id": "root",
                "leaf_count": 4,
                "parent_node": None,
                "Sibling_Divergence_P_Value": None,
            },
            {
                "node_id": "left",
                "leaf_count": 2,
                "parent_node": "root",
                "Sibling_Divergence_P_Value": 0.001,
            },
            {
                "node_id": "right",
                "leaf_count": 2,
                "parent_node": "root",
                "Sibling_Divergence_P_Value": 0.001,
            },
        ],
    )

    diagnosis = analyze_single_case(audit_path)

    assert diagnosis == {
        "mode": "ERROR",
        "reason": "Invalid audit contract; missing columns: Sibling_BH_Different",
    }


def test_analyze_single_case_uses_sibling_decision_column(tmp_path: Path) -> None:
    audit_path = tmp_path / "case_1_tbs_stats.csv"
    _write_audit(
        audit_path,
        [
            {
                "node_id": "root",
                "leaf_count": 4,
                "parent_node": None,
                "Sibling_Divergence_P_Value": 0.35,
                "Sibling_BH_Different": False,
            },
            {
                "node_id": "left",
                "leaf_count": 2,
                "parent_node": "root",
                "Sibling_Divergence_P_Value": 0.001,
                "Sibling_BH_Different": True,
            },
            {
                "node_id": "right",
                "leaf_count": 2,
                "parent_node": "root",
                "Sibling_Divergence_P_Value": 0.001,
                "Sibling_BH_Different": True,
            },
        ],
    )

    diagnosis = analyze_single_case(audit_path)

    assert diagnosis == {
        "mode": "**UNDER-SPLIT**",
        "reason": "Root split rejected (P=3.50e-01)",
    }


@pytest.mark.parametrize(
    ("root_fields", "reason"),
    [
        (
            {
                "Sibling_Divergence_Skipped": True,
                "Sibling_Divergence_Invalid": True,
                "Sibling_Test_Method": "empirical_null_no_internal_support",
                "Sibling_Divergence_P_Value_Corrected": None,
            },
            "Root sibling test skipped (empirical_null_no_internal_support)",
        ),
        (
            {"Sibling_Divergence_Skipped": True, "Sibling_Divergence_Invalid": False},
            "Root sibling test skipped",
        ),
        (
            {"Sibling_Divergence_Skipped": False, "Sibling_Divergence_Invalid": True},
            "Root sibling test invalid",
        ),
    ],
)
def test_analyze_single_case_does_not_call_a_skipped_root_test_rejected(
    tmp_path: Path, root_fields: dict[str, object], reason: str
) -> None:
    audit_path = tmp_path / "case_1_tbs_stats.csv"
    child = {"leaf_count": 2, "parent_node": "root", "Sibling_BH_Different": False}
    _write_audit(
        audit_path,
        [
            {
                "node_id": "root",
                "leaf_count": 4,
                "parent_node": None,
                # The raw p-value of a test the gate did not use must not be reported.
                "Sibling_Divergence_P_Value": 5.29e-84,
                "Sibling_BH_Different": False,
                **root_fields,
            },
            {"node_id": "left", **child},
            {"node_id": "right", **child},
        ],
    )

    assert analyze_single_case(audit_path) == {"mode": "**UNDER-SPLIT**", "reason": reason}
