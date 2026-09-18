from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from benchmarks.diagnostics.calibration.statistics.internal_support_threshold_validation import (
    evaluate_internal_support_threshold_rows,
    summarize_internal_support_threshold_rows,
)
from benchmarks.diagnostics.calibration.statistics.mixed_internal_calibration_sweep import (
    _collect_case_replicate,
    run_mixed_internal_calibration_sweep,
)
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from benchmarks.validation.sweeps.nnls_null_calibration_sweep import (
    BRANCH_SOURCE_LINKAGE,
    SPECTRAL_CONTEXT_LEAF_ONLY,
    _run_branch_source,
)


def test_mixed_sweep_reports_restored_empirical_p_values() -> None:
    case = next(
        case
        for case in get_test_cases_by_suite("binary")
        if case["name"] == "binary_low_noise_2c"
    )

    _q10_rows, threshold_rows, status = _collect_case_replicate(
        case=case,
        replicate_index=0,
        edge_alpha=0.001,
        sibling_alpha=0.01,
    )

    assert status["model_status"] == "fit"
    assert not threshold_rows.empty
    unresolved = threshold_rows[
        threshold_rows["calibration_status"].eq("internal_admissible")
    ]
    assert not unresolved.empty
    assert unresolved["calibration_p_value"].between(0.0, 1.0).all()
    assert unresolved["split_rejected_at_alpha"].notna().all()


def test_threshold_summary_excludes_unavailable_split_outcomes() -> None:
    contexts = pd.DataFrame(
        {
            "context_id": ["unavailable_null", "available_null", "available_signal"],
            "n_supported_groups": [100, 100, 100],
            "n_family_supported_groups": [100, 100, 100],
            "family_effective_sample_size": [100.0, 100.0, 100.0],
            "local_effective_sample_size": [100.0, 100.0, 100.0],
            "local_max_group_weight_share": [0.01, 0.01, 0.01],
            "leave_one_group_max_delta_log_c": [0.0, 0.0, 0.0],
            "is_null_context": [True, True, False],
            "is_signal_context": [False, False, True],
            "split_rejected_at_alpha": [np.nan, True, True],
        }
    )

    decisions = evaluate_internal_support_threshold_rows(contexts)
    summary = summarize_internal_support_threshold_rows(decisions)

    for _profile_id, group in decisions.groupby("threshold_profile_id", sort=False):
        assert group["has_outcome_labels"].tolist() == [False, True, True]
        assert pd.isna(group.iloc[0]["split_rejected_at_alpha"])
    assert summary["n_outcome_available_contexts"].eq(2).all()
    assert summary["n_admissible_null_contexts_with_outcome"].eq(1).all()
    assert summary["n_admissible_signal_contexts_with_outcome"].eq(1).all()
    assert summary["null_false_split_rate_among_admissible"].eq(1.0).all()
    assert summary["signal_retention_rate_among_admissible"].eq(1.0).all()
    assert summary["outcome_status"].eq("has_mixed_null_signal_outcomes").all()


def test_mixed_sweep_resume_rejects_stale_withheld_p_values(tmp_path: Path) -> None:
    pd.DataFrame(
        {
            "left_edge_bh_p_value": [0.5],
            "right_edge_bh_p_value": [0.5],
            "selected_hierarchy_ratio": [1.0],
            "is_null_like": [True],
            "is_edge_blocked": [False],
        }
    ).to_csv(tmp_path / "mixed_q10_sibling_records.csv", index=False)
    pd.DataFrame(
        {
            "context_id": ["stale_false_outcome"],
            "n_supported_groups": [100],
            "n_family_supported_groups": [100],
            "family_effective_sample_size": [100.0],
            "local_effective_sample_size": [100.0],
            "local_max_group_weight_share": [0.01],
            "leave_one_group_max_delta_log_c": [0.0],
            "is_null_context": [True],
            "is_signal_context": [False],
            "split_rejected_at_alpha": [False],
        }
    ).to_csv(tmp_path / "mixed_internal_support_contexts.csv", index=False)
    pd.DataFrame({"status": ["stale"]}).to_csv(
        tmp_path / "mixed_internal_calibration_status.csv",
        index=False,
    )
    (tmp_path / "manifest.json").write_text(
        json.dumps({"schema_version": "mixed_internal_calibration_sweep/v2"}) + "\n",
        encoding="utf-8",
    )

    outputs = run_mixed_internal_calibration_sweep(
        suite="binary",
        case_names=("binary_low_noise_2c",),
        n_replicates=1,
        output_dir=tmp_path,
        resume=True,
    )

    regenerated = pd.read_csv(outputs["threshold_contexts"])
    assert "stale_false_outcome" not in set(regenerated["context_id"])
    unresolved = regenerated[
        regenerated["calibration_status"].eq("internal_admissible")
    ]
    assert not unresolved.empty
    assert unresolved["split_rejected_at_alpha"].notna().all()


@pytest.fixture(scope="module")
def nnls_linkage_row() -> dict[str, object]:
    case = next(
        case
        for case in get_test_cases_by_suite("binary")
        if case["name"] == "binary_low_noise_2c"
    )
    inputs = prepare_case_inputs(case, ["tbs"])

    row, _record_rows, _eigenvalue_rows = _run_branch_source(
        data=inputs.data,
        labels=inputs.labels,
        feature_space=inputs.metadata.get("feature_space"),
        branch_source=BRANCH_SOURCE_LINKAGE,
        spectral_context=SPECTRAL_CONTEXT_LEAF_ONLY,
        edge_alpha=0.001,
        sibling_alpha=0.01,
        tree_linkage_method="average",
        pair_sample_size=None,
    )
    return row


def test_nnls_sweep_counts_restored_empirical_tests(
    nnls_linkage_row: dict[str, object],
) -> None:
    assert nnls_linkage_row["n_sibling_tested"] > 0
    assert nnls_linkage_row["n_sibling_fail_closed"] == 0


def test_nnls_sweep_marks_changed_reporting_contract_as_v5(
    nnls_linkage_row: dict[str, object],
) -> None:
    assert nnls_linkage_row["schema_version"] == "nnls_null_calibration_sweep/v5"
