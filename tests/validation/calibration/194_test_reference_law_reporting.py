from __future__ import annotations

import pytest

from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from benchmarks.validation.sweeps.nnls_null_calibration_sweep import (
    BRANCH_SOURCE_LINKAGE,
    SPECTRAL_CONTEXT_LEAF_ONLY,
    _run_branch_source,
)


@pytest.fixture(scope="module")
def nnls_linkage_row() -> dict[str, object]:
    case = next(
        case for case in get_test_cases_by_suite("binary") if case["name"] == "binary_low_noise_2c"
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
    assert nnls_linkage_row["schema_version"] == "nnls_null_calibration_sweep/v6"
