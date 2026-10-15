from __future__ import annotations

import pandas as pd
import pytest

from benchmarks.shared.benchmark_grid import benchmark_grid, benchmark_run_id
from benchmarks.shared.benchmark_runs.resume import (
    compute_resume_coverage,
    filter_methods_and_params_by_run_ids,
    run_ids_for_params,
)


def test_full_runner_resume_tracks_grid_run_ids():
    param_sets = {
        "kmeans": benchmark_grid(
            benchmark_class="canonical",
            grid_name="toy_grid",
            base_params={"n_clusters": "true", "n_init": 10},
            axes={"variant": ("a", "b")},
        )
    }
    expected_run_ids = run_ids_for_params(["kmeans"], param_sets)
    existing_results = pd.DataFrame(
        [
            {"test_case": 1, "case_id": "a", "run_id": expected_run_ids[0]},
            {"test_case": 1, "case_id": "a", "run_id": expected_run_ids[1]},
            {"test_case": 2, "case_id": "b", "run_id": expected_run_ids[0]},
        ]
    )

    completed, missing_by_case, tracked_count = compute_resume_coverage(
        existing_results,
        expected_run_ids,
    )

    assert completed == {"a"}
    assert missing_by_case == {"b": [expected_run_ids[1]]}
    assert tracked_count == 2

    methods, filtered_params = filter_methods_and_params_by_run_ids(
        methods=["kmeans"],
        param_sets=param_sets,
        run_ids={expected_run_ids[1]},
    )

    assert methods == ["kmeans"]
    assert filtered_params == {"kmeans": [param_sets["kmeans"][1]]}
    assert benchmark_run_id("kmeans", filtered_params["kmeans"][0]) == expected_run_ids[1]


def test_resume_coverage_follows_case_names_when_suite_order_changes():
    existing_results = pd.DataFrame(
        [
            {"test_case": 1, "case_id": "old_first", "run_id": "kmeans::default"},
            {"test_case": 2, "case_id": "old_second", "run_id": "kmeans::default"},
        ]
    )

    completed, missing_by_case, _ = compute_resume_coverage(
        existing_results,
        ["kmeans::default"],
    )

    # A case inserted at position 1 must not inherit old_first's rows.
    assert completed == {"old_first", "old_second"}
    assert "inserted_case" not in completed
    assert missing_by_case == {}


def test_resume_rejects_results_without_run_ids():
    existing_results = pd.DataFrame({"test_case": [1], "method": ["kmeans"], "ari": [1.0]})

    with pytest.raises(ValueError, match="predate run-id resume"):
        compute_resume_coverage(existing_results, ["kmeans::default"])
