"""Resume/grid coverage helpers for benchmark run orchestrators."""

from __future__ import annotations

import pandas as pd

from benchmarks.shared.benchmark_grid import benchmark_run_id


def compute_resume_coverage(
    existing_results: pd.DataFrame,
    expected_run_ids: list[str],
) -> tuple[set[str], dict[str, list[str]], int]:
    """Compute fully-complete and partial case coverage from an existing CSV.

    Coverage is keyed by ``case_id`` so a suite whose case order changed between
    runs does not credit one case's rows to another.
    """
    if existing_results.empty:
        return set(), {}, 0
    missing_columns = {"case_id", "run_id"} - set(existing_results.columns)
    if missing_columns:
        raise ValueError(
            f"Existing results lack {sorted(missing_columns)}; they predate run-id resume. "
            "Write this run to a new TBS_RUN_DIR."
        )

    expected_set = set(expected_run_ids)
    run_ids = existing_results["run_id"].astype(str)
    progress = pd.DataFrame(
        {
            "case_id": existing_results["case_id"],
            "run_id": run_ids.where(run_ids.isin(expected_set)),
        }
    ).dropna()
    if progress.empty:
        return set(), {}, 0

    run_ids_by_case = progress.groupby(progress["case_id"].astype(str))["run_id"].agg(set).to_dict()

    completed_cases: set[str] = set()
    missing_run_ids_by_case: dict[str, list[str]] = {}
    for case_id, seen_run_ids in run_ids_by_case.items():
        missing = [run_id for run_id in expected_run_ids if run_id not in seen_run_ids]
        if not missing:
            completed_cases.add(case_id)
        else:
            missing_run_ids_by_case[case_id] = missing

    return completed_cases, missing_run_ids_by_case, len(run_ids_by_case)


def require_result_columns(
    existing_columns: list[str],
    new_columns: list[str],
    *,
    output_path: object,
) -> None:
    """Refuse to append rows whose columns differ from the results file header."""
    if existing_columns and list(new_columns) != list(existing_columns):
        added = sorted(set(new_columns) - set(existing_columns))
        dropped = sorted(set(existing_columns) - set(new_columns))
        raise ValueError(
            f"{output_path} has a different column layout (new: {added}, missing: {dropped}, "
            "or a different order). Write this run to a new TBS_RUN_DIR."
        )


def run_ids_for_params(
    methods: list[str],
    param_sets: dict[str, list[dict[str, object]]],
) -> list[str]:
    """Return canonical benchmark run IDs for selected methods and grids."""
    return [
        benchmark_run_id(method_id, params)
        for method_id in methods
        for params in param_sets[method_id]
    ]


def filter_methods_and_params_by_run_ids(
    *,
    methods: list[str],
    param_sets: dict[str, list[dict[str, object]]],
    run_ids: set[str],
) -> tuple[list[str], dict[str, list[dict[str, object]]]]:
    """Restrict methods and grids to the requested run IDs."""
    filtered_methods: list[str] = []
    filtered_params: dict[str, list[dict[str, object]]] = {}
    for method_id in methods:
        selected_params = [
            params
            for params in param_sets[method_id]
            if benchmark_run_id(method_id, params) in run_ids
        ]
        if selected_params:
            filtered_methods.append(method_id)
            filtered_params[method_id] = selected_params
    return filtered_methods, filtered_params


__all__ = [
    "compute_resume_coverage",
    "filter_methods_and_params_by_run_ids",
    "require_result_columns",
    "run_ids_for_params",
]
