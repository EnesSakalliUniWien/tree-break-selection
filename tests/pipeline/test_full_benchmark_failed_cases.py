"""A failing case must not stop the full benchmark from running the rest."""

import pytest

import benchmarks.full.run as full_run
from benchmarks.shared.benchmark_runs.runtime import THREAD_ENV_VARS
from benchmarks.shared.cases import get_test_cases_by_suite


def test_failed_cases_are_skipped_then_reported(tmp_path, monkeypatch):
    cases = get_test_cases_by_suite("binary")[:2]
    attempted: list[str] = []

    def fail_case(*, case, **_kwargs):
        attempted.append(case["name"])
        raise RuntimeError("Case subprocess timed out after 1s")

    monkeypatch.setattr(
        full_run, "get_test_cases_by_suite", lambda _suite: [dict(c) for c in cases]
    )
    monkeypatch.setattr(full_run, "run_case_with_optional_isolation", fail_case)
    monkeypatch.setenv("TBS_RUN_DIR", str(tmp_path))
    monkeypatch.setenv("TBS_METHODS", "kmeans")
    monkeypatch.setenv("TBS_ENABLE_PLOTS", "0")
    monkeypatch.setenv("TBS_RUN_RELATIONSHIP_ANALYSIS", "0")
    # run_benchmarks setdefaults these; registering them restores the caller's values.
    for name in (*THREAD_ENV_VARS, "TBS_N_JOBS"):
        monkeypatch.setenv(name, "1")

    with pytest.raises(SystemExit, match="2 case\\(s\\) failed") as excinfo:
        full_run.run_benchmarks()

    assert attempted == [case["name"] for case in cases]
    assert all(case["name"] in str(excinfo.value) for case in cases)
    assert not any(tmp_path.glob("*_benchmark_comparison.csv"))
