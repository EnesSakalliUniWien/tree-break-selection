from __future__ import annotations

import pandas as pd

import benchmarks.shared.util.case_execution as case_execution


def _successful_benchmark():
    return pd.DataFrame([{"method": "tbs", "ari": 1.0}])


class _Queue:
    def __init__(self):
        self.payloads: list[dict[str, object]] = []

    def put(self, payload):
        self.payloads.append(payload)


def test_run_case_worker_defaults_to_file_safe_matplotlib_backend(monkeypatch):
    captured: dict[str, object] = {}

    def _fake_benchmark_fn(**kwargs):
        captured.update(kwargs)
        return _successful_benchmark()

    monkeypatch.delenv("MPLBACKEND", raising=False)
    monkeypatch.setattr(case_execution, "_get_benchmark_fn", lambda: _fake_benchmark_fn)

    q = _Queue()
    case_execution._run_case_worker(
        q,
        case={"name": "tiny_case"},
        methods_to_test=["tbs"],
        case_plot_umap=True,
        case_plot_manifold=False,
        enable_plots=True,
        pdf_path="/tmp/tiny_case.pdf",
        include_validation_page=False,
    )

    assert q.payloads and q.payloads[0].get("ok") is True
    assert captured
    assert case_execution.os.environ["MPLBACKEND"] == "Agg"


def test_run_case_with_optional_isolation_disables_cover_pages(monkeypatch):
    captured: dict[str, object] = {}

    def _fake_benchmark_fn(**kwargs):
        captured.update(kwargs)
        return pd.DataFrame()

    monkeypatch.setattr(case_execution, "_get_benchmark_fn", lambda: _fake_benchmark_fn)

    case_execution.run_case_with_optional_isolation(
        case={"name": "tiny_case"},
        methods_to_test=["tbs"],
        case_plot_umap=False,
        case_plot_manifold=False,
        enable_plots=True,
        pdf_path="/tmp/tiny_case.pdf",
        isolate_umap_cases=True,
        timeout_sec=1,
        include_validation_page=False,
    )

    assert captured.get("include_cover_pages") is False
    assert captured.get("include_validation_page") is False


def test_run_case_worker_disables_cover_pages(monkeypatch):
    captured: dict[str, object] = {}

    def _fake_benchmark_fn(**kwargs):
        captured.update(kwargs)
        return _successful_benchmark()

    monkeypatch.setattr(case_execution, "_get_benchmark_fn", lambda: _fake_benchmark_fn)

    q = _Queue()
    case_execution._run_case_worker(
        q,
        case={"name": "tiny_case"},
        methods_to_test=["tbs"],
        case_plot_umap=True,
        case_plot_manifold=False,
        enable_plots=True,
        pdf_path="/tmp/tiny_case.pdf",
        include_validation_page=False,
    )

    assert captured.get("include_cover_pages") is False
    assert captured.get("include_validation_page") is False
    assert q.payloads and q.payloads[0].get("ok") is True


def _large_payload_worker(queue, *_args) -> None:
    rows = [{f"metric_{i}": float(row * i) for i in range(40)} for row in range(300)]
    queue.put({"ok": True, "rows": rows})


def test_run_case_isolated_returns_payload_larger_than_pipe_buffer(monkeypatch):
    monkeypatch.setattr(case_execution, "_run_case_worker", _large_payload_worker)

    result = case_execution.run_case_isolated(
        case={"name": "large_case"},
        methods_to_test=["tbs"],
        method_params=None,
        case_plot_umap=True,
        case_plot_manifold=False,
        enable_plots=False,
        pdf_path=None,
        timeout_sec=60,
        include_validation_page=False,
    )

    assert result.shape == (300, 40)


def _lingering_worker(queue, *_args) -> None:
    import threading
    import time

    queue.put({"ok": True, "rows": [{"metric": 1.0}]})
    # A non-daemon thread keeps the child process alive after the result is sent.
    threading.Thread(target=time.sleep, args=(60,)).start()


def test_run_case_isolated_stops_a_child_that_lingers_after_its_result(monkeypatch):
    import time

    monkeypatch.setattr(case_execution, "_run_case_worker", _lingering_worker)

    start = time.monotonic()
    result = case_execution.run_case_isolated(
        case={"name": "lingering_case"},
        methods_to_test=["tbs"],
        method_params=None,
        case_plot_umap=True,
        case_plot_manifold=False,
        enable_plots=False,
        pdf_path=None,
        timeout_sec=5,
        include_validation_page=False,
    )

    assert result.to_dict(orient="records") == [{"metric": 1.0}]
    assert time.monotonic() - start < 30
