"""Instrumenting the spectral study must preserve native fits and retain failures."""

import gzip
import json

import numpy as np
import pytest

from benchmarks.diagnostics.calibration import subtree_research
from reports.subtree_spectral_alignment_20261109 import run
from tree_break_selection.space_separation import diffusion


@pytest.mark.parametrize("operator,link", [("symmetric", "average"), ("markov", "ward")])
def test_spectral_capture_preserves_native_pipeline_and_coordinates(operator, link):
    values, truth, _ = run.dataset("crossed", 0)
    untouched = values.copy()
    coordinate_function = diffusion.compute_diffusion_coordinates
    gate_function = subtree_research.run_gate_annotation_pipeline
    actual, detail, spec = run.capture_native(values, truth[2], operator, link)
    expected, _ = subtree_research.run_method(
        values, truth[2], operator=operator, linkage_method=link
    )
    assert actual["status"] == expected["status"] == "ok"
    assert actual["labels"] == expected["labels"]
    assert actual["root_sibling_p"] == expected["root_sibling_p"]
    assert diffusion.compute_diffusion_coordinates is coordinate_function
    assert subtree_research.run_gate_annotation_pipeline is gate_function
    np.testing.assert_array_equal(values, untouched)
    start, end = run.spectral_window(spec, values.shape[1])
    basis = spec.symmetric_vectors if operator == "symmetric" else spec.right_vectors
    coordinates = basis[:, start:end] * np.maximum(spec.eigenvalues[start:end], 0) ** 3
    np.testing.assert_allclose(
        run.pdist(coordinates), detail["baseline_diffusion_distance"], atol=1e-10
    )
    assert len(detail["native_context"]) == 3
    assert all(
        fields["principal_component_eigenvalues_by_node"] is not None
        for fields in detail["native_context"].values()
    )


def test_spectral_capture_restores_instrumentation_after_failure(monkeypatch):
    coordinate_function = diffusion.compute_diffusion_coordinates
    gate_function = subtree_research.run_gate_annotation_pipeline

    def fail(*args, **kwargs):
        raise RuntimeError("intentional native failure")

    monkeypatch.setattr(subtree_research, "run_method", fail)
    with pytest.raises(RuntimeError, match="intentional native failure"):
        run.capture_native(np.eye(4), np.zeros(4), "markov", "ward")
    assert diffusion.compute_diffusion_coordinates is coordinate_function
    assert subtree_research.run_gate_annotation_pipeline is gate_function


def test_replay_retains_each_failed_attempt_and_refuses_overwrite(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise RuntimeError("intentional replay failure")

    monkeypatch.setattr(run, "evaluate", fail)
    run.run(tmp_path, replicates=1)
    with gzip.open(tmp_path / "rows.jsonl.gz", "rt") as stream:
        records = [json.loads(line) for line in stream]
    assert len(records) == 8
    assert all(row["status"] == "error" and "replay failure" in row["reason"] for row in records)
    assert len({(row["case"], row["replicate"], row["operator"]) for row in records}) == 8
    assert json.loads((tmp_path / "verification.json").read_text())["error"] == 8
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        run.run(tmp_path, replicates=1)


def test_native_failure_keeps_original_reason_without_gate_context(monkeypatch):
    def fail(*args, **kwargs):
        return {"status": "error", "reason": "NNLS not applied: fixture"}, {}

    monkeypatch.setattr(subtree_research, "evaluate_tree", fail)
    values, truth, _ = run.dataset("crossed", 0)
    result, _ = run.evaluate(values, truth, "markov", "ward")
    assert result["native"]["status"] == "error"
    assert result["native"]["reason"] == "NNLS not applied: fixture"
    assert result["native_context"] == {}


def test_unexpected_child_refit_failure_is_error_and_does_not_hide_induced_result(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("unexpected child refit failure")

    monkeypatch.setattr(run, "capture_diffusion", fail)
    values, truth, _ = run.dataset("crossed", 0)
    result, _ = run.evaluate(values, truth, "markov", "ward")
    assert all(child["status"] == "error" for child in result["children"])
    assert all(child["induced_parent_graph"]["status"] == "ok" for child in result["children"])
    assert all(child["rebuilt_child_graph"]["status"] == "error" for child in result["children"])
