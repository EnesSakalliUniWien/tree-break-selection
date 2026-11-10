"""Native context capture must preserve decisions and structural absence."""

import json

import numpy as np
import pytest

from reports.subtree_context_statistics_20261109 import run


@pytest.fixture(scope="module")
def replay():
    values, truth, _ = run.dataset("crossed", 0)
    proposals = list(run.proposals(values))
    results = []
    for name, z, coords, basis, _ in (proposals[0], proposals[-1]):
        results.append(
            (run.topology_from_linkage(z), run.fit_tree(values, truth[2], z, coords, basis))
        )
    return values, truth[2], proposals[0], results


def test_capture_preserves_native_pipeline_and_distinguishes_basis_halves(replay):
    values, truth, (_, z, _, _, _), results = replay
    untouched = values.copy()
    data, space = run.research._frame(values)
    expected, _ = run.research.evaluate_tree(data, truth, space, z, trace=True)
    actual = results[0][1]
    assert actual["metrics"] == expected
    np.testing.assert_array_equal(values, untouched)
    root = actual["nodes"][run.key(range(len(values)))]
    assert len(root["outgoing_basis_rows"]) == 32
    assert not np.allclose(
        root["sibling_all_leaf_context"]["eigenvalues"],
        root["outgoing_edge_half_context"]["eigenvalues"],
    )
    assert sum(len(edge["test_rows"]) for edge in root["outgoing_edges"]) == 32
    for node in actual["nodes"].values():
        if node["sibling_gate_state"].startswith("active_"):
            assert np.isfinite(node["annotations"]["Sibling_Divergence_P_Value_Corrected"])


def test_same_tree_matches_context_and_root_has_no_incoming_comparison(replay):
    _, _, _, results = replay
    tree, result = results[0]
    rows = list(run.compare(tree, tree, result, result))
    assert len(rows) == 63
    assert all(row["node_exists"] and row["direct_children_preserved"] for row in rows)
    assert not any(row[f"case_{i}"] for row in rows for i in range(1, 5))
    for row in rows:
        if row["original_parent"] is None:
            assert row["same_incoming_test_rows"] is None
            assert row["reference_Child_Parent_Divergence_P_Value"] is None
        else:
            assert row["same_incoming_test_rows"] and row["same_incoming_basis_rows"]
        assert row["reference_collected_sibling_stat"] == row["alternative_collected_sibling_stat"]


def test_absent_node_has_no_alternative_test_or_traversal_decision(replay):
    _, _, _, results = replay
    (rt, r), (at, a) = results
    rows = list(run.compare(rt, at, r, a))
    absent = [row for row in rows if not row["node_exists"]]
    assert absent
    for row in absent:
        assert row["case_1"]
        for name in (
            "Sibling_Divergence_P_Value",
            "Child_Parent_Divergence_P_Value",
            "collected_sibling_stat",
            "sibling_gate_state",
            "visited",
            "final_boundary",
        ):
            assert row["alternative_" + name] is None


def test_capture_restores_functions_after_exception(monkeypatch):
    spectral = run.EDGE_MODULE.compute_child_parent_spectral_context
    gate = run.research.run_gate_annotation_pipeline

    def fail(*args, **kwargs):
        raise RuntimeError("intentional failure")

    monkeypatch.setattr(run.research, "evaluate_tree", fail)
    with pytest.raises(RuntimeError, match="intentional failure"):
        run.fit_tree(np.eye(4), np.zeros(4), None, None, None)
    assert run.EDGE_MODULE.compute_child_parent_spectral_context is spectral
    assert run.research.run_gate_annotation_pipeline is gate


def test_native_failure_without_tree_keeps_reason(monkeypatch):
    def fail(*args, **kwargs):
        return {"status": "error", "reason": "NNLS fixture"}, {}

    monkeypatch.setattr(run.research, "evaluate_tree", fail)
    result = run.fit_tree(np.eye(4), np.zeros(4), None, None, None)
    assert result["metrics"]["reason"] == "NNLS fixture"
    assert result["nodes"] == {}


def test_panel_retains_failures_and_refuses_overwrite(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise RuntimeError("intentional replay failure")

    monkeypatch.setattr(run, "fit_tree", fail)
    run.run(tmp_path, replicates=1)
    attempts = json.loads((tmp_path / "attempts.json").read_text())
    assert len(attempts) == 40
    assert all(row["status"] == "error" and "replay failure" in row["reason"] for row in attempts)
    assert json.loads((tmp_path / "verification.json").read_text())["node_comparisons"] == 0
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        run.run(tmp_path, replicates=1)


@pytest.mark.parametrize(
    "eigenvalues,blocks,input_rank,expected_rank",
    [
        ([1, 1e-14, 0, 0, 0, 0], ((0, 1), (1, 6)), 2, 1),
        ([2, 2, 1], ((0, 2), (2, 3)), 1, 2),
        ([0, 0], ((0, 2),), 2, 0),
    ],
)
def test_local_covariance_span_keeps_whole_positive_blocks(
    eigenvalues, blocks, input_rank, expected_rank
):
    eigenvalues = np.array(eigenvalues, dtype=float)
    basis = np.eye(len(eigenvalues))
    spectrum = run.align.EigenSystem(np.diag(eigenvalues), eigenvalues, basis, blocks)
    result = run.covariance_span_record(spectrum, basis[:, :input_rank], np.zeros((3, input_rank)))
    assert result["inherited_vs_rebuilt_covariance_span"]["rank_b"] == expected_rank
    if not eigenvalues.any():
        assert result["inherited_local_variance_share"] is None
