"""Research panels retain native parity, paired inputs and failed prerequisites."""

import itertools
import json

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.subtree_research import (
    NULL_CELLS,
    RECOVERY_CELLS,
    _frame,
    data_seed,
    null_data,
    recovery_data,
    run_method,
    run_panel,
    runtime_gate,
)
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS


def test_frozen_grids_and_seed_phases_are_distinct():
    assert len(RECOVERY_CELLS) == 16
    assert len(NULL_CELLS) == 24
    keys = [
        data_seed(phase=p, panel=2, cell=c, replicate=r)
        for p, c, r in itertools.product((1, 2), range(24), range(5))
    ]
    assert len(set(keys)) == len(keys)
    for cell in range(24):
        values, path, sigma, _ = null_data(cell, 0)
        n, p, cov, expected_path = NULL_CELLS[cell]
        assert values.shape == (n, p) and path == expected_path
        assert np.linalg.eigvalsh(sigma).min() > 0
        if cov == "diagonal100":
            assert np.linalg.cond(sigma) == pytest.approx(100)
    cell = RECOVERY_CELLS.index((12, True, 5.0, False))
    x, truth, _ = recovery_data(cell, 0)
    assert np.bincount(truth).tolist() == [132, 60, 36, 12]
    assert not np.array_equal(x, recovery_data(cell, 0, phase=2)[0])


def test_native_default_parity_and_truth_never_changes_fitting():
    cell = RECOVERY_CELLS.index((12, False, 5.0, False))
    x, truth, _ = recovery_data(cell, 0, n=32)
    data, space = _frame(x)
    params = dict(METHOD_SPECS["tbs_family_diffusion"].param_grid[0])
    reference = run_clustering_result(data, "tbs_family_diffusion", params, feature_space=space)
    actual, details = run_method(x, truth, trace=True)
    assert actual["status"] == reference.status.value
    if reference.labels is not None:
        np.testing.assert_array_equal(actual["labels"], reference.labels)
    assert len(details["raw_records"]) == len(x) - 1
    assert len(details["subtree_profiles"]) == len(x) - 1
    altered, _ = run_method(x, np.zeros(len(x), dtype=int))
    assert altered["status"] == actual["status"]
    if "labels" in actual:
        assert altered["labels"] == actual["labels"]


def _timings():
    return [
        dict(
            n=n,
            p=p,
            method=method,
            repeat=repeat,
            seed=7,
            data_hash="same",
            status="ok",
            elapsed_sec=1 if method == "baseline" else 4,
        )
        for n, p, repeat, method in itertools.product(
            (32, 128, 240), (12, 72), range(3), ("baseline", "candidate")
        )
    ]


def test_runtime_gate_requires_complete_genuinely_paired_successes():
    rows = _timings()
    assert runtime_gate(rows)["passed"]
    assert json.loads(json.dumps(runtime_gate(rows)))["passed"]
    assert not runtime_gate(rows[:-1])["passed"]
    for key, value in (
        ("repeat", 1),
        ("seed", 8),
        ("data_hash", "different"),
        ("status", "error"),
        ("elapsed_sec", np.nan),
    ):
        corrupted = [dict(row) for row in rows]
        corrupted[0][key] = value
        assert not runtime_gate(corrupted)["passed"], key
    slow = [dict(row) for row in rows]
    for row in slow:
        if row["method"] == "candidate":
            row["elapsed_sec"] = 6
    assert not runtime_gate(slow)["passed"]


def test_invalid_panel_does_not_create_results(tmp_path):
    target = tmp_path / "invalid"
    for stage, cells, reps in (
        ("unknown", None, 1),
        ("recovery", [-1], 1),
        ("factorial", [0, 0], 1),
        ("recovery", [0], 0),
    ):
        with pytest.raises(ValueError):
            run_panel(stage, target, cells=cells, replicates=reps)
        assert not target.exists()
