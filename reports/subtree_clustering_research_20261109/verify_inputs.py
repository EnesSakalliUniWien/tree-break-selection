"""Regenerate declared inputs, independently of fitted rows, and check their hashes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from benchmarks.diagnostics.calibration.subtree_research import (
    RECOVERY_CELLS,
    data_seed,
    recovery_data,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(value, dtype):
    return hashlib.sha256(np.ascontiguousarray(value, dtype=dtype).tobytes()).hexdigest()


def main():
    checked = {}
    distinct = set()
    for panel in (
        "recovery",
        "factorial_cell_0",
        "factorial_cell_1",
        "factorial_cell_2",
        "preflight_isolated",
    ):
        rows = [json.loads(line) for line in (HERE / panel / "runs.jsonl").read_text().splitlines()]
        inputs = {}
        for row in rows:
            cell = row["cell"]
            rep = row.get("replicate", 0)
            key = (cell, rep)
            if key not in inputs:
                if panel == "recovery":
                    x, truth, seed = recovery_data(cell, rep)
                elif panel == "preflight_isolated":
                    c = RECOVERY_CELLS.index((row["p"], False, 5.0, False))
                    x, truth, seed = recovery_data(c, 0, n=row["n"])
                elif cell == 0:
                    seed = data_seed(phase=1, panel=3, cell=0, replicate=rep)
                    x = np.random.default_rng(seed).standard_normal((64, 12))
                    truth = np.zeros(64, dtype=int)
                else:
                    c = RECOVERY_CELLS.index((12, cell == 2, 5.0, False))
                    x, truth, seed = recovery_data(c, rep)
                inputs[key] = (seed, digest(x, np.float64), digest(truth, np.int64))
            assert inputs[key] == (row["seed"], row["data_hash"], row["truth_hash"]), row["run_id"]
            assert row["status"] == "ok", row["run_id"]
            if panel != "preflight_isolated":
                distinct.add(row["data_hash"])
        before_after = json.loads((HERE / panel / "verification.json").read_text())
        assert not before_after["source_files_changed_during_run"], panel
        checked[panel] = {
            "rows": len(rows),
            "distinct_panel_inputs": len(inputs),
            "all_hashes_match": True,
        }
    initial = json.loads((HERE / "environment.json").read_text())["source_hashes"]
    production = {p: h for p, h in initial.items() if p.startswith("tree_break_selection/")}
    changed = [
        p for p, h in production.items() if hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h
    ]
    assert not changed, changed
    assert len(distinct) == 850
    result = {
        "panels": checked,
        "distinct_development_inputs": len(distinct),
        "production_files_checked": len(production),
        "production_files_changed": changed,
    }
    (HERE / "input_integrity.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
