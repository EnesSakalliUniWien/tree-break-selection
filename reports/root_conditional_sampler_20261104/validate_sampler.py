"""Frozen parity, fresh-seed and registered-case validation of the root chain API."""

import hashlib
import json
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_conditional_mc import (
    selected_diffusion_root_chain_test,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    _geometry,
    make_diffusion_root_path,
)
from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.util.case_inputs import prepare_case_inputs

HERE = Path(__file__).resolve().parent
PREVIOUS = ROOT / "reports/selective_rebuild_comparison_20261104"


def save_json(path, data):
    path.write_text(
        json.dumps(data, indent=2, default=diagnostic_json_default, allow_nan=False) + "\n"
    )


def main():
    output = HERE / "generated"
    output.mkdir(exist_ok=True)
    sources = [
        *ROOT.glob("tree_break_selection/**/*.py"),
        *ROOT.glob("benchmarks/diagnostics/calibration/selected/hierarchy/*.py"),
        *ROOT.glob("benchmarks/shared/generators/*.py"),
        ROOT / "benchmarks/shared/cases/gaussian.py",
        ROOT / "benchmarks/shared/cases/__init__.py",
        ROOT / "benchmarks/shared/util/case_inputs.py",
        Path(__file__),
        HERE / "protocol.md",
    ]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    save_json(
        HERE / "manifest.json",
        dict(
            source_sha256=hashes,
            python=sys.version,
            versions={p: version(p) for p in ["numpy", "scipy", "pandas", "pydiffmap", "pytest"]},
        ),
    )
    rows = []
    parity = 0
    memberships = 0

    def run(identifier, phase, scenario, data, sigma, chain_seed, stream, previous=None):
        nonlocal parity, memberships
        started = time.perf_counter()
        result = selected_diffusion_root_chain_test(data, sigma=sigma, seed=chain_seed)
        labels = data.index.tolist() if isinstance(data, pd.DataFrame) else None
        array = np.asarray(data, dtype=float)
        np.save(output / f"{identifier}.npy", array)
        if previous is not None:
            np.testing.assert_array_equal(result["chain_radii"], previous["chain_root_split"])
            expected = float(np.mean(previous["chain_root_split"] >= result["observed_radius"]))
            assert result["p_value"] == expected
            parity += 1
        if result["status"] == "complete":
            r = np.asarray(result["chain_radii"])
            assert len(r) == 1000 and r[result["anchor_index"]] == result["observed_radius"]
            assert result["p_value"] == np.mean(r >= result["observed_radius"])
            assert result["completed_steps"] == 999
            assert len(np.unique(r)) == result["unique_states"]
            path = make_diffusion_root_path(data, sigma=sigma)
            np.testing.assert_allclose(
                path.offset + path.observed_radius * path.slope, path.values, atol=1e-10
            )
            np.testing.assert_allclose(path.e @ path.offset, 0, atol=1e-10)
            for radius in r[[0, 500, 999]]:
                signature, _ = _geometry(path.offset + radius * path.slope, path.labels)
                groups = {i: (path.labels[i],) for i in range(len(array))}
                for step, (a, b) in enumerate(signature):
                    groups[len(array) + step] = groups[a] + groups[b]
                actual = {frozenset(groups[i]) for i in signature[-1]}
                assert actual == {frozenset(g) for g in result["root_children"]}
                memberships += 1
        record = dict(
            identifier=identifier,
            phase=phase,
            scenario=scenario,
            n=len(array),
            p=array.shape[1],
            input_sha256=hashlib.sha256(array.tobytes()).hexdigest(),
            labels=labels,
            result=result,
            seconds=time.perf_counter() - started,
        )
        stream.write(json.dumps(record, default=diagnostic_json_default, allow_nan=False) + "\n")
        stream.flush()
        rows.append(
            {k: v for k, v in record.items() if k not in {"labels", "result"}}
            | {
                k: result.get(k)
                for k in [
                    "status",
                    "decision",
                    "p_value",
                    "sampling_status",
                    "accepted_moves",
                    "unique_states",
                    "geometry_evaluations",
                    "observed_radius",
                    "sigma",
                    "seed",
                ]
            }
        )
        print(identifier, result["status"], result["decision"], result["p_value"], flush=True)

    with threadpool_limits(limits=1), (output / "records.jsonl").open("x") as stream:
        for scenario, count in [("null", 12), ("moderate", 6), ("strong", 6)]:
            for i in range(count):
                previous = np.load(PREVIOUS / f"generated/{scenario}_{i}.npz")
                run(
                    f"parity_{scenario}_{i}",
                    "parity",
                    scenario,
                    previous["input"],
                    1.0,
                    2026140401 + i,
                    stream,
                    previous,
                )
        for scenario, count, shift in [
            ("null", 40, 0.0),
            ("moderate", 10, 1.0),
            ("strong", 10, 2.0),
        ]:
            for i in range(count):
                x = np.random.default_rng(2026110701 + i).normal(size=(9, 3))
                x[:4, 0] -= shift
                x[4:, 0] += shift
                run(f"fresh_{scenario}_{i}", "fresh", scenario, x, 1.0, 2026111701 + i, stream)
        cases = {case["name"]: case for case in get_default_test_cases()}
        for i, name in enumerate(["gauss_clear_medium_continuous", "gauss_moderate_3c_continuous"]):
            case = cases[name]
            prepared = prepare_case_inputs(case, ["tbs"])
            np.testing.assert_array_equal(prepared.data.values, prepared.original_features)
            run(
                name, "registered", name, prepared.data, case["cluster_std"], 2026112701 + i, stream
            )
    table = pd.DataFrame(rows)
    table.to_csv(HERE / "results.csv", index=False)
    summary = table.groupby(["phase", "scenario"], sort=False).agg(
        cases=("status", "size"),
        complete=("status", lambda x: int((x == "complete").sum())),
        unresolved=("decision", lambda x: int((x == "unresolved").sum())),
        reject_005=("p_value", lambda x: int((x <= 0.05).sum())),
        reject_001=("p_value", lambda x: int((x <= 0.01).sum())),
        median_p=("p_value", "median"),
        minimum_moves=("accepted_moves", "min"),
    )
    summary.to_csv(HERE / "summary.csv")
    for name, digest in hashes.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    save_json(
        HERE / "verification.json",
        dict(
            parity_chains=parity,
            native_membership_checks=memberships,
            records=len(rows),
            sources_verified=len(hashes),
            failures=table[table.status != "complete"].identifier.tolist(),
        ),
    )
    print(summary.to_string(), flush=True)


if __name__ == "__main__":
    main()
