"""Run the registered cancer PC cases and default methods on six fixed seeds."""

import hashlib
import importlib.metadata
import json
import sys
import time
from pathlib import Path

from benchmarks.shared.benchmark_runs.runtime import apply_single_thread_runtime_defaults

apply_single_thread_runtime_defaults()

import numpy as np
import pandas as pd

from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.generators.generate_case_data import generate_case_data
from benchmarks.shared.result_records import benchmark_rows_to_dataframe
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.case_run import run_single_case
from benchmarks.shared.util.method_sets import DEFAULT_METHODS
from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_EDGE_ALPHA,
    DEFAULT_SIBLING_ALPHA,
)


def main():
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "generated/panel"
    output.mkdir(parents=True, exist_ok=False)
    cases = get_test_cases_by_suite("cancer_pc")
    seeds = list(range(20261104, 20261110))
    methods = list(DEFAULT_METHODS)
    params = {method: METHOD_SPECS[method].param_grid for method in methods}
    sources = sorted({
        p for root in (Path("benchmarks/shared"), Path("tree_break_selection"))
        for p in root.rglob("*.py")
    } | {Path(__file__)})
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    manifest = {
        "cases": cases, "seeds": seeds, "methods": methods, "params": params,
        "edge_alpha": DEFAULT_EDGE_ALPHA, "sibling_alpha": DEFAULT_SIBLING_ALPHA,
        "source_sha256": hashes,
        "packages": {p: importlib.metadata.version(p) for p in (
            "numpy", "scipy", "pandas", "scikit-learn", "pydiffmap", "matplotlib",
            "umap-learn", "leidenalg", "igraph", "hdbscan",
        )},
        "oracle_true_k_methods": ["kmeans", "spectral"],
        "seed_roles": {str(s): "development" if s == seeds[0] else "additional" for s in seeds},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    frames, fingerprints = [], []
    for seed in seeds:
        for recipe in cases:
            start = time.monotonic()
            index = len(frames) + 1
            case = dict(recipe, seed=seed, test_case_num=index)
            case_dir = output / str(seed) / case["name"]
            case_dir.mkdir(parents=True)
            data, truth, _, _ = generate_case_data(case)
            fingerprint = hashlib.sha256(data.to_numpy().tobytes()).hexdigest()
            rows, computed = run_single_case(
                tc=case, total_cases=len(cases) * len(seeds), selected_methods=methods,
                param_sets=params, significance_level=DEFAULT_SIBLING_ALPHA,
                edge_alpha=DEFAULT_EDGE_ALPHA, output_pdf=None,
                plots_root=case_dir / "plots", verbose=False,
            )
            frame = benchmark_rows_to_dataframe(rows)
            frame["simulation_seed"] = seed
            frame["seed_role"] = manifest["seed_roles"][str(seed)]
            assert len(frame) == sum(map(len, params.values()))
            assert set(frame.method) == set(methods)
            frame.to_csv(case_dir / "results.csv", index=False)
            frames.append(frame)
            arrays = {"values": data.to_numpy(), "truth": truth, "sample_ids": data.index.to_numpy(dtype=str)}
            for result in computed:
                np.testing.assert_array_equal(result.data.to_numpy(), arrays["values"])
                np.testing.assert_array_equal(result.y_true, truth)
                key = f"labels__{result.method}"
                assert key not in arrays, "This export expects one registered parameter setting per method."
                arrays[key] = result.labels
            np.savez_compressed(case_dir / "assignments.npz", **arrays)
            fingerprints.append({"case_id": case["name"], "simulation_seed": seed, "sha256": fingerprint})
            print(f"{index}/24 {case['name']} seed={seed}: {frame.status.value_counts().to_dict()} in {time.monotonic()-start:.2f}s", flush=True)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(output / "results.csv", index=False)
    pd.DataFrame(fingerprints).to_csv(output / "input-fingerprints.csv", index=False)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    assert len(combined) == 24 * sum(map(len, params.values()))
    print(combined.groupby(["method", "status"]).size().to_string(), flush=True)


if __name__ == "__main__":
    main()
