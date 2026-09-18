"""Resume this benchmark study through the existing per-method execution API."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--case-indices", type=int, nargs="+")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    manifest = json.loads((run / "run_manifest.json").read_text())
    root = Path(manifest["cwd"])
    os.environ.update(manifest["environment"])
    for name, expected in manifest["source_sha256"].items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Source changed since the run started: {name}")

    import numpy as np
    import pandas as pd
    from benchmarks.shared.audit_utils import export_decomposition_audit
    from benchmarks.shared.benchmark_grid import benchmark_run_id
    from benchmarks.shared.cases import get_test_cases_by_suite
    from benchmarks.shared.performance_grid import write_benchmark_performance_grid
    from benchmarks.shared.result_records import benchmark_rows_to_dataframe
    from benchmarks.shared.runners.method_registry import METHOD_SPECS
    from benchmarks.shared.util.case_inputs import prepare_case_inputs
    from benchmarks.shared.util.method_execution import run_single_method_once
    from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
        DEFAULT_EDGE_ALPHA,
        DEFAULT_SIBLING_ALPHA,
    )

    csv_path = run / "full_benchmark_comparison.csv"
    existing = pd.read_csv(csv_path) if csv_path.exists() else pd.DataFrame()
    done = set(zip(existing.get("case_id", []), existing.get("run_id", [])))
    errors_path = run / "execution_errors.jsonl"
    errors = (
        [json.loads(line) for line in errors_path.read_text().splitlines()]
        if errors_path.exists()
        else []
    )
    done.update((row["case_id"], row["run_id"]) for row in errors)
    methods = list(manifest["method_params"])
    cells = [
        (method, params, benchmark_run_id(method, params))
        for method in methods
        for params in manifest["method_params"][method]
    ]
    cases = get_test_cases_by_suite("full")
    indexed_cases = list(enumerate(cases, 1))
    if args.case_indices is not None:
        indexed_cases = [(idx, case) for idx, case in indexed_cases if idx in args.case_indices]
    assert len(indexed_cases) * len(cells) == manifest["expected_rows"]
    started = time.monotonic()
    print(f"Resuming {len(done)}/{manifest['expected_rows']} recorded outcomes", flush=True)

    for case_idx, case in indexed_cases:
        case_name = case["name"]
        pending = [
            (method, params, run_id)
            for method, params, run_id in cells
            if (case_name, run_id) not in done
        ]
        if not pending:
            continue
        inputs = prepare_case_inputs(case, methods)
        input_signature = {
            "case_id": case_name,
            "seed": case["seed"],
            "shape": list(inputs.data.shape),
            "data_sha256": hashlib.sha256(
                np.ascontiguousarray(inputs.data.values).tobytes()
            ).hexdigest(),
            "truth_sha256": hashlib.sha256(
                np.ascontiguousarray(inputs.labels).tobytes()
            ).hexdigest(),
        }
        with (run / "input_signatures.jsonl").open("a") as stream:
            stream.write(json.dumps(input_signature) + "\n")
        for method, params, run_id in pending:
            cell_start = time.monotonic()
            label = (
                method
                + "/"
                + str(params.get("tree_builder", "linkage"))
                + "/"
                + str(params["tree_linkage_method"])
            )
            print(f"[{case_idx}/{len(cases)}] {case_name}: {label}", flush=True)
            try:
                row, computed, _ = run_single_method_once(
                    method_id=method,
                    spec=METHOD_SPECS[method],
                    params=dict(params),
                    case_idx=case_idx,
                    case_name=case_name,
                    tc_seed=case["seed"],
                    significance_level=DEFAULT_SIBLING_ALPHA,
                    edge_alpha=DEFAULT_EDGE_ALPHA,
                    data_t=inputs.data,
                    y_t=inputs.labels,
                    x_original=inputs.original_features,
                    meta=inputs.metadata,
                    distance_matrix=inputs.distance_matrix,
                    distance_condensed=inputs.distance_condensed,
                    matrix_audit=False,
                )
            except Exception as exc:
                error = {
                    "test_case": case_idx,
                    "case_id": case_name,
                    "method": method,
                    "run_id": run_id,
                    "params": params,
                    "status": "execution_error",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                    "elapsed_seconds": time.monotonic() - cell_start,
                }
                errors.append(error)
                with errors_path.open("a") as stream:
                    stream.write(json.dumps(error) + "\n")
                outcome = "execution_error: " + str(exc)
            else:
                frame = benchmark_rows_to_dataframe([row])
                assert frame.iloc[0]["case_id"] == case_name
                assert frame.iloc[0]["run_id"] == run_id
                frame.to_csv(csv_path, mode="a", header=not csv_path.exists(), index=False)
                outcome = str(frame.iloc[0]["status"]) + "; ARI=" + str(frame.iloc[0]["ari"])
                if computed is not None:
                    cell_dir = (
                        run
                        / "cell_audits"
                        / f"{case_idx:03d}_{hashlib.sha256(run_id.encode()).hexdigest()[:12]}"
                    )
                    cell_dir.mkdir(parents=True, exist_ok=True)
                    (cell_dir / "identity.json").write_text(
                        json.dumps({"case_id": case_name, "run_id": run_id}) + "\n"
                    )
                    pd.DataFrame(
                        {
                            "leaf": inputs.data.index,
                            "truth": inputs.labels,
                            "cluster": computed.labels,
                        }
                    ).to_csv(cell_dir / "labels.csv", index=False)
                    export_decomposition_audit([computed], cell_dir)
                del row, computed, frame
            elapsed = time.monotonic() - cell_start
            with (run / "method_execution_times.jsonl").open("a") as stream:
                stream.write(
                    json.dumps(
                        {
                            "case_id": case_name,
                            "run_id": run_id,
                            "elapsed_seconds": elapsed,
                            "outcome": outcome,
                        }
                    )
                    + "\n"
                )
            done.add((case_name, run_id))
            print(
                f"  {outcome}; {elapsed:.2f}s; total {len(done)}/{manifest['expected_rows']}",
                flush=True,
            )
            gc.collect()
        del inputs

    results = pd.read_csv(csv_path)
    assert not results.duplicated(["case_id", "run_id"]).any()
    assert len(results) + len(errors) == manifest["expected_rows"]
    write_benchmark_performance_grid(results, run / "returned_results", source_path=csv_path)
    summaries = []
    for method, params, run_id in cells:
        rows = results[results.run_id.eq(run_id)]
        ok = rows[rows.status.eq("ok")]
        error_count = sum(error["run_id"] == run_id for error in errors)
        assert len(rows) + error_count == len(indexed_cases)
        summaries.append(
            {
                "method": method,
                "tree_builder": params.get("tree_builder", "linkage"),
                "tree_linkage_method": params["tree_linkage_method"],
                "run_id": run_id,
                "requested_cases": len(indexed_cases),
                "ok": len(ok),
                "unsupported": int(rows.status.eq("unsupported").sum()),
                "skip": int(rows.status.eq("skip").sum()),
                "execution_error": error_count,
                "mean_ari_on_ok": ok.ari.mean(),
                "median_ari_on_ok": ok.ari.median(),
                "mean_nmi_on_ok": ok.nmi.mean(),
                "exact_k_on_ok": int(ok.found_clusters.eq(ok.true_clusters).sum()),
                "over_split_on_ok": int(ok.over_split.sum()),
                "under_split_on_ok": int(ok.under_split.sum()),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(run / "configuration_summary.csv", index=False)
    status = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "exit_code": 0,
        "continuation_elapsed_seconds": time.monotonic() - started,
        "requested_outcomes": manifest["expected_rows"],
        "returned_rows": len(results),
        "execution_errors": len(errors),
        "all_configurations_attempted": True,
    }
    (run / "continuation_status.json").write_text(json.dumps(status, indent=2) + "\n")
    print(summary.drop(columns="run_id").to_string(index=False), flush=True)
    print(json.dumps(status), flush=True)


if __name__ == "__main__":
    main()
