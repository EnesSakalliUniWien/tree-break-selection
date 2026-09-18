"""Matched diagnostic ablation of the unsupported rule and existing sibling gates.

This is a study script: process-local patches never change production source.
Fixed-coordinate outcomes are clustering diagnostics, not calibrated inference.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    manifest = json.loads((run / "run_manifest.json").read_text())
    os.environ.update(manifest["environment"])
    root = Path(manifest["cwd"])
    for name, expected in manifest["source_sha256"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name

    import numpy as np
    import pandas as pd
    from benchmarks.shared.benchmark_grid import strip_benchmark_metadata
    from benchmarks.shared.cases import get_test_cases_by_suite
    from benchmarks.shared.runners import tbs_diffusion_runner, tbs_runner
    from benchmarks.shared.runners.dispatch import run_clustering_result
    from benchmarks.shared.util.case_inputs import prepare_case_inputs
    from sklearn.metrics import adjusted_rand_score

    out = run / "support_review"
    out.mkdir(exist_ok=True)
    case_names = (
        "gauss_clear_small",
        "gauss_dense_signal_highd",
        "overlap_extreme_4c",
        "gauss_null_small",
        "gauss_null_large",
        "binary_null_small",
        "binary_null_medium",
    )
    methods = manifest["method_params"]
    cells = {
        method: next(
            params
            for params in configurations
            if params.get("tree_builder", "linkage") == "linkage"
            and params["tree_linkage_method"] == "average"
        )
        for method, configurations in methods.items()
    }
    variants = [("status_rule_bypassed", 0.01)] + [
        (gate, alpha)
        for gate in ("fixed_coordinate_bh", "fixed_coordinate_by")
        for alpha in (0.005, 0.01, 0.05)
    ]
    study_manifest = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": case_names,
        "method_params": cells,
        "variants": [("production", 0.01), *variants],
        "expected_outcomes": len(case_names) * len(cells) * (1 + len(variants)),
        "edge_alpha": 0.001,
        "design": "Reuse captured diffusion distances; refit the same NNLS tree for every gate. Verify identical topology and branch lengths by SHA256. Truth used only for ARI.",
        "inference_limit": "Fixed-coordinate gates are diagnostics on a selected hierarchy, without a validated selected-tail law. Four null seeds do not establish error control.",
        "skip_policy": "A production input-contract skip applies to every downstream gate variant. Propagated cells have execution_kind=input_contract_skip and zero execution time; no gate was run for them.",
    }
    rows = []
    if (out / "results.csv").exists():
        previous = pd.read_csv(out / "results.csv")
        assert not previous.duplicated(["case_id", "method", "variant", "sibling_alpha"]).any()
        complete = previous.groupby(["case_id", "method"]).status.transform("size").eq(8)
        previous = previous[complete].copy()
        if "execution_kind" not in previous:
            previous["execution_kind"] = "method_run"
        rows = previous.to_dict("records")
    study_manifest["reused_completed_outcomes"] = len(rows)
    (out / "manifest.json").write_text(json.dumps(study_manifest, indent=2) + "\n")
    original_run = tbs_runner.run_tbs_on_distance
    signature = inspect.signature(original_run)
    cases = {case["name"]: case for case in get_test_cases_by_suite("full")}
    started = time.monotonic()

    for case_name in case_names:
        case = cases[case_name]
        inputs = prepare_case_inputs(case, list(cells))
        data_hash = hashlib.sha256(np.ascontiguousarray(inputs.data.values).tobytes()).hexdigest()
        for method, params in cells.items():
            if any(row["case_id"] == case_name and row["method"] == method for row in rows):
                continue
            captured = {}

            def capture(*positional, **keywords):
                captured.update(signature.bind(*positional, **keywords).arguments)
                return original_run(*positional, **keywords)

            def record(result, variant, alpha, elapsed, execution_kind="method_run"):
                extra = result.extra or {}
                tree = extra.get("tree")
                edge_records = (
                    sorted(
                        (repr(u), repr(v), attrs["branch_length"])
                        for u, v, attrs in tree.edges(data=True)
                    )
                    if tree is not None
                    else None
                )
                tree_hash = hashlib.sha256(json.dumps(edge_records).encode()).hexdigest()
                annotation = extra.get("annotations", pd.DataFrame())
                reason = result.unsupported_reason
                row = {
                    "case_id": case_name,
                    "seed": case["seed"],
                    "method": method,
                    "variant": variant,
                    "sibling_alpha": alpha,
                    "edge_alpha": 0.001,
                    "status": str(result.status),
                    "execution_kind": execution_kind,
                    "skip_reason": result.skip_reason,
                    "true_clusters": len(np.unique(inputs.labels)),
                    "found_clusters": result.found_clusters if result.labels is not None else None,
                    "ari": adjusted_rand_score(inputs.labels, result.labels)
                    if result.labels is not None
                    else None,
                    "unsupported_reason": reason.code.value if reason else None,
                    "admissible_support_count": int(
                        annotation["Sibling_Role_Supported"].fillna(False).astype(bool).sum()
                    )
                    if "Sibling_Role_Supported" in annotation
                    else None,
                    "data_sha256": data_hash,
                    "tree_branch_lengths_sha256": tree_hash,
                    "elapsed_seconds": elapsed,
                }
                if variant != "production":
                    baseline = next(
                        r
                        for r in rows
                        if r["case_id"] == case_name
                        and r["method"] == method
                        and r["variant"] == "production"
                    )
                    assert tree_hash == baseline["tree_branch_lengths_sha256"]
                rows.append(row)
                cell_dir = out / f"{case_name}__{method}__{variant}__{alpha}"
                cell_dir.mkdir(exist_ok=True)
                annotation.to_csv(cell_dir / "annotations.csv")
                if result.labels is not None:
                    pd.DataFrame(
                        {
                            "leaf": inputs.data.index,
                            "truth": inputs.labels,
                            "cluster": result.labels,
                        }
                    ).to_csv(cell_dir / "labels.csv", index=False)
                if reason:
                    (cell_dir / "unsupported.json").write_text(
                        json.dumps(asdict(reason), indent=2) + "\n"
                    )
                pd.DataFrame(rows).to_csv(out / "results.csv", index=False)
                print(
                    f"{len(rows)}/{study_manifest['expected_outcomes']} {case_name} {method} {variant} alpha={alpha}: {row['status']}, K={row['found_clusters']}, ARI={row['ari']}",
                    flush=True,
                )

            cell_start = time.monotonic()
            with patch.object(tbs_diffusion_runner, "run_tbs_on_distance", capture):
                result = run_clustering_result(
                    data_df=inputs.data,
                    method_id=method,
                    params=strip_benchmark_metadata(params),
                    seed=case["seed"],
                    significance_level=0.01,
                    edge_alpha=0.001,
                    distance_matrix=inputs.distance_matrix,
                    distance_condensed=None,
                    feature_space=inputs.metadata.get("feature_space"),
                )
            record(result, "production", 0.01, time.monotonic() - cell_start)
            if not captured and str(result.status) == "skip":
                for variant, alpha in variants:
                    record(result, variant, alpha, 0.0, "input_contract_skip")
                continue
            assert captured, (case_name, method, result.status)
            for variant, alpha in variants:
                keywords = dict(captured)
                keywords["sibling_significance_level"] = alpha
                cell_start = time.monotonic()
                if variant == "status_rule_bypassed":
                    with patch.object(
                        tbs_runner, "unsupported_empirical_null_reason", return_value=None
                    ):
                        result = original_run(**keywords)
                else:
                    keywords["sibling_gate_method"] = variant
                    result = original_run(**keywords)
                record(result, variant, alpha, time.monotonic() - cell_start)

    assert len(rows) == study_manifest["expected_outcomes"]
    (out / "status.json").write_text(
        json.dumps(
            {
                "exit_code": 0,
                "outcomes": len(rows),
                "elapsed_seconds": time.monotonic() - started,
                "finished_at": datetime.now(timezone.utc).isoformat(),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
