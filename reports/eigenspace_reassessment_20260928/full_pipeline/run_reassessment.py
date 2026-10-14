"""Frozen independent-dataset reassessment of the current production methods.

This report-local driver uses the repository dispatcher, feature-space builders,
traversal and reporting helpers. It does not modify production behavior.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
import platform
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.shared.generators.case_data_contracts import one_hot_encode_categorical
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from benchmarks.validation.contracts.report_contract import wilson_interval
from benchmarks.validation.statistics.calibration_support import continuous_covariance
from scipy.stats import binomtest
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.record_collection import (
    collect_sibling_pair_records,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.parent_principal_component_inputs import (
    collect_parent_principal_component_inputs_for_sibling_tests,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.projection_dimensions import (
    derive_sibling_projection_dimensions_from_child_edge_comparisons,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.feature_space import (
    bernoulli_feature_space_from_columns,
    continuous_feature_space_from_columns,
)

METHODS = (
    "tbs", "tbs_complete", "tbs_neighbor_joining", "tbs_diffusion",
    "tbs_diffusion_adaptive_nnls",
)
FAMILIES = (
    "bernoulli_iid", "bernoulli_dependent", "categorical",
    "gaussian_identity", "gaussian_ar1", "gaussian_highdim",
)
TAILS = (0.001, 0.005, 0.01, 0.05)
EDGE_ALPHA = 0.001
PRIMARY_ALPHA = 0.01


def dataset(family: str, signal: bool, seed: int):
    """Known population laws; truth never depends on fitted tree or gates."""
    rng = np.random.default_rng(seed)
    n = 40 if family == "gaussian_highdim" else 64
    minority = n // 4 if family in {"bernoulli_dependent", "gaussian_ar1"} else n // 2
    truth = np.zeros(n, dtype=int)
    if signal:
        truth[-minority:] = 1
    names = [f"S{i}" for i in range(n)]
    if family.startswith("bernoulli"):
        d = 8 if family == "bernoulli_dependent" else 24
        probabilities = np.full((n, d), 0.5)
        if signal:
            probabilities[:, : d // 2] = np.where(truth[:, None] == 0, 0.2, 0.8)
        values = (rng.uniform(size=(n, d)) < probabilities).astype(float)
        if family == "bernoulli_dependent":
            values = np.repeat(values, 3, axis=1)
        data = pd.DataFrame(values, index=names, columns=[f"F{i}" for i in range(24)])
        space = bernoulli_feature_space_from_columns(data.columns)
    elif family == "categorical":
        profiles = np.full((n, 3), 1 / 3)
        if signal:
            profiles[truth == 0] = [0.65, 0.20, 0.15]
            profiles[truth == 1] = [0.15, 0.20, 0.65]
        uniform = rng.uniform(size=(n, 8))
        values = (uniform[:, :, None] > np.cumsum(profiles, axis=1)[:, None, :]).sum(2)
        data, _, space = one_hot_encode_categorical(values, 3, names)
    else:
        d = 64 if family == "gaussian_highdim" else 8
        covariance = continuous_covariance(
            dimension=d, profile="ar1_0.6" if family == "gaussian_ar1" else "identity"
        )
        values = rng.multivariate_normal(np.zeros(d), covariance, size=n)
        if signal:
            values[:, : (8 if d == 64 else 2)] += np.where(truth[:, None] == 0, -1, 1)
        data = pd.DataFrame(values, index=names, columns=[f"X{i}" for i in range(d)])
        space = continuous_feature_space_from_columns(data.columns)
    permutation = rng.permutation(n)
    data = data.iloc[permutation].copy()
    data.index = names
    return data, truth[permutation], space


def exact_null(left: list[int], right: list[int], truth: np.ndarray) -> bool:
    """Equality of two known population-mixture parameters by integer counts."""
    return int(truth[left].sum()) * len(right) == int(truth[right].sum()) * len(left)


def number(value) -> float:
    return float(value) if value is not None and pd.notna(value) else float("nan")


def method_run(data, truth, space, key, method, alpha):
    start = perf_counter()
    params = dict(METHOD_SPECS[method].param_grid[0])
    if key["family"].startswith("gaussian") and method in METHODS[:3]:
        params["tree_distance_metric"] = "standardized_euclidean"
    row = {
        **key, "method": method, "sibling_alpha": alpha,
        "status": "runtime_error", "reason": "", "audit_error": "",
        "found_clusters": np.nan, "ari": np.nan, "false_partition": np.nan,
        "n_splits": np.nan, "n_false_splits": np.nan, "split_fdp": np.nan,
        "root_false_split": np.nan, "descendant_false_split": np.nan,
        "pass_throughs": np.nan, "root_p": np.nan, "root_true_null": np.nan,
        "root_test_rank": np.nan, "max_test_rank": np.nan,
        "n_rank_above_mp_floor": np.nan,
        "calibration_records": np.nan, "nonnull_calibration_records": np.nan,
        "calibration_base_weight": np.nan, "nonnull_calibration_base_weight": np.nan,
        "actual_config": "", "requested_params": json.dumps(params, sort_keys=True),
        "diffusion_metadata": "", "warnings": "", "elapsed_seconds": np.nan,
    }
    nodes = []
    try:
        with warnings.catch_warnings(record=True) as warning_records:
            warnings.simplefilter("always")
            result = run_clustering_result(
                data, method, params, seed=key["seed"], significance_level=alpha,
                edge_alpha=EDGE_ALPHA, feature_space=space,
            )
        row["warnings"] = json.dumps(sorted({str(w.message) for w in warning_records}))
        row["status"] = str(result.status)
        row["reason"] = str(result.unsupported_reason or result.skip_reason or "")
        if result.status == "ok":
            row.update(
                found_clusters=result.found_clusters,
                ari=adjusted_rand_score(truth, result.labels),
                false_partition=int(result.found_clusters > 1) if not key["signal"] else np.nan,
            )
        bundle = result.extra.get("gate_bundle")
        if bundle is None:
            return row, nodes
        row["actual_config"] = json.dumps(asdict(bundle.metadata), sort_keys=True)
        if bundle.metadata.edge.alpha != EDGE_ALPHA or bundle.metadata.sibling.alpha != alpha:
            raise AssertionError("Requested and executed alpha differ")
        if bundle.metadata.config.adaptive_projection_dimension_energy_fraction is not None:
            raise AssertionError("This fixed production-default study must not select adaptive rank")
        for metadata_name in ("adaptive_diffusion", "hamming_nn_diffusion"):
            if metadata_name in result.extra:
                row["diffusion_metadata"] = json.dumps(result.extra[metadata_name], sort_keys=True)
        tree = result.extra["tree"]
        root = tree.root()
        descendants = tree.compute_descendant_sets(use_labels=True)
        positions = {name: i for i, name in enumerate(data.index)}
        leaf_indices = {
            node: [positions[label] for label in labels] for node, labels in descendants.items()
        }
        trace = {}
        if result.status == "ok":
            replay = TreeDecomposition(
                tree=tree, gate_annotation_bundle=bundle, passthrough=True, trace_level="full",
            ).decompose_tree()
            replay_labels, _ = labels_and_report_from_decomposition(replay, data.index.tolist())
            if adjusted_rand_score(result.labels, replay_labels) != 1:
                raise AssertionError("Full traversal replay disagrees with actual dispatcher labels")
            trace = {r["node_id"]: r for r in replay["traversal_trace"]}
        annotations = result.extra["annotations"]
        spectral = bundle.edge_gate_result.spectral_context
        dimensions = spectral.test_projection_dimensions_by_node
        row["root_test_rank"] = dimensions.get(root, np.nan)
        row["max_test_rank"] = max(dimensions.values(), default=0)
        row["n_rank_above_mp_floor"] = sum(
            k > max(2, spectral.raw_mp_signal_counts_by_node[node])
            for node, k in dimensions.items()
        )
        # Reconstruct raw statistics for audit using the same default covariance/rank.
        sibling_dimensions = derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree, spectral_context=spectral,
        )
        projections, eigenvalues = collect_parent_principal_component_inputs_for_sibling_tests(
            sibling_dimensions, spectral_context=spectral,
        )
        config = bundle.metadata.config
        raw_records, _ = collect_sibling_pair_records(
            tree, bundle.edge_gate_result.annotated_df,
            sibling_projection_dimensions_from_edge_comparisons=sibling_dimensions,
            parent_principal_component_projections=projections,
            parent_principal_component_eigenvalues=eigenvalues,
            feature_space=space,
            continuous_covariance_policy=config.continuous_covariance_policy,
            continuous_covariance_min_child_leaf_count=config.continuous_covariance_min_child_leaf_count,
            adaptive_projection_dimension_energy_fraction=config.adaptive_projection_dimension_energy_fraction,
        )
        raw_by_parent = {r.parent: r for r in raw_records}
        for parent in tree.nodes:
            children = list(tree.successors(parent))
            if len(children) != 2:
                continue
            left, right = children
            li, ri = leaf_indices[left], leaf_indices[right]
            is_null = exact_null(li, ri, truth)
            annotation = annotations.loc[parent]
            df = number(annotation.get("Sibling_Degrees_of_Freedom"))
            weight = number(annotation.get("Sibling_Null_Weight"))
            role = annotation.get("Sibling_Role_Supported", False)
            calibration = bool(pd.notna(role) and role and df > 0 and weight > 0)
            raw = raw_by_parent.get(parent)
            action = trace.get(parent, {}).get("decision", "not_reached")
            node_row = {
                **key, "method": method, "sibling_alpha": alpha, "status": row["status"],
                "node": str(parent), "is_root": parent == root, "decision": action,
                "left_n": len(li), "right_n": len(ri),
                "left_population1": int(truth[li].sum()), "right_population1": int(truth[ri].sum()),
                "true_null": is_null, "pure_parent": len(np.unique(truth[li + ri])) == 1,
                "raw_statistic": raw.stat if raw else np.nan,
                "raw_p": raw.p_value if raw else np.nan,
                "reference_scale": raw.reference_scale if raw else np.nan,
                "adjusted_statistic": number(annotation.get("Sibling_Test_Statistic")),
                "adjusted_p": number(annotation.get("Sibling_Divergence_P_Value")),
                "bh_p": number(annotation.get("Sibling_Divergence_P_Value_Corrected")),
                "df": df, "projection_dimension": number(annotation.get("Sibling_Projection_Dimension")),
                "calibration_record": calibration, "calibration_base_weight": weight,
                "calibration_status": str(annotation.get("Sibling_Gate_P_Value_Calibration", "")),
                "calibration_role": str(annotation.get("Sibling_Gate_P_Value_Role", "")),
            }
            if raw is not None and np.isfinite(df) and df != raw.degrees_of_freedom:
                raise AssertionError("Raw-record audit disagrees with executed degrees of freedom")
            nodes.append(node_row)
        root_node = next(node for node in nodes if node["is_root"])
        row["root_p"] = root_node["adjusted_p"]
        row["root_true_null"] = root_node["true_null"]
        support = [node for node in nodes if node["calibration_record"]]
        row.update(
            calibration_records=len(support),
            nonnull_calibration_records=sum(not node["true_null"] for node in support),
            calibration_base_weight=sum(node["calibration_base_weight"] for node in support),
            nonnull_calibration_base_weight=sum(
                node["calibration_base_weight"] for node in support if not node["true_null"]
            ),
        )
        if result.status == "ok":
            splits = [node for node in nodes if node["decision"] == "split"]
            false_splits = [node for node in splits if node["true_null"]]
            row.update(
                n_splits=len(splits), n_false_splits=len(false_splits),
                split_fdp=len(false_splits) / max(1, len(splits)),
                root_false_split=int(any(node["is_root"] for node in false_splits)),
                descendant_false_split=int(any(not node["is_root"] for node in false_splits)),
                pass_throughs=sum(node["decision"] == "pass_through" for node in nodes),
            )
    except Exception as exc:
        if row["status"] != "runtime_error":
            row["audit_error"] = repr(exc)
        else:
            row["reason"] = repr(exc)
    finally:
        row["elapsed_seconds"] = perf_counter() - start
    return row, nodes


def run_dataset(job):
    family, signal, replicate, phase = job
    seed = int(np.random.SeedSequence([
        20260928, 17 if phase == "evaluation" else 91,
        FAMILIES.index(family), int(signal), replicate,
    ]).generate_state(1, dtype=np.uint64)[0])
    data, truth, space = dataset(family, signal, seed)
    key = {
        "family": family, "signal": signal, "replicate": replicate, "seed": seed,
        "n_samples": len(data), "n_features": data.shape[1],
        "data_sha256": hashlib.sha256(np.ascontiguousarray(data, dtype=np.float64).tobytes()).hexdigest(),
    }
    rows, nodes = [], []
    for method in METHODS:
        alphas = (0.005, 0.01, 0.05) if (
            family == "bernoulli_iid" and method in {"tbs", "tbs_diffusion_adaptive_nnls"}
        ) else (PRIMARY_ALPHA,)
        for alpha in alphas:
            row, node_rows = method_run(data, truth, space, key, method, alpha)
            rows.append(row)
            nodes.extend(node_rows)
    return rows, nodes


def summarize(runs: pd.DataFrame):
    summaries, tails = [], []
    for group, frame in runs.groupby(["family", "signal", "method", "sibling_alpha"]):
        family, signal, method, alpha = group
        base = dict(family=family, signal=signal, method=method, sibling_alpha=alpha)
        pipeline_ok = frame[frame.status == "ok"]
        ok = pipeline_ok[pipeline_ok.audit_error.fillna("").eq("")]
        eligible = ok[ok.root_p.notna() & ok.root_true_null.eq(True)]
        event_count = int((ok.n_false_splits > 0).sum())
        ci = wilson_interval(event_count, len(ok)) if len(ok) else (np.nan, np.nan)
        summary = {
            **base, "attempts": len(frame), "completed": len(pipeline_ok), "audited": len(ok),
            "unsupported": int((frame.status == "unsupported").sum()),
            "skipped": int((frame.status == "skip").sum()),
            "runtime_errors": int((frame.status == "runtime_error").sum()),
            "audit_errors": int(frame.audit_error.fillna("").ne("").sum()),
            "any_false_split_count": event_count,
            "any_false_split_rate": event_count / len(ok) if len(ok) else np.nan,
            "any_false_split_wilson_low": ci[0], "any_false_split_wilson_high": ci[1],
            "mean_split_fdp": ok.split_fdp.mean(),
            "mean_ari": pipeline_ok.ari.mean() if signal else np.nan,
            "mean_clusters": pipeline_ok.found_clusters.mean(),
            "mean_pass_throughs": ok.pass_throughs.mean(),
            "mean_root_test_rank": ok.root_test_rank.mean(),
            "mean_max_test_rank": ok.max_test_rank.mean(),
            "datasets_with_nonnull_calibration_records": int((ok.nonnull_calibration_records > 0).sum()),
            "calibration_records": int(ok.calibration_records.sum()),
            "nonnull_calibration_records": int(ok.nonnull_calibration_records.sum()),
            "sum_calibration_base_weight": ok.calibration_base_weight.sum(),
            "sum_nonnull_calibration_base_weight": ok.nonnull_calibration_base_weight.sum(),
        }
        if not signal and len(pipeline_ok):
            count = int(pipeline_ok.false_partition.sum())
            low, high = wilson_interval(count, len(pipeline_ok))
            summary.update(
                false_partition_count=count, false_partition_rate=count / len(pipeline_ok),
                false_partition_wilson_low=low, false_partition_wilson_high=high,
                observed_false_partition_per_attempt=count / len(frame),
                unavailable_outcome_bound_upper=(count + len(frame) - len(pipeline_ok)) / len(frame),
                binomial_inflation_p=binomtest(count, len(pipeline_ok), alpha, alternative="greater").pvalue,
            )
        summaries.append(summary)
        for cutoff in TAILS:
            count = int((eligible.root_p <= cutoff).sum())
            low, high = wilson_interval(count, len(eligible)) if len(eligible) else (np.nan, np.nan)
            tails.append({
                **base, "cutoff": cutoff, "attempts": len(frame),
                "completed": len(pipeline_ok), "audited": len(ok),
                "eligible_root_nulls": len(eligible), "rejections": count,
                "rate": count / len(eligible) if len(eligible) else np.nan,
                "wilson_low": low, "wilson_high": high,
                "binomial_inflation_p": binomtest(count, len(eligible), cutoff, alternative="greater").pvalue
                if len(eligible) else np.nan,
            })
    summary_df, tail_df = pd.DataFrame(summaries), pd.DataFrame(tails)
    # Separate prespecified families: global-null partition endpoints and root tails.
    for table in (summary_df, tail_df):
        family = ~table.signal if table is summary_df else pd.Series(True, index=table.index)
        valid = family & table.binomial_inflation_p.notna()
        adjusted = pd.Series(np.nan, index=table.index)
        adjusted.loc[family] = multipletests(
            table.loc[family, "binomial_inflation_p"].fillna(1.0), method="holm",
        )[1]
        table["holm_p"] = np.nan
        table.loc[valid, "holm_p"] = adjusted.loc[valid]
        table["holm_family_size"] = int(family.sum())
    return summary_df, tail_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["pilot", "evaluation"], required=True)
    parser.add_argument("--null-replicates", type=int, default=500)
    parser.add_argument("--signal-replicates", type=int, default=200)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    jobs = [
        (family, signal, replicate, args.phase)
        for family in FAMILIES for signal in (False, True)
        for replicate in range(args.signal_replicates if signal else args.null_replicates)
    ]
    sources = [Path(__file__), ROOT / "pyproject.toml", ROOT / "uv.lock"]
    for folder in ("tree_break_selection", "benchmarks/shared", "benchmarks/validation/statistics"):
        sources += sorted((ROOT / folder).rglob("*.py"))
    sources += [
        ROOT / "benchmarks/diagnostics/calibration/reporting.py",
        ROOT / "benchmarks/validation/contracts/report_contract.py",
        Path(__file__).with_name("protocol.md"),
    ]
    manifest = {
        "schema": "eigenspace_pipeline_reassessment/v1", "arguments": vars(args),
        "started_utc": datetime.now(UTC).isoformat(), "command": sys.argv,
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], text=True),
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "python": sys.version, "platform": platform.platform(),
        "packages": {p: importlib.metadata.version(p) for p in (
            "numpy", "scipy", "pandas", "scikit-learn", "pydiffmap", "statsmodels", "scikit-bio",
        )},
        "thread_environment": {k: os.environ.get(k) for k in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "TBS_N_JOBS",
        )},
        "protocol": "../protocol.md", "planned_datasets": len(jobs),
        "method_parameters": {m: METHOD_SPECS[m].param_grid[0] for m in METHODS},
        "primary_alpha": PRIMARY_ALPHA, "edge_alpha": EDGE_ALPHA,
        "complete": False,
    }
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    start = perf_counter()
    context = multiprocessing.get_context("spawn")
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers, mp_context=context) as pool,
        (args.output_dir / "runs.csv").open("w") as run_file,
        gzip.open(args.output_dir / "nodes.csv.gz", "wt") as node_file,
    ):
        run_writer = node_writer = None
        for completed, (result_rows, nodes) in enumerate(pool.map(run_dataset, jobs), 1):
            if run_writer is None:
                run_writer = csv.DictWriter(run_file, fieldnames=list(result_rows[0]))
                run_writer.writeheader()
            run_writer.writerows(result_rows)
            rows.extend(result_rows)
            if nodes:
                if node_writer is None:
                    node_writer = csv.DictWriter(node_file, fieldnames=list(nodes[0]))
                    node_writer.writeheader()
                node_writer.writerows(nodes)
            if completed % 50 == 0 or completed == len(jobs):
                run_file.flush()
                node_file.flush()
                print(f"{completed}/{len(jobs)} datasets; {len(rows)} method/alpha runs; {perf_counter()-start:.1f}s", flush=True)
    runs = pd.DataFrame(rows)
    summary, tails = summarize(runs)
    final_sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    if final_sources != manifest["source_sha256"]:
        raise RuntimeError("Source changed during evaluation; do not pool different implementations")
    write_diagnostic_bundle(
        output_dir=args.output_dir, tables={"summary": summary, "root_tails": tails},
        filenames={"summary": "summary.csv", "root_tails": "root_tails.csv"},
        manifest={**manifest, "complete": True, "completed_datasets": len(jobs),
                  "method_alpha_runs": len(rows), "elapsed_seconds": perf_counter() - start,
                  "raw_outputs": {"runs": "runs.csv", "nodes": "nodes.csv.gz"}},
    )
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
