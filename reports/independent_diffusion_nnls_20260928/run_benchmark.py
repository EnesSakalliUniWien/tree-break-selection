"""Discovery-only diffusion/NNLS with matched independent Gaussian testing."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import multiprocessing
import platform
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from benchmarks.validation.contracts.report_contract import wilson_interval
from scipy.stats import binomtest
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests
from tree_break_selection.hierarchy_analysis.statistics.independent_gaussian import (
    IndependentGaussianDiffusion,
    _annotate_tests,
    fit_independent_gaussian_tree,
    test_independent_gaussian_tree,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ("identity2", "ar1_8", "identity40")
GEOMETRIES = ("direct", "diffusion", "diffusion_nnls")
ALPHAS = (0.005, 0.01, 0.05)


def digest(values):
    return hashlib.sha256(np.ascontiguousarray(values, dtype=float).tobytes()).hexdigest()


def json_hash(values):
    return hashlib.sha256(json.dumps(values, sort_keys=True, default=str).encode()).hexdigest()


def frozen_hashes(model):
    return dict(
        topology_sha256=json_hash(sorted((str(u), str(v)) for u, v in model.tree.edges)),
        branches_sha256=json_hash(
            [(str(u), str(v), float(a["branch_length"])) for u, v, a in model.tree.edges(data=True)]
        ),
        pca_sha256=json_hash(
            [(str(n.parent), digest(n.projection), digest(n.eigenvalues)) for n in model.nodes]
        ),
        diffusion_sha256=None
        if model.diffusion_coordinates is None
        else digest(model.diffusion_coordinates),
    )


def source_hashes():
    paths = sorted((ROOT / "tree_break_selection").rglob("*.py")) + sorted(
        (ROOT / "benchmarks/shared").rglob("*.py")
    )
    paths += [Path(__file__), Path(__file__).with_name("protocol.md"), ROOT / "uv.lock"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run_dataset(job):
    profile, signal, replicate, namespace = job
    dimension = {"identity2": 2, "ar1_8": 8, "identity40": 40}[profile]
    sigma = (
        0.6 ** np.abs(np.subtract.outer(np.arange(dimension), np.arange(dimension)))
        if profile == "ar1_8"
        else np.eye(dimension)
    )
    seeds = np.random.SeedSequence(
        [namespace, PROFILES.index(profile), int(signal), replicate]
    ).spawn(2)
    truth = np.zeros(32, dtype=int)
    means = np.zeros((32, dimension))
    if signal:
        truth[16:] = 1
        means[:, :2] = np.where(truth[:, None] == 0, -1.0, 1.0)
    labels = [f"S{i}" for i in range(32)]
    columns = [f"X{i}" for i in range(dimension)]
    matrices = [
        pd.DataFrame(
            means + np.random.default_rng(seed).multivariate_normal(np.zeros(dimension), sigma, 32),
            index=labels,
            columns=columns,
        )
        for seed in seeds
    ]
    discovery, inference = matrices
    key = dict(
        profile=profile,
        signal=signal,
        replicate=replicate,
        discovery_sha256=digest(discovery),
        inference_sha256=digest(inference),
    )
    rows, nodes, metadata = [], [], []
    for geometry_name in GEOMETRIES:
        start = perf_counter()
        config = (
            None
            if geometry_name == "direct"
            else IndependentGaussianDiffusion(fit_nnls=geometry_name == "diffusion_nnls")
        )
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model = fit_independent_gaussian_tree(
                    discovery,
                    known_covariance=sigma,
                    observation_ids=[f"d{i}" for i in range(32)],
                    diffusion=config,
                )
                before = frozen_hashes(model)
                result = test_independent_gaussian_tree(
                    model, inference, observation_ids=[f"t{i}" for i in range(32)]
                )
            after = frozen_hashes(model)
            assert before == after, "Inference modified fitted discovery context"
            root = model.tree.root()
            true_null = {
                n.parent: int(truth[list(n.left_rows)].sum()) * len(n.right_rows)
                == int(truth[list(n.right_rows)].sum()) * len(n.left_rows)
                for n in model.nodes
            }
            metadata.append(
                {
                    **key,
                    "geometry": geometry_name,
                    "status": "ok",
                    **before,
                    "frozen_context_unchanged": True,
                    "diffusion": json.loads(model.geometry_metadata_json),
                    "nnls": json.loads(model.branch_length_metadata_json),
                    "warnings": sorted({str(w.message) for w in caught}),
                }
            )
            for item in result.tests.to_dict("records"):
                nodes.append(
                    {
                        **key,
                        "geometry": geometry_name,
                        **item,
                        "is_root": item["parent"] == root,
                        "true_null": true_null[item["parent"]],
                    }
                )
            for reference in ("conditional", "nominal"):
                tests = result.tests.copy()
                if reference == "nominal":
                    tests["p_value"] = tests.nominal_p_value
                    tests["reference_law"] = "nominal_chi_square_after_energy_selection"
                for alpha in ALPHAS:
                    if alpha == 0.01 and reference == "conditional":
                        decomposition = result.decomposition
                    else:
                        _, decomposition = _annotate_tests(
                            model, tests, alpha=alpha, edge_alpha=0.001, passthrough=True
                        )
                    predicted, _ = labels_and_report_from_decomposition(decomposition, labels)
                    split_nodes = [
                        x["node_id"]
                        for x in decomposition["traversal_trace"]
                        if x["decision"] == "split"
                    ]
                    false = sum(true_null[n] for n in split_nodes)
                    rows.append(
                        {
                            **key,
                            "geometry": geometry_name,
                            "reference": reference,
                            "alpha": alpha,
                            "status": "ok",
                            "reason": "",
                            "clusters": decomposition["num_clusters"],
                            "ari": adjusted_rand_score(truth, predicted),
                            "n_splits": len(split_nodes),
                            "n_false_splits": false,
                            "any_false_split": bool(false),
                            "root_p": float(tests.loc[root, "p_value"]),
                            "root_true_null": true_null[root],
                            "root_rank": int(tests.loc[root, "projection_dimension"]),
                            "elapsed_seconds": perf_counter() - start,
                        }
                    )
        except Exception as exc:
            reason = repr(exc)
            status = "failed_nnls" if "Discovery NNLS did not converge" in str(exc) else "error"
            metadata.append({**key, "geometry": geometry_name, "status": status, "reason": reason})
            for reference in ("conditional", "nominal"):
                for alpha in ALPHAS:
                    rows.append(
                        {
                            **key,
                            "geometry": geometry_name,
                            "reference": reference,
                            "alpha": alpha,
                            "status": status,
                            "reason": reason,
                        }
                    )
    for method in ("tbs_diffusion_adaptive", "tbs_diffusion_adaptive_nnls"):
        params = dict(METHOD_SPECS[method].param_grid[0])
        params["metric"] = "euclidean"
        start = perf_counter()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                old = run_clustering_result(
                    discovery,
                    method,
                    params,
                    seed=replicate,
                    significance_level=0.01,
                    edge_alpha=0.001,
                    feature_space=continuous_feature_space_from_columns(columns),
                )
            row = {
                **key,
                "geometry": method,
                "reference": "current_single_matrix",
                "alpha": 0.01,
                "status": old.status,
                "reason": str(old.unsupported_reason or old.skip_reason or ""),
                "elapsed_seconds": perf_counter() - start,
            }
            if old.status == "ok":
                row.update(clusters=old.found_clusters, ari=adjusted_rand_score(truth, old.labels))
            rows.append(row)
            bundle = old.extra.get("gate_bundle")
            metadata.append(
                {
                    **key,
                    "geometry": method,
                    "status": old.status,
                    "requested_parameters": params,
                    "actual_config": asdict(bundle.metadata) if bundle is not None else None,
                    "nnls": {
                        k: v
                        for k, v in old.extra.items()
                        if k.startswith("branch_length_optimization_")
                    },
                    "warnings": sorted({str(w.message) for w in caught}),
                }
            )
        except Exception as exc:
            rows.append(
                {
                    **key,
                    "geometry": method,
                    "reference": "current_single_matrix",
                    "alpha": 0.01,
                    "status": "error",
                    "reason": repr(exc),
                }
            )
    return rows, nodes, metadata


def summarize(runs):
    rows = []
    for key, g in runs.groupby(["profile", "signal", "geometry", "reference", "alpha"]):
        ok = g[g.status.eq("ok")]
        row = dict(zip(["profile", "signal", "geometry", "reference", "alpha"], key))
        row.update(
            attempts=len(g),
            completed=len(ok),
            unavailable=len(g) - len(ok),
            mean_ari=ok.ari.mean(),
            mean_clusters=ok.clusters.mean(),
        )
        if not key[1] and len(ok):
            count = int((ok.clusters > 1).sum())
            low, high = wilson_interval(count, len(ok))
            row.update(
                false_partitions=count,
                false_partition_rate=count / len(ok),
                wilson95_low=low,
                wilson95_high=high,
            )
        if len(ok) and key[3] != "current_single_matrix":
            row.update(
                any_false_split_count=int(ok.any_false_split.sum()),
                mean_false_splits=ok.n_false_splits.mean(),
            )
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["pilot", "evaluation"], required=True)
    parser.add_argument("--null-replicates", type=int, default=250)
    parser.add_argument("--signal-replicates", type=int, default=50)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    hashes = source_hashes()
    namespace = 2026092812 if args.phase == "pilot" else 2026092813
    manifest = dict(
        arguments=vars(args),
        seed_namespace=namespace,
        source_sha256=hashes,
        python=sys.version,
        platform=platform.platform(),
        complete=False,
    )
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, default=str, indent=2) + "\n"
    )
    jobs = [
        (p, s, r, namespace)
        for p in PROFILES
        for s in (False, True)
        for r in range(args.signal_replicates if s else args.null_replicates)
    ]
    rows = []
    start = perf_counter()
    with (
        ProcessPoolExecutor(
            max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")
        ) as pool,
        gzip.open(args.output_dir / "nodes.csv.gz", "wt") as node_file,
        (args.output_dir / "fitted_context.jsonl").open("w") as meta_file,
    ):
        for completed, (runs, nodes, metadata) in enumerate(pool.map(run_dataset, jobs), 1):
            rows.extend(runs)
            if nodes:
                pd.DataFrame(nodes).to_csv(node_file, index=False, header=node_file.tell() == 0)
            for item in metadata:
                meta_file.write(json.dumps(item, default=str) + "\n")
            if completed % 50 == 0 or completed == len(jobs):
                pd.DataFrame(rows).to_csv(args.output_dir / "runs.csv", index=False)
                print(
                    f"{completed}/{len(jobs)} paired datasets; {len(rows)} runs; {perf_counter() - start:.1f}s",
                    flush=True,
                )
    runs = pd.DataFrame(rows)
    summary = summarize(runs)
    tails = []
    root_runs = runs[
        runs.alpha.eq(0.01)
        & ~runs.signal
        & runs.status.eq("ok")
        & runs.reference.ne("current_single_matrix")
    ]
    for key, g in root_runs.groupby(["profile", "geometry", "reference"]):
        for cutoff in (0.001, 0.005, 0.01, 0.05):
            count = int((g.root_p <= cutoff).sum())
            low, high = wilson_interval(count, len(g))
            tails.append(
                {
                    **dict(zip(["profile", "geometry", "reference"], key)),
                    "cutoff": cutoff,
                    "n": len(g),
                    "rejections": count,
                    "rate": count / len(g),
                    "wilson95_low": low,
                    "wilson95_high": high,
                    "one_sided_binomial_p": binomtest(
                        count, len(g), cutoff, alternative="greater"
                    ).pvalue,
                }
            )
    tails = pd.DataFrame(tails)
    for reference in ("conditional", "nominal"):
        mask = tails.reference.eq(reference)
        p = tails.loc[mask, "one_sided_binomial_p"].tolist()
        tails.loc[mask, "holm_p_36"] = multipletests(p + [1.0] * (36 - len(p)), method="holm")[1][
            : len(p)
        ]
    assert hashes == source_hashes(), "Sources changed during benchmark"
    write_diagnostic_bundle(
        output_dir=args.output_dir,
        tables={"runs": runs, "summary": summary, "root_tails": tails},
        filenames={"runs": "runs.csv", "summary": "summary.csv", "root_tails": "root_tails.csv"},
        manifest={
            **manifest,
            "complete": True,
            "paired_datasets": len(jobs),
            "elapsed_seconds": perf_counter() - start,
            "source_hashes_unchanged": True,
            "status_counts": runs.status.value_counts().to_dict(),
        },
    )
    print(runs.groupby(["geometry", "reference", "status"]).size().to_string(), flush=True)


if __name__ == "__main__":
    main()
