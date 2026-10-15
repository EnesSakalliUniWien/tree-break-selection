"""Compare eigenvalue blocks, their sample geometry, trees, and TBS partitions.

Descriptive matched-case analysis; truth labels never select blocks or methods.
Uses the standard TBS runner, including its unsupported-outcome contract.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from itertools import combinations
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr, wasserstein_distance
from sklearn.metrics import adjusted_rand_score

from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from tree_break_selection.space_separation import (
    block_adaptive_diffusion_geometry,
    separate_adaptive_cosine_space,
)
from tree_break_selection.tree.construction import SUPPORTED_LINKAGE_METHODS, build_tree
from tree_break_selection.tree.topology import compute_descendant_leaf_sets


def distance_comparison(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    """Compare aligned sample-pair distances, separating scale from shape."""
    left_rms = float(np.sqrt(np.mean(left**2)))
    right_rms = float(np.sqrt(np.mean(right**2)))
    correlation = (
        float(spearmanr(left, right).statistic)
        if np.ptp(left) > 0 and np.ptp(right) > 0
        else np.nan
    )
    return {
        "distance_spearman": correlation,
        "distance_rms_ratio": left_rms / right_rms if right_rms > 0 else np.nan,
        "normalized_distance_rmse": (
            float(np.sqrt(np.mean((left / left_rms - right / right_rms) ** 2)))
            if min(left_rms, right_rms) > 0
            else np.nan
        ),
    }


def rooted_clades(tree) -> set[frozenset]:
    """Exclude trivial leaf and whole-tree clades from rooted RF comparison."""
    n = len(tree.get_leaves())
    return {s for s in compute_descendant_leaf_sets(tree).values() if 1 < len(s) < n}


def rooted_rf(left: set[frozenset], right: set[frozenset]) -> float:
    denominator = len(left) + len(right)
    return len(left ^ right) / denominator if denominator else 0.0


def compare_spectra(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    """Equal mass per eigenvalue; log distance is invariant to a common scale."""
    return {
        "boundary_ratio": float(left.min() / right.max()),
        "log_wasserstein": float(wasserstein_distance(np.log(left), np.log(right))),
        "mean_eigenvalue_ratio": float(left.mean() / right.mean()),
    }


def run_candidate(data, feature_space, distances, method, strict):
    """Retain trees even when the standard gate runner cannot supply labels."""
    options = dict(tree_linkage_method=method)
    start = perf_counter()
    tree = labels = None
    record = {"status": "error", "error": "", "n_clusters": np.nan}
    try:
        result = run_tbs_on_distance(
            data,
            distances,
            0.01,
            edge_alpha=0.001,
            feature_space=feature_space,
            enforce_internal_support_thresholds=strict,
            **options,
        )
        tree = result.extra["tree"]
        record["status"] = result.status
        if result.status == "ok":
            labels = np.asarray(result.labels)
            if labels.shape != (len(data),):
                raise ValueError("Runner labels do not cover every input sample.")
            record["n_clusters"] = len(np.unique(labels))
        annotations = result.extra.get("annotations")
        if annotations is not None:
            calibration = annotations.get("Sibling_Gate_P_Value_Calibration")
            if calibration is not None:
                record["calibration_counts"] = json.dumps(
                    calibration.dropna().astype(str).value_counts().to_dict(), sort_keys=True
                )
    except Exception as exc:
        labels = None
        record.update(status="failed_gate", error=repr(exc), n_clusters=np.nan)
        try:
            tree = build_tree(
                data,
                distances,
                linkage_method=method,
            ).tree
        except Exception as tree_exc:
            tree = None
            record.update(status="failed_tree", error=repr(tree_exc))
    record["runtime_seconds"] = perf_counter() - start
    return record, rooted_clades(tree) if tree is not None else None, labels


def analyze_case(case, args, tables):
    inputs = prepare_case_inputs(case, selected_methods=[])
    data, truth = inputs.data, inputs.labels
    base = {"case": case["name"], "seed": case["seed"]}
    tables["cases"].append(
        {
            **base,
            "generator": case["generator"],
            "n_samples": len(data),
            "n_features": data.shape[1],
            "true_clusters": len(np.unique(truth)),
            "data_sha256": hashlib.sha256(data.to_csv().encode()).hexdigest(),
        }
    )
    space = separate_adaptive_cosine_space(
        data,
        weighting="binary",
        max_rank=args.max_rank,
        min_segment_length=4,
        max_segments=8,
    )
    for index, eigenvalue in enumerate(space.eigenvalues, 1):
        tables["spectra"].append({**base, "mode": index, "eigenvalue": eigenvalue})
    eigenvalues, geometries, runs = {}, {}, {}
    for block in space.blocks:
        name = block.block_name
        values = space.eigenvalues[block.block_start - 1 : block.block_end]
        eigenvalues[name] = values
        tables["blocks"].append(
            {
                **base,
                **vars(block),
                "eigenvalue_max": values.max(),
                "eigenvalue_min": values.min(),
                "eigenvalue_mean": values.mean(),
                "width_log_ratio": float(np.log(values.max() / values.min())),
                "energy_fraction_full_operator": float(values.sum() / len(data)),
                "energy_fraction_retained": float(values.sum() / space.eigenvalues.sum()),
                "retained_rank": len(space.eigenvalues),
                "retained_energy_fraction": float(space.eigenvalues.sum() / len(data)),
            }
        )
        coordinates = space.coordinates(block)
        for geometry in args.geometries:
            key = {**base, "block": name, "geometry": geometry}
            try:
                metadata = {}
                if geometry == "raw":
                    distances = pdist(coordinates)
                else:
                    diffusion = block_adaptive_diffusion_geometry(
                        coordinates,
                        k_neighbors=15,
                        diffusion_time=3,
                        n_components=30,
                        metric="euclidean",
                        bandwidth_type="-1/(d+2)",
                        epsilon="median",
                    )
                    distances, metadata = diffusion.distance_condensed, diffusion.metadata
                if not np.isfinite(distances).all() or np.max(distances) <= 0:
                    raise ValueError("Nonfinite or identically zero block distances.")
                tables["geometries"].append({**key, "status": "ok", **metadata})
                geometries[name, geometry] = distances
            except Exception as exc:
                tables["geometries"].append(
                    {**key, "status": "failed_geometry", "error": repr(exc)}
                )
                for method in args.methods:
                    tables["runs"].append(
                        {**key, "method": method, "status": "failed_geometry", "error": repr(exc)}
                    )
                continue
            for method in args.methods:
                record, clades, labels = run_candidate(
                    data,
                    inputs.metadata.get("feature_space"),
                    distances,
                    method,
                    args.enforce_internal_support_thresholds,
                )
                record = {**key, "method": method, **record}
                record["ari_truth"] = (
                    float(adjusted_rand_score(truth, labels)) if labels is not None else np.nan
                )
                tables["runs"].append(record)
                runs[name, geometry, method] = clades, labels
                if clades is not None:
                    tables["trees"].append(
                        {
                            **key,
                            "method": method,
                            "clades": json.dumps(sorted([sorted(map(str, c)) for c in clades])),
                        }
                    )
                if labels is not None:
                    for sample, predicted, expected in zip(data.index, labels, truth, strict=True):
                        tables["assignments"].append(
                            {
                                **key,
                                "method": method,
                                "sample": sample,
                                "cluster": predicted,
                                "truth": expected,
                            }
                        )
                print(f"  {name} {geometry} {method}: {record['status']}", flush=True)
    for left, right in combinations(eigenvalues, 2):
        tables["block_pairs"].append(
            {
                **base,
                "left_block": left,
                "right_block": right,
                **compare_spectra(eigenvalues[left], eigenvalues[right]),
            }
        )
    for left, right in combinations(geometries, 2):
        if left[1] != right[1] and left[0] != right[0]:
            continue
        tables["geometry_pairs"].append(
            {
                **base,
                "left_block": left[0],
                "right_block": right[0],
                "left_geometry": left[1],
                "right_geometry": right[1],
                **distance_comparison(geometries[left], geometries[right]),
            }
        )
    for left, right in combinations(runs, 2):
        changed = [i for i in range(3) if left[i] != right[i]]
        if len(changed) != 1:
            continue
        left_clades, left_labels = runs[left]
        right_clades, right_labels = runs[right]
        tables["run_pairs"].append(
            {
                **base,
                "comparison": ("blocks", "geometries", "tree_methods")[changed[0]],
                "left_block": left[0],
                "right_block": right[0],
                "left_geometry": left[1],
                "right_geometry": right[1],
                "left_method": left[2],
                "right_method": right[2],
                "normalized_rooted_rf": (
                    rooted_rf(left_clades, right_clades)
                    if left_clades is not None and right_clades is not None
                    else np.nan
                ),
                "ari_between_partitions": (
                    float(adjusted_rand_score(left_labels, right_labels))
                    if left_labels is not None and right_labels is not None
                    else np.nan
                ),
            }
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", default="method_proof")
    parser.add_argument("--case-names", nargs="+")
    parser.add_argument("--max-rank", type=int, default=40)
    parser.add_argument(
        "--methods",
        nargs="+",
        default=list(SUPPORTED_LINKAGE_METHODS),
        choices=SUPPORTED_LINKAGE_METHODS,
    )
    parser.add_argument(
        "--geometries",
        nargs="+",
        default=["raw", "adaptive_diffusion"],
        choices=["raw", "adaptive_diffusion"],
    )
    parser.add_argument("--enforce-internal-support-thresholds", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    cases = get_test_cases_by_suite(args.suite)
    if args.case_names:
        selected = set(args.case_names)
        missing = selected - {c["name"] for c in cases}
        if missing:
            parser.error(f"Unknown cases: {sorted(missing)}")
        cases = [c for c in cases if c["name"] in selected]
    # Preserve an existing run, including incomplete runs, for audit.
    args.output_dir.mkdir(parents=True, exist_ok=False)
    tables = {
        name: []
        for name in (
            "cases",
            "spectra",
            "blocks",
            "block_pairs",
            "geometries",
            "runs",
            "trees",
            "assignments",
            "geometry_pairs",
            "run_pairs",
        )
    }
    source_paths = [Path(__file__), Path("pyproject.toml"), Path("uv.lock")]
    for root in ("tree_break_selection", "benchmarks/shared"):
        source_paths.extend(sorted(Path(root).rglob("*.py")))
    manifest = {
        "schema": "spectral_block_comparison/v1",
        "arguments": vars(args),
        "cases": cases,
        "python": platform.python_version(),
        "packages": {
            p: importlib.metadata.version(p)
            for p in (
                "numpy",
                "pandas",
                "scipy",
                "scikit-learn",
                "pydiffmap",
                "scikit-bio",
            )
        },
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], text=True),
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "weighting": "binary (raw values, including signed continuous features)",
        "edge_alpha": 0.001,
        "sibling_alpha": 0.01,
        "branch_length_variance_policy": "none",
        "limitations": [
            "Descriptive one-seed-per-case panel; pairwise rows are dependent.",
            "No truth-label block/method selection and no calibration-validity claim.",
            "NNLS refits lengths, not topology; branch-time variants are not tested.",
            "Only the retained positive cosine spectrum is compared; max_rank is a tie-aware target.",
            "Unequal block dimensions are retained; geometry differences can reflect dimension and eigenvector structure.",
        ],
    }
    for index, case in enumerate(cases, 1):
        print(f"[{index}/{len(cases)}] {case['name']}", flush=True)
        analyze_case(case, args, tables)
        write_diagnostic_bundle(
            output_dir=args.output_dir,
            tables={key: pd.DataFrame(rows) for key, rows in tables.items()},
            filenames={key: f"{key}.csv" for key in tables},
            manifest={**manifest, "completed_cases": index, "complete": index == len(cases)},
        )
    print(f"Wrote {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
