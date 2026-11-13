"""Run the frozen bounded coverage and six-leaf perturbation checks."""

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.linalg import helmert
from scipy.spatial.distance import pdist
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    select_subtree_topology,
    topology_from_linkage,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_geometry import (
    spectral_subspaces,
    ward_topology_family,
)
from tree_break_selection.tree.optimized_branch_lengths import (
    BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
)


def bounded():
    records = []
    for n, p in ((128, 12), (240, 12), (240, 72), (256, 4)):
        x = np.random.default_rng([2026110911, n, p]).normal(size=(n, p))
        current_z = linkage(x, method="ward")
        current = topology_from_linkage(current_z)
        for metric in ("euclidean", BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN):
            row = {"n": n, "p": p, "metric": metric}
            try:
                result = select_subtree_topology(
                    x, current_topology=current, current_linkage=current_z, target_metric=metric
                )
                current_fit = next(c.fit for c in result.candidates if c.topology == current)
                selected_fit = next(
                    c.fit for c in result.candidates if c.topology == result.selected_topology
                )
                row.update(
                    status="ok",
                    count=result.candidate_count,
                    spans=result.spectral_span_count,
                    pairs=result.pair_count,
                    coverage=result.coverage,
                    elapsed_sec=result.elapsed_sec,
                    changed=result.selected_topology != current,
                    current_rmse=current_fit.rmse,
                    selected_rmse=selected_fit.rmse,
                )
            except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as error:
                row.update(status="failed", error=f"{type(error).__name__}: {error}")
            records.append(row)
            print(json.dumps(row), flush=True)
    return records


def canonical_ward(values):
    distances = pdist(values)
    maximum = distances.max()
    scaled = np.round(distances / maximum, 12) if maximum > 0 else distances
    return topology_from_linkage(linkage(scaled, method="ward"))


def resampling():
    row = np.repeat([-3.0, 3.0], 3)
    column = np.tile([-2.0, 0.0, 2.0], 2)
    alternating = 0.15 * np.tile([-1.0, 1.0], 3)
    x = np.column_stack([row, row, column, column, alternating])
    family = set().union(*(ward_topology_family(s.coordinates) for s in spectral_subspaces(x)))
    mean = x.mean(axis=0)
    basis = helmert(len(x), full=False)
    centered = basis @ (x - mean)
    rms = float(np.sqrt(np.mean((x - mean) ** 2)))
    records = []
    for index, regime in enumerate(
        ("feature_bootstrap", "small_noise", "large_noise", "row_control", "gaussian_orbit")
    ):
        rng = np.random.default_rng([2026110912, index])
        counts = Counter()
        for _ in range(2000):
            if regime == "feature_bootstrap":
                data = x[:, rng.integers(x.shape[1], size=x.shape[1])]
            elif regime in {"small_noise", "large_noise"}:
                sigma = rms * (0.05 if regime == "small_noise" else 0.5)
                data = x + rng.normal(scale=sigma, size=x.shape)
            elif regime == "row_control":
                order = rng.permutation(len(x))
                data = x[order][np.argsort(order)]
            else:
                rotation, upper = np.linalg.qr(rng.normal(size=(len(basis), len(basis))))
                rotation *= np.sign(np.diag(upper))
                data = mean + basis.T @ rotation @ centered
            counts[canonical_ward(data)] += 1
        record = {
            "regime": regime,
            "replicates": sum(counts.values()),
            "distinct": len(counts),
            "outside_replicates": sum(v for k, v in counts.items() if k not in family),
            "outside_distinct": sum(k not in family for k in counts),
            "counts": [(k, v) for k, v in sorted(counts.items(), key=lambda kv: repr(kv[0]))],
        }
        records.append(record)
        print(json.dumps({k: v for k, v in record.items() if k != "counts"}), flush=True)
    return {
        "baseline_family_count": len(family),
        "baseline_family": sorted(family, key=repr),
        "baseline_full_ward": canonical_ward(x),
        "regimes": records,
    }


if __name__ == "__main__":
    command = sys.argv[1]
    with threadpool_limits(limits=1):
        result = {"coverage": bounded, "resampling": resampling}[command]()
    (HERE / f"{command}.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
