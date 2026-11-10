"""Paired subtree research panels; no production method or default is changed."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import resource
import sys
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage

from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    selected_hierarchy_null_audit as audit,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    select_subtree_topology,
    topology_from_linkage,
)
from benchmarks.shared.benchmark_grid import benchmark_grid
from benchmarks.shared.metrics import calculate_clustering_metrics
from benchmarks.shared.runners.tbs_support import unsupported_empirical_null_reason
from benchmarks.shared.statistical_profiles import build_statistical_profiles
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.space_separation.diffusion import feature_family_diffusion_geometry
from tree_break_selection.tree.construction.hierarchical import tree_from_linkage
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.optimized_branch_lengths import (
    fit_fixed_topology_nnls_branch_lengths,
)
from tree_break_selection.tree.topology import tree_root

ROOT = Path(__file__).resolve().parents[3]
RECOVERY_CELLS = tuple(itertools.product((12, 72), (False, True), (2.5, 5.0), (False, True)))
NULL_CELLS = tuple(
    itertools.product((32, 128), (4, 40), ("identity", "ar08", "diagonal100"), ((), (0,)))
)
FACTORIAL_AXES = {
    "operator": ("symmetric", "markov"),
    "linkage_method": ("average", "ward"),
    "branch_source": ("linkage", "nnls"),
    "energy": (False, True),
}


def data_seed(*, phase: int, panel: int, cell: int, replicate: int) -> int:
    """Separate data phases/panels without selecting seeds from their results."""
    if phase not in (1, 2) or min(panel, cell, replicate) < 0:
        raise ValueError("Expected phase 1/2 and nonnegative panel, cell and replicate")
    return int(
        np.random.SeedSequence([20261109, phase, panel, cell, replicate]).generate_state(1)[0]
    )


def recovery_data(cell: int, replicate: int, *, phase: int = 1, n: int = 240):
    p, imbalanced, gap, anisotropic = RECOVERY_CELLS[cell]
    if n != 240 and imbalanced:
        raise ValueError("The prespecified imbalanced counts require 240 observations")
    if n % 4:
        raise ValueError("Balanced square panels require a multiple of four observations")
    seed = data_seed(phase=phase, panel=1, cell=cell, replicate=replicate)
    rng = np.random.default_rng(seed)
    sizes = (132, 60, 36, 12) if imbalanced else (n // 4,) * 4
    truth = np.repeat(np.arange(4), sizes)
    values = rng.standard_normal((n, p))
    if anisotropic:
        values[:, :2] = (values[:, :2] * [2, 0.5]) @ (np.array([[1, -1], [1, 1]]) / np.sqrt(2))
    means = gap * np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]]) / 2
    values[:, :2] += means[truth]
    order = rng.permutation(n)
    return values[order], truth[order], seed


def null_data(cell: int, replicate: int, *, phase: int = 2):
    n, p, covariance, path = NULL_CELLS[cell]
    seed = data_seed(phase=phase, panel=2, cell=cell, replicate=replicate)
    if covariance == "ar08":
        sigma = 0.8 ** np.abs(np.arange(p)[:, None] - np.arange(p)[None, :])
    elif covariance == "diagonal100":
        sigma = np.diag(np.geomspace(1.0, 100.0, p))
    else:
        sigma = np.eye(p)
    values = np.random.default_rng(seed).standard_normal((n, p)) @ np.linalg.cholesky(sigma).T
    return values, path, sigma, seed


def _frame(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("Expected a finite matrix with at least two rows")
    # These labels have the same order as build_tree's canonical string order.
    data = pd.DataFrame(x, index=[f"S{i:06d}" for i in range(len(x))])
    data.columns = [f"F{i}" for i in range(x.shape[1])]
    return data, continuous_feature_space_from_columns(data.columns)


def evaluate_tree(data, truth, space, z, *, branch_source="nnls", energy=True, trace=False):
    """Run native gates/traversal on a supplied local topology, retaining its heights."""
    if branch_source not in ("linkage", "nnls", "zero"):
        raise ValueError("Unknown branch source")
    tree = tree_from_linkage(z, data.index.tolist())
    fit = None
    if branch_source == "nnls":
        fit = fit_fixed_topology_nnls_branch_lengths(
            tree,
            data,
            pair_sample_size=None,
            random_state=0,
            solver_tolerance=1e-5,
            max_iterations=5000,
        )
        if not fit.applied_to_tree:
            return {"status": "error", "reason": f"NNLS not applied: {fit.status}"}, {}
    elif branch_source == "zero":
        for edge in tree.edges:
            tree.edges[edge]["branch_length"] = 0.0
    tree.populate_node_distributions(data, feature_space=space)
    bundle = run_gate_annotation_pipeline(
        tree,
        tree.annotations_df,
        leaf_data=data,
        feature_space=space,
        edge_alpha=0.001,
        sibling_alpha=0.01,
        spectral_minimum_dimension=2,
        adaptive_projection_dimension_energy_fraction=0.9 if energy else None,
        continuous_covariance_policy="guarded_within_child",
        continuous_covariance_min_child_leaf_count=8,
        edge_branch_length_variance_policy="none"
        if branch_source == "zero"
        else "normalized_branch_length",
        sibling_gate_method="projected_wald_inflation",
        enforce_internal_support_thresholds=False,
    )
    annotations = bundle.annotated_df
    tree.annotations_df = annotations
    reason = unsupported_empirical_null_reason(
        annotations, sibling_gate_method="projected_wald_inflation"
    )
    row = {"status": "unsupported" if reason else "ok", "reason": str(reason) if reason else ""}
    root = tree_root(tree)
    row.update(
        root_sibling_p=float(annotations.loc[root, "Sibling_Divergence_P_Value"]),
        root_sibling_reject=bool(annotations.loc[root, "Sibling_BH_Different"]),
        mean_positive_branch=compute_mean_branch_length(tree),
        branch_source=branch_source,
        energy=bool(energy),
    )
    if fit is not None:
        row.update(nnls_status=fit.status, nnls_relative_rmse=fit.residual_rmse_to_target_mean)
    details = {}
    if reason is None:
        decomposition = TreeDecomposition(
            tree=tree,
            gate_annotation_bundle=bundle,
            passthrough=True,
            trace_level="full" if trace else "compact",
        ).decompose_tree()
        labels, report = labels_and_report_from_decomposition(decomposition, data.index.tolist())
        metrics = calculate_clustering_metrics(report, data.index, np.asarray(truth))
        row.update(
            asdict(metrics),
            found_clusters=int(len(np.unique(labels))),
            labels=np.asarray(labels).tolist(),
        )
        if trace:
            details["traversal"] = decomposition
    if trace:
        spectral = bundle.edge_gate_result.spectral_context
        dimensions = audit.derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree, spectral_context=spectral
        )
        projections, eigenvalues = (
            audit.collect_parent_principal_component_inputs_for_sibling_tests(
                dimensions, spectral_context=spectral
            )
        )
        records, _ = audit.collect_sibling_pair_records(
            tree,
            bundle.edge_gate_result.annotated_df,
            feature_space=space,
            sibling_projection_dimensions_from_edge_comparisons=dimensions,
            parent_principal_component_projections=projections,
            parent_principal_component_eigenvalues=eigenvalues,
            adaptive_projection_dimension_energy_fraction=0.9 if energy else None,
        )
        details["raw_records"] = [asdict(record) for record in records]
        keep = [
            column for column in annotations if column.startswith(("Sibling_", "Child_Parent_"))
        ]
        details["annotations"] = annotations[keep].reset_index().to_dict("records")
        details["branches"] = [
            (str(a), str(b), float(d["branch_length"])) for a, b, d in tree.edges(data=True)
        ]
        clades = tree.compute_descendant_sets()
        position = {name: i for i, name in enumerate(data.index)}
        subtrees = {
            str(node): tuple(
                tuple(position[name] for name in sorted(clades[child]))
                for child in tree.successors(node)
            )
            for node in tree
            if tree.out_degree(node) == 2
        }
        details["subtree_profiles"] = build_statistical_profiles(
            data, subtrees=subtrees
        ).subtrees.to_dict("records")
    return row, details


def run_method(
    values,
    truth,
    *,
    operator="symmetric",
    linkage_method="average",
    branch_source="nnls",
    energy=True,
    candidate=False,
    trace=False,
):
    """Complete independent geometry-to-partition timing; truth is used only for metrics."""
    start = perf_counter()
    data, space = _frame(values)
    geometry = feature_family_diffusion_geometry(data, space, coordinate_operator=operator)
    z = linkage(geometry.distance_condensed, method=linkage_method)
    extra = {"candidate_count": 1, "candidate_search_sec": 0.0}
    if candidate:
        selected = select_subtree_topology(
            data.to_numpy(),
            current_topology=topology_from_linkage(z),
            current_linkage=z,
            leaf_labels=data.index.tolist(),
            target_metric="mean_squared_standardized_euclidean",
        )
        z = selected.selected_linkage
        extra.update(
            candidate_count=selected.candidate_count,
            candidate_search_sec=selected.elapsed_sec,
            candidate_coverage=selected.coverage,
            candidate_pair_count=selected.pair_count,
        )
    row, details = evaluate_tree(
        data, truth, space, z, branch_source=branch_source, energy=energy, trace=trace
    )
    row.update(
        extra,
        elapsed_sec=perf_counter() - start,
        process_lifetime_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        * (1 if sys.platform == "darwin" else 1024),
        operator=operator,
        linkage_method=linkage_method,
        candidate=bool(candidate),
    )
    if trace:
        details.update(
            linkage=z.tolist(),
            baseline_diffusion_distance=geometry.distance_condensed.tolist(),
            geometry_metadata=geometry.metadata,
        )
        if candidate:
            details["candidate_fits"] = [asdict(item) for item in selected.candidates]
    return row, details


def source_hashes():
    paths = set((ROOT / "tree_break_selection").rglob("*.py"))
    paths |= set((ROOT / "benchmarks/shared").rglob("*.py"))
    paths |= set((ROOT / "benchmarks/diagnostics/calibration").rglob("*.py"))
    paths |= {ROOT / "uv.lock"}
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)
    }


def _json(path, value):
    path.write_text(json.dumps(value, default=diagnostic_json_default, indent=2) + "\n")


def _safe_run(values, truth, **params):
    start = perf_counter()
    try:
        return run_method(values, truth, **params)
    except (ValueError, RuntimeError, np.linalg.LinAlgError, FloatingPointError) as exc:
        return {
            "status": "error",
            "reason": f"{type(exc).__name__}: {exc}",
            "elapsed_sec": perf_counter() - start,
        }, {}


def runtime_gate(rows):
    """All cases require three successful complete paired timings and ratio <= 5."""
    frame = pd.DataFrame(rows)
    results = []
    for (n, p), group in frame.groupby(["n", "p"]):
        medians = group.groupby("method").elapsed_sec.median()
        complete = all(
            len(group[group.method == m]) == 3
            and group.loc[group.method == m, "status"].eq("ok").all()
            and set(group.loc[group.method == m, "repeat"]) == {0, 1, 2}
            for m in ("baseline", "candidate")
        )
        complete = complete and group.groupby("repeat").seed.nunique().eq(1).all()
        complete = complete and group.groupby("repeat").data_hash.nunique().eq(1).all()
        complete = (
            complete and group.elapsed_sec.gt(0).all() and np.isfinite(group.elapsed_sec).all()
        )
        ratio = float(medians.get("candidate", np.inf) / medians.get("baseline", np.nan))
        results.append(
            {
                "n": int(n),
                "p": int(p),
                "ratio": ratio,
                "complete": bool(complete),
                "passed": bool(complete and np.isfinite(ratio) and ratio <= 5),
            }
        )
    expected = {(n, p) for n in (32, 128, 240) for p in (12, 72)}
    covered = {(r["n"], r["p"]) for r in results} == expected
    return {
        "passed": covered and all(r["passed"] for r in results),
        "complete_panel": covered,
        "cases": results,
        "scope": "geometric candidate with existing empirical sibling reference; statistical validity not established",
    }


def run_panel(stage, output, *, cells=None, replicates=50):
    """Write a new, resumable-by-external-shard panel; never overwrite evidence."""
    limits = {"preflight": 6, "recovery": 16, "factorial": 3}
    if (
        stage not in limits
        or isinstance(replicates, bool)
        or not isinstance(replicates, int)
        or not 1 <= replicates <= 50
    ):
        raise ValueError("Expected a development stage and 1 through 50 replicates")
    if cells is not None and (
        not cells
        or len(set(cells)) != len(cells)
        or any(
            isinstance(c, bool) or not isinstance(c, int) or not 0 <= c < limits[stage]
            for c in cells
        )
    ):
        raise ValueError("Expected distinct cell IDs within the panel")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    before = source_hashes()
    _json(
        output / "manifest.json",
        {
            "stage": stage,
            "cells": cells,
            "replicates": 3 if stage == "preflight" else replicates,
            "source_hashes": before,
            "seed_key": [20261109],
            "phase": "development",
        },
    )
    rows = []
    with (output / "runs.jsonl").open("x") as stream:

        def record(base, values, truth, params, trace=False):
            row, detail = _safe_run(values, truth, trace=trace, **params)
            row = {
                **base,
                **row,
                "data_hash": hashlib.sha256(
                    np.ascontiguousarray(values, dtype=np.float64).tobytes()
                ).hexdigest(),
                "truth_hash": hashlib.sha256(
                    np.ascontiguousarray(truth, dtype=np.int64).tobytes()
                ).hexdigest(),
            }
            rows.append(row)
            stream.write(json.dumps(row, default=diagnostic_json_default) + "\n")
            stream.flush()
            if trace and detail:
                _json(output / f"trace_{base['run_id']}.json", detail)
            return row

        if stage == "preflight":
            for cell, (n, p) in enumerate(itertools.product((32, 128, 240), (12, 72))):
                if cells is not None and cell not in cells:
                    continue
                recovery_cell = RECOVERY_CELLS.index((p, False, 5.0, False))
                x, truth, seed = recovery_data(recovery_cell, 0, n=n)
                for repeat in range(3):
                    for method in (
                        ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
                    ):
                        record(
                            dict(
                                run_id=f"preflight_{cell}_{repeat}_{method}",
                                cell=cell,
                                n=n,
                                p=p,
                                seed=seed,
                                repeat=repeat,
                                method=method,
                            ),
                            x,
                            truth,
                            {"candidate": method == "candidate"},
                        )
                print(f"preflight cell {cell} complete", flush=True)
            _json(output / "runtime_gate.json", runtime_gate(rows))
        elif stage == "recovery":
            for cell in range(len(RECOVERY_CELLS)) if cells is None else cells:
                for rep in range(replicates):
                    x, truth, seed = recovery_data(cell, rep)
                    profiles = build_statistical_profiles(x, true_labels=truth)
                    if rep == 0:
                        profiles.groups.to_csv(output / f"profiles_{cell}.csv", index=False)
                    for method, operator, link in (
                        ("symmetric_average", "symmetric", "average"),
                        ("markov_ward", "markov", "ward"),
                    ):
                        record(
                            dict(
                                run_id=f"recovery_{cell}_{rep}_{method}",
                                cell=cell,
                                replicate=rep,
                                seed=seed,
                                n=len(x),
                                p=x.shape[1],
                                method=method,
                            ),
                            x,
                            truth,
                            {"operator": operator, "linkage_method": link},
                            trace=rep == 0,
                        )
                print(f"recovery cell {cell} complete", flush=True)
        elif stage == "factorial":
            factors = benchmark_grid(
                benchmark_class="subtree_research",
                grid_name="operator_linkage_branch_rank",
                base_params={},
                axes=FACTORIAL_AXES,
            )
            for cell in range(3) if cells is None else cells:
                for rep in range(replicates):
                    if cell == 0:
                        seed = data_seed(phase=1, panel=3, cell=cell, replicate=rep)
                        x = np.random.default_rng(seed).standard_normal((64, 12))
                        truth = np.zeros(64, dtype=int)
                    else:
                        c = RECOVERY_CELLS.index((12, cell == 2, 5.0, False))
                        x, truth, seed = recovery_data(c, rep)
                    for index, params in enumerate(factors):
                        selected = {key: params[key] for key in FACTORIAL_AXES}
                        record(
                            dict(
                                run_id=f"factorial_{cell}_{rep}_{index}",
                                cell=cell,
                                replicate=rep,
                                seed=seed,
                                n=len(x),
                                p=x.shape[1],
                                method=params["benchmark_run_id"],
                            ),
                            x,
                            truth,
                            selected,
                            trace=rep == 0,
                        )
                print(f"factorial cell {cell} complete", flush=True)
        else:
            raise ValueError("Unknown development stage")
    frame = pd.DataFrame(rows)
    frame.drop(columns=["labels"], errors="ignore").to_csv(output / "runs.csv", index=False)
    after = source_hashes()
    changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p))
    _json(
        output / "verification.json",
        {
            "attempted": len(rows),
            "source_files_changed_during_run": changed,
            "status_counts": frame.status.value_counts().to_dict(),
        },
    )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("preflight", "recovery", "factorial"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells", type=int, nargs="+")
    parser.add_argument("--replicates", type=int, default=50)
    args = parser.parse_args()
    if not 1 <= args.replicates <= 50:
        parser.error(
            "Development replicates must be 1 through 50; validation has a separate seed phase"
        )
    run_panel(args.stage, args.output, cells=args.cells, replicates=args.replicates)


if __name__ == "__main__":
    main()
