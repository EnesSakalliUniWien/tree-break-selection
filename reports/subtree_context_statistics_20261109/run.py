"""Trace matched clades through native gates under fixed subspace alternatives."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib
import json
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration import subtree_research as research
from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_spectral_alignment as align,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    topology_from_linkage,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_context import (
    compare_node_contexts,
)
from reports.subtree_spectral_alignment_20261109.alternatives.run import (
    REGIMES,
    build_ward,
    dataset,
    spans,
)
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.spectral_context import (
    split_edge_test_rows,
)
from tree_break_selection.tree.distributions import resolve_node_continuous_covariance_by_block

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EDGE_MODULE = importlib.import_module(
    "tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence."
    "child_parent_divergence_annotation.child_parent_divergence_annotation"
)


def key(rows):
    return ",".join(str(i) for i in sorted(rows))


def proposals(values):
    data, space = research._frame(values)
    for name, operator, method in (
        ("symmetric_average", "symmetric", "average"),
        ("markov_ward", "markov", "ward"),
    ):
        geometry = research.feature_family_diffusion_geometry(
            data, space, coordinate_operator=operator
        )
        yield (
            name,
            linkage(geometry.distance_condensed, method=method),
            geometry.coordinates,
            None,
            {},
        )
    cov = align.covariance_spectrum(values)
    all_spans = list(spans(values))
    leading = next(
        blocks
        for blocks, coords in all_spans
        if coords.shape[1] >= 2
        and blocks[0][0] == 0
        and all(a[1] == b[0] for a, b in zip(blocks, blocks[1:]))
    )
    for blocks, coordinates in all_spans:
        if len(blocks) != 1 and blocks != leading:
            continue
        columns = [i for start, stop in blocks for i in range(start, stop)]
        name = "cov_" + "_".join(f"{a}_{b}" for a, b in blocks)
        z, _tree = build_ward(coordinates)
        yield (
            name,
            z,
            coordinates,
            cov.eigenvectors[:, columns],
            {
                "covariance_blocks": blocks,
                "proposal_covariance_eigenvalues": cov.eigenvalues[columns],
                "proposal_covariance_trace_share": float(
                    cov.eigenvalues[columns].sum() / cov.eigenvalues.sum()
                ),
            },
        )
    z, _tree = build_ward(values)
    yield (
        "cov_full",
        z,
        values,
        np.eye(values.shape[1]),
        {
            "proposal_covariance_eigenvalues": cov.eigenvalues,
            "proposal_covariance_trace_share": 1.0,
        },
    )


def spectral_record(context, node):
    return {
        "eigenvalues": context.principal_component_eigenvalues_by_node.get(node),
        "projection": context.principal_component_projections_by_node.get(node),
        "raw_mp": context.raw_mp_signal_counts_by_node.get(node),
        "test_dimension": context.test_projection_dimensions_by_node.get(node),
    }


def covariance_span_record(raw_cov, source_basis, source_coordinates):
    positive = np.count_nonzero(
        raw_cov.eigenvalues
        > max(raw_cov.eigenvalues[0], 0)
        * np.finfo(float).eps
        * max(source_coordinates.shape[0], source_basis.shape[0])
    )
    # A block crossing the numerical zero boundary cannot supply a stable partial span.
    available = max((stop for _, stop in raw_cov.blocks if stop <= positive), default=0)
    rank = min(source_basis.shape[1], available)
    for start, stop in raw_cov.blocks:
        if start < rank <= stop:
            rank = stop
            break
    trace = float(raw_cov.eigenvalues.sum())
    return {
        "inherited_vs_rebuilt_covariance_span": asdict(
            align.compare_subspaces(source_basis, raw_cov.eigenvectors[:, :rank])
        ),
        "inherited_local_variance_share": (
            float(np.var(source_coordinates, axis=0, ddof=1).sum() / trace) if trace > 0 else None
        ),
    }


def fit_tree(values, truth, z, source_coordinates, source_basis):
    data, space = research._frame(values)
    captured = {}
    original_spectral = EDGE_MODULE.compute_child_parent_spectral_context
    original_gate = research.run_gate_annotation_pipeline

    def spectral(tree, *args, **kwargs):
        result = original_spectral(tree, *args, **kwargs)
        name = "edge_half" if kwargs.get("basis_row_indices_by_node") is not None else "all_leaf"
        captured[name] = result
        return result

    def gate(tree, *args, **kwargs):
        result = original_gate(tree, *args, **kwargs)
        captured["tree"] = tree
        return result

    with (
        patch.object(EDGE_MODULE, "compute_child_parent_spectral_context", spectral),
        patch.object(research, "run_gate_annotation_pipeline", gate),
    ):
        metrics, detail = research.evaluate_tree(data, truth, space, z, trace=True)
    if "tree" not in captured:
        return {"metrics": metrics, "nodes": {}, "trace": detail}
    tree = captured["tree"]
    position = {name: i for i, name in enumerate(data.index)}
    leaves = {
        node: tuple(sorted(position[label] for label in labels))
        for node, labels in tree.compute_descendant_sets().items()
    }
    test_rows, basis_rows = split_edge_test_rows(tree, data)
    mean_length = compute_mean_branch_length(tree)
    records = {str(row["parent"]): row for row in detail["raw_records"]}
    annotations = {str(row["node_id"]): row for row in detail["annotations"]}
    traversal = detail.get("traversal", {})
    visited = {str(row["node_id"]): row for row in traversal.get("traversal_trace", [])}
    nodes = {}
    for node in tree:
        children = list(tree.successors(node))
        if len(children) != 2:
            continue
        parents = list(tree.predecessors(node))
        parent = parents[0] if parents else None
        rows = list(leaves[node])
        raw_cov = align.covariance_spectrum(values[rows])
        covariances = resolve_node_continuous_covariance_by_block(
            tree,
            node,
            space,
            continuous_covariance_policy="guarded_within_child",
            continuous_covariance_min_child_leaf_count=8,
        )
        ann = annotations[str(node)]
        raw = records.get(str(node))
        trace = visited.get(str(node))
        eligible = (
            not bool(ann["Sibling_Divergence_Skipped"])
            and not bool(ann["Sibling_Divergence_Invalid"])
            and ann["Sibling_Gate_P_Value_Role"] == "active_traversal_sibling_gate"
            and np.isfinite(ann["Sibling_Divergence_P_Value"])
        )
        active = eligible and np.isfinite(ann["Sibling_Divergence_P_Value_Corrected"])
        record = {
            "native_node_id": str(node),
            "leaves": leaves[node],
            "parent": leaves[parent] if parent is not None else None,
            "children": sorted(leaves[child] for child in children),
            "raw_covariance_eigenvalues": raw_cov.eigenvalues,
            "guarded_covariance_eigenvalues": {
                name: np.linalg.eigvalsh(cov) if cov.ndim == 2 else cov
                for name, cov in covariances.items()
            },
            "sibling_all_leaf_context": spectral_record(captured["all_leaf"], node),
            "outgoing_edge_half_context": spectral_record(captured["edge_half"], node),
            "incoming_edge_half_context": spectral_record(captured["edge_half"], parent)
            if parent is not None
            else None,
            "incoming_test_rows": sorted(test_rows[node]) if parent is not None else None,
            "incoming_parent_test_rows": sorted(
                i for child in tree.successors(parent) for i in test_rows[child]
            )
            if parent is not None
            else None,
            "incoming_basis_rows": sorted(basis_rows.get(parent, []))
            if parent is not None
            else None,
            "outgoing_basis_rows": sorted(basis_rows.get(node, [])),
            "outgoing_edges": [
                {
                    "leaves": leaves[child],
                    "test_rows": sorted(test_rows[child]),
                    "annotations": annotations[str(child)],
                }
                for child in children
            ],
            "annotations": ann,
            "raw_sibling": raw,
            "sibling_gate_state": "active_reject"
            if active and ann["Sibling_BH_Different"]
            else "active_nonreject"
            if active
            else "eligible_not_bh_tested"
            if eligible
            else "unavailable",
            "visited": trace is not None,
            "traversal_decision": None if trace is None else trace["decision"],
            "final_boundary": None if trace is None else trace["final_boundary"],
            "mean_positive_branch": mean_length,
            "sibling_branch_factor": None
            if raw is None
            else 1.0
            if mean_length is None
            else 1 + raw["branch_length_sum"] / (2 * mean_length),
            "incoming_branch_factor": None
            if parent is None
            else 1.0
            if mean_length is None
            else 1 + tree.edges[parent, node]["branch_length"] / mean_length,
        }
        if source_basis is not None:
            record.update(covariance_span_record(raw_cov, source_basis, source_coordinates[rows]))
        nodes[key(rows)] = record
    return {"metrics": metrics, "nodes": nodes, "trace": detail}


def compare(reference_tree, alternative_tree, reference, alternative):
    for context in compare_node_contexts(reference_tree, alternative_tree):
        r = reference["nodes"][key(context.node)]
        a = alternative["nodes"].get(key(context.node))
        out = asdict(context)
        out.update(
            reference_node_id=r["native_node_id"],
            alternative_node_id=None if a is None else a["native_node_id"],
            node_size=len(context.node),
        )
        for name in (
            "sibling_gate_state",
            "visited",
            "traversal_decision",
            "final_boundary",
            "sibling_branch_factor",
        ):
            out["reference_" + name] = r[name]
            out["alternative_" + name] = None if a is None else a[name]
        fields = (
            "Sibling_Raw_Test_Statistic",
            "Sibling_Test_Statistic",
            "Sibling_Degrees_of_Freedom",
            "Sibling_Inflation_Factor",
            "Sibling_Divergence_P_Value",
            "Sibling_Divergence_P_Value_Corrected",
            "Child_Parent_Divergence_Test_Statistic",
            "Child_Parent_Divergence_df",
            "Child_Parent_Divergence_P_Value",
            "Child_Parent_Divergence_Tested",
            "Child_Parent_Divergence_Ancestor_Blocked",
            "Child_Parent_Divergence_Significant",
        )
        for name in fields:
            incoming_root = name.startswith("Child_Parent_") and context.original_parent is None
            out["reference_" + name] = None if incoming_root else r["annotations"][name]
            out["alternative_" + name] = (
                None if a is None or incoming_root else a["annotations"][name]
            )
        for name in ("stat", "degrees_of_freedom", "p_value", "reference_scale"):
            out["reference_collected_sibling_" + name] = (
                None if r["raw_sibling"] is None else r["raw_sibling"][name]
            )
            out["alternative_collected_sibling_" + name] = (
                None if a is None or a["raw_sibling"] is None else a["raw_sibling"][name]
            )
        if a is not None:
            np.testing.assert_allclose(
                r["raw_covariance_eigenvalues"], a["raw_covariance_eigenvalues"], atol=1e-12
            )
            out["same_incoming_test_rows"] = (
                None
                if context.original_parent is None
                else r["incoming_test_rows"] == a["incoming_test_rows"]
                and r["incoming_parent_test_rows"] == a["incoming_parent_test_rows"]
            )
            out["same_incoming_basis_rows"] = (
                None
                if context.original_parent is None
                else r["incoming_basis_rows"] == a["incoming_basis_rows"]
            )
            out["same_outgoing_basis_rows"] = r["outgoing_basis_rows"] == a["outgoing_basis_rows"]
        yield out


def source_hashes():
    paths = sorted((ROOT / "tree_break_selection").rglob("*.py"))
    paths += [
        Path(__file__),
        ROOT / "benchmarks/diagnostics/calibration/subtree_research.py",
        ROOT / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_context.py",
        ROOT
        / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_spectral_alignment.py",
        ROOT
        / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_candidate.py",
        ROOT / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_geometry.py",
        ROOT / "benchmarks/shared/metrics.py",
        ROOT / "benchmarks/shared/statistical_profiles.py",
        ROOT / "reports/subtree_spectral_alignment_20261109/alternatives/run.py",
    ]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run(output, replicates=4):
    if not 1 <= replicates <= 4:
        raise ValueError("Expected 1 through 4 declared replicates")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    before = source_hashes()
    (output / "manifest.json").write_text(
        json.dumps({"replicates": replicates, "source_hashes": before}, indent=2)
    )
    attempts = []
    count = 0
    with gzip.open(output / "comparisons.jsonl.gz", "xt") as comparisons:
        for case in REGIMES:
            for replicate in range(replicates):
                values, truth, seed = dataset(case, replicate)
                labels = np.zeros(len(values), dtype=int) if truth is None else truth[2]
                fitted = {}
                for name, z, coords, basis, metadata in proposals(values):
                    base = {
                        "case": case,
                        "replicate": replicate,
                        "seed": seed,
                        "method": name,
                        "data_hash": hashlib.sha256(values.tobytes()).hexdigest(),
                        **metadata,
                    }
                    try:
                        result = fit_tree(values, labels, z, coords, basis)
                        base.update(status=result["metrics"]["status"], metrics=result["metrics"])
                        if result["nodes"]:
                            fitted[name] = (topology_from_linkage(z), result)
                        with gzip.open(
                            output / f"fit_{case}_{replicate}_{name}.json.gz", "xt"
                        ) as stream:
                            json.dump(
                                {**base, "linkage": z, "source_distances": pdist(coords), **result},
                                stream,
                                default=diagnostic_json_default,
                            )
                    except Exception as exc:
                        base.update(status="error", reason=f"{type(exc).__name__}: {exc}")
                    attempts.append(base)
                for reference_name in ("symmetric_average", "markov_ward"):
                    if reference_name not in fitted:
                        continue
                    rt, r = fitted[reference_name]
                    for name, (at, a) in fitted.items():
                        if not name.startswith("cov_"):
                            continue
                        for comparison in compare(rt, at, r, a):
                            record = {
                                "case": case,
                                "replicate": replicate,
                                "reference_method": reference_name,
                                "alternative_method": name,
                                **comparison,
                            }
                            comparisons.write(
                                json.dumps(record, default=diagnostic_json_default) + "\n"
                            )
                            count += 1
                print(f"{case} {replicate}: {len(fitted)} traces", flush=True)
    (output / "attempts.json").write_text(
        json.dumps(attempts, default=diagnostic_json_default, indent=2)
    )
    after = source_hashes()
    changed = [p for p in set(before) | set(after) if before.get(p) != after.get(p)]
    verification = {
        "attempts": len(attempts),
        "status_counts": {
            s: sum(r["status"] == s for r in attempts)
            for s in sorted({r["status"] for r in attempts})
        },
        "node_comparisons": count,
        "source_files_changed": changed,
    }
    (output / "verification.json").write_text(json.dumps(verification, indent=2))
    print(json.dumps(verification), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replicates", type=int, default=4)
    args = parser.parse_args()
    run(args.output, args.replicates)
