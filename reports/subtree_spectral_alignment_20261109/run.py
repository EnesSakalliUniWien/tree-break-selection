"""Replay the native subtree pipeline and compare explicitly mapped spectral spaces."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_mutual_info_score

from benchmarks.diagnostics.calibration import subtree_research as research
from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_spectral_alignment as align,
)
from reports.subtree_spectral_alignment_20261109.alternatives.run import REGIMES, dataset
from tree_break_selection.space_separation import diffusion
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.topology import tree_root

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(values):
    return hashlib.sha256(np.ascontiguousarray(values, dtype=np.float64).tobytes()).hexdigest()


def spectral_window(spec, features):
    return diffusion._leading_diffusion_spectrum(
        spec.eigenvalues, min(features, len(spec.eigenvalues) - 2, 50)
    )


def capture_diffusion(values, operator="symmetric"):
    captured = {}
    original = diffusion.compute_diffusion_coordinates

    def observe(weights, **kwargs):
        captured["weights"] = weights.copy()
        return original(weights, **kwargs)

    data = pd.DataFrame(values)
    space = continuous_feature_space_from_columns(data.columns)
    with patch.object(diffusion, "compute_diffusion_coordinates", observe):
        geometry = diffusion.feature_family_diffusion_geometry(
            data, space, coordinate_operator=operator
        )
    return geometry, align.reversible_diffusion_spectrum(captured["weights"])


def capture_native(values, truth, operator, link):
    captured = {}
    original_coordinates = diffusion.compute_diffusion_coordinates
    original_gate = research.run_gate_annotation_pipeline

    def coordinates(weights, **kwargs):
        captured["weights"] = weights.copy()
        return original_coordinates(weights, **kwargs)

    def gate(tree, *args, **kwargs):
        bundle = original_gate(tree, *args, **kwargs)
        context = bundle.edge_gate_result.spectral_context
        nodes = [tree_root(tree)]
        nodes += list(tree.successors(nodes[0]))
        captured["native_context"] = {
            str(node): {
                field: getattr(context, field).get(node)
                for field in (
                    "test_projection_dimensions_by_node",
                    "raw_mp_signal_counts_by_node",
                    "principal_component_projections_by_node",
                    "principal_component_eigenvalues_by_node",
                    "active_feature_counts_by_node",
                    "mp_threshold_rows_by_node",
                )
            }
            for node in nodes
        }
        captured["native_nodes"] = [str(node) for node in nodes]
        return bundle

    with (
        patch.object(diffusion, "compute_diffusion_coordinates", coordinates),
        patch.object(research, "run_gate_annotation_pipeline", gate),
    ):
        row, details = research.run_method(
            values, truth, operator=operator, linkage_method=link, trace=True
        )
    spec = align.reversible_diffusion_spectrum(captured.pop("weights"))
    details.update(captured)
    return row, details, spec


def top_basis(spec, rank=2):
    end = min(rank, len(spec.eigenvalues))
    for start, stop in spec.blocks:
        if start < end <= stop:
            end = stop
            break
    return spec.eigenvectors[:, :end]


def covariance_comparison(parent, child):
    return {
        "parent_eigenvalues": parent.eigenvalues,
        "child_eigenvalues": child.eigenvalues,
        "eigenvalue_relative_difference": float(
            np.linalg.norm(parent.eigenvalues - child.eigenvalues)
            / max(np.linalg.norm(parent.eigenvalues), np.finfo(float).tiny)
        ),
        "leading_whole_block_comparison": asdict(
            align.compare_subspaces(top_basis(parent), top_basis(child))
        ),
        "whole_block_comparisons": [
            {
                "parent_block": (a, b),
                "child_block": (c, d),
                **asdict(
                    align.compare_subspaces(parent.eigenvectors[:, a:b], child.eigenvectors[:, c:d])
                ),
            }
            for a, b in parent.blocks
            for c, d in child.blocks
        ],
    }


def reflected_covariances(a, b):
    p = len(a.eigenvalues)
    reflections = []
    for axis in range(p):
        h = np.eye(p)
        h[axis, axis] = -1
        reflections.append((f"feature_{axis}_sign", h))
    if p >= 2:
        h = np.eye(p)
        h[[0, 1]] = h[[1, 0]]
        reflections.append(("swap_features_0_1", h))
    return {
        name: asdict(
            align.compare_reflection(
                top_basis(a),
                top_basis(b),
                h,
                covariance_a=a.matrix,
                covariance_b=b.matrix,
            )
        )
        for name, h in reflections
    }


def child_rows(z, n):
    groups = {i: (i,) for i in range(n)}
    for index, (a, b, *_rest) in enumerate(z):
        groups[n + index] = tuple(sorted(groups[int(a)] + groups[int(b)]))
    return [groups[int(i)] for i in z[-1, :2]]


def diffusion_comparison(parent, child, rows, features):
    ps, pe = spectral_window(parent, features)
    cs, ce = spectral_window(child, features)
    comparisons = {}
    for name in ("symmetric_vectors", "right_vectors", "left_vectors"):
        restricted = align.restrict_parent_basis(getattr(parent, name)[:, ps:pe], rows)
        comparisons[name] = {
            "parent_rank": restricted.parent_rank,
            "restricted_rank": restricted.restricted_rank,
            "rank_loss": restricted.rank_loss,
            **asdict(align.compare_subspaces(restricted.basis, getattr(child, name)[:, cs:ce])),
        }
    return {
        "parent_window": (ps, pe),
        "child_window": (cs, ce),
        "parent_selected_eigenvalues": parent.eigenvalues[ps:pe],
        "child_selected_eigenvalues": child.eigenvalues[cs:ce],
        "comparisons": comparisons,
    }


def evaluate(values, truth, operator, link):
    reference = np.zeros(len(values), dtype=int) if truth is None else truth[2]
    native, detail, parent = capture_native(values, reference, operator, link)
    result = {
        "native": native,
        "native_context": detail.get("native_context", {}),
        "native_nodes": detail.get("native_nodes", []),
    }
    if native["status"] == "ok" and truth is not None:
        result["native_factor_ami"] = {
            name: adjusted_mutual_info_score(target, native["labels"])
            for name, target in zip(("factor_a", "contrast_b", "joint"), truth, strict=True)
        }
    start, end = spectral_window(parent, values.shape[1])
    eigenvalues = np.maximum(parent.eigenvalues[start:end], 0) ** 3
    coordinates = (
        getattr(parent, "symmetric_vectors" if operator == "symmetric" else "right_vectors")[
            :, start:end
        ]
        * eigenvalues
    )
    distance_error = float(
        np.max(np.abs(pdist(coordinates) - detail["baseline_diffusion_distance"]))
    )
    assert distance_error < 1e-9
    result["diffusion"] = {
        "production_distance_max_abs_error": distance_error,
        "eigenvalues": parent.eigenvalues,
        "degrees": parent.degrees,
        "window": (start, end),
        "residuals": parent.residuals,
        "whole_blocks": [
            {
                "block": block,
                "alignments": {k: asdict(v) for k, v in parent.block_alignment(block).items()},
            }
            for block in parent.blocks
            if block[0] < end and block[1] > start
        ],
    }
    z = np.array(detail["linkage"])
    children = child_rows(z, len(values))
    assignments = np.zeros(len(values), dtype=int)
    assignments[list(children[1])] = 1
    decomposed = align.covariance_decomposition(values, assignments)
    total = align.covariance_spectrum(values)
    detail["matrix_evidence"] = {
        "parent_weights": parent.weights,
        "parent_symmetric_operator": parent.symmetric_operator,
        "parent_markov_operator": parent.markov_operator,
        "parent_raw_covariance": total.matrix,
        "children": [],
    }
    error = float(np.linalg.norm(decomposed.total - decomposed.within - decomposed.between))
    assert error < 1e-10
    result["covariance"] = {
        "eigenvalues": total.eigenvalues,
        "decomposition_error": error,
        "within_trace": float(np.trace(decomposed.within)),
        "between_trace": float(np.trace(decomposed.between)),
        "prescribed_reflections": reflected_covariances(total, total),
    }
    child_covariances = []
    result["children"] = []
    for positions in children:
        rows = np.array(positions)
        current = {"rows": positions, "n": len(rows)}
        matrices = {"rows": positions}
        detail["matrix_evidence"]["children"].append(matrices)
        if len(rows) < 3:
            current.update(status="unsupported", reason="fewer than three child rows")
            for name in ("induced_parent_graph", "rebuilt_child_graph"):
                current[name] = {"status": "unsupported", "reason": current["reason"]}
            child_covariances.append(None)
            result["children"].append(current)
            continue
        cov = align.covariance_spectrum(values[rows])
        matrices["raw_covariance"] = cov.matrix
        child_covariances.append(cov)
        current["covariance"] = covariance_comparison(total, cov)
        try:
            induced = align.reversible_diffusion_spectrum(parent.weights[np.ix_(rows, rows)])
            current["induced_parent_graph"] = {
                "status": "ok",
                **diffusion_comparison(parent, induced, rows, values.shape[1]),
            }
            matrices["induced_weights"] = induced.weights
            matrices["induced_symmetric_operator"] = induced.symmetric_operator
            matrices["induced_markov_operator"] = induced.markov_operator
        except ValueError as exc:
            current["induced_parent_graph"] = {"status": "unsupported", "reason": str(exc)}
        except Exception as exc:
            current["induced_parent_graph"] = {
                "status": "error",
                "reason": f"{type(exc).__name__}: {exc}",
            }
        try:
            _geometry, rebuilt = capture_diffusion(values[rows], operator)
            current["rebuilt_child_graph"] = {
                "status": "ok",
                **diffusion_comparison(parent, rebuilt, rows, values.shape[1]),
            }
            matrices["rebuilt_weights"] = rebuilt.weights
            matrices["rebuilt_symmetric_operator"] = rebuilt.symmetric_operator
            matrices["rebuilt_markov_operator"] = rebuilt.markov_operator
            if current["induced_parent_graph"]["status"] == "ok":
                current["induced_vs_rebuilt_operator_frobenius"] = float(
                    np.linalg.norm(induced.symmetric_operator - rebuilt.symmetric_operator)
                )
        except ValueError as exc:
            current["rebuilt_child_graph"] = {"status": "unsupported", "reason": str(exc)}
        except Exception as exc:
            current["rebuilt_child_graph"] = {
                "status": "error",
                "reason": f"{type(exc).__name__}: {exc}",
            }
        statuses = [
            current[name]["status"] for name in ("induced_parent_graph", "rebuilt_child_graph")
        ]
        current["status"] = (
            "error" if "error" in statuses else "unsupported" if "unsupported" in statuses else "ok"
        )
        result["children"].append(current)
    if all(cov is not None for cov in child_covariances):
        result["left_right_covariance"] = covariance_comparison(*child_covariances)
        result["left_right_reflections"] = reflected_covariances(*child_covariances)
    return result, detail


def source_hashes():
    paths = sorted((ROOT / "tree_break_selection").rglob("*.py"))
    paths += [Path(__file__), HERE / "alternatives/run.py"]
    paths += [
        ROOT / "benchmarks/diagnostics/calibration/subtree_research.py",
        ROOT
        / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_spectral_alignment.py",
        ROOT / "benchmarks/shared/metrics.py",
        ROOT / "benchmarks/shared/statistical_profiles.py",
    ]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run(output, replicates=20):
    if not 1 <= replicates <= 20:
        raise ValueError("Expected 1 through 20 declared replicates")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    before = source_hashes()
    (output / "manifest.json").write_text(
        json.dumps(
            {"replicates": replicates, "regimes": REGIMES, "source_hashes": before}, indent=2
        )
    )
    count, failed, child_failed = 0, 0, 0
    native_status = Counter()
    graph_status = {name: Counter() for name in ("induced_parent_graph", "rebuilt_child_graph")}
    with gzip.open(output / "rows.jsonl.gz", "xt") as stream:
        for case in REGIMES:
            for replicate in range(replicates):
                values, truth, seed = dataset(case, replicate)
                for operator, link in (("symmetric", "average"), ("markov", "ward")):
                    record = {
                        "case": case,
                        "replicate": replicate,
                        "seed": seed,
                        "operator": operator,
                        "linkage": link,
                        "data_hash": digest(values),
                    }
                    try:
                        record["result"], detail = evaluate(values, truth, operator, link)
                        record["status"] = "ok"
                        native_status[record["result"]["native"]["status"]] += 1
                        for child in record["result"]["children"]:
                            for name in graph_status:
                                graph_status[name][child[name]["status"]] += 1
                        child_failed += sum(
                            child["status"] != "ok" for child in record["result"]["children"]
                        )
                        if replicate == 0:
                            with gzip.open(
                                output / f"trace_{case}_{operator}.json.gz", "xt"
                            ) as trace:
                                json.dump(detail, trace, default=diagnostic_json_default)
                    except Exception as exc:
                        failed += 1
                        record.update(status="error", reason=f"{type(exc).__name__}: {exc}")
                    stream.write(json.dumps(record, default=diagnostic_json_default) + "\n")
                    stream.flush()
                    count += 1
            print(f"{case}: {2 * replicates} attempts", flush=True)
    after = source_hashes()
    changed = [p for p in set(before) | set(after) if before.get(p) != after.get(p)]
    verification = {
        "attempted": count,
        "error": failed,
        "non_ok_children": child_failed,
        "native_status_counts": dict(native_status),
        "child_graph_status_counts": {name: dict(counts) for name, counts in graph_status.items()},
        "source_files_changed_during_run": changed,
    }
    (output / "verification.json").write_text(json.dumps(verification, indent=2))
    print(json.dumps(verification), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=20)
    arguments = parser.parse_args()
    run(arguments.output, arguments.replicates)
