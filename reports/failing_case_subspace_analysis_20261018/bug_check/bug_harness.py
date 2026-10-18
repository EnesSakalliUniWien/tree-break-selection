"""Run canonical tbs once per failing case through the real benchmark path and check each stage.

Captures the distance, tree, annotations and labels by wrapping the dispatcher that
run_single_method_once calls. Writes bug_checks.csv next to this script; per-case
annotation tables are written only when a directory is passed as the first argument.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import chi2
from sklearn.metrics import adjusted_rand_score

import benchmarks.shared.util.method_execution as method_execution
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from benchmarks.shared.util.method_selection import resolve_selected_methods_and_param_sets
from benchmarks.shared.util.method_sets import DEFAULT_METHODS
from tree_break_selection.hierarchy_analysis.cluster_assignments import (
    build_sample_cluster_assignments,
)

warnings.filterwarnings("ignore")
EDGE_ALPHA = 0.001
SIBLING_ALPHA = 0.01


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "1.0"})


def mahalanobis_time_reference(x: np.ndarray, shrinkage: float = 0.05, ridge: float = 1e-8):
    """Independent re-implementation of the documented mahalanobis_time formula (one block)."""
    centered = x - x.mean(0)
    s = centered.T @ centered / (len(x) - 1)
    tau = np.trace(s) / s.shape[0]
    cov = (1 - shrinkage) * s + (shrinkage + ridge) * tau * np.eye(s.shape[0])
    whitened = np.linalg.solve(np.linalg.cholesky(cov), centered.T).T
    return pdist(whitened, "sqeuclidean") / s.shape[0]


def check_distance(prepared, captured: np.ndarray, metric: str) -> dict:
    x = prepared.data.to_numpy(float)
    out = {"distance_metric_used": metric, "distance_finite": bool(np.isfinite(captured).all())}
    if metric == "hamming":
        out["distance_max_abs_diff_independent"] = float(
            np.abs(captured - pdist(x, "hamming")).max()
        )
    elif metric == "standardized_euclidean":
        var = x.var(0, ddof=1)
        var[~np.isfinite(var) | (var <= 0)] = 1.0
        out["distance_max_abs_diff_independent"] = float(
            np.abs(captured - pdist(x, "seuclidean", V=var)).max()
        )
    elif metric == "mahalanobis_time":
        out["distance_max_abs_diff_independent"] = float(
            np.abs(captured - mahalanobis_time_reference(x)).max()
        )
    else:
        out["distance_max_abs_diff_independent"] = np.nan
    out["distance_cv"] = float(captured.std() / captured.mean())
    out["p_over_n"] = x.shape[1] / x.shape[0]
    return out


def check_tree(tree, data: pd.DataFrame) -> dict:
    leaves = [n for n in tree.nodes if tree.out_degree(n) == 0]
    labels = {tree.nodes[n]["label"] for n in leaves}
    root = tree.root()
    sizes = sorted(int(tree.nodes[c]["leaf_count"]) for c in tree.successors(root))
    return {
        "tree_leaves_match_index": labels == set(data.index),
        "tree_binary": all(tree.out_degree(n) in (0, 2) for n in tree.nodes),
        "root_child_sizes": "|".join(map(str, sizes)),
    }


def check_labels(result, prepared, row_ari: float) -> dict:
    if result.labels is None:
        return {"labels_returned": False}
    y = np.asarray(prepared.labels)
    labels = np.asarray(result.labels)
    out = {
        "labels_returned": True,
        "labels_length_ok": len(labels) == len(y),
        "ari_recomputed": float(adjusted_rand_score(y, labels)),
    }
    out["ari_matches_row"] = bool(np.isclose(out["ari_recomputed"], row_ari, atol=1e-9))
    # Each returned label must be exactly the leaf set of its decomposition cluster root.
    tree = result.extra["tree"]
    table = build_sample_cluster_assignments(result.extra["decomposition"])
    by_label = pd.Series(labels, index=prepared.data.index)
    mismatched = 0
    for root, group in table.groupby("cluster_root"):
        leaves = {
            tree.nodes[n]["label"]
            for n in tree.nodes
            if tree.out_degree(n) == 0
            and (n == root or n in __import__("networkx").descendants(tree, root))
        }
        rows_with_label = set(
            by_label.index[by_label.isin(by_label.loc[list(group.index)].unique())]
        )
        mismatched += int(leaves != set(group.index) or rows_with_label != set(group.index))
    out["clusters_not_matching_tree"] = mismatched
    return out


def check_edges(ann: pd.DataFrame, tree, family: str) -> dict:
    tested = truthy(ann["Child_Parent_Divergence_Tested"])
    significant = truthy(ann["Child_Parent_Divergence_Significant"])
    leaf_count = pd.to_numeric(ann["leaf_count"], errors="coerce")
    is_leaf = leaf_count == 1
    out = {
        "edges_tested": int(tested.sum()),
        "edges_rejected": int((tested & significant).sum()),
        "leaf_edges_rejected_share": float(significant[tested & is_leaf].mean())
        if (tested & is_leaf).any()
        else np.nan,
        "internal_edges_rejected_share": float(significant[tested & ~is_leaf].mean())
        if (tested & ~is_leaf).any()
        else np.nan,
    }
    if family == "continuous":
        parent_of = {child: parent for parent, child in tree.edges}
        stat = pd.to_numeric(ann["Child_Parent_Divergence_Test_Statistic"], errors="coerce")
        df = pd.to_numeric(ann["Child_Parent_Divergence_df"], errors="coerce")
        n_parent = ann["node_id"].map(
            lambda n: tree.nodes[parent_of[n]]["leaf_count"] if n in parent_of else np.nan
        )
        bound = n_parent - 1
        rows = tested & stat.notna()
        root_children = list(tree.successors(tree.root()))
        root_rows = ann["node_id"].isin(root_children) & rows
        out["root_edge_statistic"] = float(stat[root_rows].iloc[0]) if root_rows.any() else np.nan
        out["root_edge_df"] = float(df[root_rows].iloc[0]) if root_rows.any() else np.nan
        out["continuous_stat_above_n_parent_minus_1"] = int((stat[rows] > bound[rows] + 1e-6).sum())
        out["continuous_stat_over_bound_median"] = float((stat[rows] / bound[rows]).median())
        impossible = chi2.sf(bound[rows], df[rows]) > EDGE_ALPHA
        out["continuous_edges_unable_to_reject"] = int(impossible.sum())
        out["continuous_edges_tested"] = int(rows.sum())
    return out


def check_sibling(ann: pd.DataFrame) -> dict:
    skipped = (
        truthy(ann["Sibling_Divergence_Skipped"])
        if "Sibling_Divergence_Skipped" in ann
        else pd.Series(False, index=ann.index)
    )
    supported = (
        truthy(ann["Sibling_Role_Supported"])
        if "Sibling_Role_Supported" in ann
        else pd.Series(False, index=ann.index)
    )
    dof = pd.to_numeric(ann.get("Sibling_Degrees_of_Freedom"), errors="coerce")
    weight = pd.to_numeric(ann.get("Sibling_Null_Weight"), errors="coerce")
    internal = pd.to_numeric(ann["leaf_count"], errors="coerce") > 1
    return {
        "sibling_tests_run": int((~skipped & internal).sum()),
        "sibling_role_supported": int(supported.sum()),
        "sibling_admissible_recomputed": int((supported & (dof > 0) & (weight > 0)).sum()),
        "sibling_method_values": ";".join(
            sorted(
                ann.get("Sibling_Test_Method", pd.Series(dtype=str)).dropna().astype(str).unique()
            )[:3]
        ),
    }


def main(annotation_dir: Path | None) -> None:
    here = Path(__file__).parent
    names = pd.read_csv(here.parent / "failing_case_subspace_summary.csv")["case_id"].tolist()
    cases = {c["name"]: c for c in get_test_cases_by_suite("full")}
    _, param_sets = resolve_selected_methods_and_param_sets(
        methods=["tbs"],
        method_params=None,
        default_methods=DEFAULT_METHODS,
        method_specs=METHOD_SPECS,
    )
    captured: dict = {}
    original = method_execution.run_clustering_result

    def capture(**kwargs):
        result = original(**kwargs)
        captured.update(kwargs=kwargs, result=result)
        return result

    method_execution.run_clustering_result = capture
    if annotation_dir is not None:
        annotation_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in names:
        tc = cases[name]
        prepared = prepare_case_inputs(tc, ["tbs"])
        params = param_sets["tbs"][0]
        row, _ = method_execution.run_single_method_once(
            method_id="tbs",
            spec=METHOD_SPECS["tbs"],
            params=params,
            case_idx=0,
            case_name=name,
            tc_seed=tc["seed"],
            significance_level=SIBLING_ALPHA,
            edge_alpha=EDGE_ALPHA,
            data_t=prepared.data,
            y_t=prepared.labels,
            x_original=prepared.original_features,
            meta=prepared.metadata,
            distance_matrix=prepared.distance_matrix,
            distance_condensed=prepared.distance_condensed,
        )
        result = captured["result"]
        distance = np.asarray(captured["kwargs"]["distance_condensed"], float)
        metric = str(row.params_raw.get("tree_distance_metric"))
        y = np.asarray(prepared.labels)
        family = (
            prepared.metadata.get("feature_space").family_label
            if prepared.metadata.get("feature_space") is not None
            else "bernoulli"
        )
        rec = {
            "case_id": name,
            "status": str(row.status.value),
            "found": row.found_clusters,
            "ari": row.ari,
            "family": family,
            "n": len(prepared.data),
            "p": prepared.data.shape[1],
            "labels_match_meta_k": int(np.unique(y).size) == int(prepared.metadata["n_clusters"]),
            "data_finite": bool(np.isfinite(prepared.data.to_numpy(float)).all()),
            "index_unique": bool(prepared.data.index.is_unique),
        }
        rec.update(check_distance(prepared, distance, metric))
        tree = result.extra["tree"]
        rec.update(check_tree(tree, prepared.data))
        rec.update(check_labels(result, prepared, row.ari))
        ann = result.extra["annotations"].reset_index()
        if "node_id" not in ann:
            ann = ann.rename(columns={ann.columns[0]: "node_id"})
        if "leaf_count" not in ann:
            ann["leaf_count"] = ann["node_id"].map(lambda n: tree.nodes[n]["leaf_count"])
        if annotation_dir is not None:
            ann.to_csv(annotation_dir / f"{name}.csv", index=False)
        rec.update(check_edges(ann, tree, family))
        rec.update(check_sibling(ann))
        if result.unsupported_reason is not None:
            ev = result.unsupported_reason.evidence
            rec["unsupported_admissible_reported"] = ev.admissible_support_count
        rows.append(rec)
        print(name, rec["status"], rec["found"], flush=True)
    pd.DataFrame(rows).to_csv(here / "bug_checks.csv", index=False)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else None)
