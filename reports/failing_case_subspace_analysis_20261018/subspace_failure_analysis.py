"""Split-separation and subspace analysis of the benchmark cases where TBS fails.

For each failing case this measures, on the exact distance TBS used:
  * equidistance: how concentrated the pairwise distances are and how often
    nearest neighbours tie;
  * subspaces: which features separate each pair of true clusters, and how much
    the separating feature sets of different pairs overlap;
  * topology: whether the TBS tree holds the true clusters (oracle cut at the true
    K), and whether the tree built on the separating subspace or with another
    metric holds them;
  * splits: which features drive the top TBS tree splits, and whether they are
    true separating features.
"""

from __future__ import annotations

import re
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.metrics import adjusted_rand_score

from benchmarks.diagnostics.failure.debug_trace import analyze_single_case
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from tree_break_selection.tree.branch_lengths import node_id
from tree_break_selection.tree.construction.build import build_tree

warnings.filterwarnings("ignore")

RUN_DIR = Path("benchmarks/results/run_20261016_163716Z_full")
Q_SEPARATING = 0.01  # BH level for a feature to separate a pair of true clusters
TOP_FEATURES = 10
BINARY_METRICS = ("hamming", "jaccard", "cosine")
CONTINUOUS_METRICS = ("euclidean", "seuclidean", "cosine", "correlation", "cityblock")


def tbs_distance(prepared) -> tuple[np.ndarray, str]:
    """Return the condensed distance TBS used for this case, as method_execution chooses it."""
    meta = prepared.metadata
    if meta["requires_precomputed_tbs_distance"]:
        return np.asarray(meta["precomputed_distance_condensed"], float), str(
            meta.get("distance_metric", "precomputed")
        )
    return pdist(prepared.data.to_numpy(float), metric="hamming"), "hamming"


def linkage_rows(data: pd.DataFrame, distance: np.ndarray, method: str = "average"):
    """Build the tree the TBS way and return (linkage, row index of each linkage leaf)."""
    built = build_tree(data, distance, linkage_method=method)
    n = len(data)
    leaf_labels = [built.tree.nodes[node_id(i, n)]["label"] for i in range(n)]
    position = {label: row for row, label in enumerate(data.index)}
    return built.linkage_matrix, np.array([position[label] for label in leaf_labels])


def oracle_cut_ari(Z: np.ndarray, leaf_rows: np.ndarray, labels: np.ndarray, k: int) -> float:
    clusters = fcluster(Z, t=int(k), criterion="maxclust")
    row_clusters = np.empty(len(labels), dtype=int)
    row_clusters[leaf_rows] = clusters
    return float(adjusted_rand_score(labels, row_clusters))


def clades(Z: np.ndarray, leaf_rows: np.ndarray) -> list[np.ndarray]:
    """Row-index member sets of every internal node, in merge order (last = root)."""
    n = Z.shape[0] + 1
    members: list[np.ndarray] = [np.array([leaf_rows[i]]) for i in range(n)]
    for a, b, _, _ in Z:
        members.append(np.concatenate([members[int(a)], members[int(b)]]))
    return members[n:]


def clade_recall(Z, leaf_rows, labels) -> float:
    """Fraction of true clusters matched by some clade with Jaccard >= 0.9."""
    sets = [set(m.tolist()) for m in clades(Z, leaf_rows)]
    hits = 0
    for label in np.unique(labels):
        truth = set(np.flatnonzero(labels == label).tolist())
        if len(truth) == 1:
            hits += 1
            continue
        best = max(len(truth & s) / len(truth | s) for s in sets)
        hits += best >= 0.9
    return hits / len(np.unique(labels))


def welch_p(x: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-feature Welch t-test p-values between row sets a and b (constant columns get 1)."""
    if len(a) < 2 or len(b) < 2:
        return np.ones(x.shape[1])
    p = ttest_ind(x[a], x[b], axis=0, equal_var=False).pvalue
    return np.where(np.isfinite(p), p, 1.0)


def welch_t(x: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if len(a) < 2 or len(b) < 2:
        return np.abs(x[a].mean(0) - x[b].mean(0))
    t = ttest_ind(x[a], x[b], axis=0, equal_var=False).statistic
    return np.nan_to_num(np.abs(t))


def separating_sets(x: np.ndarray, labels: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    """Boolean feature mask per true-cluster pair: BH q < Q_SEPARATING."""
    groups = {g: np.flatnonzero(labels == g) for g in np.unique(labels)}
    out = {}
    for g, h in combinations(groups, 2):
        p = welch_p(x, groups[g], groups[h])
        out[(g, h)] = false_discovery_control(p, method="bh") < Q_SEPARATING
    return out


def equidistance(distance: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    square = squareform(distance)
    np.fill_diagonal(square, np.inf)
    nearest = square.min(axis=1)
    tied = np.isclose(square, nearest[:, None], rtol=0, atol=1e-12)
    same = labels[:, None] == labels[None, :]
    iu = np.triu_indices(len(labels), 1)
    within = distance[same[iu]]
    between = distance[~same[iu]]
    return {
        "distance_cv": float(distance.std() / distance.mean()) if distance.mean() > 0 else np.nan,
        "distinct_distance_fraction": float(np.unique(distance.round(12)).size / distance.size),
        "nn_tie_fraction": float((tied.sum(axis=1) > 1).mean()),
        "nn_same_label": float(((tied & same).sum(axis=1) / tied.sum(axis=1)).mean()),
        "between_over_within": float(between.mean() / within.mean())
        if within.size and between.size and within.mean() > 0
        else np.nan,
    }


def metric_trees(data, x, labels, k, binary: bool) -> dict[str, float]:
    out = {}
    for metric in BINARY_METRICS if binary else CONTINUOUS_METRICS:
        d = pdist(x, metric=metric)
        if not np.isfinite(d).all():
            continue
        Z, rows = linkage_rows(data, d)
        out[f"{metric}_average"] = oracle_cut_ari(Z, rows, labels, k)
    # Ward on Euclidean as one linkage alternative
    Z = linkage(x, method="ward")
    clusters = fcluster(Z, t=int(k), criterion="maxclust")
    out["euclidean_ward"] = float(adjusted_rand_score(labels, clusters))
    return out


def subspace_distance(x: np.ndarray, binary: bool) -> np.ndarray:
    if binary:
        return pdist(x, metric="hamming")
    variance = x.var(axis=0, ddof=1)
    variance[variance <= 0] = 1.0
    return pdist(x, metric="seuclidean", V=variance)


def describe_side(rows: np.ndarray, labels: np.ndarray) -> str:
    values, counts = np.unique(labels[rows], return_counts=True)
    order = np.argsort(-counts)
    return "+".join(f"{values[i]}:{counts[i]}" for i in order[:4]) + (
        "+..." if len(values) > 4 else ""
    )


def top_splits(Z, leaf_rows, x, labels, union_mask, pair_masks, columns, limit=4):
    """Root and the next largest internal splits: sizes, label mix and driving features."""
    n = Z.shape[0] + 1
    members = [np.array([leaf_rows[i]]) for i in range(n)]
    for a, b, _, _ in Z:
        members.append(np.concatenate([members[int(a)], members[int(b)]]))
    queue, rows = [2 * n - 2], []
    while queue and len(rows) < limit:
        queue.sort(key=lambda node: -len(members[node]))
        node = queue.pop(0)
        if node < n:
            continue
        left, right = (int(c) for c in Z[node - n, :2])
        a, b = members[left], members[right]
        score = welch_t(x, a, b)
        top = np.argsort(-score)[:TOP_FEATURES]
        pairs = [f"{g}-{h}" for (g, h), m in pair_masks.items() if m[top].any()]
        rows.append(
            {
                "split_rank": len(rows),
                "sizes": f"{len(a)}|{len(b)}",
                "min_side_fraction": min(len(a), len(b)) / (len(a) + len(b)),
                "left_labels": describe_side(a, labels),
                "right_labels": describe_side(b, labels),
                "top_features": ",".join(str(columns[i]) for i in top),
                "top_in_separating": float(union_mask[top].mean()) if union_mask.any() else 0.0,
                "true_pairs_touched": len(pairs),
            }
        )
        queue += [left, right]
    return rows


def audit_root(case_number: int) -> dict[str, object]:
    path = RUN_DIR / "audit" / f"case_{case_number}_tbs_stats.csv"
    if not path.exists():
        return {"root_diagnosis": "no audit"}
    audit = pd.read_csv(path, low_memory=False)
    diagnosis = analyze_single_case(path)

    def flag(column: str) -> int:
        if column not in audit:
            return 0
        return int(audit[column].astype(str).str.lower().eq("true").sum())

    return {
        "root_diagnosis": diagnosis["reason"],
        "edge_significant_nodes": flag("Child_Parent_Divergence_Significant"),
        "edge_tested_nodes": flag("Child_Parent_Divergence_Tested"),
        "sibling_role_supported": flag("Sibling_Role_Supported"),
        "sibling_tests_run": int(
            (~audit["Sibling_Divergence_Skipped"].astype(str).str.lower().eq("true")).sum()
        )
        if "Sibling_Divergence_Skipped" in audit
        else 0,
    }


def variable_of(column: str) -> str:
    match = re.fullmatch(r"(F\d+)_c\d+", str(column))
    return match.group(1) if match else str(column)


def analyse(case: dict, outcome: pd.Series) -> tuple[dict, list[dict]]:
    prepared = prepare_case_inputs(case, ["tbs"])
    data, labels = prepared.data, np.asarray(prepared.labels)
    x = data.to_numpy(float)
    k = int(prepared.metadata["n_clusters"])
    representation = str(outcome["feature_representation"])
    binary = representation != "continuous"
    distance, metric = tbs_distance(prepared)
    row = {
        "case_id": case["name"],
        "representation": representation,
        "tbs_metric": metric,
        "n": len(data),
        "p": x.shape[1],
        "k": k,
        "min_cluster": int(np.bincount(pd.factorize(labels)[0]).min()),
        "status_now": outcome["status_now"],
        "found_now": outcome["found_now"],
        "ari_now": outcome["ari_now"],
    }
    row.update(equidistance(distance, labels))
    Z, leaf_rows = linkage_rows(data, distance)
    row["oracle_ari_tbs_tree"] = oracle_cut_ari(Z, leaf_rows, labels, k) if k > 1 else np.nan
    row["clade_recall_tbs_tree"] = clade_recall(Z, leaf_rows, labels) if k > 1 else np.nan
    root_a, root_b = (int(c) for c in Z[-1, :2])
    sizes = [1 if c < len(data) else int(Z[c - len(data), 3]) for c in (root_a, root_b)]
    row["root_min_side"] = min(sizes)

    split_rows: list[dict] = []
    if k > 1:
        pair_masks = separating_sets(x, labels)
        stacked = np.array(list(pair_masks.values()))
        union = stacked.any(axis=0)
        per_pair = stacked.sum(axis=1)
        row["pairs"] = len(pair_masks)
        row["pairs_without_separating_feature"] = int((per_pair == 0).sum())
        row["separating_fraction_mean"] = float(per_pair.mean() / x.shape[1])
        row["separating_union_fraction"] = float(union.mean())
        variables = np.array([variable_of(c) for c in data.columns])
        row["separating_union_variables"] = int(np.unique(variables[union]).size)
        row["variables"] = int(np.unique(variables).size)
        nonempty = stacked[per_pair > 0].astype(float)
        if len(nonempty) > 1:
            inter = nonempty @ nonempty.T
            sizes_ = nonempty.sum(axis=1)
            jac = inter / (sizes_[:, None] + sizes_[None, :] - inter)
            row["pair_subspace_overlap"] = float(jac[np.triu_indices(len(jac), 1)].mean())
        else:
            row["pair_subspace_overlap"] = np.nan
        if union.any():
            sub = subspace_distance(x[:, union], binary)
            Zs, rows_s = linkage_rows(data, sub)
            row["oracle_ari_separating_subspace"] = oracle_cut_ari(Zs, rows_s, labels, k)
            row.update(
                {
                    f"subspace_{key}": value
                    for key, value in equidistance(sub, labels).items()
                    if key in ("distance_cv", "between_over_within")
                }
            )
        else:
            row["oracle_ari_separating_subspace"] = np.nan
        alternatives = metric_trees(data, x, labels, k, binary)
        best = max(alternatives, key=alternatives.get)
        row["best_other_tree"] = best
        row["oracle_ari_best_other_tree"] = alternatives[best]
        for split in top_splits(Z, leaf_rows, x, labels, union, pair_masks, list(data.columns)):
            split_rows.append({"case_id": case["name"], **split})
    row.update(audit_root(int(outcome["test_case"])))
    return row, split_rows


def classify(r: pd.Series) -> str:
    if r["k"] <= 1:
        return "null case (K=1): unsupported instead of one cluster"
    tree = r["oracle_ari_tbs_tree"]
    sub = r.get("oracle_ari_separating_subspace", np.nan)
    other = r.get("oracle_ari_best_other_tree", np.nan)
    if tree >= 0.9:
        return "tree holds the clusters; the gate does not split"
    if np.nanmax([sub, other]) >= max(0.9, tree + 0.25):
        return "tree misses the clusters; separating subspace or another metric recovers them"
    if np.nanmax([sub, other, tree]) >= 0.5:
        return "partly recoverable: every tree mixes some clusters"
    return "no tree recovers the clusters (overlap or weak signal)"


def main(out_dir: Path, rerun_csv: Path) -> None:
    old = pd.read_csv(RUN_DIR / "full_benchmark_comparison.csv")
    old = old[old.method == "tbs"].set_index("case_id")
    now = pd.read_csv(rerun_csv).set_index("case_id")
    failing = now.index
    outcomes = old.loc[failing, ["test_case", "feature_representation"]].assign(
        status_now=now["status"], found_now=now["found_clusters"], ari_now=now["ari"]
    )
    cases = {c["name"]: c for c in get_test_cases_by_suite("full")}
    rows, splits = [], []
    for name in failing:
        print("analysing", name, flush=True)
        row, split_rows = analyse(cases[name], outcomes.loc[name])
        rows.append(row)
        splits += split_rows
    table = pd.DataFrame(rows)
    table["diagnosis"] = table.apply(classify, axis=1)
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / "failing_case_subspace_summary.csv", index=False)
    pd.DataFrame(splits).to_csv(out_dir / "failing_case_top_splits.csv", index=False)
    print(table["diagnosis"].value_counts().to_string())


if __name__ == "__main__":
    here = Path(__file__).parent
    main(here, here / "tbs_rerun_failing_cases.csv")
