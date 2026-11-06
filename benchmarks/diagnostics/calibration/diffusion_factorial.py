"""Paired topology × NNLS branch time × contrast-energy rank diagnostic.

Run with ``uv run --no-sync python -m benchmarks.diagnostics.calibration.diffusion_factorial``.
No production presets or thresholds are changed. Dataset, not node, is the replicate.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata
import itertools
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import binomtest, t
from sklearn.metrics import adjusted_rand_score

from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.generators.case_data_contracts import one_hot_encode_categorical
from benchmarks.shared.runners.tbs_support import unsupported_empirical_null_reason
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.space_separation.diffusion import (
    feature_family_diffusion_geometry,
    feature_family_diffusion_unsupported_reason,
)
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.continuous_distance import standardized_euclidean_distance_condensed
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.optimized_branch_lengths import (
    fit_fixed_topology_nnls_branch_lengths,
)
from tree_break_selection.tree.topology import tree_root

SEED = 2026110203
SIGNAL_CASES = (
    "gauss_clear_medium_continuous",
    "gauss_moderate_3c_continuous",
    "binary_low_noise_4c",
    "binary_moderate_4c",
    "cat_clear_3cat_4c",
    "cat_clear_4cat_4c",
)
FACTORS = ("diffusion", "nnls", "energy")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def datasets(null_replicates, signal_replicates):
    for family_id, family in enumerate(("continuous", "bernoulli", "categorical")):
        for replicate in range(null_replicates):
            seed = SEED + family_id * 10000 + replicate
            rng = np.random.default_rng(seed)
            names = [f"S{i}" for i in range(80)]
            space = None
            if family == "categorical":
                data, _, space = one_hot_encode_categorical(
                    rng.integers(3, size=(80, 24)), 3, names
                )
            else:
                values = (
                    rng.normal(size=(80, 24))
                    if family == "continuous"
                    else rng.binomial(1, 0.5, (80, 24))
                )
                data = pd.DataFrame(values, index=names, columns=[f"F{i}" for i in range(24)])
                if family == "continuous":
                    space = continuous_feature_space_from_columns(data.columns)
            distances = (
                standardized_euclidean_distance_condensed(data.values, space)
                if family == "continuous"
                else pdist(data.values, "hamming")
            )
            yield f"iid_{family}", replicate, seed, True, data, np.zeros(80), space, distances
    cases = {case["name"]: case for case in get_default_test_cases()}
    for name in SIGNAL_CASES:
        for replicate in range(signal_replicates):
            case = {**cases[name], "seed": cases[name]["seed"] + replicate}
            inputs = prepare_case_inputs(case, ["tbs"])
            distances = inputs.metadata["precomputed_distance_condensed"]
            if distances is None:
                distances = pdist(inputs.data.values, "hamming")
            yield (
                name,
                replicate,
                case["seed"],
                False,
                inputs.data,
                inputs.labels,
                inputs.metadata.get("feature_space"),
                distances,
            )


def evaluate(tree, data, labels, space, *, nnls, energy):
    tree = copy.deepcopy(tree)
    edges_before = tuple(tree.edges())
    lengths_before = [tree.edges[e]["branch_length"] for e in edges_before]
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
        edge_branch_length_variance_policy="normalized_branch_length" if nnls else "none",
        sibling_gate_method="projected_wald_inflation",
        enforce_internal_support_thresholds=False,
    )
    annotations = bundle.annotated_df
    reason = unsupported_empirical_null_reason(
        annotations, sibling_gate_method="projected_wald_inflation"
    )
    root = tree_root(tree)
    children = list(tree.successors(root))
    mean_length = compute_mean_branch_length(tree)
    root_branch_factor = (
        1.0
        if mean_length is None
        else 1.0
        + sum(tree.edges[root, child]["branch_length"] for child in children) / (2 * mean_length)
    )
    result = {
        "status": "unsupported" if reason else "ok",
        "reason": str(reason) if reason else "",
        "root_edge_open": bool(
            annotations.loc[children, "Child_Parent_Divergence_Significant"].any()
        ),
        "root_edge_min_p": float(
            annotations.loc[children, "Child_Parent_Divergence_P_Value"].min()
        ),
        "root_edge_mean_df": float(annotations.loc[children, "Child_Parent_Divergence_df"].mean()),
        "root_sibling_p": annotations.loc[root, "Sibling_Divergence_P_Value"],
        "root_sibling_df": annotations.loc[root, "Sibling_Degrees_of_Freedom"],
        "root_sibling_rank": annotations.loc[root, "Sibling_Projection_Dimension"],
        "root_sibling_reject": bool(annotations.loc[root, "Sibling_BH_Different"]),
        "root_branch_factor": root_branch_factor,
        "edge_reject_count": int(annotations["Child_Parent_Divergence_Significant"].sum()),
        "mean_edge_df": float(annotations["Child_Parent_Divergence_df"].mean()),
        "mean_sibling_rank": float(annotations["Sibling_Projection_Dimension"].mean()),
        "invalid_sibling_count": int(annotations["Sibling_Divergence_Invalid"].sum()),
    }
    if reason is None:
        decomposer = TreeDecomposition(
            tree=tree, gate_annotation_bundle=bundle, passthrough=True, trace_level="compact"
        )
        partition, _ = labels_and_report_from_decomposition(
            decomposer.decompose_tree(), data.index.tolist()
        )
        result.update(
            clusters=len(np.unique(partition)),
            ari=adjusted_rand_score(labels, partition),
            labels=json.dumps(np.asarray(partition).tolist()),
        )
    else:
        result.update(clusters=np.nan, ari=np.nan, labels="")
    assert tuple(tree.edges()) == edges_before
    assert [tree.edges[e]["branch_length"] for e in edges_before] == lengths_before
    if not nnls:
        assert mean_length is None and root_branch_factor == 1.0
    return result, annotations.loc[[root, *children]].copy()


def summarize(rows, output):
    frame = pd.DataFrame(rows)
    summaries = []
    for keys, group in frame.groupby(["case", *FACTORS]):
        valid = group[group.status == "ok"]
        n, splits = len(valid), int((valid.clusters > 1).sum())
        interval = binomtest(splits, n).proportion_ci() if n else None
        summaries.append(
            dict(zip(["case", *FACTORS], keys))
            | {
                "n_total": len(group),
                "n_supported": n,
                "split_count": splits,
                "split_fraction": splits / n if n else np.nan,
                "split_ci_low": interval.low if interval else np.nan,
                "split_ci_high": interval.high if interval else np.nan,
                "mean_ari": valid.ari.mean(),
                "mean_clusters": valid.clusters.mean(),
                "root_edge_open_count": group.root_edge_open.sum(),
                "mean_root_sibling_rank": group.root_sibling_rank.mean(),
            }
        )
    pd.DataFrame(summaries).to_csv(output / "summary.csv", index=False)
    contrasts = []
    for factor in FACTORS:
        controls = [f for f in FACTORS if f != factor]
        for keys, group in frame.groupby(["case", *controls]):
            paired = group[group[factor] == 0].merge(
                group[group[factor] == 1],
                on="replicate",
                suffixes=("_off", "_on"),
                validate="one_to_one",
            )
            valid = paired[(paired.status_off == "ok") & (paired.status_on == "ok")]
            for metric in (
                "ari",
                "clusters",
                "root_edge_open",
                "root_sibling_rank",
                "mean_edge_df",
            ):
                used = (
                    valid
                    if metric in ("ari", "clusters")
                    else paired.dropna(subset=[metric + "_on", metric + "_off"])
                )
                differences = used[metric + "_on"].astype(float) - used[metric + "_off"].astype(
                    float
                )
                n = len(differences)
                mean = differences.mean()
                half = (
                    t.ppf(0.975, n - 1) * differences.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
                )
                contrasts.append(
                    dict(zip(["case", *controls], keys))
                    | {
                        "factor": factor,
                        "metric": metric,
                        "n_pairs": n,
                        "mean_on_minus_off": mean,
                        "ci_low": mean - half,
                        "ci_high": mean + half,
                        "increased": int((differences > 1e-12).sum()),
                        "decreased": int((differences < -1e-12).sum()),
                    }
                )
    pd.DataFrame(contrasts).to_csv(output / "paired_contrasts.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--null-replicates", type=int, default=50)
    parser.add_argument("--signal-replicates", type=int, default=5)
    parser.add_argument(
        "--output", type=Path, default=Path("raw/assets/diffusion-factorial-20261102")
    )
    args = parser.parse_args()
    if min(args.null_replicates, args.signal_replicates) < 0 or not (
        args.null_replicates + args.signal_replicates
    ):
        parser.error("replicate counts must be nonnegative and at least one must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    source_paths = sorted(
        set(Path("tree_break_selection").rglob("*.py"))
        | set(Path("benchmarks/shared").rglob("*.py"))
        | {Path(__file__), Path("uv.lock")}
    )
    manifest = {
        "seed": SEED,
        "null_replicates": args.null_replicates,
        "signal_replicates": args.signal_replicates,
        "signal_cases": SIGNAL_CASES,
        "python": platform.python_version(),
        "versions": {
            p: importlib.metadata.version(p)
            for p in ("numpy", "scipy", "pandas", "scikit-learn", "pydiffmap", "networkx")
        },
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "source_hashes": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "edge_alpha": 0.001,
        "sibling_alpha": 0.01,
        "energy_fraction": 0.9,
        "continuous_covariance_policy": "guarded_within_child",
        "minimum_child_size": 8,
        "branch_off": "all branch lengths zero; edge policy none; sibling factor one",
        "nnls": "all pairs; original feature mean squared standardized Euclidean; tol1e-5; maxiter1000; no nonconverged fits",
        "inference": "dataset-paired descriptive effects; unadjusted t intervals; no selected-null calibration claim",
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows, root_rows = [], []
    for case, replicate, seed, is_null, data, labels, space, raw_distances in datasets(
        args.null_replicates, args.signal_replicates
    ):
        base = dict(
            case=case,
            replicate=replicate,
            seed=seed,
            is_null=is_null,
            n=len(data),
            p=data.shape[1],
            data_hash=digest(data.values.tolist()),
        )
        reason = feature_family_diffusion_unsupported_reason(data, space)
        distances = {0: raw_distances}
        if reason is None:
            geometry = feature_family_diffusion_geometry(data, space)
            distances[1] = geometry.distance_condensed
        block = []
        for diffusion in (0, 1):
            if diffusion not in distances:
                for nnls, energy in itertools.product((0, 1), repeat=2):
                    block.append(
                        base
                        | dict(
                            diffusion=diffusion,
                            nnls=nnls,
                            energy=energy,
                            status="unsupported",
                            reason=reason,
                        )
                    )
                continue
            build = build_tree(data, distances[diffusion], linkage_method="average")
            topology_hash = digest(build.linkage_matrix[:, :2].tolist())
            for nnls in (0, 1):
                tree = copy.deepcopy(build.tree)
                fit = None
                if nnls:
                    fit = fit_fixed_topology_nnls_branch_lengths(
                        tree,
                        data,
                        pair_sample_size=50000,
                        random_state=0,
                        solver_tolerance=1e-5,
                        max_iterations=1000,
                    )
                    if not fit.applied_to_tree:
                        raise RuntimeError(
                            f"NNLS did not converge for {case}/{replicate}/{diffusion}: {fit}"
                        )
                else:
                    for edge in tree.edges:
                        tree.edges[edge]["branch_length"] = 0.0
                assert tuple(tree.edges()) == tuple(build.tree.edges())
                branch_hash = digest(
                    [(str(a), str(b), d["branch_length"]) for a, b, d in tree.edges(data=True)]
                )
                for energy in (0, 1):
                    treatment = base | dict(
                        diffusion=diffusion,
                        nnls=nnls,
                        energy=energy,
                        topology_hash=topology_hash,
                        branch_hash=branch_hash,
                        nnls_status=fit.status if fit else "off",
                        nnls_relative_rmse=fit.residual_rmse_to_target_mean if fit else np.nan,
                    )
                    result, root_annotations = evaluate(
                        tree, data, labels, space, nnls=nnls, energy=energy
                    )
                    block.append(treatment | result)
                    root_annotations["node"] = root_annotations.index
                    for key in ("case", "replicate", *FACTORS):
                        root_annotations[key] = treatment[key]
                    root_rows.append(root_annotations)
        assert len(block) == 8
        supported_geometry = pd.DataFrame(block).dropna(subset=["topology_hash"])
        assert supported_geometry.groupby("diffusion").topology_hash.nunique().eq(1).all()
        assert supported_geometry.groupby(["diffusion", "nnls"]).branch_hash.nunique().eq(1).all()
        rows.extend(block)
        pd.DataFrame(rows).to_csv(args.output / "runs.csv", index=False)
        print(f"{case} replicate={replicate} complete ({len(rows)} runs)", flush=True)
    pd.concat(root_rows, ignore_index=True).to_csv(
        args.output / "root_annotations.csv.gz", index=False
    )
    summarize(rows, args.output)
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
