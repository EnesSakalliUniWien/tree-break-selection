"""Known-branch Gaussian diagnostic for locating sibling-test skew.

The fixed tree has three named internal nodes (R, A, B) and four observed
leaf groups (A0, A1, B0, B1). The diagnostic crosses known versus sample
pooled covariance with a fixed full basis versus same-data PCA/MP projection,
then evaluates the selected root separately and sends exact known-covariance
local p-values through production sibling BH and ``TreeDecomposition``.
Nothing in this module changes production defaults or calibration methods.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shlex
import subprocess
import sys
from dataclasses import asdict, dataclass
from fractions import Fraction
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.diagnostics.calibration.statistics.statistic_distribution_shape_panel import (
    expected_chi_square_mixture_skew,
)
from benchmarks.validation.contracts.report_contract import wilson_interval
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.contrast_covariance import (
    build_null_whitened_tangent_matrix,
    compute_whitened_wald_contrast,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projection_dimension_estimation.projection_dimension_estimators import (
    estimate_marchenko_pastur_dimension,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflated_projected_wald_annotation.fdr_annotation import (
    apply_traversal_aligned_sibling_bh_results,
    init_sibling_annotation_df,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns
from tree_break_selection.tree.poset_tree import PosetTree

SCHEMA_VERSION = "known_branch_skew_diagnostic/v1"
GENERATED_BY = "benchmarks.diagnostics.calibration.sibling.known_branch_skew_diagnostic"
PATTERNS = (
    "global_null",
    "root_only",
    "A_deep_only",
    "separate_deep_branches",
    "strong_root_weak_A",
)
PROFILES = ("identity", "ar1_0.6")
INTERNAL_CHILDREN = {"R": ("A", "B"), "A": ("A0", "A1"), "B": ("B0", "B1")}
LEAF_NAMES = ("A0", "A1", "B0", "B1")
VARIANTS = (
    ("oracle_covariance", "fixed_full"),
    ("pooled_sample_covariance", "fixed_full"),
    ("oracle_covariance", "same_data_pca_mp"),
    ("pooled_sample_covariance", "same_data_pca_mp"),
)
NODE_POPULATION_CHILD_LEAVES = {
    "R": (("A0", "A1"), ("B0", "B1")),
    "A": (("A0",), ("A1",)),
    "B": (("B0",), ("B1",)),
}


@dataclass(frozen=True)
class KnownBranchConfig:
    output_dir: Path
    replicates: int = 500
    base_seed: int = 20260928
    n_per_group: int = 32
    dimension: int = 8
    alpha: float = 0.01
    ridge: float = 0.0
    patterns: tuple[str, ...] = PATTERNS
    covariance_profiles: tuple[str, ...] = PROFILES

    def __post_init__(self) -> None:
        if self.replicates < 1 or self.n_per_group < 4 or self.dimension < 2:
            raise ValueError(
                "replicates, n_per_group, or dimension is below its supported minimum."
            )
        if not 0.0 < self.alpha < 1.0:
            raise ValueError("alpha must lie in (0, 1).")
        if not math.isfinite(self.ridge) or self.ridge < 0.0:
            raise ValueError("ridge must be finite and non-negative.")
        if not self.patterns or set(self.patterns) - set(PATTERNS):
            raise ValueError(f"patterns must be selected from {PATTERNS!r}.")
        if not self.covariance_profiles or set(self.covariance_profiles) - set(PROFILES):
            raise ValueError(f"covariance_profiles must be selected from {PROFILES!r}.")


def covariance_matrix(profile: str, dimension: int) -> np.ndarray:
    """Return the declared Gaussian population covariance."""
    if profile == "identity":
        return np.eye(dimension)
    if profile == "ar1_0.6":
        indexes = np.arange(dimension)
        return 0.6 ** np.abs(indexes[:, None] - indexes[None, :])
    raise ValueError(f"Unknown covariance profile: {profile!r}.")


def population_means(pattern: str, dimension: int) -> dict[str, np.ndarray]:
    """Return exact leaf-group population means for a named branch pattern."""
    zero = np.zeros(dimension)
    mu = np.zeros(dimension)
    nu = np.zeros(dimension)
    mu[0] = 1.25
    nu[1] = 1.0
    if pattern == "global_null":
        return {leaf: zero.copy() for leaf in LEAF_NAMES}
    if pattern == "root_only":
        return {"A0": mu.copy(), "A1": mu.copy(), "B0": -mu.copy(), "B1": -mu.copy()}
    if pattern == "A_deep_only":
        return {"A0": mu.copy(), "A1": -mu.copy(), "B0": zero.copy(), "B1": zero.copy()}
    if pattern == "separate_deep_branches":
        return {"A0": mu.copy(), "A1": -mu.copy(), "B0": nu.copy(), "B1": -nu.copy()}
    if pattern == "strong_root_weak_A":
        # A0/A1 average to +mu; their within-A contrast is 0.5*nu.
        return {
            "A0": mu + 0.25 * nu,
            "A1": mu - 0.25 * nu,
            "B0": -mu.copy(),
            "B1": -mu.copy(),
        }
    raise ValueError(f"Unknown branch pattern: {pattern!r}.")


def expected_node_contrasts(
    means: dict[str, np.ndarray],
) -> dict[str, tuple[np.ndarray, bool]]:
    """Compute named population sibling contrasts and exact equality labels."""
    result: dict[str, tuple[np.ndarray, bool]] = {}
    for node, (left, right) in INTERNAL_CHILDREN.items():
        left_leaves, right_leaves = NODE_POPULATION_CHILD_LEAVES[node]
        left_mean, left_exact = _population_mean(left_leaves, means)
        right_mean, right_exact = _population_mean(right_leaves, means)
        contrast = left_mean - right_mean
        result[node] = (contrast, left_exact == right_exact)
    return result


def _population_mean(
    leaves: Sequence[str], means: dict[str, np.ndarray]
) -> tuple[np.ndarray, tuple[Fraction, ...]]:
    exact = tuple(
        sum((Fraction(str(float(means[leaf][j]))) for leaf in leaves), Fraction()) / len(leaves)
        for j in range(len(means[leaves[0]]))
    )
    return np.asarray([float(value) for value in exact]), exact


def _fixed_tree() -> PosetTree:
    tree = PosetTree()
    tree.add_edges_from(
        (("R", "A"), ("R", "B"), ("A", "A0"), ("A", "A1"), ("B", "B0"), ("B", "B1"))
    )
    tree.graph["root"] = "R"
    for leaf in LEAF_NAMES:
        tree.nodes[leaf].update(label=leaf, is_leaf=True)
    for node in INTERNAL_CHILDREN:
        tree.nodes[node].update(label=node, is_leaf=False)
    return tree


def _local_test(
    left: np.ndarray,
    right: np.ndarray,
    *,
    population_covariance: np.ndarray,
    covariance_estimator: str,
    projection_kind: str,
    ridge: float,
) -> dict[str, object]:
    pooled = np.vstack((left, right))
    left_mean, right_mean = left.mean(axis=0), right.mean(axis=0)
    parent_mean = pooled.mean(axis=0)
    covariance = (
        population_covariance
        if covariance_estimator == "oracle_covariance"
        else np.cov(pooled, rowvar=False, ddof=1)
    )
    feature_space = continuous_feature_space_from_columns(
        tuple(f"x{i}" for i in range(left.shape[1]))
    )
    covariance_by_block = {"continuous": covariance}
    z = compute_whitened_wald_contrast(
        left_mean,
        right_mean,
        float(len(left)),
        float(len(right)),
        comparison="sibling",
        feature_space=feature_space,
        continuous_covariance_by_block=covariance_by_block,
        ridge=ridge,
    )
    rank = z.size
    raw_mp_count = 0
    if projection_kind == "same_data_pca_mp":
        # All raw rows enter the same coordinates used for z; using the four
        # group means here would violate the production descendant-row contract.
        tangent_rows = build_null_whitened_tangent_matrix(
            pooled,
            parent_mean,
            feature_space=feature_space,
            continuous_covariance_by_block=covariance_by_block,
            ridge=ridge,
        )
        eigen = eigendecompose_covariance(tangent_rows, compute_eigenvectors=True)
        if eigen is None:
            return {"status": "unsupported", "failure_reason": "zero_rank_tangent"}
        dimension = estimate_marchenko_pastur_dimension(
            eigen.eigenvalues,
            n_samples=tangent_rows.shape[0],
            n_features=eigen.active_feature_count,
            effective_independent_rows=len(pooled),
            mp_threshold_rows=len(pooled),
            minimum_projection_dimension=1,
        )
        rank = int(dimension.test_projection_dimension)
        raw_mp_count = int(dimension.raw_mp_signal_count)
        projection, eigenvalues = build_pca_projection(
            eigen,
            projection_dimension=rank,
            n_features_total=z.size,
        )
        result = run_projected_wald_kernel(
            z,
            spectral_k=rank,
            pca_projection=projection,
            pca_eigenvalues=eigenvalues,
        )
        statistic, df, pvalue = result.statistic, result.degrees_of_freedom, result.p_value
    else:
        statistic = float(np.dot(z, z))
        df = rank
        pvalue = float(chi2.sf(statistic, df))
    return {
        "status": "ok",
        "failure_reason": "",
        "statistic": float(statistic),
        "degrees_of_freedom": float(df),
        "p_value": float(pvalue),
        "statistic_per_df": float(statistic / df) if df else math.nan,
        "projection_dimension": int(rank),
        "raw_mp_signal_count": int(raw_mp_count),
        "covariance_condition": float(np.linalg.cond(covariance)),
    }


def _selected_root_test(
    data: np.ndarray,
    labels: np.ndarray,
    means: dict[str, np.ndarray],
    covariance: np.ndarray,
) -> dict[str, object]:
    names = [f"sample_{index:04d}" for index in range(data.shape[0])]
    z = data.copy()
    # Average-linkage Euclidean hierarchy is selected from these same rows.
    from scipy.spatial.distance import pdist

    tree = build_tree(
        pd.DataFrame(z, index=names), pdist(z, metric="euclidean"), linkage_method="average"
    ).tree
    root = tree.root()
    children = list(tree.successors(root))
    if len(children) != 2:
        return {"status": "unsupported", "failure_reason": "selected_root_nonbinary"}
    descendants = tree.compute_descendant_sets(use_labels=True)
    name_to_index = {name: i for i, name in enumerate(names)}
    node_data = [z[[name_to_index[name] for name in descendants[child]]] for child in children]
    node_labels = [
        labels[[name_to_index[name] for name in descendants[child]]] for child in children
    ]
    child_mean_pairs = [
        _population_mean([str(label) for label in subset], means) for subset in node_labels
    ]
    child_means = [pair[0] for pair in child_mean_pairs]
    population_contrast = child_means[0] - child_means[1]
    test = _local_test(
        node_data[0],
        node_data[1],
        population_covariance=covariance,
        covariance_estimator="oracle_covariance",
        projection_kind="fixed_full",
        ridge=0.0,
    )
    return {
        **test,
        "status": str(test["status"]),
        "focal_node": "selected_root",
        "left_sample_size": int(len(node_data[0])),
        "right_sample_size": int(len(node_data[1])),
        "left_leaf_group_counts": json.dumps(
            {leaf: int(np.count_nonzero(node_labels[0] == leaf)) for leaf in LEAF_NAMES},
            sort_keys=True,
        ),
        "right_leaf_group_counts": json.dumps(
            {leaf: int(np.count_nonzero(node_labels[1] == leaf)) for leaf in LEAF_NAMES},
            sort_keys=True,
        ),
        "population_contrast_norm": float(np.linalg.norm(population_contrast)),
        "population_null": child_mean_pairs[0][1] == child_mean_pairs[1][1],
        "left_population_mean": json.dumps(child_means[0].tolist()),
        "right_population_mean": json.dumps(child_means[1].tolist()),
        "population_contrast": json.dumps(population_contrast.tolist()),
    }


def _traversal_row(
    tree: PosetTree,
    tests: dict[str, dict[str, object]],
    truths: dict[str, tuple[np.ndarray, bool]],
    *,
    alpha: float,
    passthrough: bool,
) -> dict[str, object]:
    annotations = init_sibling_annotation_df(
        pd.DataFrame(
            {"Child_Parent_Divergence_Significant": pd.Series(True, index=list(tree.nodes))}
        )
    )
    parents = list(INTERNAL_CHILDREN)
    results = [
        (
            float(tests[node]["statistic"]),
            float(tests[node]["degrees_of_freedom"]),
            float(tests[node]["p_value"]),
        )
        for node in parents
    ]
    annotations = apply_traversal_aligned_sibling_bh_results(
        tree, annotations, parents, results, alpha=alpha
    )
    result = TreeDecomposition(
        tree=tree,
        annotations_df=annotations,
        passthrough=passthrough,
        trace_level="full",
    ).decompose_tree()
    trace = result["traversal_trace"]
    internal_rows = [row for row in trace if row["node_id"] in truths]
    reached_splits = [row for row in internal_rows if row["decision"] == "split"]
    false_splits = [row for row in reached_splits if truths[row["node_id"]][1]]
    return {
        "status": "ok",
        "passthrough": bool(passthrough),
        "num_clusters": int(result["num_clusters"]),
        "n_reached_internal": len(internal_rows),
        "n_reached_splits": len(reached_splits),
        "n_reached_false_splits": len(false_splits),
        "reached_false_split_fdp": len(false_splits) / max(1, len(reached_splits)),
        "false_split_node_ids": ",".join(str(row["node_id"]) for row in false_splits),
        "injected_local_tests": json.dumps(
            {
                node: {
                    "statistic": float(tests[node]["statistic"]),
                    "degrees_of_freedom": float(tests[node]["degrees_of_freedom"]),
                    "p_value": float(tests[node]["p_value"]),
                }
                for node in parents
            },
            sort_keys=True,
        ),
        "traversal_decisions": json.dumps(
            {str(row["node_id"]): str(row["decision"]) for row in internal_rows}, sort_keys=True
        ),
    }


def run_diagnostic(config: KnownBranchConfig) -> dict[str, pd.DataFrame]:
    """Run the matched Gaussian design and return replicate-level tables."""
    local_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []
    traversal_rows: list[dict[str, object]] = []
    seed_sequence = np.random.SeedSequence(config.base_seed)
    child_seeds = iter(
        seed_sequence.spawn(
            len(config.patterns) * len(config.covariance_profiles) * config.replicates
        )
    )
    for profile in config.covariance_profiles:
        sigma = covariance_matrix(profile, config.dimension)
        tree = _fixed_tree()
        for pattern in config.patterns:
            means = population_means(pattern, config.dimension)
            truths = expected_node_contrasts(means)
            for replicate in range(config.replicates):
                dataset_seed_sequence = next(child_seeds)
                rng = np.random.default_rng(dataset_seed_sequence)
                dataset_seed = int(dataset_seed_sequence.generate_state(1, dtype=np.uint32)[0])
                labels = np.repeat(np.asarray(LEAF_NAMES), config.n_per_group)
                data = np.vstack(
                    [
                        rng.multivariate_normal(means[leaf], sigma, size=config.n_per_group)
                        for leaf in LEAF_NAMES
                    ]
                )
                dataset_hash = hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()
                group_data = {
                    leaf: data[index * config.n_per_group : (index + 1) * config.n_per_group]
                    for index, leaf in enumerate(LEAF_NAMES)
                }
                node_data = {
                    "A0": group_data["A0"],
                    "A1": group_data["A1"],
                    "B0": group_data["B0"],
                    "B1": group_data["B1"],
                    "A": np.vstack((group_data["A0"], group_data["A1"])),
                    "B": np.vstack((group_data["B0"], group_data["B1"])),
                    "R": data,
                }
                exact_tests: dict[str, dict[str, object]] = {}
                for node, (left_name, right_name) in INTERNAL_CHILDREN.items():
                    left_key, right_key = (left_name, right_name)
                    left_data, right_data = node_data[left_key], node_data[right_key]
                    population_contrast, is_null = truths[node]
                    for covariance_estimator, projection_kind in VARIANTS:
                        try:
                            result = _local_test(
                                left_data,
                                right_data,
                                population_covariance=sigma,
                                covariance_estimator=covariance_estimator,
                                projection_kind=projection_kind,
                                ridge=config.ridge,
                            )
                            status, reason = (
                                str(result["status"]),
                                str(result.get("failure_reason", "")),
                            )
                        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                            result, status, reason = {}, "runtime_error", str(exc)
                        local_rows.append(
                            {
                                "pattern": pattern,
                                "covariance_profile": profile,
                                "replicate": replicate,
                                "dataset_seed": dataset_seed,
                                "dataset_sha256": dataset_hash,
                                "seed_spawn_key": ".".join(
                                    map(str, rng.bit_generator._seed_seq.spawn_key)
                                ),
                                "focal_node": node,
                                "left_child": left_name,
                                "right_child": right_name,
                                "covariance_estimator": covariance_estimator,
                                "projection_kind": projection_kind,
                                "covariance_label": (
                                    "known_population_covariance"
                                    if covariance_estimator == "oracle_covariance"
                                    else "sample_pooled_covariance_including_between_sibling_means"
                                ),
                                "status": status,
                                "failure_reason": reason,
                                "population_null": is_null,
                                "left_population_mean": json.dumps(
                                    _population_mean(NODE_POPULATION_CHILD_LEAVES[node][0], means)[
                                        0
                                    ].tolist()
                                ),
                                "right_population_mean": json.dumps(
                                    _population_mean(NODE_POPULATION_CHILD_LEAVES[node][1], means)[
                                        0
                                    ].tolist()
                                ),
                                "population_contrast": json.dumps(population_contrast.tolist()),
                                "population_contrast_norm": float(
                                    np.linalg.norm(population_contrast)
                                ),
                                "degrees_of_freedom": result.get("degrees_of_freedom", math.nan),
                                "statistic": result.get("statistic", math.nan),
                                "statistic_per_df": result.get("statistic_per_df", math.nan),
                                "p_value": result.get("p_value", math.nan),
                                "projection_dimension": result.get(
                                    "projection_dimension", math.nan
                                ),
                                "raw_mp_signal_count": result.get("raw_mp_signal_count", math.nan),
                                "covariance_condition": result.get(
                                    "covariance_condition", math.nan
                                ),
                            }
                        )
                    exact_tests[node] = _local_test(
                        left_data,
                        right_data,
                        population_covariance=sigma,
                        covariance_estimator="oracle_covariance",
                        projection_kind="fixed_full",
                        ridge=config.ridge,
                    )
                for passthrough in (False, True):
                    try:
                        outcome = _traversal_row(
                            tree,
                            exact_tests,
                            truths,
                            alpha=config.alpha,
                            passthrough=passthrough,
                        )
                        status, reason = str(outcome["status"]), ""
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                        outcome, status, reason = {}, "runtime_error", str(exc)
                    traversal_rows.append(
                        {
                            "pattern": pattern,
                            "covariance_profile": profile,
                            "replicate": replicate,
                            "dataset_seed": dataset_seed,
                            "dataset_sha256": dataset_hash,
                            "seed_spawn_key": ".".join(
                                map(str, rng.bit_generator._seed_seq.spawn_key)
                            ),
                            "injected_p_value_source": "oracle_covariance_fixed_full_basis_exact_gaussian_null",
                            "edge_prerequisite_regime": "all_edges_deliberately_open",
                            "alpha": config.alpha,
                            "status": status,
                            "failure_reason": reason,
                            **outcome,
                        }
                    )
                try:
                    selected = _selected_root_test(data, labels, means, sigma)
                except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                    selected = {"status": "runtime_error", "failure_reason": str(exc)}
                selected_rows.append(
                    {
                        "pattern": pattern,
                        "covariance_profile": profile,
                        "replicate": replicate,
                        "dataset_seed": dataset_seed,
                        "dataset_sha256": dataset_hash,
                        "selection_method": "same_data_average_linkage_euclidean",
                        "status": selected["status"],
                        "failure_reason": selected.get("failure_reason", ""),
                        "population_null": selected.get("population_null", math.nan),
                        "left_sample_size": selected.get("left_sample_size", math.nan),
                        "right_sample_size": selected.get("right_sample_size", math.nan),
                        "left_leaf_group_counts": selected.get("left_leaf_group_counts", ""),
                        "right_leaf_group_counts": selected.get("right_leaf_group_counts", ""),
                        "left_population_mean": selected.get("left_population_mean", ""),
                        "right_population_mean": selected.get("right_population_mean", ""),
                        "population_contrast": selected.get("population_contrast", ""),
                        "population_contrast_norm": selected.get(
                            "population_contrast_norm", math.nan
                        ),
                        "degrees_of_freedom": selected.get("degrees_of_freedom", math.nan),
                        "statistic": selected.get("statistic", math.nan),
                        "statistic_per_df": selected.get("statistic_per_df", math.nan),
                        "p_value": selected.get("p_value", math.nan),
                    }
                )
    return {
        "local": pd.DataFrame.from_records(local_rows),
        "selected_root": pd.DataFrame.from_records(selected_rows),
        "traversal": pd.DataFrame.from_records(traversal_rows),
    }


def summarize_tables(tables: dict[str, pd.DataFrame], *, alpha: float) -> pd.DataFrame:
    """Summarize valid outcomes while retaining status counts."""
    rows: list[dict[str, object]] = []
    local = tables["local"]
    group_columns = [
        "pattern",
        "covariance_profile",
        "focal_node",
        "covariance_estimator",
        "projection_kind",
    ]
    for key, frame in local.groupby(group_columns, dropna=False, sort=True):
        ok = frame[frame["status"] == "ok"]
        null = ok[ok["population_null"].astype(bool)]
        signal = ok[~ok["population_null"].astype(bool)]
        null_rejections = int((null["p_value"] <= alpha).sum())
        signal_rejections = int((signal["p_value"] <= alpha).sum())
        null_ci = wilson_interval(null_rejections, len(null)) if len(null) else (math.nan, math.nan)
        signal_ci = (
            wilson_interval(signal_rejections, len(signal)) if len(signal) else (math.nan, math.nan)
        )
        rows.append(
            {
                **dict(zip(group_columns, key)),
                "table": "local",
                "n_total": len(frame),
                "n_ok": len(ok),
                "n_unavailable": len(frame) - len(ok),
                "null_n": len(null),
                "null_rejections": null_rejections,
                "null_rejection_rate": null_rejections / len(null) if len(null) else math.nan,
                "null_rejection_ci95_low": null_ci[0],
                "null_rejection_ci95_high": null_ci[1],
                "signal_n": len(signal),
                "signal_rejections": signal_rejections,
                "signal_rejection_rate": signal_rejections / len(signal)
                if len(signal)
                else math.nan,
                "signal_rejection_ci95_low": signal_ci[0],
                "signal_rejection_ci95_high": signal_ci[1],
                "mean_statistic_per_df": float(ok["statistic_per_df"].mean())
                if len(ok)
                else math.nan,
                "p_value_q05": float(ok["p_value"].quantile(0.05)) if len(ok) else math.nan,
                "p_value_median": float(ok["p_value"].median()) if len(ok) else math.nan,
                "p_value_q95": float(ok["p_value"].quantile(0.95)) if len(ok) else math.nan,
                "null_statistic_skew": float(null["statistic"].skew())
                if len(null) >= 3
                else math.nan,
                "expected_null_chi_square_mixture_skew": (
                    expected_chi_square_mixture_skew(null["degrees_of_freedom"])
                    if len(null)
                    else math.nan
                ),
                "mean_projection_dimension": float(ok["projection_dimension"].mean())
                if len(ok)
                else math.nan,
                "mean_raw_mp_signal_count": float(ok["raw_mp_signal_count"].mean())
                if len(ok)
                else math.nan,
            }
        )
    for table_name in ("selected_root", "traversal"):
        frame = tables[table_name]
        by = ["pattern", "covariance_profile"] + (
            ["passthrough"] if table_name == "traversal" else []
        )
        for key, group in frame.groupby(by, dropna=False, sort=True):
            valid = group[group["status"] == "ok"]
            row = dict(zip(by, key if isinstance(key, tuple) else (key,)))
            row.update(
                {
                    "table": table_name,
                    "n_total": len(group),
                    "n_ok": len(valid),
                    "n_unavailable": len(group) - len(valid),
                    "selected_root_null_rate": (
                        float(valid["population_null"].astype(bool).mean())
                        if table_name == "selected_root" and len(valid)
                        else math.nan
                    ),
                    "selected_root_rejection_rate": (
                        float((valid["p_value"] <= alpha).mean())
                        if table_name == "selected_root" and len(valid)
                        else math.nan
                    ),
                    "selected_root_null_n": (
                        int(valid["population_null"].astype(bool).sum())
                        if table_name == "selected_root" and len(valid)
                        else 0
                    ),
                    "selected_root_null_rejections": (
                        int(
                            (
                                valid.loc[valid["population_null"].astype(bool), "p_value"] <= alpha
                            ).sum()
                        )
                        if table_name == "selected_root" and len(valid)
                        else 0
                    ),
                    "reached_false_split_rate": (
                        float((valid["n_reached_false_splits"] > 0).mean())
                        if table_name == "traversal" and len(valid)
                        else math.nan
                    ),
                    "mean_reached_false_split_fdp": (
                        float(valid["reached_false_split_fdp"].mean())
                        if table_name == "traversal" and len(valid)
                        else math.nan
                    ),
                    "false_split_datasets": (
                        int((valid["n_reached_false_splits"] > 0).sum())
                        if table_name == "traversal" and len(valid)
                        else 0
                    ),
                    "reached_false_split_ci95_low": (
                        wilson_interval(
                            int((valid["n_reached_false_splits"] > 0).sum()), len(valid)
                        )[0]
                        if table_name == "traversal" and len(valid)
                        else math.nan
                    ),
                    "reached_false_split_ci95_high": (
                        wilson_interval(
                            int((valid["n_reached_false_splits"] > 0).sum()), len(valid)
                        )[1]
                        if table_name == "traversal" and len(valid)
                        else math.nan
                    ),
                }
            )
            rows.append(row)
    return pd.DataFrame.from_records(rows)


def _write_protocol_report(
    path: Path,
    *,
    config: KnownBranchConfig,
    summaries: pd.DataFrame,
    git_commit: str,
) -> None:
    pattern_rows = [
        ("global_null", "R, A, B null", "all four leaf means are zero"),
        ("root_only", "R signal; A and B null", "A0=A1=+mu; B0=B1=-mu"),
        ("A_deep_only", "R and B null; A signal", "A0=+mu; A1=-mu; B0=B1=0"),
        ("separate_deep_branches", "R null; A and B signal", "A0=+mu; A1=-mu; B0=+nu; B1=-nu"),
        (
            "strong_root_weak_A",
            "R strong; A weak; B null",
            "A0/A1 average +mu with contrast 0.5 nu; B0=B1=-mu",
        ),
    ]
    pattern_table = "| Pattern | True sibling nodes | Leaf population means |\n|---|---|---|\n"
    pattern_table += "\n".join(f"| `{a}` | {b} | {c} |" for a, b, c in pattern_rows)
    global_null = summaries[
        (summaries.get("table") == "traversal") & (summaries.get("pattern") == "global_null")
    ]
    traversal_lines = "\n".join(
        f"- `{row.covariance_profile}`, passthrough={bool(row.passthrough)}: "
        f"observed any reached false split {row.reached_false_split_rate:.4f} "
        f"(n={int(row.n_ok)}); global-null ideal reference "
        f"{config.alpha if not row.passthrough else 1.0 - (1.0 - config.alpha) ** 2:.4f}."
        for row in global_null.itertuples()
        if pd.notna(row.reached_false_split_rate)
    )
    text = f"""# Known-branch skew diagnostic

Generated for diagnostic attribution only. It does not establish that the selected production hierarchy is calibrated and does not change production methods or defaults.

## Protocol

Each replicate is one independent dataset with four leaf groups, each containing {config.n_per_group} observations in {config.dimension} dimensions. Population covariance profiles are identity and AR(1) with correlation 0.6. Each replicate is reused across matched local variants. For each fixed sibling node, the saved left and right population means and their vector contrast are computed directly from the named leaf-group means; `population_null` means exact equality of those population vectors. The JSON vector columns preserve the actual contrast, including nodes whose descendant groups differ while their parent means coincide.

| Stage | Fixed component | Changed component | Interpretation |
|---|---|---|---|
| Covariance | fixed node membership; full whitened coordinates | known population covariance vs sample pooled covariance | The sample pooled covariance includes between-sibling mean variation; it is labeled as such and is not silently called a within-child estimator. |
| Projection/MP | fixed node membership; covariance held within matched rows | fixed full basis vs PCA and production MP rank from all raw descendant observations | Production Cholesky whitening puts PCA rows and the tested contrast in the same coordinates. Eigenvalue ties are retained by production rank logic; this controlled diagnostic uses MP floor 1. |
| Selection | same dataset; oracle covariance and full basis | fixed root membership vs same-data average-linkage root membership | The selected root is a separate node estimand. Its exact composition-weighted population means are saved; selected roots are not pooled with fixed-node null calibration. |
| FDR/traversal | known tree; exact Gaussian local p-values | production traversal-aligned sibling BH and `TreeDecomposition`, with passthrough disabled/enabled | Every edge prerequisite is deliberately open. Outcomes are empirical for this controlled regime and count only reached `split` decisions. |

## Planted branch patterns

{pattern_table}

`mu = 1.25 e1` and `nu = e2`; all remaining coordinates are zero. Every internal node has a stored population contrast vector, not only a purity or label summary.

## Outcome

{traversal_lines or "No valid global-null traversal rows were produced."}

For the ideal global-null traversal reference, known-covariance full-basis sibling contrasts at R, A, and B are independent in this balanced Gaussian design. The root family is a singleton and the deeper family is BH over A/B. With passthrough disabled, only a root split can be reached; with passthrough enabled, a deeper BH rejection can be reached through a root pass-through. This analytic reference relies on exact Gaussian known-covariance p-values and deliberately open edges; it is not evidence about empirical production p-values or selected hierarchies.

## Run record

- Replicates per pattern/profile: {config.replicates}
- Seed: {config.base_seed}
- Alpha: {config.alpha}
- Source commit: `{git_commit}`
- Rows include `ok`, `unsupported`, and `runtime_error` status plus failure reason.
- Selected-tree root p-values are descriptive after same-data hierarchy selection.
- Raw-statistic skew is reported beside the expected chi-square skew at each observed degree of freedom; p-value rejection rates and quantiles are reported separately.
"""
    path.write_text(text, encoding="utf-8")


def write_diagnostic(
    config: KnownBranchConfig, tables: dict[str, pd.DataFrame], *, run_command: str
) -> dict[str, Path]:
    if config.output_dir.exists() and any(config.output_dir.iterdir()):
        raise FileExistsError(
            f"Refusing to overwrite populated output directory: {config.output_dir}"
        )
    summaries = summarize_tables(tables, alpha=config.alpha)
    try:
        from benchmarks.validation.contracts.report_contract import (
            read_git_commit,
            read_git_worktree_status,
        )

        commit = read_git_commit()
        status = read_git_worktree_status()
    except (OSError, subprocess.SubprocessError, ValueError):
        commit, status = "unavailable", ["unavailable"]
    repository_root = Path(__file__).resolve().parents[4]
    source_paths = (
        Path(__file__).resolve(),
        repository_root
        / "tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py",
        repository_root / "tree_break_selection/hierarchy_analysis/tree_decomposition.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/decomposition.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py",
        repository_root
        / "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py",
        repository_root / "tree_break_selection/eigenvalue_groups.py",
        repository_root / "tree_break_selection/tree/feature_space.py",
        repository_root / "tree_break_selection/tree/construction/hierarchical.py",
        repository_root / "tree_break_selection/tree/poset_tree.py",
        repository_root / "benchmarks/diagnostics/calibration/reporting.py",
        repository_root
        / "benchmarks/diagnostics/calibration/statistics/statistic_distribution_shape_panel.py",
        repository_root / "benchmarks/validation/contracts/report_contract.py",
    )
    protocol_path = repository_root / "reports/known_branch_skew_20260928/protocol.md"
    if protocol_path.exists():
        source_paths = (*source_paths, protocol_path)
    source_hashes = {
        str(path.relative_to(repository_root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in source_paths
    }
    output = write_diagnostic_bundle(
        output_dir=config.output_dir,
        tables={**tables, "summary": summaries},
        filenames={
            "local": "known_branch_local_rows.csv",
            "selected_root": "known_branch_selected_root_rows.csv",
            "traversal": "known_branch_traversal_rows.csv",
            "summary": "known_branch_summary.csv",
        },
        manifest={
            "schema_version": SCHEMA_VERSION,
            "generated_by": GENERATED_BY,
            "study_role": "diagnostic_only_known_branch_stage_attribution",
            "run_command": run_command,
            "source_commit": commit,
            "git_worktree_status": status,
            "runtime": {
                "python": platform.python_version(),
                "numpy": version("numpy"),
                "scipy": version("scipy"),
                "pandas": version("pandas"),
                "networkx": version("networkx"),
            },
            "source_sha256": source_hashes,
            "random_seed_policy": "SeedSequence(base_seed), one spawned generator per dataset replicate, reused across matched stages",
            "config": asdict(config),
            "design": {
                "fixed_tree": ["R", "A", "B", "A0", "A1", "B0", "B1"],
                "patterns": list(config.patterns),
                "covariance_profiles": list(config.covariance_profiles),
                "matched_local_variants": [list(item) for item in VARIANTS],
                "sample_covariance_label": "sample_pooled_covariance_including_between_sibling_means",
                "selected_root_estimand": "same_data_average_linkage_euclidean_root_with_known_composition_weighted_population_means",
                "traversal_edge_regime": "all_edges_deliberately_open",
                "traversal_counts": "actual visited decision=split from TreeDecomposition full trace",
            },
            "limitations": [
                "Gaussian-only bounded attribution design; categorical covariance validation is outside this run.",
                "Selected-root p-values condition on a hierarchy selected from the same observations and are descriptive.",
                "Traversal evidence is empirical under independent known-covariance Gaussian local p-values and open edge prerequisites.",
                "No production method, default, dependency, or existing report is changed.",
            ],
        },
    )
    _write_protocol_report(
        config.output_dir / "report.md", config=config, summaries=summaries, git_commit=commit
    )
    output["report"] = config.output_dir / "report.md"
    return output


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--base-seed", type=int, default=20260928)
    parser.add_argument("--n-per-group", type=int, default=32)
    parser.add_argument("--dimension", type=int, default=8)
    parser.add_argument("--alpha", type=float, default=0.01)
    args = parser.parse_args(argv)
    config = KnownBranchConfig(
        output_dir=args.output_dir,
        replicates=args.replicates,
        base_seed=args.base_seed,
        n_per_group=args.n_per_group,
        dimension=args.dimension,
        alpha=args.alpha,
    )
    tables = run_diagnostic(config)
    outputs = write_diagnostic(
        config,
        tables,
        run_command=" ".join(
            shlex.quote(part) for part in [sys.executable, "-m", GENERATED_BY, *sys.argv[1:]]
        ),
    )
    print(json.dumps({key: str(path) for key, path in outputs.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
