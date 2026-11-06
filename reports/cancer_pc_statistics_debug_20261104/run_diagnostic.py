"""Frozen-tree factor checks; oracle calculations are diagnostics, not calibrated tests."""

import hashlib
import importlib
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from benchmarks.shared.benchmark_runs.runtime import apply_single_thread_runtime_defaults

apply_single_thread_runtime_defaults()

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.metrics import adjusted_rand_score

from benchmarks.diagnostics.calibration.diffusion_factorial import evaluate
from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    selected_hierarchy_null_audit as audit,
)
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.generators.generate_case_data import generate_case_data
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.spectral_context import (
    split_edge_test_rows,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_projected_wald.child_parent_projected_wald_test import (
    run_child_parent_projected_wald_test,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    decide_empirical_null_calibration,
    fit_empirical_null_inflation_model,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.distributions import resolve_node_continuous_covariance_by_block
from tree_break_selection.tree.topology import tree_root


def save(path, value):
    # Support diagnostics can contain infinity; preserve it as an explicit JSON string.
    value = json.loads(json.dumps(value, default=diagnostic_json_default), parse_constant=str)
    path.write_text(
        json.dumps(value, indent=2, default=diagnostic_json_default, allow_nan=False) + "\n"
    )


def oracle_covariance(left, right, row_variance, *, nested):
    weights = np.zeros(len(row_variance))
    weights[left] += 1 / len(left)
    weights[right] -= 1 / len(right)
    base = 1 / len(left) + (-1 if nested else 1) / len(right)
    assert base > 0 and np.isclose(weights.sum(), 0)
    variance = float((weights**2) @ row_variance)
    np.testing.assert_allclose((weights**2).sum(), base)
    return variance / base, variance


def run_case(recipe, output, previous):
    data, truth, values, meta = generate_case_data(recipe)
    space = meta["feature_space"]
    tag = {"case": recipe["name"], "seed": recipe["seed"]}
    module = importlib.import_module(
        "tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.child_parent_divergence_annotation"
    )
    original = module.run_child_parent_tests_across_tree
    captured = {}

    def capture(**kwargs):
        assert not captured
        captured.update(kwargs)
        return original(**kwargs)

    with patch.object(module, "run_child_parent_tests_across_tree", capture):
        result = run_clustering_result(
            data,
            "tbs_family_diffusion",
            dict(METHOD_SPECS["tbs_family_diffusion"].param_grid[0]),
            feature_space=space,
        )
    extra = result.extra
    tree, annotations = extra["tree"], extra["annotations"]
    root = tree_root(tree)
    children = list(tree.successors(root))
    members = {}
    for node in reversed(list(nx.topological_sort(tree))):
        members[node] = (
            [data.index.get_loc(tree.nodes[node]["label"])]
            if tree.out_degree(node) == 0
            else [i for child in tree.successors(node) for i in members[child]]
        )
    means = np.asarray(meta["population_means"])[truth]
    row_variance = np.ones(len(data))
    row_variance[0] = recipe["hub_scale"] ** 2
    mean_branch = compute_mean_branch_length(tree)
    test_rows, _ = split_edge_test_rows(tree, data)
    bundle = extra["gate_bundle"]
    spectral = bundle.edge_gate_result.spectral_context
    dimensions = audit.derive_sibling_projection_dimensions_from_child_edge_comparisons(
        tree, spectral_context=spectral
    )
    projections, eigenvalues = audit.collect_parent_principal_component_inputs_for_sibling_tests(
        dimensions, spectral_context=spectral
    )
    records, _ = audit.collect_sibling_pair_records(
        tree,
        bundle.edge_gate_result.annotated_df,
        feature_space=space,
        sibling_projection_dimensions_from_edge_comparisons=dimensions,
        parent_principal_component_projections=projections,
        parent_principal_component_eigenvalues=eigenvalues,
        adaptive_projection_dimension_energy_fraction=0.9,
    )
    model = fit_empirical_null_inflation_model(records)
    root_record = next(r for r in records if r.parent == root)
    decision = decide_empirical_null_calibration(model, root_record)
    if pd.notna(annotations.loc[root, "Sibling_Divergence_P_Value"]):
        np.testing.assert_allclose(
            decision.p_value, annotations.loc[root, "Sibling_Divergence_P_Value"], rtol=1e-10
        )
    expected = previous[(recipe["name"], recipe["seed"])]
    assert hashlib.sha256(values.tobytes()).hexdigest() == expected["matrix_sha256"]
    assert result.found_clusters == expected["clusters"]
    np.testing.assert_allclose(
        adjusted_rand_score(truth, result.labels), expected["ari"], atol=1e-12
    )

    def covariance(parent):
        return resolve_node_continuous_covariance_by_block(tree, parent, space)

    def sibling_inputs(record):
        left, right = members[record.left], members[record.right]
        return dict(
            left_distribution=values[left].mean(0),
            right_distribution=values[right].mean(0),
            left_sample_size=len(left),
            right_sample_size=len(right),
            branch_length_sum=record.branch_length_sum,
            mean_branch_length=mean_branch,
            projection_dimension_from_edge_comparisons=dimensions[record.parent],
            parent_principal_component_projection=projections[record.parent],
            parent_principal_component_eigenvalues=eigenvalues[record.parent],
            feature_space=space,
            continuous_covariance_by_block=covariance(record.parent),
            adaptive_projection_dimension_energy_fraction=0.9,
        )

    sibling_rows, oracle_records, pair_rows = [], [], []
    for record in records:
        left, right = members[record.left], members[record.right]
        true_norm = float(np.linalg.norm(means[left].mean(0) - means[right].mean(0)))
        inputs = sibling_inputs(record)
        q, scale, df, raw_p = sibling_divergence_test(**inputs)
        np.testing.assert_allclose(
            [q, scale, df],
            [record.stat, record.reference_scale, record.degrees_of_freedom],
            rtol=1e-8,
            atol=1e-8,
        )
        cov = next(iter(inputs["continuous_covariance_by_block"].values()))
        eig = np.linalg.eigvalsh(cov) if cov.ndim == 2 else cov
        noise, _ = oracle_covariance(left, right, row_variance, nested=False)
        oracle = dict(
            inputs,
            continuous_covariance_by_block={
                key: np.eye(30) * noise for key in inputs["continuous_covariance_by_block"]
            },
        )
        oq, os, odf, op = sibling_divergence_test(**oracle)
        if record.n_parent == 2 and true_norm < 1e-10:
            delta = values[left].mean(0) - values[right].mean(0)
            branch_factor = 1 + record.branch_length_sum / (2 * mean_branch)
            expected_q = float((projections[record.parent][0] @ delta) ** 2) / (
                2 * noise * branch_factor
            )
            np.testing.assert_allclose(oq, expected_q, rtol=1e-8, atol=1e-8)
            fixed = dict(
                oracle,
                parent_principal_component_projection=np.eye(30)[:1],
                parent_principal_component_eigenvalues=np.ones(1),
                projection_dimension_from_edge_comparisons=1,
                adaptive_projection_dimension_energy_fraction=None,
            )
            fq, _, fdf, fp = sibling_divergence_test(**fixed)
            assert fdf == odf == 1
            pair_rows.append(
                tag
                | {
                    "parent": record.parent,
                    "q_native": q,
                    "p_native": raw_p,
                    "q_known_noise_pca": oq,
                    "p_known_noise_pca": op,
                    "q_known_noise_fixed_axis": fq,
                    "p_known_noise_fixed_axis": fp,
                    "p_known_noise_full_space": chi2.sf(
                        delta @ delta / (2 * noise * branch_factor), 30
                    ),
                    "branch_factor": branch_factor,
                    "native_q_times_branch_factor": q * branch_factor,
                }
            )
        oracle_records.append(
            replace(record, stat=oq, reference_scale=os, degrees_of_freedom=odf, p_value=op)
        )
        sibling_rows.append(
            tag
            | asdict(record)
            | {
                "true_mean_norm": true_norm,
                "true_mean_null": true_norm < 1e-10,
                "left_n": len(left),
                "right_n": len(right),
                "eligible": pd.notna(
                    annotations.loc[record.parent, "Sibling_Divergence_P_Value_Corrected"]
                ),
                "reject": bool(annotations.loc[record.parent, "Sibling_BH_Different"]),
                "supported": record.has_internal_empirical_null_support,
                "weighted_stat": record.sibling_null_weight * record.stat
                if record.has_internal_empirical_null_support
                else 0,
                "cov_min_eigenvalue": float(eig.min()),
                "cov_max_eigenvalue": float(eig.max()),
                "cov_rank": int(np.count_nonzero(eig > 1e-8)),
                "oracle_stat": oq,
                "oracle_df": odf,
            }
        )
    siblings = pd.DataFrame(sibling_rows)
    supported = siblings[siblings.supported]
    np.testing.assert_allclose(
        (
            supported.weighted_stat.sum()
            / (supported.sibling_null_weight * supported.degrees_of_freedom).sum()
        ),
        model.baseline_scale_estimate,
    )
    siblings.to_csv(output / f"{recipe['name']}-{recipe['seed']}-siblings.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(
        output / f"{recipe['name']}-{recipe['seed']}-pairs.csv", index=False
    )
    traced = TreeDecomposition(
        tree=tree, gate_annotation_bundle=bundle, passthrough=True, trace_level="full"
    ).decompose_tree()
    traced_labels, _ = labels_and_report_from_decomposition(traced, data.index.tolist())
    assert adjusted_rand_score(traced_labels, result.labels) == 1
    traversal = pd.DataFrame(traced["traversal_trace"]).merge(
        siblings[["parent", "true_mean_norm"]], left_on="node_id", right_on="parent", how="left"
    )
    traversal.to_csv(output / f"{recipe['name']}-{recipe['seed']}-traversal.csv", index=False)

    edge_rows, factor_rows = [], []
    for idx, child in enumerate(captured["child_ids"]):
        parent = captured["parent_ids"][idx]
        left = test_rows[child]
        right = [i for other in tree.successors(parent) for i in test_rows[other]]
        true_norm = float(np.linalg.norm(means[left].mean(0) - means[right].mean(0)))
        edge_rows.append(
            tag
            | {
                "parent": parent,
                "child": child,
                "true_mean_norm": true_norm,
                "true_mean_null": true_norm < 1e-10,
                "n_child_test": len(left),
                "n_parent_test": len(right),
                **{
                    k: annotations.loc[child, k]
                    for k in annotations
                    if k.startswith("Child_Parent_")
                },
            }
        )
        if parent != root:
            continue
        inputs = dict(
            child_dist=captured["child_distributions"][idx],
            parent_dist=captured["parent_distributions"][idx],
            n_child=int(captured["child_leaf_counts"][idx]),
            n_parent=int(captured["parent_leaf_counts"][idx]),
            spectral_k=captured["spectral_dims"][root],
            pca_projection=captured["pca_projections"][root],
            pca_eigenvalues=captured["pca_eigenvalues"][root],
            feature_space=space,
            continuous_covariance_by_block=covariance(root),
            branch_length=tree.edges[parent, child]["branch_length"],
            mean_branch_length=mean_branch,
            adaptive_projection_dimension_energy_fraction=0.9,
        )
        np.testing.assert_allclose(inputs["child_dist"], values[left].mean(0), atol=1e-10)
        np.testing.assert_allclose(inputs["parent_dist"], values[right].mean(0), atol=1e-10)
        noise, variance = oracle_covariance(left, right, row_variance, nested=True)
        variants = {
            "native": inputs,
            "no_branch_time": dict(inputs, branch_length=0.0),
            "common_identity_covariance": dict(
                inputs,
                continuous_covariance_by_block={
                    key: np.eye(30) for key in inputs["continuous_covariance_by_block"]
                },
            ),
            "known_noise_covariance": dict(
                inputs,
                continuous_covariance_by_block={
                    key: np.eye(30) * noise for key in inputs["continuous_covariance_by_block"]
                },
            ),
            "full_projection": dict(
                inputs, spectral_k=30, adaptive_projection_dimension_energy_fraction=None
            ),
        }
        native_q = annotations.loc[child, "Child_Parent_Divergence_Test_Statistic"]
        branch_factor = 1 + inputs["branch_length"] / mean_branch
        for variant, args in variants.items():
            q, df, p = run_child_parent_projected_wald_test(**args)
            if variant == "native":
                np.testing.assert_allclose(q, native_q, rtol=1e-10)
            if variant == "no_branch_time":
                np.testing.assert_allclose(q, native_q * branch_factor, rtol=1e-8, atol=1e-8)
            factor_rows.append(
                tag
                | {
                    "gate": "edge",
                    "node": child,
                    "variant": variant,
                    "q": q,
                    "df": df,
                    "c_hat": 1.0,
                    "p": p,
                    "branch_factor": branch_factor,
                    "known_variance": variance,
                }
            )
    edges = pd.DataFrame(edge_rows)
    edges.to_csv(output / f"{recipe['name']}-{recipe['seed']}-edges.csv", index=False)

    inputs = sibling_inputs(root_record)
    noise, variance = oracle_covariance(
        members[root_record.left], members[root_record.right], row_variance, nested=False
    )
    variants = {
        "native": inputs,
        "no_branch_time": dict(inputs, branch_length_sum=0.0),
        "common_identity_covariance": dict(
            inputs,
            continuous_covariance_by_block={
                key: np.eye(30) for key in inputs["continuous_covariance_by_block"]
            },
        ),
        "known_noise_covariance": dict(
            inputs,
            continuous_covariance_by_block={
                key: np.eye(30) * noise for key in inputs["continuous_covariance_by_block"]
            },
        ),
        "full_projection": dict(
            inputs,
            projection_dimension_from_edge_comparisons=30,
            adaptive_projection_dimension_energy_fraction=None,
        ),
    }
    for variant, args in variants.items():
        q, _, df, _ = sibling_divergence_test(**args)
        if variant == "no_branch_time":
            np.testing.assert_allclose(
                q,
                root_record.stat * (1 + root_record.branch_length_sum / (2 * mean_branch)),
                rtol=1e-8,
                atol=1e-8,
            )
        factor_rows.append(
            tag
            | {
                "gate": "sibling",
                "node": root,
                "variant": variant,
                "q": q,
                "df": df,
                "c_hat": decision.c_hat,
                "p": chi2.sf(q / decision.c_hat, df),
                "branch_factor": 1 + root_record.branch_length_sum / (2 * mean_branch),
                "known_variance": variance,
            }
        )
    null_ids = set(siblings.loc[siblings.true_mean_null, "parent"])
    calibration = []
    selections = {
        "native": records,
        "known_mean_null_support_only": [r for r in records if r.parent in null_ids],
        "tested_null_role_only": [r for r in records if not r.is_edge_blocked],
        "known_noise_in_calibration_only": oracle_records,
    }
    for variant, selection in selections.items():
        try:
            fitted = fit_empirical_null_inflation_model(selection)
            d = decide_empirical_null_calibration(fitted, root_record)
            calibration.append(
                tag
                | {
                    "variant": variant,
                    "c_hat": d.c_hat,
                    "p": d.p_value,
                    "status": d.status,
                    "estimator": d.estimator,
                    "family_c": max(1, fitted.baseline_scale_estimate),
                }
            )
        except ValueError as exc:
            calibration.append(
                tag | {"variant": variant, "status": "undefined", "reason": str(exc)}
            )
    calibration.append(
        tag
        | {
            "variant": "c_equals_one",
            "c_hat": 1,
            "p": chi2.sf(
                root_record.stat / root_record.reference_scale, root_record.degrees_of_freedom
            ),
            "status": "diagnostic",
        }
    )
    null_support = supported[supported.true_mean_null]
    null_family_c = max(
        1,
        null_support.weighted_stat.sum()
        / (null_support.sibling_null_weight * null_support.degrees_of_freedom).sum(),
    )
    calibration.append(
        tag
        | {
            "variant": "known_mean_null_support_family_only",
            "c_hat": null_family_c,
            "p": chi2.sf(root_record.stat / null_family_c, root_record.degrees_of_freedom),
            "status": "diagnostic",
        }
    )
    enforced = decide_empirical_null_calibration(
        model, root_record, enforce_support_thresholds=True
    )
    save(
        output / f"{recipe['name']}-{recipe['seed']}-root-calibration.json",
        {"native": decision, "enforced": enforced},
    )
    summary = tag | {
        "clusters": result.found_clusters,
        "ari": expected["ari"],
        "root_sizes": str(sorted(len(members[c]) for c in children)),
        "root_edge_open": bool(
            annotations.loc[children, "Child_Parent_Divergence_Significant"].any()
        ),
        "root_sibling_reject": bool(annotations.loc[root, "Sibling_BH_Different"]),
        "root_sibling_p": decision.p_value,
        "root_c": decision.c_hat,
        "root_q": root_record.stat,
        "support_enforced_status": enforced.status,
        "support_records": len(supported),
        "support_true_nonnull": int((~supported.true_mean_null).sum()),
        "nonnull_stat_share": supported.loc[~supported.true_mean_null, "weighted_stat"].sum()
        / supported.weighted_stat.sum(),
        "edge_eligible": int(edges.Child_Parent_Divergence_Tested.sum()),
        "edge_rejections": int(edges.Child_Parent_Divergence_Significant.sum()),
        "edge_true_null_rejections": int(
            (edges.true_mean_null & edges.Child_Parent_Divergence_Significant).sum()
        ),
        "sibling_eligible": int(siblings.eligible.sum()),
        "sibling_rejections": int(siblings.reject.sum()),
        "sibling_true_null_rejections": int((siblings.true_mean_null & siblings.reject).sum()),
    }
    summary.update(
        edge_true_null_raw_below_alpha=int(
            (edges.true_mean_null & (edges.Child_Parent_Divergence_P_Value <= 0.001)).sum()
        ),
        sibling_true_null_raw_below_alpha=int(
            (siblings.true_mean_null & (siblings.p_value <= 0.01)).sum()
        ),
        true_null_siblings=int(siblings.true_mean_null.sum()),
    )
    print(summary, flush=True)
    return summary, factor_rows, calibration, (tree, data, truth, space)


def main():
    output = (
        Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "generated/complete"
    )
    output.mkdir(parents=True, exist_ok=False)
    old = json.loads(Path("raw/assets/cancer-pc-surrogate-20261104-final/results.json").read_text())
    previous = {(r["case"], r["seed"]): r for r in old}
    recipes = get_test_cases_by_suite("cancer_pc")
    sources = sorted(
        {
            p
            for folder in (
                "tree_break_selection",
                "benchmarks/shared",
                "benchmarks/diagnostics/calibration",
            )
            for p in Path(folder).rglob("*.py")
        }
        | {Path(__file__)}
    )
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    summary, factors, calibration, swaps = [], [], [], []
    for seed in (20261104, 20261106):
        states = {}
        for recipe in recipes:
            row, f, c, state = run_case(dict(recipe, seed=seed), output, previous)
            summary.append(row)
            factors.extend(f)
            calibration.extend(c)
            states[recipe["name"]] = state
        for source in ("cancer_pc_hub_disconnected", "cancer_pc_nohub_disconnected"):
            for target in ("cancer_pc_hub_disconnected", "cancer_pc_nohub_disconnected"):
                tree = states[source][0]
                _, data, truth, space = states[target]
                row, _ = evaluate(tree, data, truth, space, nnls=True, energy=True)
                row.pop("labels")
                if source == target:
                    assert row["clusters"] == previous[(target, seed)]["clusters"]
                    np.testing.assert_allclose(
                        row["ari"], previous[(target, seed)]["ari"], atol=1e-12
                    )
                swaps.append(
                    {"seed": seed, "tree_and_branch_source": source, "data_source": target, **row}
                )
    for name, rows in (
        ("summary", summary),
        ("root_factors", factors),
        ("calibration_factors", calibration),
        ("frozen_tree_swaps", swaps),
    ):
        pd.DataFrame(rows).to_csv(output / f"{name}.csv", index=False)
    z = np.random.default_rng(2026110401).normal(size=(10000, 30))
    q = (z**2).sum(axis=1)
    for vector, expected in zip(z[:10], q[:10], strict=True):
        check = run_projected_wald_kernel(
            vector,
            spectral_k=1,
            pca_projection=(vector / np.linalg.norm(vector))[None, :],
            pca_eigenvalues=np.ones(1),
            adaptive_dimension_energy_fraction=0.9,
        )
        np.testing.assert_allclose(check.statistic, expected)
        assert check.degrees_of_freedom == 1
    save(
        output / "fixed_pair_reference_check.json",
        {
            "replicates": len(z),
            "seed": 2026110401,
            "features": 30,
            "alpha": 0.01,
            "same_pair_pca_rank_one_rejections": int((chi2.sf(q, 1) <= 0.01).sum()),
            "fixed_first_axis_rejections": int((chi2.sf(z[:, 0] ** 2, 1) <= 0.01).sum()),
            "full_space_rejections": int((chi2.sf(q, 30) <= 0.01).sum()),
            "exact_same_pair_pca_rejection_probability": chi2.sf(chi2.isf(0.01, 1), 30),
            "scope": "Independent Gaussian fixed-pair control, known identity covariance, no tree selection, no branch time, no empirical inflation.",
        },
    )
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    save(
        output / "verification.json",
        {
            "native_runs": len(summary),
            "prior_input_and_metric_parity": True,
            "all_sibling_replays": 495 * len(summary),
            "root_edge_replays": 2 * len(summary),
            "branch_rescaling_identity": True,
            "contrast_weight_variance_identity": True,
            "matched_frozen_tree_replays": 4,
            "full_trace_partition_replays": len(summary),
            "known_noise_two_leaf_rank_one_identity": True,
            "fixed_pair_native_kernel_replays": 10,
            "sources_unchanged_during_run": True,
            "scope": "Observed selected-tree mechanism diagnosis; not a Type-I-error or FDR estimate.",
        },
    )
    save(
        output / "manifest.json",
        {
            "recipes": recipes,
            "seeds": [20261104, 20261106],
            "method": "tbs_family_diffusion",
            "params": METHOD_SPECS["tbs_family_diffusion"].param_grid[0],
            "source_sha256": hashes,
        },
    )


if __name__ == "__main__":
    main()
