"""Bounded Gaussian reference-law controls and a saved selected-tree replay."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import beta, chi2, f, kstest

from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.spectral_context import (
    split_edge_test_rows,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_projected_wald.child_parent_projected_wald_test import (
    run_child_parent_projected_wald_test,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.tree.distributions import (
    populate_distributions,
    resolve_node_continuous_covariance_by_block,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

ROOT = Path(__file__).resolve().parent
DRAWS = 20_000
RNG = np.random.default_rng(20261104)


def check(statistic, nominal_df, reference_p):
    nominal_p = chi2.sf(statistic, nominal_df)
    return {
        "mean_q": float(np.mean(statistic)),
        "nominal_rejections_at_0.05": int(np.sum(nominal_p <= 0.05)),
        "reference_rejections_at_0.05": int(np.sum(reference_p <= 0.05)),
        "nominal_rejections_at_0.01": int(np.sum(nominal_p <= 0.01)),
        "reference_rejections_at_0.01": int(np.sum(reference_p <= 0.01)),
        "reference_uniform_ks_distance": float(kstest(reference_p, "uniform").statistic),
    }


def covariance_controls():
    n, p, k = 12, 6, 3
    data = RNG.normal(size=(DRAWS, n, p))
    contrasts = (data[:, :6].mean(axis=1) - data[:, 6:].mean(axis=1)) * np.sqrt(3.0)
    fixed_q = np.sum(contrasts[:, :k] ** 2, axis=1)

    parent_residual = data - data.mean(axis=1, keepdims=True)
    parent_factor = np.linalg.cholesky(
        parent_residual.transpose(0, 2, 1) @ parent_residual / (n - 1)
    )
    total_z = np.linalg.solve(parent_factor, contrasts[..., None])[..., 0]
    parent_q = np.sum(total_z[:, :k] ** 2, axis=1)

    groups = data.reshape(DRAWS, 2, 6, p)
    within_residual = (groups - groups.mean(axis=2, keepdims=True)).reshape(DRAWS, n, p)
    within_factor = np.linalg.cholesky(
        within_residual.transpose(0, 2, 1) @ within_residual / (n - 2)
    )
    within_z = np.linalg.solve(within_factor, contrasts[..., None])[..., 0]
    within_q = np.sum(within_z**2, axis=1)
    heldout = data[:, ::2, 0]
    heldout_z = (heldout[:, :3].mean(axis=1) - heldout[:, 3:].mean(axis=1)) * np.sqrt(1.5)
    heldout_q = (heldout_z / within_factor[:, 0, 0]) ** 2
    residual_q = np.sum(within_residual[:, :, 0] ** 2, axis=1)
    scalar = data[..., 0]
    pooled_coeff = np.r_[np.ones(6), -np.ones(6)] / np.sqrt(12)
    heldout_coeff = np.zeros(n)
    heldout_coeff[::2] = np.r_[np.ones(3), -np.ones(3)] / np.sqrt(6)
    residual_coeff = heldout_coeff - (heldout_coeff @ pooled_coeff) * pooled_coeff
    overlap = residual_coeff @ residual_coeff
    np.testing.assert_allclose(overlap, 0.5)
    np.testing.assert_allclose(
        heldout_z, np.sqrt(1 - overlap) * (scalar @ pooled_coeff) + scalar @ residual_coeff
    )

    return {
        "known_covariance_fixed_projection": check(fixed_q, k, chi2.sf(fixed_q, k)),
        "same_row_parent_total_fixed_projection": check(
            parent_q, k, beta.sf(parent_q / (n - 1), k / 2, (n - 1 - k) / 2)
        ),
        "within_child_full_means_full_projection": check(
            within_q, p, f.sf(within_q * (n - 2 - p + 1) / ((n - 2) * p), p, n - 2 - p + 1)
        ),
        "within_child_heldout_means_scalar": {
            "residual_contrast_overlap_squared": float(overlap),
            "contrast_square_vs_residual_sum_squares_covariance": float(
                np.cov(heldout_z**2, residual_q)[0, 1]
            ),
            "expected_covariance": float(2 * overlap),
            "chi_square_rejections_at_0.05": int(np.sum(chi2.sf(heldout_q, 1) <= 0.05)),
            "inapplicable_F_rejections_at_0.05": int(np.sum(f.sf(heldout_q, 1, n - 2) <= 0.05)),
        },
    }


def guarded_covariance_replay():
    rng = np.random.default_rng(20261105)
    leaves_per_child, dimension, native_replays = 8, 6, 25
    residual_df = 2 * (leaves_per_child - 1)
    groups = rng.normal(size=(DRAWS, 2, leaves_per_child, dimension))
    means = groups.mean(axis=2)
    residuals = groups - means[:, :, None, :]
    pooled = np.einsum("bgri,bgrj->bij", residuals, residuals) / residual_df
    full_delta = means[:, 0] - means[:, 1]
    test_means = groups[:, :, ::2].mean(axis=2)
    edge_delta = test_means[:, 0] - test_means.mean(axis=1)
    full_variance_factor = 2 / leaves_per_child
    edge_variance_factor = 1 / (leaves_per_child // 2) - 1 / leaves_per_child

    def quadratic(delta, variance_factor):
        return (
            np.sum(delta * np.linalg.solve(pooled, delta[..., None])[..., 0], axis=1)
            / variance_factor
        )

    full_q = quadratic(full_delta, full_variance_factor)
    edge_q = quadratic(edge_delta, edge_variance_factor)
    denominator_df = residual_df - dimension + 1
    f_factor = denominator_df / (residual_df * dimension)

    tree = nx.DiGraph()
    tree.add_nodes_from(["root", "left", "right"], is_leaf=False)
    tree.add_edges_from([("root", "left"), ("root", "right")])
    labels = [f"leaf_{index:02}" for index in range(2 * leaves_per_child)]
    for index, label in enumerate(labels):
        tree.add_node(label, is_leaf=True, label=label)
        tree.add_edge("left" if index < leaves_per_child else "right", label)
    columns = tuple(f"feature_{index}" for index in range(dimension))
    feature_space = continuous_feature_space_from_columns(columns)
    block_name = feature_space.continuous_blocks[0].name
    basis = np.eye(dimension)
    eigenvalues = np.arange(dimension, 0, -1, dtype=float)
    max_error = 0.0
    for index in range(native_replays):
        frame = pd.DataFrame(groups[index].reshape(-1, dimension), index=labels, columns=columns)
        populate_distributions(tree, frame, feature_space=feature_space)
        covariance = resolve_node_continuous_covariance_by_block(tree, "root", feature_space)
        assert covariance is not None and covariance[block_name].shape == (dimension, dimension)
        np.testing.assert_allclose(covariance[block_name], pooled[index], atol=1e-14)
        sibling_q, _, sibling_df, _ = sibling_divergence_test(
            tree.nodes["left"]["distribution"],
            tree.nodes["right"]["distribution"],
            leaves_per_child,
            leaves_per_child,
            branch_length_sum=0.0,
            mean_branch_length=None,
            projection_dimension_from_edge_comparisons=dimension,
            parent_principal_component_projection=basis,
            parent_principal_component_eigenvalues=eigenvalues,
            feature_space=feature_space,
            continuous_covariance_by_block=covariance,
        )
        test_rows, _ = split_edge_test_rows(tree, frame)
        left_rows = test_rows["left"]
        parent_rows = left_rows + test_rows["right"]
        edge_native_q, edge_df, _ = run_child_parent_projected_wald_test(
            frame.iloc[left_rows].to_numpy().mean(axis=0),
            frame.iloc[parent_rows].to_numpy().mean(axis=0),
            len(left_rows),
            len(parent_rows),
            spectral_k=dimension,
            pca_projection=basis,
            pca_eigenvalues=eigenvalues,
            feature_space=feature_space,
            continuous_covariance_by_block=covariance,
        )
        assert sibling_df == edge_df == dimension
        error = max(abs(sibling_q - full_q[index]), abs(edge_native_q - edge_q[index]))
        max_error = max(max_error, error)
    assert max_error < 1e-7

    def f_p(q):
        return f.sf(q * f_factor, dimension, denominator_df)

    full = check(full_q, dimension, f_p(full_q))
    assert full["reference_uniform_ks_distance"] < 0.02
    return {
        "seed": 20261105,
        "draws": DRAWS,
        "native_replays_per_contrast": native_replays,
        "native_max_abs_q_error_vs_unregularized_formula": max_error,
        "full_means_hotelling": full,
        "alternating_test_half_means": {
            "mean_q": float(edge_q.mean()),
            "nominal_rejections_at_0.05": int(np.sum(chi2.sf(edge_q, dimension) <= 0.05)),
            "inapplicable_hotelling_rejections_at_0.05": int(np.sum(f_p(edge_q) <= 0.05)),
        },
    }


def projection_controls():
    z = RNG.normal(size=(DRAWS, 30))
    q = np.sum(z**2, axis=1)
    for row, expected in zip(z[:25], q[:25], strict=True):
        native = run_projected_wald_kernel(
            row,
            spectral_k=1,
            pca_projection=(row / np.linalg.norm(row))[None, :],
            pca_eigenvalues=np.ones(1),
        )
        assert native.degrees_of_freedom == 1
        np.testing.assert_allclose(native.statistic, expected)
    adaptive_z = RNG.normal(size=(DRAWS, 8))
    basis = np.eye(8)
    eigenvalues = np.arange(8, 0, -1)
    adaptive_results = [
        run_projected_wald_kernel(
            row,
            spectral_k=2,
            pca_projection=basis,
            pca_eigenvalues=eigenvalues,
            adaptive_dimension_energy_fraction=0.9,
        )
        for row in adaptive_z
    ]
    adaptive = np.asarray(
        [(result.degrees_of_freedom, result.p_value) for result in adaptive_results]
    )
    return {
        "known_covariance_fitted_pair_projection": check(q, 1, chi2.sf(q, 30)),
        "fixed_basis_adaptive_rank_0.9": {
            str(k): {
                "draws": int(np.sum(adaptive[:, 0] == k)),
                "nominal_rejections_at_0.05": int(
                    np.sum((adaptive[:, 0] == k) & (adaptive[:, 1] <= 0.05))
                ),
                "nominal_rejections_at_0.01": int(
                    np.sum((adaptive[:, 0] == k) & (adaptive[:, 1] <= 0.01))
                ),
            }
            for k in range(1, 9)
        },
    }


def saved_pairs():
    path = (
        ROOT.parent
        / "cancer_pc_statistics_debug_20261104/generated/complete"
        / "cancer_pc_hub_disconnected-20261106-pairs.csv"
    )
    if not path.exists():
        return {"unavailable": str(path.relative_to(ROOT.parent.parent))}
    with path.open() as source:
        pairs = list(csv.DictReader(source))
    products = [float(row["native_q_times_branch_factor"]) for row in pairs]
    assert len(pairs) == 183
    np.testing.assert_allclose(products, 30, atol=4e-5, rtol=0)
    return {
        "pairs": len(pairs),
        "native_p_at_most_0.01": sum(float(row["p_native"]) <= 0.01 for row in pairs),
        "q_times_branch_factor_min": min(products),
        "q_times_branch_factor_max": max(products),
    }


def main():
    covariance = covariance_controls()
    projection = projection_controls()
    for result in [
        *list(covariance.values())[:3],
        projection["known_covariance_fitted_pair_projection"],
    ]:
        assert result["reference_uniform_ks_distance"] < 0.02
    assert (
        sum(row["draws"] for row in projection["fixed_basis_adaptive_rank_0.9"].values()) == DRAWS
    )
    output = {
        "seed": 20261104,
        "draws_per_control": DRAWS,
        **covariance,
        **projection,
        "saved_selected_pairs": saved_pairs(),
        "native_guard_8_per_child_p6": guarded_covariance_replay(),
    }
    destination = ROOT / "generated/reference/reference_probe.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
