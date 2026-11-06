"""Attribute Gaussian root null behavior to tree selection versus covariance.

Run with ``uv run --no-sync python -m benchmarks.diagnostics.calibration.diffusion_attribution``.
Trees use the preceding factorial's full diffusion geometry. Fresh independent
observations populate each fixed tree without rebuilding or relabeling it.
Full-space identity projections, zero branch lengths, and c_hat=1 isolate a 2x2
comparison. The root sibling test is primary; full-row child-parent comparisons
are algebra diagnostics, different from production's split-row edge annotation.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, f

from benchmarks.diagnostics.calibration.diffusion_factorial import datasets, digest
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_projected_wald.child_parent_projected_wald_test import (
    run_child_parent_projected_wald_test,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.types.sibling_pair_record import (
    SiblingPairRecord,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.space_separation.diffusion import feature_family_diffusion_geometry
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.distributions import (
    CONTINUOUS_COVARIANCE_BY_BLOCK,
    resolve_node_continuous_covariance_by_block,
)
from tree_break_selection.tree.topology import tree_root

INDEPENDENT_SEED = 2026110303
COVARIANCE_MODES = ("known_identity", "production_guarded")


@dataclass(frozen=True)
class TreeDataPair:
    replicate: int
    selected_seed: int
    independent_seed: int
    selected_data: pd.DataFrame
    independent_data: pd.DataFrame
    feature_space: object
    tree: object
    topology_hash: str
    data_hash: str


def gaussian_tree_pairs(replicates=50):
    """Reproduce continuous factorial datasets and freeze their diffusion trees."""
    for case, replicate, seed, _, selected, _, space, _ in itertools.islice(
        datasets(replicates, 0), replicates
    ):
        assert case == "iid_continuous"
        geometry = feature_family_diffusion_geometry(selected, space)
        build = build_tree(selected, geometry.distance_condensed, linkage_method="average")
        for edge in build.tree.edges:
            build.tree.edges[edge]["branch_length"] = 0.0
        independent_seed = INDEPENDENT_SEED + replicate
        independent = pd.DataFrame(
            np.random.default_rng(independent_seed).normal(size=selected.shape),
            index=selected.index.copy(),
            columns=selected.columns.copy(),
        )
        yield TreeDataPair(
            replicate,
            seed,
            independent_seed,
            selected,
            independent,
            space,
            build.tree,
            digest(build.linkage_matrix[:, :2].tolist()),
            digest(selected.values.tolist()),
        )


def covariance_for_node(tree, node, space, mode):
    if mode == "known_identity":
        return {block.name: np.eye(block.raw_dimension) for block in space.continuous_blocks}
    if mode != "production_guarded":
        raise ValueError(f"Unsupported covariance mode: {mode!r}.")
    return resolve_node_continuous_covariance_by_block(
        tree,
        node,
        space,
        continuous_covariance_policy="guarded_within_child",
        continuous_covariance_min_child_leaf_count=8,
    )


def full_space_sibling_records(tree, data, space, covariance_mode):
    """Return a populated copy and full-dimensional records for every binary node.

    Null-role metadata is deliberately absent. Restoring empirical inflation
    requires assigning roles and weights from its declared edge annotations.
    The input topology, annotations, and branch lengths are preserved.
    """
    populated = copy.deepcopy(tree)
    for edge in populated.edges:
        populated.edges[edge]["branch_length"] = 0.0
    populated.populate_node_distributions(data, feature_space=space)
    p = space.contrast_dimension
    projection, eigenvalues = np.eye(p), np.ones(p)
    records = []
    for parent in populated.nodes:
        children = list(populated.successors(parent))
        if len(children) != 2:
            continue
        left, right = children
        first, second = populated.nodes[left], populated.nodes[right]
        statistic, reference_scale, degrees, p_value = sibling_divergence_test(
            first["distribution"],
            second["distribution"],
            first["leaf_count"],
            second["leaf_count"],
            branch_length_sum=0.0,
            mean_branch_length=None,
            projection_dimension_from_edge_comparisons=p,
            parent_principal_component_projection=projection,
            parent_principal_component_eigenvalues=eigenvalues,
            feature_space=space,
            continuous_covariance_by_block=covariance_for_node(
                populated, parent, space, covariance_mode
            ),
        )
        records.append(
            SiblingPairRecord(
                parent=parent,
                left=left,
                right=right,
                stat=statistic,
                reference_scale=reference_scale,
                degrees_of_freedom=degrees,
                p_value=p_value,
                branch_length_sum=0.0,
                n_parent=populated.nodes[parent]["leaf_count"],
                is_null_like=False,
                sibling_projection_dimension=degrees,
                parent_positive_eigenvalue_count=float(p),
                feature_family="continuous",
            )
        )
    return populated, records


def root_covariance_diagnostics(tree, root, space, mode):
    blocks = covariance_for_node(tree, root, space, mode)
    block = next(iter(blocks.values()))
    eigenvalues = np.linalg.eigvalsh(block) if block.ndim == 2 else block
    children = list(tree.successors(root))
    counts = [tree.nodes[child]["leaf_count"] for child in children]
    source = "known_identity"
    if mode == "production_guarded":
        source = (
            "parent_total_fallback"
            if blocks is tree.nodes[root][CONTINUOUS_COVARIANCE_BY_BLOCK]
            else "pooled_within_child"
        )
    fallback = ""
    if source == "parent_total_fallback":
        fallback = "child_count_below_8" if min(counts) < 8 else "residual_df_below_raw_dimension"
    matrix = block if block.ndim == 2 else np.diag(block)
    return {
        "covariance_source": source,
        "fallback_reason": fallback,
        "covariance_representation": "dense" if block.ndim == 2 else "diagonal",
        "covariance_rank": int(np.linalg.matrix_rank(matrix)),
        "covariance_condition": float(np.linalg.cond(matrix)),
        "covariance_min_eigenvalue": float(np.min(eigenvalues)),
        "covariance_max_eigenvalue": float(np.max(eigenvalues)),
        "covariance_trace": float(np.sum(eigenvalues)),
        "pooled_residual_df": int(sum(counts) - 2),
    }


def baseline_rows(pair):
    root_rows = []
    for data_mode, data in (
        ("selected", pair.selected_data),
        ("independent", pair.independent_data),
    ):
        for mode in COVARIANCE_MODES:
            tree, records = full_space_sibling_records(pair.tree, data, pair.feature_space, mode)
            assert tuple(tree.edges()) == tuple(pair.tree.edges())
            root = tree_root(tree)
            record = next(record for record in records if record.parent == root)
            covariance = covariance_for_node(tree, root, pair.feature_space, mode)
            covariance_metadata = root_covariance_diagnostics(tree, root, pair.feature_space, mode)
            children = (record.left, record.right)
            p = data.shape[1]
            edge_results = []
            for child in children:
                edge_results.append(
                    run_child_parent_projected_wald_test(
                        tree.nodes[child]["distribution"],
                        tree.nodes[root]["distribution"],
                        tree.nodes[child]["leaf_count"],
                        tree.nodes[root]["leaf_count"],
                        spectral_k=p,
                        pca_projection=np.eye(p),
                        pca_eigenvalues=np.ones(p),
                        feature_space=pair.feature_space,
                        continuous_covariance_by_block=covariance,
                    )
                )
            # Both edges are the same standardized contrast with opposite signs.
            np.testing.assert_allclose(
                [result[0] for result in edge_results], record.stat, rtol=1e-8, atol=1e-10
            )
            exact_f_statistic = exact_f_p = ridge_statistic_difference = np.nan
            if (
                covariance_metadata["covariance_source"] == "pooled_within_child"
                and covariance_metadata["covariance_representation"] == "dense"
            ):
                # Only this actual matrix is the ordinary pooled full covariance.
                block = next(iter(covariance.values()))
                difference = (
                    tree.nodes[record.left]["distribution"]
                    - tree.nodes[record.right]["distribution"]
                )
                variance = sum(1 / tree.nodes[child]["leaf_count"] for child in children)
                unregularized = float(difference @ np.linalg.solve(block, difference) / variance)
                nu = len(data) - 2
                exact_f_statistic = (nu - p + 1) * unregularized / (p * nu)
                exact_f_p = float(f.sf(exact_f_statistic, p, nu - p + 1))
                ridge_statistic_difference = record.stat - unregularized
            root_rows.append(
                {
                    "replicate": pair.replicate,
                    "selected_seed": pair.selected_seed,
                    "independent_seed": pair.independent_seed,
                    "data_mode": data_mode,
                    "covariance_mode": mode,
                    "topology_hash": pair.topology_hash,
                    "selected_data_hash": pair.data_hash,
                    "test_data_hash": digest(data.values.tolist()),
                    "root": root,
                    "left": record.left,
                    "right": record.right,
                    "n_left": tree.nodes[record.left]["leaf_count"],
                    "n_right": tree.nodes[record.right]["leaf_count"],
                    "n": len(data),
                    "p": p,
                    "statistic": record.stat,
                    "df": record.degrees_of_freedom,
                    "p_value": record.p_value,
                    "reject_sibling_alpha_0_01": record.p_value < 0.01,
                    "edge_left_statistic": edge_results[0][0],
                    "edge_right_statistic": edge_results[1][0],
                    "edge_left_p": edge_results[0][2],
                    "edge_right_p": edge_results[1][2],
                    "reject_edge_alpha_0_001": min(result[2] for result in edge_results) < 0.001,
                    "dense_pooled_unregularized_f_statistic": exact_f_statistic,
                    "dense_pooled_unregularized_f_p": exact_f_p,
                    "production_ridge_statistic_difference": ridge_statistic_difference,
                    **covariance_metadata,
                }
            )
    assert (
        len({(row["left"], row["right"], row["n_left"], row["n_right"]) for row in root_rows}) == 1
    )
    return root_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument(
        "--output", type=Path, default=Path("raw/assets/diffusion-attribution-20261102")
    )
    args = parser.parse_args()
    if not 1 <= args.replicates <= 100:
        parser.error("replicates must be between 1 and 100 to keep seed ranges disjoint")
    args.output.mkdir(parents=True, exist_ok=True)
    previous_path = Path("raw/assets/diffusion-factorial-20261102/runs.csv")
    previous = pd.read_csv(previous_path)
    reference = (
        previous[(previous["case"] == "iid_continuous") & (previous.diffusion == 1)]
        .groupby("replicate")
        .first()
    )
    rows = []
    for pair in gaussian_tree_pairs(args.replicates):
        if pair.replicate in reference.index:
            assert pair.topology_hash == reference.loc[pair.replicate, "topology_hash"]
            assert pair.data_hash == reference.loc[pair.replicate, "data_hash"]
        rows.extend(baseline_rows(pair))
        pd.DataFrame(rows).to_csv(args.output / "base_runs.csv", index=False)
        print(f"baseline replicate={pair.replicate} complete", flush=True)
    frame = pd.DataFrame(rows)
    summary = []
    for (data_mode, mode), group in frame.groupby(["data_mode", "covariance_mode"]):
        rejected = int(group.reject_sibling_alpha_0_01.sum())
        ci = binomtest(rejected, len(group)).proportion_ci()
        summary.append(
            {
                "data_mode": data_mode,
                "covariance_mode": mode,
                "replicates": len(group),
                "sibling_rejections": rejected,
                "sibling_rejection_fraction": rejected / len(group),
                "sibling_fraction_ci_low": ci.low,
                "sibling_fraction_ci_high": ci.high,
                "full_row_edge_rejections": int(group.reject_edge_alpha_0_001.sum()),
                "mean_statistic": float(group.statistic.mean()),
                "mean_p_value": float(group.p_value.mean()),
                "parent_total_fallbacks": int(
                    (group.covariance_source == "parent_total_fallback").sum()
                ),
                "diagonal_covariances": int((group.covariance_representation == "diagonal").sum()),
            }
        )
    pd.DataFrame(summary).to_csv(args.output / "base_summary.csv", index=False)
    strata = (
        frame.groupby(
            ["data_mode", "covariance_mode", "covariance_source", "covariance_representation"]
        )
        .agg(
            replicates=("replicate", "size"),
            sibling_rejections=("reject_sibling_alpha_0_01", "sum"),
            mean_statistic=("statistic", "mean"),
            mean_p_value=("p_value", "mean"),
            minimum_left_count=("n_left", "min"),
            minimum_right_count=("n_right", "min"),
        )
        .reset_index()
    )
    strata.to_csv(args.output / "base_covariance_strata.csv", index=False)
    sources = [
        Path(__file__),
        Path("tree_break_selection/tree/distributions.py"),
        Path("tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py"),
    ]
    manifest = {
        "replicates": args.replicates,
        "selected_seed": 2026110203,
        "independent_seed": INDEPENDENT_SEED,
        "n": 80,
        "p": 24,
        "covariance_modes": COVARIANCE_MODES,
        "branch_lengths": 0,
        "c_hat": 1,
        "projection": "identity full 24-dimensional space, no PCA or adaptive rank",
        "primary": "root sibling p < .01; no FDR; dataset is replicate",
        "secondary": "full-row child-parent means p < .001, not split-row production annotation",
        "known_identity_null": "independent test data conditional on fixed tree gives chi-square24 before numerical ridge",
        "production_covariance_reference": "exact production guarded estimator and chi-square plug-in reference; not assumed calibrated",
        "dense_pooled_f": "unregularized ordinary pooled covariance only when actual representation is dense; F(p,n-p-1) reference exact for independent Gaussian rows, descriptive for selected rows; production adds 1e-12 contrast covariance jitter",
        "prior_runs": str(previous_path),
        "prior_matching_topologies": args.replicates <= len(reference),
        "source_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources
        },
        "limitation": "50 dataset replicates are attribution evidence, not a tail probability proof",
    }
    (args.output / "base_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
