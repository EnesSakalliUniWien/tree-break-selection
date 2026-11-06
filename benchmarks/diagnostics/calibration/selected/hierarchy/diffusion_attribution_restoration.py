"""Direct projection/inflation attribution on fixed diffusion trees; diagnostic only.

Edge roles are recomputed per data matrix once using production guarded covariance,
then frozen across the crossed interventions. Native PCA is refitted per matrix and fixed across covariance arms;
contrast-energy rank is off, but MP rank and tie expansion remain data dependent.
All branch lengths are zero. No traversal or calibrated error-control claim.
"""

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration.diffusion_attribution import (
    full_space_sibling_records,
    gaussian_tree_pairs,
    root_covariance_diagnostics,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence import (
    annotate_child_parent_divergence,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    decide_empirical_null_calibration,
    fit_empirical_null_inflation_model,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.pair_observations import (
    extract_sibling_pair_observations,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.record_collection import (
    collect_sibling_pair_records,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.projection_dimensions import (
    derive_sibling_projection_dimensions_from_child_edge_comparisons,
)

ROLE_FIELDS = (
    "is_null_like",
    "is_edge_blocked",
    "sibling_null_weight",
    "calibration_dependency_group",
    "calibration_dependency_weight",
)
EDGE_ALPHA = 0.001
SIBLING_ALPHA = 0.01


def _projection_records(tree, annotations, space, spectral):
    ranks = derive_sibling_projection_dimensions_from_child_edge_comparisons(
        tree, spectral_context=spectral
    )
    return collect_sibling_pair_records(
        tree,
        annotations,
        sibling_projection_dimensions_from_edge_comparisons=ranks,
        parent_principal_component_projections=spectral.principal_component_projections_by_node,
        parent_principal_component_eigenvalues=spectral.principal_component_eigenvalues_by_node,
        feature_space=space,
        adaptive_projection_dimension_energy_fraction=None,
    )[0], ranks


def compare_matrix(tree, data, space):
    """Return all crossed arms, preserving edge-role metadata across arms."""
    guarded_tree, guarded_full = full_space_sibling_records(
        tree, data, space, covariance_mode="production_guarded"
    )
    annotations, guarded_spectral = annotate_child_parent_divergence(
        guarded_tree,
        pd.DataFrame(index=list(guarded_tree.nodes)),
        leaf_data=data,
        feature_space=space,
        significance_level_alpha=EDGE_ALPHA,
        adaptive_projection_dimension_energy_fraction=None,
    )
    guarded_projected, guarded_ranks = _projection_records(
        guarded_tree, annotations, space, guarded_spectral
    )
    roles = {r.parent: {key: getattr(r, key) for key in ROLE_FIELDS} for r in guarded_projected}
    root = next(node for node in guarded_tree if guarded_tree.in_degree(node) == 0)
    output = []
    for covariance_mode in ("known_identity", "production_guarded"):
        current_tree, full_records = (
            (guarded_tree, guarded_full)
            if covariance_mode == "production_guarded"
            else full_space_sibling_records(tree, data, space, covariance_mode=covariance_mode)
        )
        if covariance_mode == "production_guarded":
            projected, ranks = guarded_projected, guarded_ranks
        else:
            # Hold the exact native all-leaf basis/rank fixed across covariance arms.
            spectral = guarded_spectral
            projected, ranks = guarded_projected, guarded_ranks
            updated = []
            for record in projected:
                left, right, nl, nr, _, _ = extract_sibling_pair_observations(
                    current_tree, record.parent, record.left, record.right
                )
                stat, scale, df, p = sibling_divergence_test(
                    left,
                    right,
                    nl,
                    nr,
                    branch_length_sum=0.0,
                    mean_branch_length=None,
                    projection_dimension_from_edge_comparisons=ranks[record.parent],
                    parent_principal_component_projection=spectral.principal_component_projections_by_node[
                        record.parent
                    ],
                    parent_principal_component_eigenvalues=spectral.principal_component_eigenvalues_by_node[
                        record.parent
                    ],
                    feature_space=space,
                    continuous_covariance_by_block={
                        block.name: np.eye(block.raw_dimension) for block in space.continuous_blocks
                    },
                    adaptive_projection_dimension_energy_fraction=None,
                )
                updated.append(
                    replace(
                        record,
                        stat=stat,
                        reference_scale=scale,
                        degrees_of_freedom=df,
                        p_value=p,
                        sibling_projection_dimension=df,
                    )
                )
            projected = updated
        for projection, records in [("full_feature", full_records), ("production", projected)]:
            records = [replace(r, **roles[r.parent]) for r in records]
            model = (
                fit_empirical_null_inflation_model(records)
                if any(r.has_internal_empirical_null_support for r in records)
                else None
            )
            for record in records:
                covariance_metadata = root_covariance_diagnostics(
                    current_tree, record.parent, space, covariance_mode
                )
                decision = (
                    decide_empirical_null_calibration(model, record) if model is not None else None
                )
                for inflation in ("none", "empirical"):
                    available = (
                        inflation == "none" or decision is not None and decision.p_value is not None
                    )
                    p = (
                        record.p_value
                        if inflation == "none"
                        else decision.p_value
                        if decision
                        else None
                    )
                    c = 1.0 if inflation == "none" else decision.c_hat if decision else None
                    output.append(
                        dict(
                            parent=record.parent,
                            is_root=record.parent == root,
                            left_leaf_count=current_tree.nodes[record.left]["leaf_count"],
                            right_leaf_count=current_tree.nodes[record.right]["leaf_count"],
                            spectral_covariance_chart=("native_parent_total"),
                            full_span=(record.degrees_of_freedom == data.shape[1]),
                            **covariance_metadata,
                            left=record.left,
                            right=record.right,
                            n_parent=record.n_parent,
                            covariance=covariance_mode,
                            projection=projection,
                            inflation=inflation,
                            raw_statistic=record.stat,
                            raw_df=record.degrees_of_freedom,
                            raw_p=record.p_value,
                            c_hat=c,
                            adjusted_p=p,
                            nominal_projection_dimension=(
                                data.shape[1]
                                if projection == "full_feature"
                                else ranks[record.parent]
                            ),
                            available=available,
                            calibration_status=(
                                "unadjusted_diagnostic"
                                if inflation == "none"
                                else decision.status
                                if decision
                                else "undefined_no_internal_support"
                            ),
                            estimator=(
                                "none"
                                if inflation == "none"
                                else decision.estimator
                                if decision
                                else None
                            ),
                            support=json.dumps(
                                decision.support if decision else {}, sort_keys=True
                            ),
                            gate_eligible=not record.is_null_like,
                            descriptive_reject_at_0_01=(
                                bool(p < SIBLING_ALPHA) if p is not None else None
                            ),
                            production_gate_p=(
                                p if inflation == "empirical" and not record.is_null_like else None
                            ),
                            gate_status=(
                                "skipped_null_like"
                                if record.is_null_like
                                else "unadjusted_diagnostic"
                                if inflation == "none"
                                else decision.status
                                if decision
                                else "undefined_no_internal_support"
                            ),
                            branch_length_sum=record.branch_length_sum,
                            edge_roles_recomputed=False,
                            edge_role_source="per_matrix_production_guarded_edge_alpha_0.001",
                            **roles[record.parent],
                        )
                    )
    return pd.DataFrame(output)


def run_panel(replicates=50):
    frames = []
    for pair in gaussian_tree_pairs(replicates=replicates):
        for source, data in [
            ("selected", pair.selected_data),
            ("independent", pair.independent_data),
        ]:
            rows = compare_matrix(pair.tree, data, pair.feature_space)
            rows["replicate"] = pair.replicate
            rows["data_source"] = source
            rows["selected_seed"] = pair.selected_seed
            rows["independent_seed"] = pair.independent_seed
            rows["topology_hash"] = pair.topology_hash
            rows["data_hash"] = hashlib.sha256(
                np.ascontiguousarray(data.to_numpy()).tobytes()
            ).hexdigest()
            frames.append(rows)
        print(f"completed restoration replicate {pair.replicate}", flush=True)
    return pd.concat(frames, ignore_index=True)


def write_panel(rows, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(output_dir / "restoration_nodes.csv.gz", index=False)
    rows[rows.is_root].to_csv(output_dir / "restoration_roots.csv", index=False)
    keys = ["data_source", "covariance", "projection", "inflation"]
    summary = (
        rows[rows.is_root]
        .groupby(keys)
        .agg(
            independent_dataset_count=("replicate", "nunique"),
            available_count=("available", "sum"),
            raw_rejections=("raw_p", lambda p: int((p < SIBLING_ALPHA).sum())),
            descriptive_adjusted_rejections=(
                "adjusted_p",
                lambda p: int((p < SIBLING_ALPHA).sum()),
            ),
            gate_eligible_count=("gate_eligible", "sum"),
            production_gate_available=("production_gate_p", lambda p: int(p.notna().sum())),
            production_gate_rejections=(
                "production_gate_p",
                lambda p: int((p < SIBLING_ALPHA).sum()),
            ),
            full_span_count=("full_span", "sum"),
            median_raw_statistic=("raw_statistic", "median"),
            median_df=("raw_df", "median"),
            median_c_hat=("c_hat", "median"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / "restoration_root_summary.csv", index=False)
    (output_dir / "restoration_contract.json").write_text(
        json.dumps(
            dict(
                role="diagnostic_direct_attribution_not_production_calibration",
                edge_alpha=EDGE_ALPHA,
                sibling_alpha=SIBLING_ALPHA,
                branch_lengths=0,
                contrast_energy_rank=False,
                edge_roles="recomputed per selected/independent matrix with guarded covariance, frozen across arms",
                projection="same native production all-leaf parent-total chart in both covariance arms; MP rank and tie expansion remain adaptive",
                covariance_modes=["known_identity", "production_guarded"],
                inflation="existing fit/decision with default enforce_support_thresholds=False; predictions on null-like rows descriptive only and withheld from production_gate_p",
                conditional_selected_law="unresolved; no traversal/BH/error-control claim",
                node_count_unit="dependent within dataset; use root independent dataset count",
            ),
            indent=2,
        )
        + "\n"
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(write_panel(run_panel(args.replicates), args.output_dir).to_string(index=False))
