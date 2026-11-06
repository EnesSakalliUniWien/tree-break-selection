"""Audit calibration aggregation, correction arithmetic and traversal boundaries."""

import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from tree_break_selection.hierarchy_analysis.decomposition.gates.gate_evaluator import (
    GateEvaluator,
    TraversalDecision,
)
from tree_break_selection.hierarchy_analysis.statistics.multiple_testing.tree_bh import (
    apply_tree_bh_correction,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflated_projected_wald_annotation.fdr_annotation import (
    apply_traversal_aligned_sibling_bh_results,
    init_sibling_annotation_df,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    fit_empirical_null_inflation_model,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.types.sibling_pair_record import (
    SiblingPairRecord,
)

tree = nx.DiGraph(
    [("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2")]
)
ids = ["A", "B", "A1", "A2", "B1", "B2"]
adjustment = apply_tree_bh_correction(
    tree, np.array([0.001, 0.5, 0.012, 0.03, 1e-9, 1e-9]), ids, 0.05
)
np.testing.assert_allclose(
    adjustment.child_parent_edge_corrected_p_values_by_tree_bh, [0.002, 0.5, 0.024, 0.03, 1, 1]
)
assert adjustment.child_parent_edge_null_rejected_by_tree_bh.tolist() == [
    True,
    False,
    True,
    False,
    False,
    False,
]

edge_flags = {node: node != "root" for node in tree}
annotations = init_sibling_annotation_df(
    pd.DataFrame({"Child_Parent_Divergence_Significant": pd.Series(edge_flags)})
)
out = apply_traversal_aligned_sibling_bh_results(
    tree, annotations, ["root", "A", "B"], [(1, 1, 0.9), (10, 1, 0.001), (1, 1, 0.9)], 0.01
)
assert not out.loc["root", "Sibling_BH_Different"] and out.loc["A", "Sibling_BH_Different"]
children = {node: list(tree.successors(node)) for node in tree}
different = out.Sibling_BH_Different.to_dict()
skipped = {node: False for node in tree}
plain = GateEvaluator(tree, edge_flags, different, skipped, children)
passthrough = GateEvaluator(tree, edge_flags, different, skipped, children, passthrough=True)
assert plain.decision("root") is TraversalDecision.BOUNDARY
assert passthrough.decision("root") is TraversalDecision.PASS_THROUGH
closed = GateEvaluator(
    tree, edge_flags | {"A": False, "B": False}, different, skipped, children, passthrough=True
)
assert closed.decision("root") is TraversalDecision.BOUNDARY


def support_record(parent, group, statistic):
    return SiblingPairRecord(
        parent=parent,
        left=parent + "L",
        right=parent + "R",
        stat=statistic,
        reference_scale=1,
        degrees_of_freedom=30,
        p_value=0.5,
        branch_length_sum=0,
        n_parent=16,
        is_null_like=True,
        is_edge_blocked=True,
        sibling_null_weight=1,
        calibration_dependency_group=group,
        calibration_dependency_weight=1,
        sibling_projection_dimension=30,
        parent_positive_eigenvalue_count=30,
        feature_family="continuous",
    )


two = fit_empirical_null_inflation_model(
    [support_record("p1", "g1", 30), support_record("p2", "g2", 300)]
)
eleven = fit_empirical_null_inflation_model(
    [support_record(f"p{i}", "g1", 30) for i in range(10)] + [support_record("p10", "g2", 300)]
)
np.testing.assert_allclose(two.baseline_scale_estimate, 5.5)
np.testing.assert_allclose(eleven.baseline_scale_estimate, 20 / 11)
np.testing.assert_allclose([two.effective_sample_size, eleven.effective_sample_size], [2, 2])

# Mean equality at an ancestor need not be the intersection of descendant equalities.
means = {"A1": -1, "A2": 1, "B1": 0, "B2": 0}
for node in ("A", "B", "root"):
    means[node] = np.mean([means[child] for child in children[node]])
assert means["A"] == means["root"] and means["A1"] != means["A"] and means["A2"] != means["A"]

# Equal, invalid root p-values U² pass the correction iff U² <= alpha.
root_only = nx.DiGraph([("root", "L"), ("root", "R")])
for u in np.linspace(0.01, 0.99, 99):
    result = apply_tree_bh_correction(root_only, np.array([u * u, u * u]), ["L", "R"], 0.05)
    assert result.child_parent_edge_null_rejected_by_tree_bh.all() == (u * u <= 0.05)

report = {
    "tree_bh_family_arithmetic": "passed",
    "sibling_frontier_reaches_A_despite_root_sibling_p_0_9": True,
    "root_without_passthrough": plain.decision("root").value,
    "root_with_passthrough": passthrough.decision("root").value,
    "closed_edge_blocks_passthrough": closed.decision("root").value,
    "calibration_two_records_c": two.baseline_scale_estimate,
    "calibration_eleven_records_c": eleven.baseline_scale_estimate,
    "dependency_group_ess_both_models": [two.effective_sample_size, eleven.effective_sample_size],
    "mean_null_is_not_descendant_null_intersection": means,
    "invalid_p_uniform_squared_root_false_rejection_probability": float(np.sqrt(0.05)),
    "scope": "Deterministic code-boundary and algebra checks; not a full-pipeline error estimate.",
}
destination = Path(__file__).parent / "generated/downstream"
destination.mkdir(parents=True, exist_ok=True)
(destination / "probe.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
