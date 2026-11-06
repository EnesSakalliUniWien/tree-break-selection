"""Replay the first Bernoulli NNLS false-split reversal from the factorial panel."""

import copy
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration import diffusion_factorial as factorial
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    selected_hierarchy_null_audit as audit,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    decide_empirical_null_calibration,
    fit_empirical_null_inflation_model,
)


def main():
    output = Path("raw/assets/diffusion-factorial-20261102")
    case = next(
        item for item in factorial.datasets(4, 0) if item[0] == "iid_bernoulli" and item[1] == 3
    )
    _, _, _, _, data, labels, space, _ = case
    distances = factorial.feature_family_diffusion_geometry(data, space).distance_condensed
    original = factorial.build_tree(data, distances, linkage_method="average").tree
    results, raw_rows = [], []
    for nnls in (0, 1):
        tree = copy.deepcopy(original)
        if nnls:
            fit = factorial.fit_fixed_topology_nnls_branch_lengths(
                tree,
                data,
                pair_sample_size=50000,
                random_state=0,
                solver_tolerance=1e-5,
                max_iterations=1000,
            )
            assert fit.applied_to_tree
        else:
            for edge in tree.edges:
                tree.edges[edge]["branch_length"] = 0.0
        tree.populate_node_distributions(data, feature_space=space)
        bundle = factorial.run_gate_annotation_pipeline(
            tree,
            tree.annotations_df,
            leaf_data=data,
            feature_space=space,
            edge_alpha=0.001,
            sibling_alpha=0.01,
            adaptive_projection_dimension_energy_fraction=0.9,
            edge_branch_length_variance_policy="normalized_branch_length" if nnls else "none",
        )
        edges = bundle.edge_gate_result.annotated_df
        spectral = bundle.edge_gate_result.spectral_context
        dimensions = audit.derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree, spectral_context=spectral
        )
        projections, eigenvalues = (
            audit.collect_parent_principal_component_inputs_for_sibling_tests(
                dimensions, spectral_context=spectral
            )
        )
        records, _ = audit.collect_sibling_pair_records(
            tree,
            edges,
            feature_space=space,
            sibling_projection_dimensions_from_edge_comparisons=dimensions,
            parent_principal_component_projections=projections,
            parent_principal_component_eigenvalues=eigenvalues,
            adaptive_projection_dimension_energy_fraction=0.9,
        )
        model = fit_empirical_null_inflation_model(records)
        root = factorial.tree_root(tree)
        record = next(r for r in records if r.parent == root)
        decision = decide_empirical_null_calibration(model, record)
        result, _ = factorial.evaluate(tree, data, labels, space, nnls=nnls, energy=1)
        np.testing.assert_allclose(decision.p_value, result["root_sibling_p"], rtol=1e-12)
        results.append(
            {
                "nnls": nnls,
                "root_record": asdict(record),
                "decision": asdict(decision),
                "run": result,
            }
        )
        raw_rows.extend({"nnls": nnls, **asdict(r)} for r in records)
    off, on = results
    np.testing.assert_allclose(
        off["root_record"]["stat"],
        on["root_record"]["stat"] * on["run"]["root_branch_factor"],
        rtol=1e-10,
    )
    assert off["run"]["clusters"] == 1 and on["run"]["clusters"] == 2
    (output / "nnls_reversal_trace.json").write_text(json.dumps(results, indent=2) + "\n")
    pd.DataFrame(raw_rows).to_csv(output / "nnls_reversal_records.csv", index=False)
    for item in results:
        print(
            item["nnls"],
            "stat",
            item["root_record"]["stat"],
            "c_hat",
            item["decision"]["c_hat"],
            "p",
            item["decision"]["p_value"],
        )


if __name__ == "__main__":
    main()
