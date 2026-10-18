"""Is edge rejection on null data caused by fitting the PCA basis to the tested leaves?"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

import tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence.child_parent_divergence_annotation.child_parent_divergence_annotation as edge_module
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.tree.construction.build import build_tree

warnings.filterwarnings("ignore")

original = edge_module.compute_child_parent_spectral_context
state = {"basis_data": None}


def patched(tree, leaf_data, **kwargs):
    basis = leaf_data if state["basis_data"] is None else state["basis_data"]
    return original(tree, basis, **kwargs)


edge_module.compute_child_parent_spectral_context = patched
original_tests = edge_module.run_child_parent_tests_across_tree
captured_tests = {}


def capture_tests(**kwargs):
    stats, dfs, ps = original_tests(**kwargs)
    captured_tests.update(children=list(kwargs["child_ids"]), stats=stats, dfs=dfs, ps=ps)
    return stats, dfs, ps


edge_module.run_child_parent_tests_across_tree = capture_tests


def edge_rates(data, tree_distance, basis_data):
    tree = build_tree(data, tree_distance, linkage_method="average").tree
    tree.populate_node_distributions(data)
    state["basis_data"] = basis_data
    bundle = run_gate_annotation_pipeline(
        tree, tree.annotations_df.copy(), leaf_data=data, edge_alpha=0.001, sibling_alpha=0.01
    )
    a = bundle.annotated_df
    tested = a["Child_Parent_Divergence_Tested"].fillna(False).astype(bool)
    sig = a["Child_Parent_Divergence_Significant"].fillna(False).astype(bool)
    raw = pd.to_numeric(a["Child_Parent_Divergence_P_Value"], errors="coerce") < 0.001
    leaf = (
        pd.Series({n: tree.out_degree(n) == 0 for n in tree.nodes})
        .reindex(a.index)
        .fillna(False)
        .astype(bool)
    )
    stat = pd.to_numeric(a["Child_Parent_Divergence_Test_Statistic"], errors="coerce")
    df = pd.to_numeric(a["Child_Parent_Divergence_df"], errors="coerce")
    all_p = np.asarray(captured_tests["ps"])
    all_ratio = np.asarray(captured_tests["stats"]) / np.asarray(captured_tests["dfs"])
    is_leaf_edge = np.array([tree.out_degree(c) == 0 for c in captured_tests["children"]])
    unblocked = {
        "all_edges": len(all_p),
        "all_raw_reject": float((all_p < 0.001).mean()),
        "all_leaf_raw_reject": float((all_p[is_leaf_edge] < 0.001).mean()),
        "all_internal_raw_reject": float((all_p[~is_leaf_edge] < 0.001).mean()),
        "all_mean_stat_over_df": float(np.nanmean(all_ratio)),
    }
    return {
        **unblocked,
        "tested": int(tested.sum()),
        "rejected_bh": float(sig[tested].mean()),
        "raw_p_below_0.001": float(raw[tested].mean()),
        "leaf_raw_reject": float(raw[tested & leaf].mean()),
        "internal_raw_reject": float(raw[tested & ~leaf].mean()),
        "mean_stat_over_df": float((stat[tested] / df[tested]).mean()),
    }


rng = np.random.default_rng(7)
case = {c["name"]: c for c in get_test_cases_by_suite("full")}["gauss_null_small"]
real = prepare_case_inputs(case, ["tbs"]).data.astype(float)
datasets = {"gauss_null_small": real}
for r in range(2):
    datasets[f"bernoulli_null_{r}"] = pd.DataFrame(
        (rng.random(real.shape) < 0.5).astype(float), index=real.index, columns=real.columns
    )
rows = []
for name, data in datasets.items():
    independent = pd.DataFrame(
        (rng.random(data.shape) < data.mean().to_numpy()).astype(float),
        index=data.index,
        columns=data.columns,
    )
    random_tree_distance = pdist(rng.normal(size=(len(data), 5)))
    for tree_kind, dist in (
        ("data_tree", pdist(data.to_numpy(), "hamming")),
        ("random_tree", random_tree_distance),
    ):
        for basis_kind, basis in (("same_data_basis", None), ("independent_basis", independent)):
            rows.append(
                {
                    "data": name,
                    "tree": tree_kind,
                    "basis": basis_kind,
                    **edge_rates(data, dist, basis),
                }
            )
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(Path(__file__).parent / "edge_basis_experiment.csv", index=False)
