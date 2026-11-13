"""Rescore saved CHOIR and four-seed imbalance assignments by descriptive profiles."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
import sklearn

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.shared.metrics import matched_group_recovery
from benchmarks.shared.statistical_profiles import build_statistical_profiles

ANALYSES = ROOT / "results/analyses"
FOLLOWUP = ANALYSES / "root_validation_followup_20261109/metrics"
AUDIT = ANALYSES / "three_tests_audit_20261109"
CHOIR = ANALYSES / "choir_cell_lines_tbs_lbfgsb_20261106"
RULE_OUT = ANALYSES / "ceiling_rule_real_20261109"
TIES = ANALYSES / "subspace_distance_ties_20261108/real_splits.jsonl"
METHODS = {
    "tbs_family_diffusion": "symmetric_average",
    "tbs_family_diffusion_markov_ward": "markov_ward",
}
ARMS = ("A", "B", "C", "rule_out")
INPUTS = {
    Path(__file__),
    ROOT / "benchmarks/shared/metrics.py",
    ROOT / "benchmarks/shared/statistical_profiles.py",
}


def read(path, **kwargs):
    INPUTS.add(path)
    return pd.read_csv(path, **kwargs)


def check_previous_match(current, old, keys):
    fields = ["precision", "recall", "f1", "rare", "matched_90_90", "legacy_90_90"]
    check = current[keys + fields].merge(old[keys + fields], on=keys, validate="one_to_one")
    assert len(check) == len(current) == len(old)
    for field in fields:
        np.testing.assert_allclose(check[f"{field}_x"], check[f"{field}_y"], atol=1e-12)


def synthetic():
    source = FOLLOWUP / "metrics.py"
    INPUTS.add(source)
    spec = importlib.util.spec_from_file_location("previous_imbalance_profile", source)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    assigned = read(FOLLOWUP / "fresh_assignments.csv.gz").query(
        "regime == 'imbalanced' and not perturbed"
    )
    old = read(FOLLOWUP / "fresh_groups.csv").query("regime == 'imbalanced' and not perturbed")
    profiles, pairs, recovery = [], [], []
    for seed, panel in assigned.groupby("seed", sort=True):
        x, truth = audit.simulate(*audit.REGIMES["imbalanced"], seed)
        profile = build_statistical_profiles(x, true_labels=truth)
        profiles.append(profile.groups.assign(seed=seed))
        pairs.append(profile.pairs.assign(seed=seed))
        for (method, arm), fit in panel.groupby(["method", "arm"], sort=True):
            assert np.array_equal(fit.cell, np.arange(len(x)))
            assert np.array_equal(fit.truth, truth)
            recovery.append(
                matched_group_recovery(truth, fit.cluster, rare_labels=[3]).assign(
                    seed=seed, method=method, arm=arm
                )
            )
    groups = pd.concat(profiles, ignore_index=True)
    matched = pd.concat(recovery, ignore_index=True)
    check_previous_match(matched, old, ["seed", "method", "arm", "group"])
    pd.concat(pairs, ignore_index=True).to_csv(HERE / "synthetic_pair_profiles.csv", index=False)
    joined = matched.merge(groups, on=["seed", "group"], validate="many_to_one")
    assert len(joined) == 128
    joined.to_csv(HERE / "synthetic_profile_recovery.csv", index=False)
    joined.groupby(["method", "arm", "group"], sort=True).agg(
        seeds=("seed", "size"),
        size=("size", "first"),
        mean_f1=("f1", "mean"),
        mean_precision=("precision", "mean"),
        mean_recall=("recall", "mean"),
        mean_effective_rank=("effective_rank", "mean"),
        mean_nearest_centroid_distance=("nearest_centroid_distance", "mean"),
    ).to_csv(HERE / "synthetic_profile_summary.csv")


def children_leaves(linkage, node):
    n = len(linkage) + 1

    def leaves(start):
        pending, out = [start], []
        while pending:
            i = pending.pop()
            if i < n:
                out.append(i)
            else:
                pending.extend(linkage[i - n, :2].astype(int))
        return out

    a, b = linkage[node - n, :2].astype(int)
    return leaves(a), leaves(b)


def choir():
    meta = read(
        CHOIR / "cell_metadata.csv",
        usecols=["cell_id", "Cell_line"],
        index_col="cell_id",
        dtype=str,
    )
    points_path = CHOIR / "pca.npy"
    INPUTS.add(points_path)
    x = np.load(points_path, mmap_mode="r")
    assert x.shape == (48_879, 30) and meta.index.is_unique and len(meta) == len(x)
    truth = meta.Cell_line.to_numpy()
    groups = build_statistical_profiles(x, true_labels=truth).groups
    matched = []
    subtrees, geometry, subtree_index = {}, {}, []
    prior_nodes = read(FOLLOWUP / "choir_split_relationships_annotated.csv")
    INPUTS.add(TIES)
    ties = pd.read_json(TIES, lines=True)
    selected = prior_nodes.loc[
        prior_nodes.test_tie.notna(), ["tree", "node", "n", "n_small", "B_pass", "C_pass"]
    ]
    selected = selected.merge(
        ties[["tree", "node", "n", "n_small", "tree_test_disagree"]],
        on=["tree", "node", "n", "n_small"],
        validate="one_to_one",
    )
    assert len(selected) == 209
    for method, tree in METHODS.items():
        for arm in ARMS:
            base = RULE_OUT if arm == "rule_out" else AUDIT
            fit = read(base / method / f"assignments_{arm}.csv.gz")
            assert np.array_equal(fit.cell, meta.index)
            assert np.array_equal(fit.Cell_line, truth)
            matched.append(
                matched_group_recovery(truth, fit.cluster).assign(method=method, arm=arm)
            )
        linkage_path = CHOIR / method / "linkage.npy"
        INPUTS.add(linkage_path)
        linkage = np.load(linkage_path, mmap_mode="r")
        assert len(linkage) == len(x) - 1
        for node in selected.loc[selected.tree == tree].itertuples(index=False):
            name = f"{tree}:{node.node}"
            a, b = children_leaves(linkage, node.node)
            assert len(a) + len(b) == node.n and min(len(a), len(b)) == node.n_small
            subtrees[name] = (a, b)
            geometry[name] = node.tree_test_disagree
            subtree_index.append(
                {
                    "subtree": name,
                    "tree": tree,
                    "node": node.node,
                    "B_pass": node.B_pass,
                    "C_pass": node.C_pass,
                }
            )
    recovery = pd.concat(matched, ignore_index=True)
    check_previous_match(recovery, read(FOLLOWUP / "choir_groups.csv"), ["method", "arm", "group"])
    result = recovery.merge(groups, on="group", validate="many_to_one")
    assert len(groups) == 190 and len(result) == 1520
    groups.to_csv(HERE / "choir_group_profiles.csv", index=False)
    result.drop(columns="covariance_eigenvalues").to_csv(
        HERE / "choir_profile_recovery.csv", index=False
    )

    geometry_rows = build_statistical_profiles(
        x, subtrees=subtrees, geometry_disagreement=geometry
    ).subtrees.merge(pd.DataFrame(subtree_index), on="subtree", validate="one_to_one")
    geometry_rows.to_csv(HERE / "choir_subtree_profiles.csv", index=False)
    rows = []
    for feature in [
        "size",
        "nearest_mahalanobis_distance",
        "effective_rank",
        "nearest_relative_covariance_difference",
        "nearest_centroid_overlap_share",
    ]:
        cuts = pd.qcut(groups[feature], q=4, duplicates="drop")
        strata = pd.Series(cuts.cat.codes.to_numpy() + 1, index=groups.group)
        strata[strata == 0] = -1
        view = result.assign(stratum=result.group.map(strata))
        summary = view.groupby(["method", "arm", "stratum"], sort=True).agg(
            groups=("group", "size"),
            mean_f1=("f1", "mean"),
            mean_precision=("precision", "mean"),
            mean_recall=("recall", "mean"),
            min_profile=(feature, "min"),
            max_profile=(feature, "max"),
        )
        rows.append(summary.reset_index().assign(profile=feature))
    pd.concat(rows, ignore_index=True).to_csv(HERE / "choir_profile_strata.csv", index=False)


def main():
    synthetic()
    choir()
    (HERE / "manifest.json").write_text(
        json.dumps(
            {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
                "sklearn": sklearn.__version__,
                "sha256": {
                    str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in sorted(INPUTS)
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
