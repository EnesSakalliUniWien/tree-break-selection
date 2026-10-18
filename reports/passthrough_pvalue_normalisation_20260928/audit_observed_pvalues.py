"""Separate unchanged node p-values from changed traversal selection."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
INPUT = BASE.parent / "passthrough_benchmark_20260928/evaluation/known_branch_traversal_rows.csv"


def main():
    rows = pd.read_csv(INPUT)
    keys = ["pattern", "covariance_profile", "replicate"]
    assert len(rows) == 20_000 and rows.status.eq("ok").all()
    assert not rows.duplicated(keys + ["passthrough"]).any()
    pairs = rows.groupby(keys)
    assert pairs.size().eq(2).all()
    assert pairs[["dataset_sha256", "injected_local_tests"]].nunique().eq(1).all().all()
    records = []
    for row in rows.itertuples():
        tests = json.loads(row.injected_local_tests)
        trace = json.loads(row.traversal_decisions)
        names = ["R", "A", "B"]
        p = np.array([tests[node]["p_value"] for node in names])
        order = np.argsort(p[1:]) + 1
        q = p.copy()
        q[order[0]] = min(2 * p[order[0]], p[order[1]])
        for index, node in enumerate(names):
            records.append(
                {
                    **{key: getattr(row, key) for key in keys},
                    "passthrough": row.passthrough,
                    "node": node,
                    "raw_p": p[index],
                    "bh_p": q[index],
                    "reached": node in trace,
                    "decision": trace.get(node, "unreached"),
                }
            )
    nodes = pd.DataFrame(records)
    paired_nodes = nodes.groupby(keys + ["node"])
    assert paired_nodes[["raw_p", "bh_p"]].nunique().eq(1).all().all()
    summaries = []
    for key, group in nodes.groupby(["pattern", "covariance_profile", "passthrough", "node"]):
        reached = group[group.reached]
        summaries.append(
            {
                **dict(
                    zip(["pattern", "covariance_profile", "passthrough", "node"], key, strict=True)
                ),
                "datasets": len(group),
                "reached_count": len(reached),
                "split_count": int(group.decision.eq("split").sum()),
                "all_raw_p_mean": group.raw_p.mean(),
                "reached_raw_p_mean": reached.raw_p.mean(),
                "reached_raw_p_median": reached.raw_p.median(),
                "reached_bh_p_mean": reached.bh_p.mean(),
            }
        )
    pd.DataFrame(summaries).to_csv(BASE / "node_reach_pvalue_summary.csv", index=False)
    changes = []
    for (pattern, profile), group in rows.groupby(keys[:2]):
        wide = group.pivot(index="replicate", columns="passthrough", values="num_clusters")
        reach = group.pivot(index="replicate", columns="passthrough", values="n_reached_internal")
        assert (wide[True] >= wide[False]).all()
        assert (reach[True] >= reach[False]).all()
        changes.append(
            {
                "pattern": pattern,
                "covariance_profile": profile,
                "matched_datasets": len(wide),
                "changed_partition_count": int((wide[True] != wide[False]).sum()),
                "newly_reached_node_occurrences": int((reach[True] - reach[False]).sum()),
                "mean_clusters_disabled": wide[False].mean(),
                "mean_clusters_enabled": wide[True].mean(),
            }
        )
    pd.DataFrame(changes).to_csv(BASE / "paired_partition_changes.csv", index=False)
    local_input = INPUT.parent / "known_branch_local_rows.csv"
    selected_input = INPUT.parent / "known_branch_selected_root_rows.csv"
    local = pd.read_csv(local_input)
    selected = pd.read_csv(selected_input)
    assert local.status.eq("ok").all() and selected.status.eq("ok").all()
    null_root = local[local.pattern.eq("global_null") & local.focal_node.eq("R")]
    factors = []
    for profile, group in null_root.groupby("covariance_profile"):
        baseline = (
            group[
                group.covariance_estimator.eq("oracle_covariance")
                & group.projection_kind.eq("fixed_full")
            ]
            .set_index("replicate")
            .sort_index()
        )
        alternatives = list(group.groupby(["covariance_estimator", "projection_kind"]))
        alternatives.append(
            (
                ("oracle_covariance", "selected_root_full"),
                selected[
                    selected.pattern.eq("global_null") & selected.covariance_profile.eq(profile)
                ],
            )
        )
        for (covariance, projection), other in alternatives:
            other = other.set_index("replicate").sort_index()
            assert baseline.index.equals(other.index)
            assert baseline.dataset_sha256.equals(other.dataset_sha256)
            before = baseline.p_value <= 0.01
            after = other.p_value <= 0.01
            factors.append(
                {
                    "covariance_profile": profile,
                    "covariance_estimator": covariance,
                    "projection_kind": projection,
                    "matched_datasets": len(other),
                    "baseline_rejections_at_0.01": int(before.sum()),
                    "rejections_at_0.01": int(after.sum()),
                    "new_rejections": int((~before & after).sum()),
                    "lost_rejections": int((before & ~after).sum()),
                    "median_p": other.p_value.median(),
                    "median_p_difference_from_oracle_full": (
                        other.p_value - baseline.p_value
                    ).median(),
                    "min_df": other.degrees_of_freedom.min(),
                    "max_df": other.degrees_of_freedom.max(),
                }
            )
    pd.DataFrame(factors).to_csv(BASE / "upstream_null_root_factor_summary.csv", index=False)
    result = {
        "passed": True,
        "matched_datasets": len(pairs),
        "matched_node_tests": len(paired_nodes),
        "changed_raw_p_values": 0,
        "changed_reconstructed_bh_p_values": 0,
        "changed_partition_count": sum(row["changed_partition_count"] for row in changes),
        "newly_reached_node_occurrences": sum(
            row["newly_reached_node_occurrences"] for row in changes
        ),
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "upstream_inputs_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (local_input, selected_input)
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Cached known-covariance Gaussian fixed-tree tests with open edges; p equality is by construction. BH is reconstructed independently for families {R} and {A,B}. Reached-p summaries are descriptive selections, not new independent calibration trials. Partition changes here are count changes; fixed three-node topology makes every additional descent refine the partition.",
    }
    (BASE / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
