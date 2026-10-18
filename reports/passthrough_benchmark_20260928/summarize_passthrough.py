"""Summarize matched passthrough traversal rows for the 2026-09-28 report."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from benchmarks.validation.contracts.report_contract import wilson_interval

ROOT = Path(__file__).resolve().parents[2]
EVAL = ROOT / "reports/passthrough_benchmark_20260928/evaluation"
traversal = pd.read_csv(EVAL / "known_branch_traversal_rows.csv")
local = pd.read_csv(EVAL / "known_branch_local_rows.csv")
expected_patterns = {
    "global_null": set(),
    "root_only": {"R"},
    "A_deep_only": {"A"},
    "separate_deep_branches": {"A", "B"},
    "strong_root_weak_A": {"R", "A"},
}
assert len(traversal) == 20_000 and traversal["status"].eq("ok").all()
assert len(local) == 120_000 and local["status"].eq("ok").all()
keys = ["pattern", "covariance_profile", "replicate"]
paired = traversal.groupby(keys)[["dataset_seed", "dataset_sha256"]].nunique()
assert paired.eq(1).all().all()
assert not traversal.duplicated(keys + ["passthrough"]).any()

traversal["any_false_split"] = traversal["n_reached_false_splits"].gt(0)
pattern_rows = []
node_rows = []
partition_rows = []
for (pattern, profile), frame in traversal.groupby(["pattern", "covariance_profile"], sort=True):
    wide_false = frame.pivot(index="replicate", columns="passthrough", values="any_false_split")
    upward = int((~wide_false[False] & wide_false[True]).sum())
    downward = int((wide_false[False] & ~wide_false[True]).sum())
    n = len(frame[frame["passthrough"].eq(False)])
    for passthrough, group in frame.groupby("passthrough", sort=True):
        k = int(group["any_false_split"].sum())
        lo, hi = wilson_interval(k, len(group))
        pattern_rows.append(
            {
                "pattern": pattern,
                "covariance_profile": profile,
                "passthrough": bool(passthrough),
                "n": len(group),
                "any_reached_false_split_count": k,
                "any_reached_false_split_rate": k / len(group),
                "false_split_wilson95_low": lo,
                "false_split_wilson95_high": hi,
                "mean_reached_false_split_fdp": float(group["reached_false_split_fdp"].mean()),
                "mean_final_partition_count": float(group["num_clusters"].mean()),
                "pass_through_R_count": int(
                    group["traversal_decisions"]
                    .map(lambda s: json.loads(s).get("R") == "pass_through")
                    .sum()
                ),
                "pass_through_A_count": int(
                    group["traversal_decisions"]
                    .map(lambda s: json.loads(s).get("A") == "pass_through")
                    .sum()
                ),
                "pass_through_B_count": int(
                    group["traversal_decisions"]
                    .map(lambda s: json.loads(s).get("B") == "pass_through")
                    .sum()
                ),
                "paired_false_split_up_when_enabled": upward if passthrough else 0,
                "paired_false_split_down_when_enabled": downward if passthrough else 0,
                "paired_rate_difference_enabled_minus_disabled": (upward - downward) / n
                if passthrough
                else 0.0,
            }
        )
        for node in expected_patterns[pattern]:
            decisions = group["traversal_decisions"].map(json.loads)
            reached = decisions.map(lambda d: node in d)
            split = decisions.map(lambda d: d.get(node) == "split")
            count = int(split.sum())
            nlo, nhi = wilson_interval(count, len(group))
            node_rows.append(
                {
                    "pattern": pattern,
                    "covariance_profile": profile,
                    "passthrough": bool(passthrough),
                    "signal_node": node,
                    "n_datasets": len(group),
                    "reached_count": int(reached.sum()),
                    "reach_rate": float(reached.mean()),
                    "split_count": count,
                    "unconditional_split_power": count / len(group),
                    "split_power_wilson95_low": nlo,
                    "split_power_wilson95_high": nhi,
                    "conditional_split_rate_given_reached": count / int(reached.sum())
                    if reached.any()
                    else float("nan"),
                    "pass_through_count": int(
                        decisions.map(lambda d: d.get(node) == "pass_through").sum()
                    ),
                }
            )
        for k_value in sorted(group["num_clusters"].unique()):
            partition_rows.append(
                {
                    "pattern": pattern,
                    "covariance_profile": profile,
                    "passthrough": bool(passthrough),
                    "final_partition_count": int(k_value),
                    "dataset_count": int(group["num_clusters"].eq(k_value).sum()),
                    "dataset_rate": float(group["num_clusters"].eq(k_value).mean()),
                }
            )

pd.DataFrame(pattern_rows).to_csv(EVAL / "passthrough_pattern_summary.csv", index=False)
pd.DataFrame(node_rows).to_csv(EVAL / "passthrough_signal_node_summary.csv", index=False)
pd.DataFrame(partition_rows).to_csv(EVAL / "passthrough_partition_summary.csv", index=False)
print(
    f"wrote pattern={len(pattern_rows)} signal_node={len(node_rows)} partition={len(partition_rows)} summary rows"
)
