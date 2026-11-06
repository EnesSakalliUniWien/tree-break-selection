"""Compare this recorded full-suite rerun with its unchanged baseline."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_dir", type=Path)
run_dir = parser.parse_args().run_dir.resolve()
repo = Path(__file__).resolve().parents[2]
manifest = json.loads((run_dir / "matched-run-manifest.json").read_text())
completion = json.loads((run_dir / "completion.json").read_text())
assert completion["exit_code"] == 0, completion
baseline_path = Path(manifest["baseline"])
assert hashlib.sha256(baseline_path.read_bytes()).hexdigest() == manifest["baseline_sha256"]
for name, digest in manifest["source_sha256"].items():
    assert hashlib.sha256((repo / name).read_bytes()).hexdigest() == digest, name
baseline = pd.read_csv(baseline_path)
current_path = run_dir / "full_benchmark_comparison.csv"
current = pd.read_csv(current_path)
keys = ["case_id", "run_id"]
assert len(baseline) == len(current) == manifest["expected_rows"]
paired = baseline.merge(
    current,
    on=keys,
    suffixes=("_before", "_after"),
    validate="one_to_one",
    how="outer",
    indicator=True,
)
assert paired["_merge"].eq("both").all()
matched_columns = [
    "method",
    "test_case",
    "params",
    "samples",
    "features",
    "true_clusters",
    "case_category",
    "benchmark_class",
    "benchmark_grid",
    "benchmark_repeat",
]
for column in matched_columns:
    assert paired[f"{column}_before"].fillna("").eq(paired[f"{column}_after"].fillna("")).all(), (
        column
    )
paired["ari_delta"] = paired["ari_after"] - paired["ari_before"]
paired["status_changed"] = paired["status_before"] != paired["status_after"]
paired["k_changed"] = ~np.isclose(
    paired["found_clusters_before"], paired["found_clusters_after"], rtol=0, atol=0, equal_nan=True
)
paired["ari_changed"] = ~np.isclose(
    paired["ari_before"], paired["ari_after"], rtol=0, atol=1e-12, equal_nan=True
)
selected_columns = keys + [
    "method_before",
    "case_category_before",
    "true_clusters_before",
    "status_before",
    "status_after",
    "found_clusters_before",
    "found_clusters_after",
    "ari_before",
    "ari_after",
    "ari_delta",
    "nmi_before",
    "nmi_after",
    "status_changed",
    "k_changed",
    "ari_changed",
]
paired[selected_columns].to_csv(run_dir / "matched_case_results.csv", index=False)
changed = paired.loc[
    paired[["status_changed", "k_changed", "ari_changed"]].any(axis=1), selected_columns
]
changed = changed.sort_values("ari_delta", ascending=False, na_position="last")
changed.to_csv(run_dir / "changed_case_results.csv", index=False)
summaries = []
for method, group in paired.groupby("method_before", sort=False):
    row = {"method": method, "total_cases": len(group)}
    for suffix in ("before", "after"):
        ok = group[f"status_{suffix}"].eq("ok")
        null = group["true_clusters_before"].eq(1) & ok
        row.update(
            {
                f"ok_{suffix}": int(ok.sum()),
                f"skip_{suffix}": int(group[f"status_{suffix}"].eq("skip").sum()),
                f"unsupported_{suffix}": int(group[f"status_{suffix}"].eq("unsupported").sum()),
                f"mean_ari_{suffix}": group.loc[ok, f"ari_{suffix}"].mean(),
                f"mean_nmi_{suffix}": group.loc[ok, f"nmi_{suffix}"].mean(),
                f"exact_k_{suffix}": int(
                    (ok & group[f"found_clusters_{suffix}"].eq(group["true_clusters_before"])).sum()
                ),
                f"null_cases_ok_{suffix}": int(null.sum()),
                f"null_split_{suffix}": int((null & group[f"found_clusters_{suffix}"].gt(1)).sum()),
            }
        )
    row["mean_ari_delta"] = row["mean_ari_after"] - row["mean_ari_before"]
    row["ari_improved"] = int(group["ari_delta"].gt(1e-12).sum())
    row["ari_degraded"] = int(group["ari_delta"].lt(-1e-12).sum())
    row["ari_changed"] = int(group["ari_changed"].sum())
    row["k_changed"] = int(group["k_changed"].sum())
    row["status_changed"] = int(group["status_changed"].sum())
    summaries.append(row)
summary = pd.DataFrame(summaries)
summary.to_csv(run_dir / "matched_method_summary.csv", index=False)
checks = {
    "paired_rows": len(paired),
    "matched_cases": current["case_id"].nunique(),
    "matched_run_cells": current["run_id"].nunique(),
    "matched_columns": matched_columns,
    "baseline_hash_unchanged": True,
    "source_hashes_unchanged": len(manifest["source_sha256"]),
    "changed_result_rows": len(changed),
    "ari_change_tolerance": 1e-12,
    "current_results_sha256": hashlib.sha256(current_path.read_bytes()).hexdigest(),
}
(run_dir / "comparison_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
columns = [
    "method",
    "ok_before",
    "ok_after",
    "mean_ari_before",
    "mean_ari_after",
    "mean_ari_delta",
    "exact_k_before",
    "exact_k_after",
    "ari_improved",
    "ari_degraded",
    "null_split_before",
    "null_split_after",
]
report = "# Matched full benchmark rerun — 2026-11-03\n\n"
report += "Baseline: `run_20261019_064234Z_full`, commit `8e0093d3`. All 1,342 case/run pairs and recorded parameters, seeds, sample/feature counts, true K and benchmark metadata match. Case/generator source and lockfile match the baseline commit; original-run matrix hashes were not recorded. New reconstructed input hashes are saved.\n\n"
report += summary[columns].to_markdown(index=False, floatfmt=".6f") + "\n\n"
report += "## Changed results\n\n"
report += (
    changed[
        [
            "case_id",
            "method_before",
            "ari_before",
            "ari_after",
            "found_clusters_before",
            "found_clusters_after",
        ]
    ].to_markdown(index=False, floatfmt=".6f")
    + "\n\n"
)
report += "## Interpretation limits\n\nThe current working tree includes both the self-distance correction and the prior empirical-inflation repair. This is an integrated implementation comparison, not an isolated causal estimate of the self-distance change. K-means and spectral clustering receive true K. Skipped cases are excluded from quality means and counted separately. Six suite null cases are not a calibrated Type-I-error study. Stage timings omit parts of diffusion/NNLS and do not provide comparable end-to-end method runtimes.\n"
(run_dir / "matched_comparison.md").write_text(report)
print(json.dumps(checks))
print(summary[columns].to_string(index=False))
