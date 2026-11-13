"""Reconstruct paired spectral comparisons and verify shared input identities."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reports.subtree_spectral_alignment_20261109.alternatives.run import dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    panel = HERE / "parent_child"
    with gzip.open(panel / "rows.jsonl.gz", "rt") as stream:
        rows = [json.loads(line) for line in stream]
    alternatives = pd.read_csv(HERE / "alternatives/rows.csv.gz")
    assert len(rows) == 160 and len(alternatives) == 5120
    assert len({(r["case"], r["replicate"], r["operator"]) for r in rows}) == 160
    assert all(row["status"] == "ok" for row in rows)
    a, children = [], []
    for row in rows:
        x, _truth, seed = dataset(row["case"], row["replicate"])
        h = hashlib.sha256(x.tobytes()).hexdigest()
        assert h == row["data_hash"] and seed == row["seed"]
        other = alternatives[
            (alternatives.regime == row["case"]) & (alternatives.replicate == row["replicate"])
        ]
        assert len(other) == 64 and other.data_sha256.eq(h).all() and other.seed.eq(seed).all()
        r = row["result"]
        core = {key: row[key] for key in ("case", "replicate", "operator")}
        d, blocks = r["diffusion"], r["diffusion"]["whole_blocks"]
        a.append(
            core
            | {
                "native_ami": r["native"]["ami"],
                "groups": r["native"]["found_clusters"],
                "root_p": r["native"]["root_sibling_p"],
                "root_reject": r["native"]["root_sibling_reject"],
                "raw_left_right_gap": max(
                    b["alignments"]["right_vs_left"]["projector_frobenius"] for b in blocks
                ),
                "transported_gap": max(
                    b["alignments"]["symmetric_vs_degree_right"]["projector_frobenius"]
                    for b in blocks
                ),
                "production_distance_error": d["production_distance_max_abs_error"],
                "covariance_error": r["covariance"]["decomposition_error"],
                "between_fraction": r["covariance"]["between_trace"]
                / (r["covariance"]["within_trace"] + r["covariance"]["between_trace"]),
                "raw_mp": r["native_context"][r["native_nodes"][0]]["raw_mp_signal_counts_by_node"],
            }
        )
        for index, child in enumerate(r["children"]):
            assert child["status"] == "ok"
            comp = child["rebuilt_child_graph"]["comparisons"]
            children.append(
                core
                | {
                    "child": index,
                    "n": child["n"],
                    "raw_covariance_gap": child["covariance"]["leading_whole_block_comparison"][
                        "projector_frobenius"
                    ],
                    "raw_eigenvalue_gap": child["covariance"]["eigenvalue_relative_difference"],
                    "right_gap": comp["right_vectors"]["projector_frobenius"],
                    "left_gap": comp["left_vectors"]["projector_frobenius"],
                    "symmetric_gap": comp["symmetric_vectors"]["projector_frobenius"],
                    "rank_loss": comp["right_vectors"]["rank_loss"],
                    "operator_gap": child["induced_vs_rebuilt_operator_frobenius"],
                }
            )
    roots, child = pd.DataFrame(a), pd.DataFrame(children)
    roots.to_csv(panel / "root_summary.csv", index=False)
    child.to_csv(panel / "child_summary.csv", index=False)
    assert len(child) == 320
    sources = json.loads((panel / "manifest.json").read_text())["source_hashes"]
    changed = [
        path
        for path, value in sources.items()
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != value
    ]
    assert not changed, changed
    result = {
        "native_fits": len(rows),
        "child_comparisons": len(child),
        "alternative_rows": len(alternatives),
        "unique_shared_datasets": len({r["data_hash"] for r in rows}),
        "all_input_hashes_regenerated_and_paired": True,
        "source_files_changed": changed,
        "maximum_errors": roots[
            ["transported_gap", "production_distance_error", "covariance_error"]
        ]
        .max()
        .to_dict(),
        "child_rank_loss_count": int(child.rank_loss.gt(0).sum()),
        "induced_rebuilt_gap_above_1e_10": int(child.operator_gap.gt(1e-10).sum()),
        "operator_gap_range": [float(child.operator_gap.min()), float(child.operator_gap.max())],
        "native_grouped": roots.groupby(["case", "operator"])
        .agg(
            mean_ami=("native_ami", "mean"),
            mean_groups=("groups", "mean"),
            root_rejections=("root_reject", "sum"),
            median_raw_left_right_gap=("raw_left_right_gap", "median"),
            median_between_fraction=("between_fraction", "median"),
        )
        .reset_index()
        .to_dict("records"),
        "child_grouped": child.groupby(["case", "operator"])
        .agg(
            median_raw_covariance_gap=("raw_covariance_gap", "median"),
            median_raw_eigenvalue_gap=("raw_eigenvalue_gap", "median"),
            median_right_gap=("right_gap", "median"),
            median_operator_gap=("operator_gap", "median"),
            positive_rank_loss=("rank_loss", lambda x: int(x.gt(0).sum())),
        )
        .reset_index()
        .to_dict("records"),
    }
    assert all(np.isfinite(value) for value in result["maximum_errors"].values())
    (HERE / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in result.items() if not key.endswith("grouped")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
