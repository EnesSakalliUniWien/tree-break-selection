"""Reconstruct the research conclusions from completed, paired development rows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from benchmarks.diagnostics.calibration.subtree_research import FACTORIAL_AXES, RECOVERY_CELLS
from benchmarks.validation.sweeps.subtree_research_validation import paired_recovery_summary

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def interval(values):
    x = np.asarray(values, dtype=float)
    if len(x) < 2:
        return [None, None]
    mean = x.mean()
    half = t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return [float(mean - half), float(mean + half)]


def main():
    factorial = pd.concat(
        [pd.read_csv(HERE / f"factorial_cell_{cell}/runs.csv") for cell in range(3)],
        ignore_index=True,
    )
    recovery = pd.read_csv(HERE / "recovery/runs.csv")
    assert len(factorial) == 2400 and len(recovery) == 1600
    assert not factorial.duplicated(["cell", "replicate", "method"]).any()
    assert not recovery.duplicated(["cell", "replicate", "method"]).any()
    for frame in (factorial, recovery):
        assert frame.groupby(["cell", "replicate"]).data_hash.nunique().eq(1).all()
        assert frame.groupby(["cell", "replicate"]).truth_hash.nunique().eq(1).all()

    summary = (
        factorial.groupby(["cell", *FACTORIAL_AXES])
        .agg(
            attempted=("status", "size"),
            successful=("status", lambda x: int(x.eq("ok").sum())),
            root_rejections=("root_sibling_reject", "sum"),
            mean_ami=("ami", "mean"),
            mean_ari=("ari", "mean"),
            mean_matched_f1=("matched_macro_f1", "mean"),
            mean_groups=("found_clusters", "mean"),
            mean_split_share=("pair_split_share", "mean"),
            mean_merge_share=("pair_merge_share", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(HERE / "factorial_summary.csv", index=False)
    effects = []
    for cell, group in factorial.groupby("cell"):
        for factor, levels in FACTORIAL_AXES.items():
            for metric in ("ami", "matched_macro_f1", "found_clusters", "root_sibling_reject"):
                pivot = group.assign(
                    value=pd.to_numeric(group[metric], errors="coerce")
                ).pivot_table(index="replicate", columns=factor, values="value", aggfunc="mean")
                supported = group.groupby("replicate").apply(
                    lambda g: len(g) == 16 and g.status.eq("ok").all(), include_groups=False
                )
                delta = (pivot[levels[1]] - pivot[levels[0]]).loc[supported].dropna()
                effects.append(
                    dict(
                        cell=int(cell),
                        factor=factor,
                        level0=str(levels[0]),
                        level1=str(levels[1]),
                        metric=metric,
                        complete_datasets=len(delta),
                        mean_difference=float(delta.mean()),
                        interval=interval(delta),
                    )
                )
    (HERE / "factorial_main_effects.json").write_text(json.dumps(effects, indent=2) + "\n")

    recovery["phase"] = 1
    comparison = paired_recovery_summary(
        recovery, "symmetric_average", "markov_ward", expected_replicates=50
    )
    (HERE / "recovery_comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    grouped = (
        recovery.groupby(["cell", "method"])
        .agg(
            attempted=("status", "size"),
            successful=("status", lambda x: int(x.eq("ok").sum())),
            ami=("ami", "mean"),
            ari=("ari", "mean"),
            matched_macro_f1=("matched_macro_f1", "mean"),
            matched_rare_f1=("matched_rare_f1", "mean"),
            pair_split_share=("pair_split_share", "mean"),
            pair_merge_share=("pair_merge_share", "mean"),
            clusters=("found_clusters", "mean"),
        )
        .reset_index()
    )
    for index, name in enumerate(("features", "imbalanced", "gap", "anisotropic")):
        grouped[name] = grouped.cell.map(lambda c: RECOVERY_CELLS[c][index])
    grouped.to_csv(HERE / "recovery_summary.csv", index=False)

    initial = json.loads((HERE / "environment.json").read_text())["source_hashes"]
    production_changed = [
        p
        for p, h in initial.items()
        if p.startswith("tree_break_selection/")
        and hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h
    ]
    runtime = json.loads((HERE / "preflight_isolated/runtime_gate.json").read_text())
    decision = {
        "status": "not_promoted",
        "candidate": "bounded spectral Ward alternatives, common-target minimum NNLS residual, native empirical sibling reference",
        "runtime": runtime,
        "model_specific_candidate_reference": "not_established: common-event Ward family reference does not cover the new diffusion/minimum-residual procedure or its empirical eligibility/traversal",
        "family_power": "no true family rejections in 50 mixed-root development datasets at 999 steps",
        "frozen_calibration": {
            "status": "not_run_prerequisites_failed",
            "planned_datasets": 120000,
            "cells": 24,
            "datasets_per_cell": 5000,
        },
        "frozen_recovery": {
            "status": "not_run_prerequisites_failed",
            "planned_datasets": 3200,
            "cells": 16,
            "datasets_per_cell": 200,
        },
        "candidate_full_recovery": {
            "status": "not_run_runtime_preflight_failed",
            "baseline_recovery_fits_completed": 1600,
        },
        "production_files_changed": production_changed,
        "development": {
            "factorial_fits": len(factorial),
            "recovery_fits": len(recovery),
            "distinct_input_matrices": len(set(factorial.data_hash) | set(recovery.data_hash)),
            "factorial_status": factorial.status.value_counts().to_dict(),
            "recovery_status": recovery.status.value_counts().to_dict(),
        },
    }
    assert not production_changed
    (HERE / "promotion_gate.json").write_text(json.dumps(decision, indent=2) + "\n")
    print(
        json.dumps(
            {
                "decision": decision["status"],
                "development": decision["development"],
                "recovery_comparison": comparison["aggregate"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
