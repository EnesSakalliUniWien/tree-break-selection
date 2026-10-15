"""Describe the frozen reassessment outputs without selecting a method or rule."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
METHOD_NAMES = {
    "tbs": "Average linkage",
    "tbs_complete": "Complete linkage",
    "tbs_neighbor_joining": "Neighbor joining",
    "tbs_diffusion": "Hamming diffusion",
    "tbs_diffusion_adaptive_nnls": "Adaptive diffusion + NNLS",
}
FAMILY_NAMES = {
    "bernoulli_iid": "Independent Bernoulli",
    "bernoulli_dependent": "Dependent Bernoulli",
    "categorical": "Categorical, three categories",
    "gaussian_identity": "Gaussian, identity covariance",
    "gaussian_ar1": "Gaussian, AR1 covariance",
    "gaussian_highdim": "Gaussian, 64 features / 40 samples",
}


def bootstrap_mean_ci(values, rng):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.nan, np.nan
    means = np.mean(values[rng.integers(0, len(values), size=(2000, len(values)))], axis=1)
    return tuple(np.quantile(means, [0.025, 0.975]))


def summarize_full_pipeline():
    path = HERE / "full_pipeline/evaluation"
    manifest = json.loads((path / "manifest.json").read_text())
    assert manifest["complete"]
    runs = pd.read_csv(path / "runs.csv", dtype={"seed": str}, low_memory=False)
    summaries = pd.read_csv(path / "summary.csv")
    assert len(runs) == 23800
    keys = ["family", "signal", "replicate", "method", "sibling_alpha"]
    assert not runs.duplicated(keys).any()
    dataset_keys = ["family", "signal", "replicate"]
    datasets = runs.groupby(dataset_keys)
    assert len(datasets) == 4200
    assert datasets.data_sha256.nunique().eq(1).all()
    assert datasets.seed.nunique().eq(1).all()
    assert runs.groupby("seed").ngroups == 4200
    null = summaries[(~summaries.signal) & summaries.sibling_alpha.eq(0.01)].copy()
    null.to_csv(HERE / "null_partition_summary.csv", index=False)
    # Descriptive bootstrap intervals resample independent datasets within a cell.
    # The same dataset enters paired method differences only once.
    rng = np.random.default_rng(2026092897)
    quality, matched = [], []
    primary = runs[runs.sibling_alpha.eq(0.01)]
    for family, frame in primary[primary.signal].groupby("family"):
        baseline = frame[frame.method.eq("tbs") & frame.status.eq("ok")]
        for method, values in frame.groupby("method"):
            ok = values[values.status.eq("ok")]
            audited = ok[ok.audit_error.fillna("").eq("")]
            low, high = bootstrap_mean_ci(ok.ari, rng)
            fdp_low, fdp_high = bootstrap_mean_ci(audited.split_fdp, rng)
            quality.append(
                {
                    "family": family,
                    "method": method,
                    "attempts": len(values),
                    "completed": len(ok),
                    "audited": len(audited),
                    "mean_ari": ok.ari.mean(),
                    "ari_bootstrap95_low": low,
                    "ari_bootstrap95_high": high,
                    "mean_clusters": ok.found_clusters.mean(),
                    "undercluster_count": int((ok.found_clusters < 2).sum()),
                    "overcluster_count": int((ok.found_clusters > 2).sum()),
                    "exact_partition_count": int(np.isclose(ok.ari, 1, atol=1e-12).sum()),
                    "mean_split_fdp": audited.split_fdp.mean(),
                    "split_fdp_bootstrap95_low": fdp_low,
                    "split_fdp_bootstrap95_high": fdp_high,
                    "mean_false_splits": audited.n_false_splits.mean(),
                    "mean_true_splits": (audited.n_splits - audited.n_false_splits).mean(),
                    "mean_pass_throughs": audited.pass_throughs.mean(),
                }
            )
            if method == "tbs":
                continue
            pair = baseline.merge(ok, on=dataset_keys, suffixes=("_baseline", "_method"))
            assert pair.data_sha256_baseline.eq(pair.data_sha256_method).all()
            differences = pair.ari_method - pair.ari_baseline
            low, high = bootstrap_mean_ci(differences, rng)
            matched.append(
                {
                    "family": family,
                    "baseline": "tbs",
                    "method": method,
                    "jointly_completed": len(pair),
                    "mean_ari_difference": differences.mean(),
                    "bootstrap95_low": low,
                    "bootstrap95_high": high,
                    "mean_cluster_difference": (
                        pair.found_clusters_method - pair.found_clusters_baseline
                    ).mean(),
                }
            )
    quality = pd.DataFrame(quality)
    quality.to_csv(HERE / "signal_recovery_summary.csv", index=False)
    pd.DataFrame(matched).to_csv(HERE / "paired_signal_comparisons.csv", index=False)
    alpha = summaries[
        summaries.family.eq("bernoulli_iid")
        & summaries.method.isin(
            [
                "tbs",
                "tbs_diffusion_adaptive_nnls",
            ]
        )
    ]
    alpha.to_csv(HERE / "alpha_sensitivity.csv", index=False)
    overview = {
        "independent_datasets": len(datasets),
        "method_alpha_attempts": len(runs),
        "statuses": runs.status.value_counts().to_dict(),
        "audit_errors": int(runs.audit_error.fillna("").ne("").sum()),
        "unavailable_reasons": runs.loc[runs.status.ne("ok"), "reason"].value_counts().to_dict(),
        "primary_null_cells": len(null),
        "primary_null_cells_with_results": int((null.completed > 0).sum()),
        "primary_null_cells_holm_excess_partitions": int((null.holm_p < 0.05).sum()),
        "kernel_adaptive_fraction_in_full_pipeline": None,
        "bootstrap_intervals": "Descriptive percentile 95%, 2000 dataset draws, seed2026092897; no simultaneous or selection claim",
    }
    (HERE / "analysis_summary.json").write_text(json.dumps(overview, indent=2) + "\n")
    return null, quality


def plot_null_partitions(null):
    colors = ["#35618d", "#5597a5", "#8781a8", "#c39145", "#aa5044"]
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.5), sharex=True, sharey=True)
    methods = list(METHOD_NAMES)
    for ax, (family, title) in zip(axes.flat, FAMILY_NAMES.items(), strict=True):
        frame = null[null.family.eq(family)].set_index("method")
        ax.axvline(0.01, color="#8f2935", linestyle="--", linewidth=1)
        for i, method in enumerate(methods):
            row = frame.loc[method]
            if row.completed == 0:
                ax.text(
                    0.1, i, "unsupported for this input", color="#777777", va="center", fontsize=9
                )
                continue
            value, low, high = (
                row.false_partition_rate,
                row.false_partition_wilson_low,
                row.false_partition_wilson_high,
            )
            ax.errorbar(
                value,
                i,
                xerr=[[max(0, value - low)], [max(0, high - value)]],
                fmt="o",
                color=colors[i],
                markersize=6,
                capsize=3,
                linewidth=1.4,
            )
            ax.text(
                high + 0.025 if value < 0.8 else low - 0.025,
                i,
                f"{int(row.false_partition_count)}/{int(row.completed)}",
                va="center",
                ha="left" if value < 0.8 else "right",
                fontsize=9,
            )
        ax.set_title(title, fontsize=11, loc="left", fontweight="bold")
        ax.set_yticks(range(len(methods)), [METHOD_NAMES[m] for m in methods], fontsize=9)
        ax.set_xlim(-0.03, 1.035)
        ax.set_ylim(4.5, -0.5)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
        ax.grid(axis="x", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Complete pipeline: false partitions under a one-population null",
        fontsize=16,
        x=0.02,
        ha="left",
    )
    fig.supxlabel("Fraction of completed datasets with final K > 1", fontsize=11, y=0.085)
    fig.text(
        0.02,
        0.025,
        "500 independent datasets per cell; points and Wilson 95% intervals. Dashed line: 0.01 diagnostic target.\n"
        "Every dataset refits the method. Missing/unsupported outcomes are not one-cluster results. Labels give events/completed datasets.",
        fontsize=9,
    )
    fig.tight_layout(rect=[0.015, 0.13, 0.99, 0.93], w_pad=2, h_pad=2)
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"null_partition_rates.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_wald():
    counts = pd.read_csv(HERE / "wald_rank/counts.csv")
    frame = counts[
        (counts.role == "held_out")
        & (counts.version == "current_complete_groups")
        & counts.selection.eq("adaptive_90")
        & counts.alpha.eq(0.01)
    ]
    # Fail visibly if the upstream schema differs; never silently plot no data.
    assert len(frame) == 6, (
        len(frame),
        counts.role.unique(),
        counts.version.unique(),
        counts.selection.unique(),
    )
    fig, ax = plt.subplots(figsize=(8, 4.8))
    labels = {
        "fully_tied_2": "(1, 1)\nfully tied",
        "distinct_2": "(2, 1)\ndistinct",
        "partial_tie_4": "(4, 2, 2, 1)\npartially tied",
    }
    ax.axhline(0.01, color="#8f2935", linestyle="--", label="Nominal alpha = 0.01")
    for seed_idx, (seed, group) in enumerate(frame.groupby("seed")):
        x = np.arange(3) + (-0.08 if seed_idx == 0 else 0.08)
        group = group.set_index("spectrum").loc[list(labels)]
        y = group.rejection_rate.to_numpy()
        ax.errorbar(
            x,
            y,
            yerr=[y - group.cp95_low.to_numpy(), group.cp95_high.to_numpy() - y],
            fmt="o",
            capsize=4,
            color=["#35618d", "#c39145"][seed_idx],
            label=f"Evaluation seed {seed}",
        )
    ax.set_xticks(range(3), list(labels.values()))
    ax.set_ylim(0.0088, 0.017)
    ax.set_ylabel("Null rejection probability")
    ax.set_title(
        "Complete-group adaptive Wald rank remains anti-conservative", loc="left", fontsize=13
    )
    ax.grid(axis="y", alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.text(
        0.1,
        0.015,
        "One million independent Gaussian contrasts per point; exact pointwise binomial 95% intervals.\n"
        "Fixed basis, energy fraction 0.90. No tree or empirical inflation is included.",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.1, 1, 1])
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"wald_rank_calibration.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    null, _quality = summarize_full_pipeline()
    plot_null_partitions(null)
    plot_wald()
