"""Plot saved benchmark rows and assignments without refitting methods."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

from benchmarks.shared.util.method_sets import DEFAULT_METHODS

CASES = [
    "cancer_pc_hub_disconnected", "cancer_pc_nohub_disconnected",
    "cancer_pc_hub_connected", "cancer_pc_nohub_connected",
]
CASE_LABELS = ["Dense row\nDisconnected", "Ordinary row\nDisconnected", "Dense row\nConnected", "Ordinary row\nConnected"]
METHODS = ["tbs_family_diffusion", "tbs", "spectral", "kmeans", "leiden", "louvain", "dbscan", "optics", "hdbscan", "tbs_diffusion", "tbs_diffusion_adaptive_nnls"]
LABELS = ["Family-diffusion TBS + NNLS", "TBS (standardized Euclidean)", "Spectral (known K)*", "K-means (known K)*", "Leiden", "Louvain", "DBSCAN", "OPTICS", "HDBSCAN", "Hamming-diffusion TBS", "Adaptive Hamming TBS + NNLS"]
COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]


def save_figure(fig, output, name):
    fig.savefig(output / f"{name}.png", dpi=170, bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent
    panel = output / "generated/panel"
    rows = pd.read_csv(panel / "results.csv")
    assert set(rows.method) == set(DEFAULT_METHODS) == set(METHODS)
    assert len(rows) == 264 and not rows.duplicated(["case_id", "simulation_seed", "run_id"]).any()
    additional = rows[rows.seed_role.eq("additional")]
    summaries = []
    ari = np.full((len(METHODS), len(CASES)), np.nan)
    k_error = ari.copy()
    k_text = {}
    for i, method in enumerate(METHODS):
        for j, case in enumerate(CASES):
            group = additional[(additional.method == method) & (additional.case_id == case)]
            assert len(group) == 5
            good = group[group.status.eq("ok")]
            item = {"method": method, "case_id": case, "ok": len(good), "skip": int(group.status.eq("skip").sum())}
            if not good.empty:
                assert len(good) == 5
                ari[i, j] = good.ari.mean()
                k_error[i, j] = np.abs(good.found_clusters - 6).mean()
                low, high = good.found_clusters.min(), good.found_clusters.max()
                k_text[i, j] = str(int(low)) if low == high else f"{int(low)}–{int(high)}"
                item.update(mean_ari=good.ari.mean(), min_ari=good.ari.min(), max_ari=good.ari.max(), min_k=low, max_k=high)
            summaries.append(item)
    pd.DataFrame(summaries).to_csv(output / "additional_seed_summary.csv", index=False)
    rows[rows.seed_role.eq("development")].to_csv(output / "development_results.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"width_ratios": [1, 1]})
    fig.subplots_adjust(left=0.23, right=0.96, top=0.80, bottom=0.19, wspace=0.16)
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#ECEFF1")
    kmap = plt.get_cmap("YlOrBr").copy()
    kmap.set_bad("#ECEFF1")
    left = axes[0].imshow(np.ma.masked_invalid(ari), vmin=-0.05, vmax=1, cmap=cmap, aspect="auto")
    right = axes[1].imshow(np.ma.masked_invalid(k_error), vmin=0, vmax=6, cmap=kmap, aspect="auto")
    for j, ax in enumerate(axes):
        ax.set_xticks(range(4), CASE_LABELS, fontsize=10)
        ax.xaxis.tick_top()
        ax.tick_params(length=0, pad=8)
        ax.set_yticks(range(len(METHODS)), LABELS if j == 0 else [""] * len(METHODS), fontsize=10)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(METHODS), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.5)
        ax.tick_params(which="minor", length=0)
        for row in range(len(METHODS)):
            for col in range(4):
                missing = np.isnan(ari[row, col])
                text = "Skipped" if missing else (f"{ari[row, col]:.3f}" if j == 0 else k_text[row, col])
                dark = not missing and (ari[row, col] > 0.6 if j == 0 else k_error[row, col] > 4)
                ax.text(col, row, text, ha="center", va="center", fontsize=11, color="white" if dark else "#17252F")
    axes[0].set_title("Mean adjusted Rand index", pad=53, fontsize=13, weight="bold")
    axes[1].set_title("Returned cluster count: range across seeds", pad=53, fontsize=13, weight="bold")
    fig.colorbar(left, ax=axes[0], orientation="horizontal", fraction=0.045, pad=0.035, label="ARI: larger means closer to the six generating populations")
    fig.colorbar(right, ax=axes[1], orientation="horizontal", fraction=0.045, pad=0.035, label="Color: mean absolute difference from six clusters")
    fig.suptitle("Cancer-inspired benchmark: five additional seeds", x=0.57, y=0.965, fontsize=19, weight="bold")
    fig.text(0.23, 0.065, "496 observations · 30 continuous features · six generating populations · unchanged registered parameters\n* Known K=6 is supplied to spectral clustering and K-means. Other methods infer K.\nCluster counts exclude noise (K=0 means all noise). Gray cells are unsupported Hamming inputs, not failures or ARI=0.", fontsize=10, linespacing=1.7)
    save_figure(fig, output, "method_comparison")

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.17, top=0.84, hspace=0.44, wspace=0.13)
    selected = [None, "tbs_family_diffusion", "spectral", "kmeans"]
    titles = ["Generating populations", "Family-diffusion TBS + NNLS", "Spectral (known K=6)", "K-means (known K=6)"]
    for i, case in enumerate(CASES):
        arrays = np.load(panel / "20261104" / case / "assignments.npz")
        x, truth = arrays["values"], arrays["truth"]
        for j, method in enumerate(selected):
            labels = truth if method is None else arrays[f"labels__{method}"]
            unique = np.unique(labels)
            counts = np.array([[np.sum((truth == t) & (labels == label)) for label in unique] for t in range(6)])
            true_ids, cluster_ids = linear_sum_assignment(-counts)
            color_map = {unique[c]: t for t, c in zip(true_ids, cluster_ids, strict=True)}
            assert len(color_map) == len(unique)
            color = np.array([color_map[label] for label in labels])
            ax = axes[i, j]
            ax.scatter(x[:, 0], x[:, 1], c=color, cmap=ListedColormap(COLORS), vmin=0, vmax=5, s=10, alpha=0.65, linewidths=0, rasterized=True)
            ax.scatter(x[0, 0], x[0, 1], s=100, facecolors="none", edgecolors="black", linewidths=1.5)
            score = adjusted_rand_score(truth, labels)
            ax.set_title(f"K={len(unique)}, ARI={score:.3f}", fontsize=10)
            if i == 0:
                position = ax.get_position()
                fig.text((position.x0 + position.x1) / 2, 0.885, titles[j], ha="center", fontsize=11, weight="bold")
            if j == 0:
                ax.set_ylabel(CASE_LABELS[i] + "\nSecond coordinate", fontsize=10)
            if i == 3:
                ax.set_xlabel("First coordinate")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(alpha=0.12)
    fig.suptitle("Where the populations go: development-seed assignments", fontsize=18, weight="bold", y=0.98)
    fig.text(0.10, 0.045, "Only the first two of 30 features are shown; these contain all population-mean signal. All 30 features enter clustering.\nColors are aligned to generating populations by maximum label overlap; the black ring marks the designated observation.\nConnected controls move the two isolated population means closer, changing both graph connectivity and population separation.", fontsize=10, linespacing=1.7)
    fig.legend([Line2D([], [], marker="o", linestyle="", color=c) for c in COLORS], [f"Population {i+1}" for i in range(6)], loc="upper center", bbox_to_anchor=(0.55, 0.945), ncol=6, frameon=False)
    save_figure(fig, output, "population_recovery")
    (output / "plot_manifest.json").write_text(json.dumps({"summary_seeds": list(range(20261105, 20261110)), "assignment_seed": 20261104, "methods": METHODS, "cases": CASES, "source": str(panel / "results.csv")}, indent=2) + "\n")


if __name__ == "__main__":
    main()
