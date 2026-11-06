"""Plot saved diagnostic outputs; no clustering fits."""

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path(__file__).parent
SOURCE = ROOT / "generated/complete"
factors = pd.read_csv(SOURCE / "root_factors.csv")
calibration = pd.read_csv(SOURCE / "calibration_factors.csv")
swaps = pd.read_csv(SOURCE / "frozen_tree_swaps.csv")
pairs = pd.read_csv(SOURCE / "cancer_pc_hub_disconnected-20261106-pairs.csv")
plt.rcParams.update(
    {"font.size": 10, "axes.titlesize": 11, "axes.spines.top": False, "axes.spines.right": False}
)
fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout="constrained")
fig.suptitle(
    "Cancer PC surrogate: different mechanisms produce rejection and stopping", fontsize=15
)
colors = ["#665E59", "#B68A62", "#BAA790", "#4B817E", "#225C62"]

ax = axes[0, 0]
order = [
    "native",
    "no_branch_time",
    "full_projection",
    "common_identity_covariance",
    "known_noise_covariance",
]
labels = [
    "Native",
    "Remove branch scaling",
    "All 30 projection dimensions",
    "Common identity covariance",
    "Known row noise covariance",
]
edge = (
    factors[
        (factors.case == "cancer_pc_hub_disconnected")
        & (factors.seed == 20261104)
        & (factors.node == "L0")
    ]
    .set_index("variant")
    .loc[order]
)
ratio = edge.q / chi2.isf(0.001, edge.df)
ax.barh(labels, ratio, color=colors)
ax.axvline(1, color="black", linestyle="--", linewidth=1)
ax.set_xscale("log")
ax.set_xlim(0.003, 30)
ax.invert_yaxis()
ax.set_xlabel("Q / raw chi-square cutoff at alpha = 0.001 (log scale)")
ax.set_title("A. Singleton root: covariance dominates\nHub/disconnected, seed 20261104", loc="left")
for i, value in enumerate(ratio):
    ax.text(value * 1.15, i, f"{value:.3g}", va="center", fontsize=9)

ax = axes[0, 1]
root = factors[
    (factors.case == "cancer_pc_hub_disconnected")
    & (factors.seed == 20261106)
    & (factors.gate == "sibling")
].set_index("variant")
clean = calibration[
    (calibration.case == "cancer_pc_hub_disconnected")
    & (calibration.seed == 20261106)
    & (calibration.variant == "known_mean_null_support_family_only")
].iloc[0]
native_c = root.loc["native", "c_hat"]
q = root.loc["native", "q"]
df = root.loc["native", "df"]
critical_c = q / chi2.isf(0.01, df)
ax.barh(
    ["Native family calibration", "Keep known mean-null support"],
    [native_c, clean.c_hat],
    color=[colors[0], colors[-1]],
)
ax.axvline(critical_c, color="black", linestyle="--", linewidth=1)
ax.invert_yaxis()
ax.set_xlim(0, 18)
ax.set_xlabel("Empirical-null scale c; root statistic held fixed")
ax.set_title(
    "B. Calibration closes the 16/480 root gate\n11 mean-different records supply 0.708\nof the calibration's weighted statistic",
    loc="left",
)
for i, (c, p) in enumerate([(native_c, root.loc["native", "p"]), (clean.c_hat, clean.p)]):
    ax.text(c + 0.2, i, f"c={c:.2f}\np={p:.3g}", va="center", fontsize=9)
ax.text(critical_c + 0.2, 0.47, f"c cutoff {critical_c:.2f}", fontsize=9)

ax = axes[1, 0]
columns = ["p_native", "p_known_noise_pca", "p_known_noise_fixed_axis", "p_known_noise_full_space"]
counts = [(pairs[c] <= 0.01).sum() for c in columns]
assert len(pairs) == 183 and counts[0] == 182
ax.barh(
    [
        "Native covariance + PCA",
        "Known noise, frozen native PCA",
        "Known noise, fixed first axis",
        "Known noise, full 30 dimensions",
    ],
    counts,
    color=colors[:4],
)
ax.invert_yaxis()
ax.set_xlim(0, 205)
ax.set_xlabel("Raw p <= 0.01 among 183 selected, same-mean leaf pairs")
ax.set_title(
    "C. Raw p-values at same-mean leaf pairs\nZero native sibling rejections among these pairs",
    loc="left",
)
for i, count in enumerate(counts):
    ax.text(count + 3, i, str(count), va="center")

ax = axes[1, 1]
names = ["cancer_pc_hub_disconnected", "cancer_pc_nohub_disconnected"]
matrix = (
    swaps[swaps.seed == 20261106]
    .pivot(index="tree_and_branch_source", columns="data_source", values="clusters")
    .loc[names, names]
    .to_numpy()
)
np.testing.assert_array_equal(matrix, [[1, 1], [6, 6]])
ax.imshow(matrix, cmap="Blues", vmin=0, vmax=7, aspect="auto")
ax.set_xticks([0, 1], ["Hub data", "Ordinary-row data"])
ax.set_yticks([0, 1], ["Hub tree + branch lengths", "Ordinary-row tree + branch lengths"])
for row in range(2):
    for col in range(2):
        ax.text(
            col,
            row,
            f"{matrix[row, col]} cluster" + ("s" if matrix[row, col] != 1 else ""),
            ha="center",
            va="center",
            color="white" if row else "black",
            fontsize=15,
        )
ax.set_title(
    "D. Frozen-tree swaps preserve the outcome\nSame cluster counts in both tested seeds",
    loc="left",
)
fig.supxlabel(
    "Diagnostic interventions use generator truth and fixed observed inputs; they are not validated replacement tests.",
    fontsize=10,
)
fig.savefig(ROOT / "statistics_diagnosis.png", dpi=180)
plt.close(fig)
paths = [
    Path(__file__),
    *[
        SOURCE / name
        for name in (
            "root_factors.csv",
            "calibration_factors.csv",
            "frozen_tree_swaps.csv",
            "cancer_pc_hub_disconnected-20261106-pairs.csv",
        )
    ],
]
(ROOT / "plot_manifest.json").write_text(
    json.dumps(
        {
            "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            "output": "statistics_diagnosis.png",
        },
        indent=2,
    )
    + "\n"
)
