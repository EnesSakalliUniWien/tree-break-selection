"""Render the saved comparison without rerunning inference: python <this file>."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent
spectra = pd.read_csv(ROOT / "spectra.csv")
blocks = pd.read_csv(ROOT / "blocks.csv")
runs = pd.read_csv(ROOT / "runs.csv")
cases = pd.read_csv(ROOT / "cases.csv")
METHODS = [
    "average",
    "complete",
    "weighted",
    "single",
    "centroid",
    "median",
    "ward",
    "neighbor_joining",
]
COLORS = [
    "#3b638c",
    "#d0873b",
    "#559b84",
    "#ac6484",
    "#8775aa",
    "#91a851",
    "#79674b",
    "#3ba4ad",
    "#b36a61",
]


def spectrum(ax, case):
    values = spectra[spectra.case.eq(case)]
    selected = blocks[blocks.case.eq(case)]
    for i, row in enumerate(selected.itertuples()):
        part = values[values["mode"].between(row.block_start, row.block_end)]
        ax.plot(
            part["mode"],
            part.eigenvalue,
            ".-",
            color=COLORS[i],
            lw=1.5,
            label=f"Modes {row.block_start}–{row.block_end}",
        )
        ax.axvspan(row.block_start - 0.4, row.block_end + 0.4, color=COLORS[i], alpha=0.08)
    ax.set_yscale("log")
    ax.set_xlabel("Ordered eigenvalue number")
    ax.set_ylabel("Cosine eigenvalue (log scale)")
    ax.grid(axis="y", alpha=0.18)
    ax.spines[["top", "right"]].set_visible(False)


fig, axes = plt.subplots(4, 3, figsize=(15, 13), constrained_layout=True)
for ax, case in zip(axes.flat, cases.case):
    spectrum(ax, case)
    ax.set_title(case.replace("_", " "), fontsize=10, wrap=True)
for ax in list(axes.flat)[len(cases) :]:
    ax.axis("off")
fig.suptitle(
    "Eigenvalue blocks in all 11 method_proof benchmark examples\n"
    "Colors mark the existing log-spectrum segmentation; only retained modes are shown.",
    fontsize=15,
)
fig.savefig(ROOT / "eigenvalue_blocks.png", dpi=160)
plt.close(fig)


with PdfPages(ROOT / "benchmark_clusters.pdf") as pdf:
    for case in cases.case:
        selected = blocks[blocks.case.eq(case)]
        names = selected.block_name.tolist()
        labels = [f"Modes {r.block_start}–{r.block_end}" for r in selected.itertuples()]
        fig = plt.figure(figsize=(14, 10), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.4])
        ax = fig.add_subplot(grid[0, :])
        spectrum(ax, case)
        ax.legend(ncol=min(8, len(names)), loc="upper right", frameon=False)
        for col, geometry in enumerate(["raw", "adaptive_diffusion"]):
            ax = fig.add_subplot(grid[1, col])
            subset = runs[runs.case.eq(case) & runs.geometry.eq(geometry)]
            values = np.full((len(names), len(METHODS)), np.nan)
            annotations = {}
            for i, name in enumerate(names):
                for j, method in enumerate(METHODS):
                    record = subset[subset.block.eq(name) & subset.method.eq(method)]
                    if record.empty:
                        annotations[i, j] = "Missing"
                        continue
                    row = record.iloc[0]
                    if row.status == "ok":
                        values[i, j] = row.ari_truth
                        count = int(row.n_clusters)
                        annotations[i, j] = (
                            f"{row.ari_truth:.2f}\n{count} cluster{'s' if count != 1 else ''}"
                        )
                    else:
                        annotations[i, j] = {
                            "failed_tree": "Tree\nfailed",
                            "unsupported": "Unsupported",
                            "failed_gate": "Gate\nfailed",
                            "failed_geometry": "Geometry\nfailed",
                        }.get(row.status, row.status)
            cmap = plt.colormaps["cividis"].copy()
            cmap.set_bad("#e6e6e6")
            im = ax.imshow(values, aspect="auto", vmin=-1, vmax=1, cmap=cmap)
            for (i, j), text in annotations.items():
                color = "white" if np.isfinite(values[i, j]) and values[i, j] < -0.1 else "#191919"
                ax.text(j, i, text, ha="center", va="center", fontsize=8, color=color)
            ax.set_xticks(
                range(len(METHODS)),
                [x.replace("_", " ").title() for x in METHODS],
                rotation=45,
                ha="right",
            )
            ax.set_yticks(range(len(names)), labels)
            ax.set_title(
                "Raw block coordinates"
                if geometry == "raw"
                else "Adaptive diffusion of block coordinates"
            )
            fig.colorbar(
                im, ax=ax, shrink=0.65, label="Adjusted Rand index against generator labels"
            )
        truth_k = int(cases.loc[cases.case.eq(case), "true_clusters"].iloc[0])
        fig.suptitle(
            f"{case.replace('_', ' ')}\nGenerator clusters: {truth_k}; each cell shows agreement and TBS cluster count",
            fontsize=15,
        )
        fig.supxlabel(
            "One fixed seed per example. Failed and unsupported cells are unavailable. Support thresholds are not enforced; results are descriptive.",
            fontsize=10,
        )
        pdf.savefig(fig)
        if case == cases.case.iloc[0]:
            fig.savefig(ROOT / "binary_example.png", dpi=160)
        plt.close(fig)
print("Rendered eigenvalue_blocks.png, binary_example.png, and benchmark_clusters.pdf")
