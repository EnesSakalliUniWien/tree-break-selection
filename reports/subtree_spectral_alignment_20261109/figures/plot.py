"""Render crossed partitions and finite-panel topology coverage."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "alternatives" / "summary.csv"
STEM = "crossed_partitions_and_topology_coverage"
BLUE = "#18648a"
ORANGE = "#a64324"
INK = "#162337"
MUTED = "#43536a"
POINTS = {"A": (-2, -1), "B": (-2, 1), "C": (2, -1), "D": (2, 1)}
REGIMES = ("crossed", "tied_population_crossed", "compatible")


def coverage():
    data = pd.read_csv(SOURCE)
    assert not data.duplicated(["regime", "replicate"]).any()
    grouped = data.groupby("regime")
    assert grouped.size().eq(20).all()
    assert set(grouped.groups) == {*REGIMES, "common_mean_null"}
    null = grouped.get_group("common_mean_null")
    assert null.unique_topologies.eq(63).all() and null.union_has_both.isna().all()
    rows = data[data.regime.isin(REGIMES)]
    columns = ["union_has_both", "any_tree_has_both"]
    assert rows[columns].notna().all().all()
    assert set(rows[columns].to_numpy().ravel()) == {False, True}
    return rows.groupby("regime")[columns].sum().reindex(REGIMES).astype(int)


def draw_partition(ax, groups, axis):
    boundary = ax.axvline if axis == "x" else ax.axhline
    boundary(0, color="#9faab8", linewidth=1, ls="--")
    for group, color in zip(groups, (BLUE, ORANGE), strict=True):
        x, y = np.array([POINTS[leaf] for leaf in group]).T
        ax.plot(x, y, color=color, alpha=0.30, linewidth=7, solid_capstyle="round")
        ax.scatter(x, y, color=color, s=55, zorder=3)
        for name, xx, yy in zip(group, x, y, strict=True):
            ax.annotate(name, (xx, yy), xytext=(8, 4), textcoords="offset points", color=INK)
    ax.set(xlim=(-2.7, 2.9), ylim=(-1.55, 1.6), xticks=[-2, 0, 2], yticks=[-1, 0, 1])
    ax.set_title(f"${axis}$ split   {'|'.join(groups)}", fontsize=11, fontweight="bold")
    ax.set_xlabel("x", labelpad=0)
    ax.set_ylabel("y", rotation=0, labelpad=5)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(length=2)


def draw_tree(ax, groups):
    for j, (group, color) in enumerate(zip(groups, (BLUE, ORANGE), strict=True)):
        x = 2 * j
        ax.plot([x, x, x + 1, x + 1], [0, 1, 1, 0], color=color, linewidth=2)
        ax.plot([x + 0.5, x + 0.5], [1, 1.65], color=INK, linewidth=1.5)
        for i, leaf in enumerate(group):
            ax.text(x + i, -0.13, leaf, ha="center", va="top", color=color, fontweight="bold")
    ax.plot([0.5, 2.5], [1.65, 1.65], color=INK, linewidth=1.5)
    ax.plot([1.5, 1.5], [1.65, 1.92], color=INK, linewidth=1.5)
    ax.set(xlim=(-0.45, 3.45), ylim=(-0.5, 2.05))
    ax.set_axis_off()


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(11.7, 5.3), facecolor="white")
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(1, 1, 1.9),
        height_ratios=(1.18, 0.7),
        left=0.07,
        right=0.97,
        top=0.76,
        bottom=0.24,
        wspace=0.36,
        hspace=0.37,
    )
    fig.text(0.06, 0.91, "When two partitions cross", fontsize=22, fontweight="bold", color=INK)
    fig.text(0.06, 0.827, "A  Two factors on the same four leaves", fontweight="bold", color=INK)
    for j, (groups, axis) in enumerate([(("AB", "CD"), "x"), (("AC", "BD"), "y")]):
        draw_partition(fig.add_subplot(grid[0, j]), groups, axis)
        draw_tree(fig.add_subplot(grid[1, j]), groups)

    ax = fig.add_subplot(grid[:, 2])
    ax.set_title("B  Both contrasts represented", loc="left", fontweight="bold", pad=14, color=INK)
    positions = np.arange(3)
    colors = (BLUE, ORANGE)
    for j, (column, label, color) in enumerate(
        zip(
            coverage().to_numpy().T,
            ("Across the spectral family", "Together in one family tree"),
            colors,
            strict=True,
        )
    ):
        bars = ax.bar(positions + (j - 0.5) * 0.35, column, width=0.32, color=color, label=label)
        ax.bar_label(
            bars, labels=[f"{value}/20" for value in column], padding=3, fontsize=9, color=color
        )
    ax.set_xticks(positions, ["Crossed", "Tied-population\ncrossed", "Compatible\n(nested B)"])
    ax.set(ylim=(0, 29), yticks=[0, 5, 10, 15, 20], ylabel="Datasets out of 20")
    ax.grid(axis="y", color="#e5e9ef")
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    fig.text(
        0.06,
        0.13,
        "A  Both partitions describe real factors in the example. Their clades cross, so a single rooted tree cannot contain both.\n"
        "B  Descriptive 20-seed panels. In the compatible case, contrast B is a nested group versus its complement.\n"
        "    The common-mean null still produced 63 distinct candidates per dataset. Counts describe geometric sensitivity.",
        color=MUTED,
        fontsize=9,
        linespacing=1.48,
        va="top",
    )
    for extension in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{STEM}.{extension}", dpi=200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
