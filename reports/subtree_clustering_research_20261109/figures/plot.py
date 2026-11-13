"""Plot completed Gaussian factorial counts and isolated runtime preflight."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGURE_NAME = "gaussian_null_and_runtime_preflight"
ROW_ORDER = [
    (operator, method) for operator in ("symmetric", "markov") for method in ("average", "ward")
]
COLUMN_ORDER = [(source, energy) for source in ("linkage", "nnls") for energy in (False, True)]


def null_counts():
    path = ROOT / "factorial_cell_0" / "runs.csv"
    rows = pd.read_csv(path)
    assert len(rows) == 800 and rows.status.eq("ok").all()
    assert rows.cell.eq(0).all() and rows.n.eq(64).all() and rows.p.eq(12).all()
    assert set(rows.root_sibling_reject.unique()) == {True, False}
    keys = ["operator", "linkage_method", "branch_source", "energy"]
    assert not rows.duplicated(keys + ["replicate"]).any()
    groups = rows.groupby(keys).root_sibling_reject.agg(["size", "sum"])
    assert len(groups) == 16 and groups["size"].eq(50).all()
    return np.array(
        [
            [
                groups.loc[(operator, method, source, energy), "sum"]
                for source, energy in COLUMN_ORDER
            ]
            for operator, method in ROW_ORDER
        ]
    )


def runtime_ratios():
    path = ROOT / "preflight_isolated" / "runtime_gate.json"
    gate = json.loads(path.read_text())
    assert gate["complete_panel"] and not gate["passed"]
    cases = pd.DataFrame(gate["cases"])
    assert len(cases) == 6 and cases.complete.all()
    assert cases.passed.eq(cases.ratio.le(5)).all()
    return cases.pivot(index="n", columns="p", values="ratio").loc[[32, 128, 240], [12, 72]]


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, (null, runtime) = plt.subplots(
        1, 2, figsize=(12.6, 5.5), gridspec_kw={"width_ratios": [1.08, 1]}
    )
    fig.subplots_adjust(left=0.16, right=0.96, bottom=0.34, top=0.75, wspace=0.34)
    fig.text(
        0.055,
        0.94,
        "Subtree clustering diagnostics",
        fontsize=22,
        fontweight="bold",
        color="#162337",
    )

    counts = null_counts()
    null.imshow(counts, cmap="Blues", vmin=0, vmax=50, aspect="auto")
    null.set_title("A   Selected-root native sibling test", loc="left", fontweight="bold", pad=26)
    null.text(
        0,
        1.045,
        "Gaussian null, n = 64, p = 12, nominal α = 0.01",
        transform=null.transAxes,
        color="#43536a",
        fontsize=10,
    )
    null.set_yticks(
        range(4), ["Symmetric / average", "Symmetric / Ward", "Markov / average", "Markov / Ward"]
    )
    null.set_xticks(
        range(4),
        [
            "Linkage\nNo energy\ntruncation",
            "Linkage\n0.90 energy",
            "NNLS\nNo energy\ntruncation",
            "NNLS\n0.90 energy",
        ],
    )
    null.tick_params(length=0, pad=8)
    null.set_xticks(np.arange(5) - 0.5, minor=True)
    null.set_yticks(np.arange(5) - 0.5, minor=True)
    null.grid(which="minor", color="white", linewidth=1.5)
    null.tick_params(which="minor", bottom=False, left=False)
    null.spines[["left", "bottom"]].set_visible(False)
    for (i, j), count in np.ndenumerate(counts):
        null.text(
            j,
            i,
            f"{count}/50",
            va="center",
            ha="center",
            fontweight="bold",
            color="white" if count >= 31 else "#142536",
        )
    null.set_xlabel("Branch source and contrast-energy rule", labelpad=10)

    ratios = runtime_ratios()
    runtime.set_title("B   Isolated runtime preflight", loc="left", fontweight="bold", pad=26)
    runtime.text(
        0,
        1.045,
        "Candidate / baseline median complete-pipeline time",
        transform=runtime.transAxes,
        color="#43536a",
        fontsize=10,
    )
    runtime.axhline(5, color="#6a6b72", linewidth=1.3, linestyle="--")
    runtime.text(
        0.015, 5.3, "5× limit", transform=runtime.get_yaxis_transform(), fontsize=9, color="#454951"
    )
    for p, color, marker, direction in [(12, "#18648a", "o", -1), (72, "#a64324", "D", 1)]:
        runtime.plot(
            ratios.index, ratios[p], label=f"p = {p}", color=color, marker=marker, lw=2, ms=6
        )
        for n, ratio in ratios[p].items():
            runtime.annotate(
                f"{ratio:.2f}×",
                (n, ratio),
                xytext=(0, 8 * direction),
                textcoords="offset points",
                ha="center",
                va="bottom" if direction == 1 else "top",
                color=color,
                fontsize=9,
            )
    runtime.set_yscale("log")
    runtime.set_yticks([1, 2, 5, 10, 25], ["1×", "2×", "5×", "10×", "25×"])
    runtime.set_ylim(0.78, 37)
    runtime.set_xlim(15, 257)
    runtime.set_xticks(ratios.index)
    runtime.set_xlabel("Subtree size n", labelpad=8)
    runtime.grid(axis="y", which="major", color="#e5e9ef")
    runtime.set_axisbelow(True)
    runtime.legend(loc="upper left", frameon=False, fontsize=10)

    fig.text(
        0.055,
        0.065,
        "A  Root BH rejection counts on 50 shared Gaussian datasets per configuration. The base native projection remains in both arms.\n"
        "B  Three complete runs per method per case. Both p = 72 cases at n ≥ 128 exceeded the runtime limit.",
        fontsize=9.5,
        color="#43536a",
        linespacing=1.55,
    )
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{FIGURE_NAME}.{suffix}", dpi=200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
