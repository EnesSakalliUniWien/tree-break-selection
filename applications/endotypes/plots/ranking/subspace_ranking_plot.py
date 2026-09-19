"""Render ranking inputs on aligned rows in specificity-aware order."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def create_subspace_ranking(ranking: pd.DataFrame) -> plt.Figure:
    top = ranking[ranking["status"].eq("ok")].sort_values("display_rank").head(20)
    fig, axes = plt.subplots(
        1, 3, figsize=(18, max(5.5, 0.55 * len(top) + 2.5)), sharey=True,
        gridspec_kw={"width_ratios": [1.1, 1.0, 1.4]},
    )
    y = np.arange(len(top))
    labels = [
        f"{int(row.display_rank):02d}  {row.weighting} / {row.block_name}"
        for row in top.itertuples(index=False)
    ]
    tier_ax, specificity_ax, bic_ax = axes
    tier_ax.set_yticks(y, labels, fontsize=11)
    tier_ax.set_xlim(0, 1)
    tier_ax.set_xticks([])
    tier_ax.set_title("Quality tier (lower first)", fontsize=12)
    for pos, row in enumerate(top.itertuples(index=False)):
        tier_ax.text(
            0.04, pos, f"{int(row.quality_tier)}  {row.quality_tier_label}",
            va="center", fontsize=11,
        )
    for ax, column, title, color in [
        (specificity_ax, "specificity_score", "Specificity score (higher first)", "#278577"),
        (bic_ax, "go_bic_active_per_gene", "GO-BIC active / gene (lower is better)", "#4c78a8"),
    ]:
        bars = ax.barh(y, top[column], height=0.6, color=color)
        ax.bar_label(bars, fmt="%.3f" if column == "specificity_score" else "%.2f", padding=5)
        ax.margins(x=0.22)
        ax.set_title(title, fontsize=12)
        ax.grid(axis="x", alpha=0.18)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", left=False)
    specificity_ax.set_xlim(0, 1.2)
    tier_ax.set_ylim(len(top) - 0.5, -0.5)
    for ax in axes:
        for pos in range(1, len(top)):
            if top.iloc[pos].quality_tier != top.iloc[pos - 1].quality_tier:
                ax.axhline(pos - 0.5, color="#777777", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Adaptive diffusion subspaces: specificity-aware ranking", fontsize=15)
    fig.text(
        0.02, 0.025,
        "Order: quality tier → specificity score → specific-cluster fraction → weighted specificity delta → GO-BIC/gene.\n"
        "GO-BIC values are comparable within quality tier; lower GO-BIC alone does not determine rank.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.1, 1, 0.96))
    return fig
