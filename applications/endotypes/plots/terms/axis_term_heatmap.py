"""Render signed GO-term heatmaps with a shared PNG/PDF drawing function."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_terms import (
    select_axis_term_heatmap,
)
from applications.endotypes.plots.shared._plot_text import wrap_labels


def draw_axis_term_heatmap(ax: plt.Axes, pivot: pd.DataFrame, *, label_width: int) -> object:
    limit = float(np.nanmax(np.abs(pivot.to_numpy())))
    if not math.isfinite(limit) or limit <= 0.0:
        limit = 1.0
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(pivot.columns)), [f"{int(axis):02d}" for axis in pivot.columns])
    ax.set_yticks(
        np.arange(len(pivot.index)), wrap_labels(pd.Series(pivot.index), width=label_width)
    )
    return ax.figure.colorbar(image, ax=ax, label="feature loading", fraction=0.028, pad=0.02)


def create_axis_term_heatmap(
    axis_frame: pd.DataFrame,
    title: str = "Combined GO-term loadings across subspace axes",
    *,
    terms_per_axis: int = 8,
    max_terms: int = 60,
    report: bool = False,
) -> plt.Figure | None:
    pivot = select_axis_term_heatmap(axis_frame, terms_per_axis=terms_per_axis, max_terms=max_terms)
    if pivot.empty and not report:
        return None
    if report:
        fig, ax = plt.subplots(figsize=(18, 10.5))
    else:
        height = max(9.0, min(28.0, 0.32 * len(pivot) + 3.0))
        width = max(12.0, min(28.0, 0.42 * len(pivot.columns) + 8.5))
        fig, ax = plt.subplots(figsize=(width, height))
    if pivot.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.5, "No GO-term loading data available", ha="center", va="center", fontsize=13
        )
        return fig
    cbar = draw_axis_term_heatmap(ax, pivot, label_width=62 if report else 58)
    if report:
        ax.set_title(title, fontsize=15, pad=12)
        ax.set_xlabel("cosine mode", fontsize=12)
        ax.tick_params(axis="x", labelsize=11)
        ax.tick_params(axis="y", labelsize=8.8)
        cbar.ax.tick_params(labelsize=10)
        fig.subplots_adjust(left=0.37, right=0.91, top=0.9, bottom=0.08)
    else:
        ax.set_xlabel("cosine mode")
        ax.set_title(title)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=8)
        fig.tight_layout()
    return fig
