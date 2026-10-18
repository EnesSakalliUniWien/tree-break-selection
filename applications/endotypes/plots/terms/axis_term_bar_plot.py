"""Render signed GO-term loading bars for PNG and PDF reports."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from applications.endotypes.plots.shared._plot_text import wrap_labels


def draw_axis_term_bars(ax: plt.Axes, frame: pd.DataFrame, title: str, *, top_n: int) -> None:
    top = frame.sort_values("abs_loading", ascending=False).head(top_n).sort_values("loading")
    colors = np.where(top["loading"] >= 0, "#4c78a8", "#e45756")
    ax.barh(wrap_labels(top["go_term"], width=42), top["loading"], color=colors)
    ax.axvline(0.0, color="#333333", linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel("feature loading")
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(axis="x", alpha=0.18)


def create_axis_term_bars(
    frame: pd.DataFrame, title: str, *, top_n: int, figsize: tuple[float, float] | None = None
) -> plt.Figure:
    if figsize is None:
        height = max(5.5, min(10.5, 0.34 * min(len(frame), top_n) + 2.0))
        figsize = (11.5, height)
    fig, ax = plt.subplots(figsize=figsize)
    draw_axis_term_bars(ax, frame, title, top_n=top_n)
    fig.tight_layout()
    return fig
