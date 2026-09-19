"""Compose full-width embedding panels and native-size GO-term annotations."""

import matplotlib.pyplot as plt
import numpy as np

from applications.endotypes.plots.shared._plot_text import wrap_annotation_lines


def create_embedding_comparison(
    subspace_image: np.ndarray | None,
    diffusion_image: np.ndarray | None,
    title: str,
    axis_term_lines: list[str],
) -> plt.Figure:
    split = (len(axis_term_lines) + 1) // 2
    columns = [
        wrap_annotation_lines(lines, width=78)
        for lines in (axis_term_lines[:split], axis_term_lines[split:])
    ]
    text_height = max(0.9, 0.25 * max(len(text.splitlines()) for text in columns) + 0.8)
    # Preserve the scale of the original eight-inch-wide plots, including legends.
    image_height = max(
        (8.0 * img.shape[0] / img.shape[1] for img in (subspace_image, diffusion_image)
         if img is not None), default=6.5,
    )
    fig = plt.figure(figsize=(18, image_height + text_height + 1.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[image_height, text_height])
    for col, (img, label) in enumerate([
        (subspace_image, "Subspace embedding"), (diffusion_image, "Adaptive diffusion embedding")
    ]):
        ax = fig.add_subplot(grid[0, col])
        ax.axis("off")
        ax.set_title(label, fontsize=13, loc="left")
        if img is None:
            ax.text(0.5, 0.5, "not available", ha="center", va="center", fontsize=12)
        else:
            ax.imshow(img)
        text_ax = fig.add_subplot(grid[1, col])
        text_ax.axis("off")
        text_ax.text(
            0, 0.98, columns[col], va="top", fontsize=12, linespacing=1.2,
            transform=text_ax.transAxes,
        )
        if col == 0:
            text_ax.set_title("Highest GO-term loadings by cosine eigenmode", loc="left", fontsize=12)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig
