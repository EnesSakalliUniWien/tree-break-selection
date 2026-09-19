"""GO-term metadata and cosine feature loadings for subspace analysis."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize


def parse_go_term(column: str) -> tuple[str, str]:
    import re

    match = re.search(r"(GO:\d{7})", column)
    go_id = match.group(1) if match else ""
    term = column.replace(f"({go_id})", "").strip() if go_id else column
    return term, go_id


def binary_entropy_bits(probabilities: np.ndarray) -> np.ndarray:
    p = np.asarray(probabilities, dtype=float)
    out = np.zeros_like(p)
    valid = (p > 0.0) & (p < 1.0)
    q = p[valid]
    out[valid] = -(q * np.log2(q) + (1.0 - q) * np.log2(1.0 - q))
    return out


def component_feature_loadings(
    values: np.ndarray, eigvals: np.ndarray, eigvecs: np.ndarray
) -> np.ndarray:
    """Return feature loadings for the sample-cosine eigenvectors.

    The cosine operator is X_norm X_norm^T. For positive eigenvalue lambda and
    sample eigenvector u, the corresponding feature loading is
    X_norm^T u / sqrt(lambda).
    """

    row_normed = normalize(values, norm="l2", axis=1)
    loadings = np.zeros((values.shape[1], len(eigvals)), dtype=float)
    for idx, eigval in enumerate(eigvals):
        if eigval <= 1e-12:
            continue
        loadings[:, idx] = row_normed.T @ eigvecs[:, idx] / math.sqrt(float(eigval))
    return np.nan_to_num(loadings)


def build_term_metadata(data: pd.DataFrame) -> pd.DataFrame:
    supports = data.sum(axis=0).to_numpy(dtype=int)
    prevalence = supports / max(len(data), 1)
    parsed = [parse_go_term(col) for col in data.columns]
    return pd.DataFrame(
        {
            "go_term": [term for term, _go_id in parsed],
            "go_id": [go_id for _term, go_id in parsed],
            "column": data.columns,
            "support": supports,
            "prevalence": prevalence,
            "entropy_bits": binary_entropy_bits(prevalence),
        }
    )


def build_axis_loadings(
    *,
    term_metadata: pd.DataFrame,
    loadings: np.ndarray,
    block_start: int,
    block_end: int,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    top_rows = []
    for component in range(block_start, block_end + 1):
        component_values = loadings[:, component - 1]
        frame = term_metadata.copy()
        frame.insert(0, "axis", int(component))
        frame.insert(1, "axis_name", f"mode_{component:02d}")
        frame["loading"] = component_values
        frame["abs_loading"] = np.abs(component_values)
        frame["loading_sign"] = np.where(frame["loading"] >= 0.0, "positive", "negative")
        frame["abs_rank"] = frame["abs_loading"].rank(method="first", ascending=False).astype(int)
        rows.append(frame)

        pos = frame.sort_values("loading", ascending=False).head(top_n).copy()
        pos.insert(2, "selection", "top_positive")
        neg = frame.sort_values("loading", ascending=True).head(top_n).copy()
        neg.insert(2, "selection", "top_negative")
        abs_top = frame.sort_values("abs_loading", ascending=False).head(top_n).copy()
        abs_top.insert(2, "selection", "top_absolute")
        top_rows.extend([pos, neg, abs_top])

    all_loadings = pd.concat(rows, ignore_index=True)
    top_loadings = pd.concat(top_rows, ignore_index=True)
    return all_loadings, top_loadings


def axis_term_summary_lines(top_terms: pd.DataFrame) -> list[str]:
    """Summarize the strongest term loadings for compact embedding annotations."""

    absolute = top_terms[top_terms["selection"].eq("top_absolute")].copy()
    if absolute.empty:
        absolute = top_terms.copy()
    n_axes = int(absolute["axis"].nunique())
    terms_per_axis = 2 if n_axes <= 6 else 1
    lines: list[str] = []
    for axis, frame in absolute.groupby("axis", sort=True):
        top = frame.sort_values("abs_loading", ascending=False).head(terms_per_axis)
        labels = []
        for row in top.itertuples(index=False):
            sign = "+" if float(row.loading) >= 0 else "-"
            labels.append(f"{row.go_term} ({sign}{abs(float(row.loading)):.3f})")
        lines.append(f"mode {int(axis):02d}: " + "; ".join(labels))
    return lines


def select_axis_term_heatmap(
    axis_frame: pd.DataFrame, *, terms_per_axis: int, max_terms: int
) -> pd.DataFrame:
    selected_frames = []
    for _axis, frame in axis_frame.groupby("axis", sort=True):
        selected_frames.append(
            frame.sort_values("abs_loading", ascending=False).head(terms_per_axis)
        )
    if not selected_frames:
        return pd.DataFrame()
    selected = pd.concat(selected_frames, ignore_index=True)
    term_order = (
        selected.groupby(["go_term", "go_id"], as_index=False)
        .agg(max_abs_loading=("abs_loading", "max"))
        .sort_values("max_abs_loading", ascending=False)
        .head(max_terms)
    )
    selected_terms = set(zip(term_order["go_term"], term_order["go_id"], strict=False))
    frame = axis_frame[
        axis_frame.apply(lambda row: (row["go_term"], row["go_id"]) in selected_terms, axis=1)
    ].copy()
    frame["term_label"] = frame.apply(
        lambda row: f"{row['go_term']} ({row['go_id']})" if row["go_id"] else str(row["go_term"]),
        axis=1,
    )
    order = [
        f"{row.go_term} ({row.go_id})" if row.go_id else str(row.go_term)
        for row in term_order.itertuples(index=False)
    ]
    return (
        frame.pivot_table(index="term_label", columns="axis", values="loading", aggfunc="first")
        .reindex(order)
        .fillna(0.0)
    )
