"""Completeness and statistical summaries for the frozen subtree research panels."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binomtest, t

_STATUSES = ("ok", "unresolved", "unsupported", "error", "failed", "skip")
_AUXILIARY = ("ari", "matched_macro_f1", "matched_rare_f1", "pair_split_share", "pair_merge_share")


def _panel(frame, expected_cells, expected_replicates, columns, keys):
    for value in (expected_cells, expected_replicates):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError("Expected positive integer cell and replicate counts.")
    required = {"cell", "replicate", "status", *columns}
    if not required <= set(frame) or not frame.columns.is_unique:
        raise ValueError(f"Missing or duplicate panel columns: {sorted(required - set(frame))}.")
    frame = frame.reset_index(drop=True).copy()
    for column, stop in (("cell", expected_cells), ("replicate", expected_replicates)):
        values = pd.to_numeric(frame[column], errors="coerce")
        if (
            frame[column].map(lambda v: isinstance(v, (bool, np.bool_))).any()
            or not (values.ge(0) & values.lt(stop) & values.mod(1).eq(0)).fillna(False).all()
        ):
            raise ValueError(f"Unknown {column} ID.")
        frame[column] = values.astype(int)
    if not frame.status.isin(_STATUSES).all():
        raise ValueError("Unknown panel status.")
    if frame.duplicated(keys).any():
        raise ValueError("Panel contains duplicate attempt keys.")
    if "phase" in frame and (not frame.phase.isin([1, 2]).all() or frame.phase.nunique() > 1):
        raise ValueError("Expected one known data phase per panel.")
    return frame


def _metric(series, low, high, *, required=None):
    values = pd.to_numeric(series, errors="coerce").astype(float)
    valid = values.between(low, high)
    boolean = series.map(lambda value: isinstance(value, (bool, np.bool_))).astype(bool)
    present = series.notna() if required is None else (series.notna() | required)
    if (present & (~valid | boolean)).any():
        raise ValueError(f"Malformed {series.name} value.")
    return values


def _holm(values):
    """Keep absent hypotheses in the declared family size."""
    order = sorted((value, i) for i, value in enumerate(values) if value is not None)
    adjusted = [None] * len(values)
    running = 0.0
    for rank, (value, index) in enumerate(order):
        running = max(running, min(1.0, (len(values) - rank) * value))
        adjusted[index] = float(running)
    return adjusted


def _status_counts(status, expected):
    counts = status.value_counts()
    return {key: int(counts.get(key, 0)) for key in _STATUSES} | {
        "missing": int(expected - status.notna().sum())
    }


def _frozen_shape(frame, actual, prescribed):
    return (
        actual == prescribed
        and not frame.empty
        and "phase" in frame
        and bool(frame.phase.eq(2).all())
    )


def calibration_summary(
    frame: pd.DataFrame,
    expected_cells: int = 24,
    expected_replicates: int = 5000,
    alphas: tuple[float, ...] = (0.01, 0.05),
) -> dict:
    """Summarize per-dataset calibration with exact intervals and Holm excess tests.

    Attempted unresolved/failed rows remain in the denominator, and only
    completed p-values contribute rejections. Missing attempts are reported
    separately. The source decision is retained as a count, while threshold
    comparisons use the completed p-values. An unresolved decision overrides
    an apparently successful attempt status. A pass is local to the
    prescribed phase-2 panel and requires resolved, nondegenerate outputs.
    """
    alphas = tuple(float(value) for value in alphas)
    if not alphas or len(set(alphas)) != len(alphas) or not all(0 < a < 1 for a in alphas):
        raise ValueError("Alphas must be distinct probabilities strictly between zero and one.")
    frame = _panel(
        frame,
        expected_cells,
        expected_replicates,
        ("p_value", "decision"),
        ["cell", "replicate"],
    )
    expected_cells, expected_replicates = int(expected_cells), int(expected_replicates)
    decisions = frame.decision.astype("string").str.strip().str.lower()
    resolved = decisions.isin(("reject", "do_not_reject"))
    unresolved = decisions.eq("unresolved").fillna(False)
    if (decisions.notna() & ~(resolved | unresolved)).any():
        raise ValueError("Unknown calibration decision.")
    if (frame.status.eq("ok") & decisions.isna()).any() or (
        ~frame.status.eq("ok") & resolved
    ).any():
        raise ValueError("Calibration decision conflicts with attempt status.")
    frame.loc[frame.status.eq("ok") & unresolved, "status"] = "unresolved"
    frame["p_value"] = _metric(frame.p_value, 0, 1, required=frame.status.eq("ok"))
    expected = expected_cells * expected_replicates
    counts = _status_counts(frame.status, expected)
    frozen = _frozen_shape(
        frame, (expected_cells, expected_replicates, alphas), (24, 5000, (0.01, 0.05))
    )
    comparisons, degenerate = [], []
    for cell in range(expected_cells):
        group = frame.loc[frame.cell == cell]
        completed = group.loc[group.status.eq("ok"), "p_value"]
        n = len(group)
        status = _status_counts(group.status, expected_replicates)
        zeros, ones = int(completed.eq(0).sum()), int(completed.eq(1).sum())
        if len(completed) > 1 and completed.nunique() == 1:
            degenerate.append(cell)
        for alpha in alphas:
            rejections = int(completed.le(alpha).sum())
            interval, excess = None, None
            if n:
                test = binomtest(rejections, n, p=alpha, alternative="greater")
                exact = binomtest(rejections, n).proportion_ci(method="exact")
                interval, excess = [float(exact.low), float(exact.high)], float(test.pvalue)
            comparisons.append(
                {
                    "cell": cell,
                    "alpha": alpha,
                    "attempts": n,
                    "rejections": rejections,
                    "status_counts": status,
                    "zero_p_values": zeros,
                    "one_p_values": ones,
                    "rate": rejections / n if n else None,
                    "interval": interval,
                    "excess_p": excess,
                }
            )
    adjusted = _holm([row["excess_p"] for row in comparisons])
    for row, value in zip(comparisons, adjusted, strict=True):
        row["excess_p_holm"] = value
    excess_detected = any(value is not None and value <= 0.05 for value in adjusted)
    all_successful = counts["ok"] == expected
    return {
        "design": "frozen_shape" if frozen else "development_or_custom",
        "complete": counts["missing"] == 0,
        "all_successful": all_successful,
        "passed": frozen and all_successful and not degenerate and not excess_detected,
        "expected_attempts": expected,
        "attempts": len(frame),
        "counts": counts,
        "decision_counts": {
            str(key): int(v) for key, v in frame.decision.value_counts(dropna=False).items()
        },
        "degenerate_cells": degenerate,
        "excess_detected": excess_detected,
        "comparisons": comparisons,
    }


def _difference_stats(differences):
    values = np.asarray(differences, dtype=float)
    n = len(values)
    if not n:
        return None, None, None
    mean = float(values.mean())
    if n < 2:
        return mean, None, None
    if np.ptp(values) == 0:
        return mean, [mean, mean], float(mean >= 0)
    stderr = float(values.std(ddof=1) / np.sqrt(n))
    margin = float(t.ppf(0.975, n - 1)) * stderr
    return mean, [mean - margin, mean + margin], float(t.cdf(mean / stderr, n - 1))


def paired_recovery_summary(
    frame: pd.DataFrame,
    baseline: str,
    candidate: str,
    expected_cells: int = 16,
    expected_replicates: int = 200,
) -> dict:
    """Compare matched datasets, weighting all cells equally within a replicate block.

    The primary interval is Student t across blocks with a successful pair
    in every prescribed cell. Per-cell loss tests use paired differences and
    Holm over every prescribed cell. Incomplete pairs and blocks are counted
    explicitly and preclude a frozen pass.
    """
    if not baseline or not candidate or baseline == candidate:
        raise ValueError("Supply two distinct methods.")
    frame = _panel(
        frame,
        expected_cells,
        expected_replicates,
        ("method", "ami", "data_hash", "seed"),
        ["cell", "replicate", "method"],
    )
    expected_cells, expected_replicates = int(expected_cells), int(expected_replicates)
    if not frame.method.isin([baseline, candidate]).all():
        raise ValueError("Unknown method ID in paired frame.")
    if (
        not frame.data_hash.map(lambda value: isinstance(value, str)).all()
        or not frame.data_hash.astype(str).str.fullmatch(r"[0-9a-fA-F]{64}").all()
    ):
        raise ValueError("Expected a SHA-256 data_hash for each attempt.")
    frame["data_hash"] = frame.data_hash.astype(str).str.lower()
    seeds = pd.to_numeric(frame.seed, errors="coerce")
    if (
        frame.seed.map(lambda v: isinstance(v, (bool, np.bool_))).any()
        or not (np.isfinite(seeds) & seeds.ge(0) & seeds.mod(1).eq(0)).fillna(False).all()
    ):
        raise ValueError("Expected a nonnegative integer seed for each attempt.")
    frame["seed"] = seeds
    frame["ami"] = _metric(frame.ami, -1, 1, required=frame.status.eq("ok"))
    auxiliary = [metric for metric in _AUXILIARY if metric in frame]
    for metric in auxiliary:
        frame[metric] = _metric(frame[metric], -1 if metric == "ari" else 0, 1)
    expected = expected_cells * expected_replicates
    frozen = _frozen_shape(frame, (expected_cells, expected_replicates), (16, 200))
    index = pd.MultiIndex.from_product(
        [range(expected_cells), range(expected_replicates)], names=["cell", "replicate"]
    )
    base, trial = (
        frame.loc[frame.method == method].set_index(["cell", "replicate"]).reindex(index)
        for method in (baseline, candidate)
    )
    both_attempted = base.status.notna() & trial.status.notna()
    dataset_rows = frame.drop_duplicates(["cell", "replicate"])
    for key in ("seed", "data_hash"):
        if (both_attempted & base[key].ne(trial[key])).any():
            raise ValueError(f"Paired dataset {key} mismatch.")
        if dataset_rows[key].duplicated().any():
            raise ValueError(f"Dataset {key} is reused across distinct attempt keys.")
    successful = base.status.eq("ok") & trial.status.eq("ok")
    ami_delta = (trial.ami - base.ami).where(successful)
    cells = []
    for cell in range(expected_cells):
        valid = successful.loc[cell]
        mean, interval, loss = _difference_stats(ami_delta.loc[cell].dropna())
        other = {}
        for metric in auxiliary:
            delta = (trial[metric].loc[cell] - base[metric].loc[cell]).where(valid).dropna()
            other[metric] = {
                "pairs": len(delta),
                "mean_difference": float(delta.mean()) if len(delta) else None,
            }
        cells.append(
            {
                "cell": cell,
                "paired_attempted": int(both_attempted.loc[cell].sum()),
                "paired_successful": int(valid.sum()),
                "mean_ami_difference": mean,
                "interval": interval,
                "loss_p": loss,
                "auxiliary": other,
                "baseline_status_counts": _status_counts(
                    base.status.loc[cell], expected_replicates
                ),
                "candidate_status_counts": _status_counts(
                    trial.status.loc[cell], expected_replicates
                ),
            }
        )
    adjusted = _holm([row["loss_p"] for row in cells])
    for row, value in zip(cells, adjusted, strict=True):
        row["loss_p_holm"] = value
    loss_detected = any(value is not None and value <= 0.05 for value in adjusted)
    blocks = ami_delta.unstack("cell").dropna().mean(axis=1)
    mean, interval, _ = _difference_stats(blocks)
    all_successful = bool(successful.all())
    return {
        "design": "frozen_shape" if frozen else "development_or_custom",
        "baseline": baseline,
        "candidate": candidate,
        "complete": bool(both_attempted.all()),
        "all_successful": all_successful,
        "passed": frozen
        and all_successful
        and interval is not None
        and interval[0] > 0
        and not loss_detected,
        "loss_detected": loss_detected,
        "coverage": {
            "expected_pairs": expected,
            "paired_attempted": int(both_attempted.sum()),
            "paired_successful": int(successful.sum()),
            "paired_share": float(successful.sum() / expected),
            "baseline_status_counts": _status_counts(base.status, expected),
            "candidate_status_counts": _status_counts(trial.status, expected),
        },
        "aggregate": {
            "blocks": len(blocks),
            "expected_blocks": expected_replicates,
            "mean_ami_difference": mean,
            "interval": interval,
        },
        "cells": cells,
    }
