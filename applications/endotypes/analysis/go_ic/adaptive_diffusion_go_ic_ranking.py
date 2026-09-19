"""GO-IC quality tiers and specificity-aware subspace ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd

from applications.endotypes.plots.go_ic_tree_summary_plots import (
    quality_tier,
)


def add_specificity_aware_rank(
    ranking: pd.DataFrame, coherence_long: pd.DataFrame | None
) -> pd.DataFrame:
    out = ranking.copy()
    out["old_display_rank"] = out.get("display_rank", pd.Series(index=out.index, dtype=float))
    out["specific_cluster_count"] = 0
    out["specific_cluster_fraction"] = 0.0
    out["weighted_specificity_delta"] = 0.0
    out["median_specificity_delta"] = 0.0
    out["mean_specificity_delta"] = 0.0
    out["specificity_score"] = 0.0
    if coherence_long is not None and not coherence_long.empty:
        metrics: dict[str, dict[str, float]] = {}
        for run_id, frame in coherence_long.groupby("run_id"):
            deltas = pd.to_numeric(frame["top_term_prevalence_delta"], errors="coerce").fillna(0.0)
            sizes = pd.to_numeric(frame["cluster_size"], errors="coerce").fillna(0.0)
            specific = frame["coherent_by_rule"].astype(bool) & deltas.gt(0.0)
            specific_deltas = deltas[specific]
            weight_denom = float(sizes[specific].sum())
            weighted = (
                float((deltas[specific] * sizes[specific]).sum() / weight_denom)
                if weight_denom > 0
                else 0.0
            )
            fraction = float(specific.sum() / max(len(frame), 1))
            metrics[str(run_id)] = {
                "specific_cluster_count": float(specific.sum()),
                "specific_cluster_fraction": fraction,
                "weighted_specificity_delta": weighted,
                "median_specificity_delta": float(specific_deltas.median())
                if len(specific_deltas)
                else 0.0,
                "mean_specificity_delta": float(specific_deltas.mean())
                if len(specific_deltas)
                else 0.0,
                "specificity_score": 0.5 * fraction + 0.5 * weighted,
            }
        for column in [
            "specific_cluster_count",
            "specific_cluster_fraction",
            "weighted_specificity_delta",
            "median_specificity_delta",
            "mean_specificity_delta",
            "specificity_score",
        ]:
            out[column] = out["run_id"].map(
                lambda run_id: metrics.get(str(run_id), {}).get(column, 0.0)
            )
        out["specific_cluster_count"] = out["specific_cluster_count"].astype(int)

    ok_mask = out["status"].eq("ok")
    if ok_mask.any():
        order = out.loc[ok_mask].sort_values(
            [
                "quality_tier",
                "specificity_score",
                "specific_cluster_fraction",
                "weighted_specificity_delta",
                "go_bic_active_per_gene",
            ],
            ascending=[True, False, False, False, True],
        )
        ranks = dict(zip(order["run_id"], np.arange(1, len(order) + 1), strict=False))
        out.loc[ok_mask, "specificity_aware_rank"] = (
            out.loc[ok_mask, "run_id"].map(ranks).astype(int)
        )
        out.loc[ok_mask, "display_rank"] = out.loc[ok_mask, "specificity_aware_rank"].astype(int)
    return out.sort_values(["status", "display_rank"], ascending=[False, True]).reset_index(
        drop=True
    )


def rank_go_ic_results(summary_rows: list[dict[str, object]]) -> pd.DataFrame:
    ranking = pd.DataFrame.from_records(summary_rows)
    ok_mask = ranking["status"].eq("ok")
    if ok_mask.any():
        raw_order = ranking.loc[ok_mask].sort_values(
            [
                "go_bic_active",
                "coherent_cluster_fraction",
                "weighted_mean_within_tfidf_cosine",
                "singleton_gene_fraction",
            ],
            ascending=[True, False, False, True],
        )
        raw_rank = dict(zip(raw_order["run_id"], np.arange(1, len(raw_order) + 1), strict=False))
        tiers = ranking.loc[ok_mask].apply(quality_tier, axis=1, result_type="expand")
        ranking.loc[ok_mask, "quality_tier"] = tiers[0].astype(int).to_numpy()
        ranking.loc[ok_mask, "quality_tier_label"] = tiers[1].astype(str).to_numpy()
        ranking.loc[ok_mask, "raw_go_ic_rank"] = (
            ranking.loc[ok_mask, "run_id"].map(raw_rank).astype(int)
        )
        ok_sorted = ranking.loc[ok_mask].sort_values(
            [
                "quality_tier",
                "go_bic_active",
                "coherent_cluster_fraction",
                "weighted_mean_within_tfidf_cosine",
                "singleton_gene_fraction",
            ],
            ascending=[True, True, False, False, True],
        )
        display_rank = dict(
            zip(ok_sorted["run_id"], np.arange(1, len(ok_sorted) + 1), strict=False)
        )
        ranking.loc[ok_mask, "display_rank"] = (
            ranking.loc[ok_mask, "run_id"].map(display_rank).astype(int)
        )
        ranking = ranking.sort_values(
            ["status", "display_rank"], ascending=[False, True]
        ).reset_index(drop=True)
    return ranking
