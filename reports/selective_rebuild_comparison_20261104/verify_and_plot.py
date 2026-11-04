"""Recompute saved ranks and inspect the bounded replacement comparison."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from run_comparison import (
    ROOT,
    importance_summary,
    make_diffusion_root_path,
    root_groups,
    split_score,
)

HERE = Path(__file__).resolve().parent


def main():
    table = pd.read_csv(HERE / "results.csv", keep_default_na=False, na_values=[""])
    assert len(table) == 192 and set(table.scenario) == {"null", "moderate", "strong"}
    assert table.groupby("scenario").replicate.nunique().to_dict() == {
        "null": 12,
        "moderate": 6,
        "strong": 6,
    }
    reference = np.load(HERE / "generated/known_null_scores.npy")
    assert len(reference) == 999 and np.isfinite(reference).all()
    checked = 0
    membership_checks = 0
    for (scenario, replicate), group in table.groupby(["scenario", "replicate"]):
        arrays = np.load(HERE / f"generated/{scenario}_{replicate}.npz")
        path = make_diffusion_root_path(arrays["input"], sigma=1)
        np.testing.assert_allclose(
            path.offset + path.observed_radius * path.slope, arrays["input"], atol=1e-12
        )
        np.testing.assert_allclose(path.e @ path.offset, 0, atol=1e-12)
        for key in ["naive", "importance"]:
            flags = arrays[f"{key}_membership"]
            assert np.all(~flags[:, 0] | flags[:, 1])
        for row in group.itertuples():
            index = 0 if row.event == "full_history" else 1
            if row.method == "conditional_rejection":
                flags = arrays["naive_membership"][:, index]
                r = arrays["naive_radii"][flags]
                p = (1 + np.count_nonzero(r >= path.observed_radius)) / (1 + len(r))
                assert len(r) == row.accepted
            elif row.method == "conditional_importance":
                result = importance_summary(
                    arrays["importance_radii"],
                    arrays["importance_membership"][:, index],
                    path.observed_radius,
                    3,
                )
                p = result["p_value"]
                np.testing.assert_allclose(result["ess"], row.ess)
            elif row.method == "conditional_chain":
                r = arrays[f"chain_{row.event}"]
                assert len(r) == 1000 and r[int(row.anchor)] == path.observed_radius
                p = np.mean(r >= path.observed_radius)
                assert len(np.unique(r)) == row.unique_states
                for sample in r[[0, len(r) // 2, -1]]:
                    point = path.evaluate(float(sample))
                    assert point.signature
                    if index == 0:
                        assert point.signature == path.observed_signature
                    else:
                        assert root_groups(point.signature, 9) == root_groups(
                            path.observed_signature, 9
                        )
                    membership_checks += 1
            else:
                score = split_score(arrays["input"], path.observed_signature)
                np.testing.assert_allclose(score, row.score)
                null = (
                    reference
                    if row.method == "rebuild_known_covariance"
                    else arrays["fitted_null_scores"]
                )
                p = (1 + np.count_nonzero(null >= score)) / (1 + len(null))
            if p is None:
                assert pd.isna(row.p_value)
            else:
                np.testing.assert_allclose(p, row.p_value, rtol=1e-12, atol=1e-15)
            checked += 1
    manifest = json.loads((HERE / "manifest.json").read_text())
    import hashlib

    for name, digest in manifest["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    verification = dict(
        saved_statistics_recomputed=checked,
        sampled_chain_memberships_replayed=membership_checks,
        source_hashes_verified=len(manifest["source_sha256"]),
        baseline_tests=16,
        scope="Known-I Gaussian roots only; numerical membership, no exact-real certificate or full-pipeline error validation.",
    )
    (HERE / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    methods = [
        ("conditional_rejection", "full_history", "Ordinary\nfull history"),
        ("conditional_chain", "full_history", "Chain\nfull history"),
        ("conditional_chain", "root_split", "Chain\nroot split"),
        ("conditional_importance", "root_split", "Weighted*\nroot split"),
        ("rebuild_known_covariance", "one_gaussian_population", "Rebuild\nknown covariance"),
        ("rebuild_fitted_covariance", "one_gaussian_population", "Rebuild*\nfitted covariance"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True, layout="constrained")
    for ax, scenario, color in zip(
        axes, ["null", "moderate", "strong"], ["#286c8e", "#bd6b22", "#338064"]
    ):
        for i, (method, event, label) in enumerate(methods):
            rows = table[
                (table.scenario == scenario) & (table.method == method) & (table.event == event)
            ].sort_values("replicate")
            p = rows.p_value.to_numpy()
            assert np.isfinite(p).all() and np.all((p >= -1e-12) & (p <= 1 + 1e-12))
            p = np.clip(p, 0, 1)
            ax.scatter(
                i + np.linspace(-0.16, 0.16, len(p)),
                np.maximum(p, 1e-6),
                color=color,
                s=31,
                alpha=0.8,
                edgecolor="white",
                linewidth=0.4,
            )
        ax.set_title(
            f"{scenario.capitalize()} ({len(table[table.scenario == scenario].replicate.unique())} datasets)"
        )
        ax.set_xticks(
            range(len(methods)), [m[2] for m in methods], rotation=42, ha="right", fontsize=9
        )
        ax.set_yscale("log")
        ax.set_ylim(8e-7, 1.4)
        ax.axhline(0.05, color="#555555", ls="--", lw=0.8)
        ax.axhline(0.01, color="#888888", ls=":", lw=0.8)
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Root p-value / approximate tail estimate (log scale)")
    fig.suptitle("Candidate replacement rules on the same diffusion trees", fontsize=15)
    fig.supxlabel(
        "Known-identity Gaussian controls, 9 rows × 3 features. Lines: 0.05 and 0.01.\n*Approximation; numerical event membership throughout. Conditional and rebuilding tests ask different questions.",
        fontsize=10,
    )
    fig.savefig(HERE / "comparison.png", dpi=170)
    fig.savefig(HERE / "comparison.pdf")
    plt.close(fig)
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
