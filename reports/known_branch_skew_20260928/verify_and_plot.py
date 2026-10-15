"""Independent arithmetic/traversal audit and figures for the saved evaluation."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.stats import chi2

BASE = Path(__file__).resolve().parent
EVALUATION = BASE / "evaluation"
REPOSITORY = BASE.parents[1]


def regenerate(manifest, local, selected):
    """Rebuild inputs and oracle statistics without importing the runner."""
    config = manifest["config"]
    p, n = config["dimension"], config["n_per_group"]
    assert config["ridge"] == 0
    means_by_pattern = {
        "global_null": [(0, 0)] * 4,
        "root_only": [(1.25, 0)] * 2 + [(-1.25, 0)] * 2,
        "A_deep_only": [(1.25, 0), (-1.25, 0), (0, 0), (0, 0)],
        "separate_deep_branches": [(1.25, 0), (-1.25, 0), (0, 1), (0, -1)],
        "strong_root_weak_A": [(1.25, 0.25), (1.25, -0.25), (-1.25, 0), (-1.25, 0)],
    }
    rows = local[
        local.covariance_estimator.eq("oracle_covariance") & local.projection_kind.eq("fixed_full")
    ].set_index(["covariance_profile", "pattern", "replicate", "focal_node"])
    selection = selected.set_index(["covariance_profile", "pattern", "replicate"])
    seeds = iter(np.random.SeedSequence(config["base_seed"]).spawn(len(selection)))
    checked = 0
    for profile in config["covariance_profiles"]:
        sigma = (
            np.eye(p)
            if profile == "identity"
            else 0.6 ** np.abs(np.arange(p)[:, None] - np.arange(p))
        )
        for pattern in config["patterns"]:
            means = np.zeros((4, p))
            means[:, :2] = means_by_pattern[pattern]
            for rep in range(config["replicates"]):
                seed = next(seeds)
                rng = np.random.default_rng(seed)
                data = np.vstack([rng.multivariate_normal(mu, sigma, size=n) for mu in means])
                digest = hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()
                for node, left, right in (
                    ("R", data[: 2 * n], data[2 * n :]),
                    ("A", data[:n], data[n : 2 * n]),
                    ("B", data[2 * n : 3 * n], data[3 * n :]),
                ):
                    row = rows.loc[(profile, pattern, rep, node)]
                    assert row.dataset_sha256 == digest
                    assert int(row.dataset_seed) == int(seed.generate_state(1, dtype=np.uint32)[0])
                    d = left.mean(axis=0) - right.mean(axis=0)
                    t = float(d @ np.linalg.solve(sigma, d) / (1 / len(left) + 1 / len(right)))
                    np.testing.assert_allclose(row.statistic, t, atol=2e-10, rtol=2e-11)
                    assert row.degrees_of_freedom == p
                # Independently extract the selected root directly from scipy's tree.
                tree = to_tree(linkage(data, method="average", metric="euclidean"))
                indices = [child.pre_order() for child in (tree.left, tree.right)]
                left, right = (data[idx] for idx in indices)
                d = left.mean(axis=0) - right.mean(axis=0)
                t = float(d @ np.linalg.solve(sigma, d) / (1 / len(left) + 1 / len(right)))
                row = selection.loc[(profile, pattern, rep)]
                assert row.dataset_sha256 == digest
                np.testing.assert_allclose(row.statistic, t, atol=2e-10, rtol=2e-11)
                counts = [np.bincount(np.asarray(idx) // n, minlength=4) for idx in indices]
                saved_counts = [
                    json.loads(row.left_leaf_group_counts),
                    json.loads(row.right_leaf_group_counts),
                ]
                names = ("A0", "A1", "B0", "B1")
                saved_vectors = [tuple(c[name] for name in names) for c in saved_counts]
                assert sorted(saved_vectors) == sorted(tuple(map(int, c)) for c in counts)
                exact = [
                    tuple(
                        sum(Fraction(str(means[g, j])) * int(c[g]) for g in range(4)) / int(c.sum())
                        for j in range(p)
                    )
                    for c in counts
                ]
                assert row.population_null == (exact[0] == exact[1])
                checked += 1
    return checked


def paired_results(local, alpha):
    keys = ["pattern", "covariance_profile", "focal_node", "replicate"]
    matched = local.pivot(
        index=keys, columns=["covariance_estimator", "projection_kind"], values="p_value"
    )
    contrasts = [
        (
            "covariance_at_full",
            ("oracle_covariance", "fixed_full"),
            ("pooled_sample_covariance", "fixed_full"),
        ),
        (
            "projection_at_oracle",
            ("oracle_covariance", "fixed_full"),
            ("oracle_covariance", "same_data_pca_mp"),
        ),
        (
            "projection_at_sample",
            ("pooled_sample_covariance", "fixed_full"),
            ("pooled_sample_covariance", "same_data_pca_mp"),
        ),
        (
            "covariance_at_projection",
            ("oracle_covariance", "same_data_pca_mp"),
            ("pooled_sample_covariance", "same_data_pca_mp"),
        ),
    ]
    summaries, all_pairs = [], []
    for label, before, after in contrasts:
        pairs = matched[[before, after]].copy()
        pairs.columns = ["before_p", "after_p"]
        pairs["comparison"] = label
        pairs["available"] = pairs.before_p.notna() & pairs.after_p.notna()
        pairs["reject_change"] = np.where(
            pairs.available,
            pairs.after_p.le(alpha).astype(int) - pairs.before_p.le(alpha).astype(int),
            np.nan,
        )
        all_pairs.append(pairs.reset_index())
        for key, group in pairs.groupby(level=keys[:3]):
            valid = group[group.available]
            delta = valid.reject_change
            summaries.append(
                dict(zip(keys[:3], key))
                | {
                    "comparison": label,
                    "n_pairs": len(valid),
                    "n_unavailable": len(group) - len(valid),
                    "gained_rejections": int(delta.eq(1).sum()),
                    "lost_rejections": int(delta.eq(-1).sum()),
                    "paired_rate_change": delta.mean(),
                    "paired_standard_error": delta.std(ddof=1) / np.sqrt(len(valid)),
                }
            )
    pd.concat(all_pairs, ignore_index=True).to_csv(BASE / "paired_stage_rows.csv", index=False)
    pd.DataFrame(summaries).to_csv(BASE / "paired_stage_summary.csv", index=False)


def wilson(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, center - half), min(1.0, center + half)


def main() -> None:
    local = pd.read_csv(EVALUATION / "known_branch_local_rows.csv")
    selected = pd.read_csv(EVALUATION / "known_branch_selected_root_rows.csv")
    traversal = pd.read_csv(EVALUATION / "known_branch_traversal_rows.csv")
    manifest = json.loads((EVALUATION / "manifest.json").read_text())
    alpha = manifest["config"]["alpha"]
    keys = ["pattern", "covariance_profile", "replicate"]
    local_keys = keys + ["focal_node", "covariance_estimator", "projection_kind"]
    assert not local.duplicated(local_keys).any()
    assert not selected.duplicated(keys).any()
    assert not traversal.duplicated(keys + ["passthrough"]).any()
    assert (local.groupby(keys).size() == 12).all()
    assert (traversal.groupby(keys).size() == 2).all()
    assert len(local) == 12 * len(selected)
    assert len(traversal) == 2 * len(selected)
    assert local.groupby(keys).dataset_sha256.nunique().eq(1).all()
    assert traversal.groupby(keys).dataset_sha256.nunique().eq(1).all()
    assert (
        local.groupby(keys)
        .dataset_sha256.first()
        .equals(traversal.groupby(keys).dataset_sha256.first())
    )
    checks = {"datasets": len(selected), "local_rows": len(local), "traversal_rows": len(traversal)}
    for path, digest in manifest["source_sha256"].items():
        assert hashlib.sha256((REPOSITORY / path).read_bytes()).hexdigest() == digest, path
    checks["verified_source_hashes"] = len(manifest["source_sha256"])
    checks["regenerated_datasets_and_oracle_statistics"] = regenerate(manifest, local, selected)
    paired_results(local, alpha)
    for name, table in (("local", local), ("selected_root", selected)):
        valid = table[table.status.eq("ok")]
        np.testing.assert_allclose(
            valid.p_value,
            chi2.sf(valid.statistic, valid.degrees_of_freedom),
            atol=2e-14,
            rtol=2e-12,
        )
        np.testing.assert_allclose(
            valid.statistic_per_df,
            valid.statistic / valid.degrees_of_freedom,
            atol=2e-13,
            rtol=2e-12,
        )
        checks[f"{name}_arithmetic_rows"] = len(valid)
        checks[f"{name}_statuses"] = table.status.value_counts().to_dict()
    expected = {
        "global_null": {"R": (0.0, 0.0), "A": (0.0, 0.0), "B": (0.0, 0.0)},
        "root_only": {"R": (2.5, 0.0), "A": (0.0, 0.0), "B": (0.0, 0.0)},
        "A_deep_only": {"R": (0.0, 0.0), "A": (2.5, 0.0), "B": (0.0, 0.0)},
        "separate_deep_branches": {"R": (0.0, 0.0), "A": (2.5, 0.0), "B": (0.0, 2.0)},
        "strong_root_weak_A": {"R": (2.5, 0.0), "A": (0.0, 0.5), "B": (0.0, 0.0)},
    }
    for (pattern, node), frame in local.groupby(["pattern", "focal_node"]):
        want = np.zeros(manifest["config"]["dimension"])
        want[:2] = expected[pattern][node]
        contrasts = np.array(frame.population_contrast.map(json.loads).tolist())
        np.testing.assert_allclose(
            contrasts, np.broadcast_to(want, contrasts.shape), atol=0, rtol=0
        )
        assert frame.population_null.eq(not np.any(want)).all()
        left = np.array(frame.left_population_mean.map(json.loads).tolist())
        right = np.array(frame.right_population_mean.map(json.loads).tolist())
        np.testing.assert_allclose(left - right, contrasts, atol=0, rtol=0)
    oracle = local[
        local.covariance_estimator.eq("oracle_covariance")
        & local.projection_kind.eq("fixed_full")
        & local.status.eq("ok")
    ]
    p = oracle.pivot(index=keys, columns="focal_node", values="p_value")
    truth = oracle.pivot(index=keys, columns="focal_node", values="population_null")
    audited = 0
    recovery = []
    for row in traversal[traversal.status.eq("ok")].itertuples():
        key = (row.pattern, row.covariance_profile, row.replicate)
        values = p.loc[key]
        root_open = values.R <= alpha
        ab = np.array([values.A, values.B])
        # Independent reimplementation of BH for exactly two depth-one tests.
        child_open = np.ones(2, dtype=bool) if ab.max() <= alpha else ab <= alpha / 2
        descend = bool(root_open or (row.passthrough and child_open.any()))
        flags = np.array([root_open, descend and child_open[0], descend and child_open[1]])
        nulls = truth.loc[key, ["R", "A", "B"]].to_numpy(dtype=bool)
        nfalse = int(np.sum(flags & nulls))
        assert row.n_reached_splits == flags.sum()
        assert row.n_reached_false_splits == nfalse
        assert np.isclose(row.reached_false_split_fdp, nfalse / max(1, flags.sum()))
        assert row.num_clusters == (2 + int(child_open.sum()) if descend else 1)
        trace = json.loads(row.traversal_decisions)
        assert trace["R"] == ("split" if root_open else "pass_through" if descend else "boundary")
        for i, node in enumerate(("A", "B")):
            assert (node in trace) == descend
            if descend:
                assert trace[node] == ("split" if child_open[i] else "boundary")
        for i, node in enumerate(("R", "A", "B")):
            recovery.append(
                {
                    "pattern": row.pattern,
                    "covariance_profile": row.covariance_profile,
                    "passthrough": row.passthrough,
                    "replicate": row.replicate,
                    "node": node,
                    "population_null": bool(nulls[i]),
                    "reached": node in trace,
                    "split": bool(flags[i]),
                    "pass_through": trace.get(node) == "pass_through",
                    "num_clusters": row.num_clusters,
                }
            )
        audited += 1
    checks["independently_reconstructed_traversals"] = audited
    checks["traversal_statuses"] = traversal.status.value_counts().to_dict()
    checks["scope"] = (
        "Regenerated inputs, direct oracle quadratic forms and selected-root composition truth, saved arithmetic, source hashes, and independent three-node BH/traversal reconstruction; not a general calibration proof."
    )
    recovery = pd.DataFrame(recovery)
    recovery.to_csv(BASE / "node_recovery_rows.csv", index=False)
    recovery_summary = []
    for key, group in recovery.groupby(
        ["pattern", "covariance_profile", "passthrough", "node", "population_null"]
    ):
        k, n = int(group.split.sum()), len(group)
        lo, hi = wilson(k, n)
        recovery_summary.append(
            dict(
                zip(
                    ["pattern", "covariance_profile", "passthrough", "node", "population_null"], key
                )
            )
            | {
                "n": n,
                "reached": int(group.reached.sum()),
                "split_count": k,
                "split_rate": k / n,
                "split_wilson_low": lo,
                "split_wilson_high": hi,
                "pass_through_count": int(group.pass_through.sum()),
                "mean_num_clusters": group.num_clusters.mean(),
            }
        )
    pd.DataFrame(recovery_summary).to_csv(BASE / "node_recovery_summary.csv", index=False)
    checks["passed"] = True
    (BASE / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
    profiles = list(manifest["config"]["covariance_profiles"])
    labels = [
        "Known covariance / full",
        "Sample pooled covariance / full",
        "Known covariance / PCA–MP",
        "Sample pooled covariance / PCA–MP",
        "Selected root / known covariance",
    ]
    variants = [
        ("oracle_covariance", "fixed_full"),
        ("pooled_sample_covariance", "fixed_full"),
        ("oracle_covariance", "same_data_pca_mp"),
        ("pooled_sample_covariance", "same_data_pca_mp"),
    ]
    fig, axes = plt.subplots(1, len(profiles), figsize=(12, 4.5), sharey=True, squeeze=False)
    summary = []
    for ax, profile in zip(axes[0], profiles):
        candidates = [
            local[
                local.pattern.eq("global_null")
                & local.covariance_profile.eq(profile)
                & local.focal_node.eq("R")
                & local.covariance_estimator.eq(cov)
                & local.projection_kind.eq(proj)
            ]
            for cov, proj in variants
        ]
        candidates.append(
            selected[selected.pattern.eq("global_null") & selected.covariance_profile.eq(profile)]
        )
        for y, (label, frame) in enumerate(zip(labels, candidates)):
            valid = frame[frame.status.eq("ok")]
            if valid.empty:
                continue
            n = len(valid)
            k = int(valid.p_value.le(alpha).sum())
            lo, hi = wilson(k, n)
            value = k / n
            ax.errorbar(
                value, y, xerr=[[value - lo], [hi - value]], fmt="o", color="#28647c", capsize=4
            )
            ax.annotate(
                f"{k}/{n}",
                (value, y),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )
            summary.append(
                {
                    "covariance_profile": profile,
                    "variant": label,
                    "rejects": k,
                    "n_ok": n,
                    "n_unavailable": len(frame) - n,
                    "rate": value,
                    "wilson_low": lo,
                    "wilson_high": hi,
                }
            )
        ax.axvline(alpha, color="#ad462c", linestyle="--", label=f"Nominal alpha {alpha:g}")
        ax.set_title(profile)
        ax.set_ylim(-0.6, len(labels) - 0.65)
        ax.set_xscale("symlog", linthresh=0.002)
        ax.set_xlim(0, 1.6)
        ax.set_xticks([0, 0.01, 0.1, 1], ["0", "1%", "10%", "100%"])
        ax.set_xlabel("Null rejection rate (pointwise Wilson 95% interval)")
        ax.grid(axis="x", alpha=0.2)
    axes[0, 0].set_yticks(range(len(labels)), labels)
    axes[0, 0].invert_yaxis()
    fig.suptitle(
        f"Known-root controls and same-data selected roots\nDashed line: nominal alpha {alpha:g}",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(BASE / "null_stage_comparison.png", dpi=180, bbox_inches="tight")
    fig.savefig(BASE / "null_stage_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(summary).to_csv(BASE / "null_root_comparison.csv", index=False)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
