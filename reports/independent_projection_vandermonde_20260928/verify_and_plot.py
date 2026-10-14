"""Independent arithmetic, stream and NumPy projection audit of the control."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation", type=Path, default=BASE / "control_evaluation")
    args = parser.parse_args()
    out = args.evaluation
    frozen = json.loads((out / "frozen_manifest.json").read_text())
    completion = json.loads((out / "completion.json").read_text())
    config = frozen["config"]
    assert completion["source_hashes_unchanged"]
    for name, digest in frozen["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    for name, digest in completion["output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest, name
    rows = pd.read_csv(out / "control_rows.csv.gz")
    datasets = pd.read_csv(out / "datasets.csv")
    keys = ["covariance_profile", "pattern", "replicate"]
    assert not rows.duplicated(keys + ["focal_node", "regime", "basis_source", "rank_source"]).any()
    assert not datasets.duplicated(keys).any()
    assert len(rows) == config["replicates"] * 486
    valid = rows[rows.status.eq("ok")]
    np.testing.assert_allclose(valid.p_value, chi2.sf(valid.statistic, valid.actual_rank), atol=2e-14, rtol=2e-11)
    original = pd.read_csv(ROOT / "reports/known_branch_skew_20260928/evaluation/known_branch_local_rows.csv")
    original = original[original.covariance_estimator.eq("oracle_covariance") & original.projection_kind.eq("same_data_pca_mp")]
    original = original.set_index(keys + ["focal_node"])
    same = rows[rows.regime.eq("canonical_mp") & rows.basis_source.eq("evaluation") & rows.rank_source.eq("evaluation_mp")].set_index(keys + ["focal_node"])
    original = original.loc[same.index]
    assert same.dataset_sha256.equals(original.dataset_sha256)
    np.testing.assert_allclose(same.p_value, original.p_value, atol=2e-14, rtol=2e-11)
    np.testing.assert_array_equal(same.actual_rank, original.projection_dimension)
    literal = {"global_null": [(0, 0)]*4, "root_only": [(1.25, 0)]*2 + [(-1.25, 0)]*2,
               "A_deep_only": [(1.25, 0), (-1.25, 0), (0, 0), (0, 0)],
               "separate_deep_branches": [(1.25, 0), (-1.25, 0), (0, 1), (0, -1)],
               "strong_root_weak_A": [(1.25, .25), (1.25, -.25), (-1.25, 0), (-1.25, 0)]}
    selected = rows[rows.replicate.isin([0, 17, 499, 999])].set_index(keys)
    checked, regenerated = 0, 0
    p, n = config["dimension"], config["n_per_group"]
    for d in datasets.itertuples():
        covariance = np.eye(p) if d.covariance_profile == "identity" else .6 ** np.abs(np.arange(p)[:, None]-np.arange(p))
        means = np.zeros((4, p))
        means[:, :2] = literal[d.pattern]
        sources = {}
        for source in ("evaluation", "basis_training", "rank_training"):
            key = tuple(json.loads(getattr(d, source + "_spawn_key")))
            if source != "evaluation":
                assert key == (17021983, 1 if source == "basis_training" else 2, *json.loads(d.evaluation_spawn_key))
            rng = np.random.default_rng(np.random.SeedSequence(config["base_seed"], spawn_key=key))
            data = np.vstack([rng.multivariate_normal(mean, covariance, size=n) for mean in means])
            assert hashlib.sha256(data.tobytes()).hexdigest() == getattr(d, source + "_sha256")
            sources[source] = data
            regenerated += 1
        if d.replicate not in (0, 17, 499, 999):
            continue
        subset = selected.loc[(d.covariance_profile, d.pattern, d.replicate)]
        for node, start, middle, stop in (("R", 0, 2*n, 4*n), ("A", 0, n, 2*n), ("B", 2*n, 3*n, 4*n)):
            chol = np.linalg.cholesky(covariance)
            data = sources["evaluation"]
            delta = data[start:middle].mean(axis=0) - data[middle:stop].mean(axis=0)
            z = np.linalg.solve(chol, delta) / np.sqrt(1/(middle-start) + 1/(stop-middle))
            eigen = {}
            for source, data in sources.items():
                part = data[start:stop]
                w = np.linalg.solve(chol, (part-part.mean(axis=0)).T).T
                values, vectors = np.linalg.eigh(w.T @ w / len(w))
                eigen[source] = (values[::-1], vectors[:, ::-1])
                assert np.min(np.abs(np.diff(values))) > 1e-10*np.max(np.abs(values)), "Audit sample requires explicit tie handling"
            for row in subset[subset.focal_node.eq(node)].itertuples():
                assert row.status == "ok"
                if row.rank_source.startswith("fixed_"):
                    rank = int(row.rank_source.split("_")[1])
                else:
                    source = "evaluation" if row.rank_source == "evaluation_mp" else "rank_training"
                    count = int(np.sum(eigen[source][0] > (1+np.sqrt(p/row.mp_threshold_rows))**2))
                    assert count == row.raw_mp_signal_count
                    rank = max(1, count)
                assert rank == row.requested_rank == row.actual_rank
                basis = np.eye(p) if row.basis_source == "fixed_coordinates" else eigen["evaluation" if row.basis_source == "evaluation" else "basis_training"][1]
                projected = basis[:, :rank].T @ z
                np.testing.assert_allclose(row.statistic, projected @ projected, atol=2e-9, rtol=1e-9)
                checked += 1
    summaries = pd.read_csv(out / "summary.csv")
    groupkeys = ["regime", "pattern", "covariance_profile", "focal_node", "basis_source", "rank_source"]
    expected = valid.assign(reject=valid.p_value.le(config["alpha"])).groupby(groupkeys).agg(rejections=("reject", "sum"), n_ok=("reject", "size")).sort_index()
    actual = summaries.set_index(groupkeys).loc[expected.index]
    np.testing.assert_array_equal(actual.rejections, expected.rejections)
    np.testing.assert_array_equal(actual.n_ok, expected.n_ok)
    results = {"passed": True, "regenerated_input_streams": regenerated, "matched_original_node_tests": len(same), "chi_square_rows": len(valid), "independent_numpy_projection_checks": checked, "summary_cells": len(expected), "frozen_source_hashes": len(frozen["source_sha256"]), "status_counts": rows.status.value_counts().to_dict()}
    (BASE / "verification.json").write_text(json.dumps(results, indent=2) + "\n")
    # The same basis/rank cross is shown separately for canonical and stress rows.
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=True)
    bases = ("evaluation", "independent", "fixed_coordinates")
    for i, regime in enumerate(("canonical_mp", "threshold_rows_stress")):
        for j, profile in enumerate(("identity", "ar1_0.6")):
            ax = axes[i, j]
            frame = summaries[(summaries.regime == regime) & (summaries.covariance_profile == profile) & (summaries.pattern == "global_null") & (summaries.focal_node == "R")]
            for offset, rank_source, color in ((-.08, "evaluation_mp", "#a14339"), (.08, "independent_mp", "#176d81")):
                for k, basis in enumerate(bases):
                    r = frame[(frame.basis_source == basis) & (frame.rank_source == rank_source)].iloc[0]
                    ax.errorbar(k+offset, r.rejection_rate, yerr=[[r.rejection_rate-r.wilson95_low], [r.wilson95_high-r.rejection_rate]], fmt="o", color=color, capsize=3, label=(rank_source.replace("_mp", " MP rank") if k == 0 else None))
            ax.axhline(.01, color="gray", linestyle="--", linewidth=1)
            ax.set_xticks(range(3), ["Same-data PCA", "Independent PCA", "Fixed coordinates"], rotation=10)
            ax.set_title(profile + (" / canonical MP" if i == 0 else " / changed-threshold stress"), fontsize=10)
            ax.set_ylabel("Null rejection rate")
            ax.grid(axis="y", alpha=.2)
    axes[0, 1].legend(fontsize=8)
    fig.suptitle("Separating basis data from rank data: fixed global-null roots\nPointwise Wilson 95% intervals; dashed line = 1%")
    fig.tight_layout()
    fig.savefig(BASE / "projection_controls.png", dpi=180)
    fig.savefig(BASE / "projection_controls.pdf")
    plt.close(fig)
    atoms = pd.read_csv(BASE / "vandermonde/atoms.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True)
    for ax, profile in zip(axes, ("identity", "ar1_0.6")):
        frame = atoms[atoms.case.eq("population_" + profile)]
        for anchor, color in ((0., "#a14339"), (.137, "#176d81")):
            selected_atoms = frame[frame.anchor.eq(anchor)]
            ax.vlines(selected_atoms.frequency, 0, selected_atoms.positive_weight, color=color, alpha=.6)
            ax.scatter(selected_atoms.frequency, selected_atoms.positive_weight, color=color, s=25, label=f"Anchor {anchor:g}")
        ax.set_title(profile)
        ax.set_xlabel("Frequency (cycles per feature index)")
        ax.set_xlim(-.02, 1)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Positive atom weight")
    fig.suptitle("Different Vandermonde atoms reconstruct the same full-rank covariance")
    fig.tight_layout()
    fig.savefig(BASE / "vandermonde_nonuniqueness.png", dpi=180)
    fig.savefig(BASE / "vandermonde_nonuniqueness.pdf")
    plt.close(fig)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
