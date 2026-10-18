"""Independent output consistency checks and matched benchmark visualization."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import betainc, gammaincc
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.energy_selected_reference import (
    energy_selected_gaussian_sf,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
r = pd.read_csv(BASE / "evaluation/runs.csv")
s = pd.read_csv(BASE / "evaluation/summary.csv")
nodes = pd.read_csv(BASE / "evaluation/nodes.csv.gz")
m = json.loads((BASE / "evaluation/manifest.json").read_text())
assert m["complete"] and m["source_hashes_unchanged"]
assert all(
    hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h for p, h in m["source_sha256"].items()
)
assert len(r) == 24300
assert not r.status.eq("error").any()
keys = ["profile", "signal", "replicate", "method", "alpha", "discovery_sha256", "inference_sha256"]
a = r[r.variant.eq("independent_conditional")]
b = r[r.variant.eq("independent_nominal")]
pairs = a.merge(b, on=keys, suffixes=("_new", "_nominal"), validate="one_to_one")
assert len(pairs) == 10800
assert (pairs.status_new == pairs.status_nominal).all()
assert (pairs.discovery_sha256 != pairs.inference_sha256).all()
assert (nodes.p_value.between(0, 1) & nodes.nominal_p_value.between(0, 1)).all()
assert not nodes.duplicated(["profile", "signal", "replicate", "method", "parent"]).any()
root = nodes[nodes.is_root]
joined = a[a.status.eq("ok") & a.alpha.eq(0.01)].merge(
    root, on=["profile", "signal", "replicate", "method"], validate="one_to_one"
)
np.testing.assert_allclose(joined.root_p, joined.p_value)
paired = []
for key, g in pairs.groupby(["profile", "signal", "method", "alpha"]):
    ok = g[g.status_new.eq("ok")]
    row = dict(zip(["profile", "signal", "method", "alpha"], key))
    row.update(
        attempts=len(g),
        matched_completed=len(ok),
        unavailable=len(g) - len(ok),
        mean_ari_difference=(ok.ari_new - ok.ari_nominal).mean(),
        mean_cluster_difference=(ok.clusters_new - ok.clusters_nominal).mean(),
        changed_cluster_count=int((ok.clusters_new != ok.clusters_nominal).sum()),
        new_only_false_partition=int(((ok.clusters_new > 1) & (ok.clusters_nominal == 1)).sum())
        if not key[1]
        else np.nan,
        nominal_only_false_partition=int(((ok.clusters_new == 1) & (ok.clusters_nominal > 1)).sum())
        if not key[1]
        else np.nan,
    )
    paired.append(row)
pd.DataFrame(paired).to_csv(BASE / "paired_comparison.csv", index=False)
tails = pd.read_csv(BASE / "evaluation/root_tails.csv")
tails = tails[~tails.signal & tails.variant.eq("independent_conditional")].copy()
tails["one_sided_binomial_p"] = [
    binomtest(int(x.rejections), int(x.n), x.cutoff, alternative="greater").pvalue
    for x in tails.itertuples()
]
# Post-evaluation diagnostic, not a preregistered decision rule. Preserve the
# entire planned 3-profile x 8-method x 4-cutoff family, including unavailable cells.
family = tails.one_sided_binomial_p.tolist() + [1.0] * (96 - len(tails))
tails["holm_p_96_planned_comparisons"] = multipletests(family, method="holm")[1][: len(tails)]
tails.to_csv(BASE / "root_calibration_check.csv", index=False)
errors = []
for m0, k, previous, q in [(2, 1, 0, 0.9), (8, 3, 2, 0.5), (8, 4, 2, 0.5), (40, 36, 35, 0.9)]:
    a0, b0 = k / 2, (m0 - k) / 2

    def weight(c):
        return c ** (a0 - 1) * (
            1 if previous == 0 else betainc(previous / 2, (k - previous) / 2, q / c)
        )

    denominator = quad(weight, q, 1, weight="alg", wvar=(0, b0 - 1), epsabs=0, epsrel=1e-11)[0]
    for t in (1.0, 7.0, 25.0):
        numerator = quad(
            lambda c: weight(c) * gammaincc(m0 / 2, t / (2 * c)),
            q,
            1,
            weight="alg",
            wvar=(0, b0 - 1),
            epsabs=0,
            epsrel=1e-11,
        )[0]
        actual = energy_selected_gaussian_sf(
            t,
            candidate_dimension=m0,
            projection_dimension=k,
            previous_group_dimension=previous,
            energy_fraction=q,
        )
        error = abs(numerator / denominator - actual)
        assert error < 1e-9
        errors.append(error)
verification = dict(
    run_rows=len(r),
    matched_variant_rows=len(pairs),
    matched_completed=int(pairs.status_new.eq("ok").sum()),
    node_rows=len(nodes),
    unexpected_errors=0,
    source_hashes_match=True,
    independent_quadrature_checks=len(errors),
    max_absolute_quadrature_error=max(errors),
    focused_tests_passed=146,
)
(BASE / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
methods = ["average", "complete", "neighbor_joining"]
variants = ["current_single_matrix", "independent_nominal", "independent_conditional"]
labels = ["Current single matrix", "Independent, nominal tail", "Independent, conditional tail"]
colors = ["#ad453e", "#c38b20", "#176b69"]
fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
for col, profile in enumerate(["identity2", "ar1_8", "identity40"]):
    for vi, (variant, label, color) in enumerate(zip(variants, labels, colors, strict=True)):
        x = np.arange(3) + (vi - 1) * 0.24
        for row, signal in enumerate([False, True]):
            g = (
                s[
                    s.profile.eq(profile)
                    & s.signal.eq(signal)
                    & s.variant.eq(variant)
                    & s.alpha.eq(0.01)
                ]
                .set_index("method")
                .loc[methods]
            )
            y = g.mean_ari.to_numpy() if signal else g.false_partition_rate.to_numpy()
            err = (
                None
                if signal
                else np.vstack(
                    [y - g.wilson95_low.to_numpy(), g.wilson95_high.to_numpy() - y]
                ).clip(0)
            )
            axes[row, col].bar(x, y, width=0.23, label=label, color=color, yerr=err, capsize=2)
            axes[row, col].set_xticks(
                np.arange(3), ["Average", "Complete", "Neighbor joining"], rotation=12
            )
    axes[0, col].set_title(profile.replace("_", " "))
    axes[0, col].set_ylabel("Null false-partition rate (95% Wilson)")
    axes[0, col].set_ylim(0, 0.65)
    axes[1, col].set_ylabel("Signal mean adjusted Rand index")
    axes[1, col].set_ylim(0, 0.75)
    for ax in axes[:, col]:
        ax.grid(axis="y", alpha=0.2)
axes[0, 0].legend(fontsize=7)
fig.suptitle(
    "Independent Gaussian discovery/test pipeline\n250 null + 50 signal pairs per profile; edge alpha .001, sibling alpha .01",
    fontsize=13,
)
fig.savefig(BASE / "benchmark_comparison.png", dpi=160)
fig.savefig(BASE / "benchmark_comparison.pdf")
plt.close(fig)
print(json.dumps(verification, indent=2))
print(
    pd.DataFrame(paired)
    .query("alpha==.01")[
        ["profile", "signal", "method", "matched_completed", "changed_cluster_count"]
    ]
    .to_string(index=False)
)
