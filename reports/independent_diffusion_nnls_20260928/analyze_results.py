"""Verify matched benchmark records and render primary-alpha comparisons."""

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "evaluation"
manifest = json.loads((DATA / "manifest.json").read_text())
assert manifest["complete"] and manifest["source_hashes_unchanged"]
assert all(
    hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == value
    for path, value in manifest["source_sha256"].items()
)
runs = pd.read_csv(DATA / "runs.csv")
nodes = pd.read_csv(DATA / "nodes.csv.gz")
meta = pd.DataFrame(json.loads(x) for x in (DATA / "fitted_context.jsonl").read_text().splitlines())
keys = ["profile", "signal", "replicate"]
assert len(runs) == 18_000 and len(nodes) == 83_700 and len(meta) == 4_500
assert runs.status.eq("ok").all(), runs.status.value_counts().to_dict()
assert not runs.duplicated(keys + ["geometry", "reference", "alpha"]).any()
assert runs.groupby(keys).size().eq(20).all()
for name in ["discovery_sha256", "inference_sha256"]:
    assert runs.groupby(keys)[name].nunique().eq(1).all()
assert runs.discovery_sha256.ne(runs.inference_sha256).all()
assert nodes[["p_value", "nominal_p_value"]].apply(lambda s: s.between(0, 1).all()).all()
ind = meta[meta.geometry.isin(["direct", "diffusion", "diffusion_nnls"])]
assert ind.frozen_context_unchanged.eq(True).all()


def ordered(table, geometry, extra=()):
    return table[table.geometry.eq(geometry)].set_index(keys + list(extra)).sort_index()


a = ordered(nodes, "diffusion", ["parent"])
b = ordered(nodes, "diffusion_nnls", ["parent"])
pd.testing.assert_frame_equal(
    a.drop(columns="geometry"), b.drop(columns="geometry"), check_exact=True
)
am = ordered(meta, "diffusion")
bm = ordered(meta, "diffusion_nnls")
for name in ["topology_sha256", "pca_sha256", "diffusion_sha256"]:
    assert am[name].equals(bm[name])
branch_changes = int(am.branches_sha256.ne(bm.branches_sha256).sum())
assert branch_changes == 900
ar = ordered(runs, "diffusion", ["reference", "alpha"])
br = ordered(runs, "diffusion_nnls", ["reference", "alpha"])
pd.testing.assert_frame_equal(
    ar.drop(columns=["geometry", "elapsed_seconds"]),
    br.drop(columns=["geometry", "elapsed_seconds"]),
    check_exact=True,
)
direct = ordered(meta, "direct")
topology_changes = int(direct.topology_sha256.ne(am.topology_sha256).sum())
paired = runs[runs.reference.eq("conditional")].merge(
    runs[runs.reference.eq("nominal")],
    on=keys + ["geometry", "alpha"],
    suffixes=("_conditional", "_nominal"),
    validate="one_to_one",
)
assert len(paired) == 8100
paired["cluster_count_changed"] = paired.clusters_conditional.ne(paired.clusters_nominal)
paired["ari_change"] = paired.ari_conditional - paired.ari_nominal
equation = (
    paired.groupby(["profile", "signal", "geometry", "alpha"])
    .agg(
        pairs=("cluster_count_changed", "size"),
        cluster_count_changes=("cluster_count_changed", "sum"),
        mean_ari_change=("ari_change", "mean"),
    )
    .reset_index()
)
equation.to_csv(HERE / "equation_effects.csv", index=False)
nnls = pd.DataFrame(bm.nnls.tolist())
assert nnls.status.eq("ok").all() and nnls.applied_to_tree.all()
solver = nnls[
    ["optimality", "iterations", "residual_rmse_to_target_mean", "branch_length_min"]
].describe()
solver.to_csv(HERE / "nnls_diagnostics.csv")
warnings = sorted({w for items in meta.warnings for w in items})
verification = dict(
    source_hashes_match=True,
    paired_datasets=900,
    run_rows=len(runs),
    node_rows=len(nodes),
    status_counts=runs.status.value_counts().to_dict(),
    independent_frozen_contexts=len(ind),
    unchanged_diffusion_nnls_node_records=len(a),
    unchanged_diffusion_nnls_run_pairs=len(ar),
    changed_nnls_branch_fits=branch_changes,
    changed_direct_vs_diffusion_topologies=topology_changes,
    matched_equation_pairs=len(paired),
    primary_alpha_equation_cluster_count_changes=int(
        paired.loc[paired.alpha.eq(0.01), "cluster_count_changed"].sum()
    ),
    warnings=warnings,
    nnls_solver_messages=nnls.solver_message.value_counts().to_dict(),
)
(HERE / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
summary = pd.read_csv(DATA / "summary.csv")
primary = summary[summary.alpha.eq(0.01)]
primary.to_csv(HERE / "primary_comparison.csv", index=False)
profiles = ["identity2", "ar1_8", "identity40"]
variants = [
    ("direct", "conditional", "Independent direct"),
    ("diffusion", "conditional", "Independent diffusion"),
    ("diffusion_nnls", "conditional", "Independent diffusion + NNLS"),
    ("tbs_diffusion_adaptive", "current_single_matrix", "Current diffusion"),
    ("tbs_diffusion_adaptive_nnls", "current_single_matrix", "Current diffusion + NNLS"),
]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
colors = ["#0072B2", "#009E73", "#CC79A7", "#E69F00", "#D55E00"]
for j, (geometry, reference, label) in enumerate(variants):
    group = primary[primary.geometry.eq(geometry) & primary.reference.eq(reference)]
    for ax, signal, value in zip(axes, [False, True], ["false_partition_rate", "mean_ari"]):
        g = group[group.signal.eq(signal)].set_index("profile").loc[profiles]
        x = np.arange(3) + (j - 2) * 0.16
        errors = (
            None
            if signal
            else np.stack(
                [
                    g[value] - g.wilson95_low,
                    g.wilson95_high - g[value],
                ]
            )
        )
        # CSV roundoff can put a Wilson endpoint infinitesimally above zero.
        if errors is not None:
            errors = np.maximum(errors, 0)
        ax.bar(x, g[value], width=0.15, color=colors[j], label=label, yerr=errors, capsize=2)
for ax, title, ylabel in zip(
    axes,
    ["Global null: 250 pairs per profile", "Signal: 50 pairs per profile"],
    ["Fraction with more than one cluster (95% Wilson)", "Mean adjusted Rand index"],
):
    ax.set(title=title, ylabel=ylabel, xticks=np.arange(3), xticklabels=profiles)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.2)
axes[0].set_ylim(0, 0.06)
axes[1].set_ylim(-0.03, 1)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside lower center", ncols=3, frameon=False, fontsize=9)
fig.suptitle("Diffusion NNLS: sibling alpha 0.01, edge alpha 0.001, passthrough enabled")
fig.savefig(HERE / "comparison.png", dpi=180)
fig.savefig(HERE / "comparison.pdf")
print(json.dumps(verification, indent=2))
print(primary.to_string(index=False))
