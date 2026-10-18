from pathlib import Path

import pandas as pd
from scipy.stats import chi2

HERE = Path(__file__).parent
b = pd.read_csv(HERE / "bug_checks.csv").set_index("case_id")
groups = {
    "1 tree holds clusters": [
        "binary_balanced_low_noise__2",
        "binary_many_features",
        "binary_noise_feat_30i_500n",
        "binary_noise_feat_80i_400n",
        "phylo_divergent_8taxa",
        "phylo_dna_16taxa_low_mut",
        "phylo_protein_12taxa",
        "phylo_large_32taxa",
        "phylo_large_64taxa",
        "gauss_dense_signal_highd",
        "gauss_dense_signal_highd_continuous",
        "gauss_extreme_noise_many",
        "binary_perfect_2c",
        "binary_perfect_4c",
        "binary_perfect_8c",
        "gauss_clear_medium_continuous",
        "gauss_moderate_3c_continuous",
        "gauss_single_outlier_4c_continuous",
        "gauss_outlier_cluster_4c_continuous",
        "cat_mod_4cat_6c",
        "binary_low_noise_2c",
    ],
    "2 shared subspace drowned": [
        "dim_consolidated_4c_24f",
        "dim_consolidated_4c_72f",
        "dim_consolidated_4c_272f",
        "dim_consolidated_4c_72f_continuous",
    ],
    "3 pair-specific subspaces": [
        "dim_diffuse_6c_36f",
        "dim_diffuse_6c_136f",
        "dim_diffuse_6c_536f",
        "dim_diffuse_6c_136f_continuous",
    ],
    "4 metric or linkage": ["cat_mod_3cat_4c", "phylo_brownian_null_16taxa"],
    "5 spread over weak features": [
        "bar_binary_balanced_4c",
        "bar_binary_unbalanced_95_5",
        "traversal_deep_signal_under_same_parent",
        "cat_highcard_10cat_4c",
        "cat_highcard_20cat_4c",
    ],
    "6 no location signal": [
        "mp_spike_above_bbp_continuous",
        "cont_lowrank_pggn_shrinkage",
        "cat_simplex_face_rare_20cat",
        "cat_dirichlet_overdispersed_20cat",
        "sbm_moderate",
        "sbm_hard",
        "overlap_heavy_4c_small_feat",
        "overlap_heavy_4c_med_feat",
        "overlap_heavy_8c_large_feat",
        "overlap_extreme_4c",
        "gauss_sparse_signal_highd_noise",
    ],
    "7 null case": ["gauss_null_small", "gauss_null_large"],
}
group_of = {c: g for g, cs in groups.items() for c in cs}
assert set(group_of) == set(b.index), set(b.index) ^ set(group_of)
plumbing_cols = [
    "labels_match_meta_k",
    "data_finite",
    "index_unique",
    "distance_finite",
    "tree_leaves_match_index",
    "tree_binary",
]
rows = []
for case, r in b.iterrows():
    plumbing = bool(r[plumbing_cols].all()) and (
        r.distance_max_abs_diff_independent < 1e-9
        if pd.notna(r.distance_max_abs_diff_independent)
        else True
    )
    if r.labels_returned:
        plumbing &= bool(r.ari_matches_row) and r.clusters_not_matching_tree == 0
    if r.status == "unsupported":
        plumbing &= r.unsupported_admissible_reported == r.sibling_admissible_recomputed
    over = r.internal_edges_rejected_share >= 0.9
    capped = False
    if r.family == "continuous":
        # The root edge cannot reach the 0.999 cutoff whatever the data.
        capped = chi2.ppf(0.999, r.root_edge_df) > r.n - 1
    mahal = r.distance_metric_used == "mahalanobis_time"
    possible = 0.8 <= r.internal_edges_rejected_share < 0.9
    bugs = [
        name
        for name, on in (
            ("edge basis fitted to tested leaves", over),
            ("continuous root edge cannot reject (statistic capped at n-1)", capped),
            ("mahalanobis_time erases separation", mahal),
        )
        if on
    ]
    if group_of[case].startswith("1") and bugs:
        role = "bug is the cause: the tree already holds the clusters"
    elif group_of[case].startswith("7") and bugs:
        role = "bug is the cause: null data should give one cluster"
    elif bugs and group_of[case][0] in "23456":
        role = "bug present; the data or geometry would still limit recovery"
    elif possible:
        role = "edge over-rejection possibly contributes"
    else:
        role = "no bug implicated"
    rows.append(
        {
            "case_id": case,
            "group": group_of[case],
            "status": r.status,
            "found": r.found,
            "plumbing_checks_pass": plumbing,
            "internal_edges_rejected_share": round(r.internal_edges_rejected_share, 3),
            "bugs_implicated": "; ".join(bugs) or "none",
            "role": role,
        }
    )
v = pd.DataFrame(rows).sort_values(["group", "case_id"])
v.to_csv(HERE / "bug_verdicts.csv", index=False)
pd.set_option("display.width", 260)
pd.set_option("display.max_rows", 60)
pd.set_option("display.max_colwidth", 70)
print(v.to_string(index=False))
print(v.role.value_counts().to_string())
print("plumbing all pass:", v.plumbing_checks_pass.all())
