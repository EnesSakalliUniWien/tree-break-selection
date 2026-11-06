"""Independent checks over saved records, source provenance and reported examples."""

import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

root = Path(__file__).parent
source = root / "generated/complete"
manifest = json.loads((root / "manifest.json").read_text())
assert all(
    hashlib.sha256(Path(p).read_bytes()).hexdigest() == h
    for p, h in manifest["source_sha256"].items()
)
summary = pd.read_csv(root / "summary.csv")
factors = pd.read_csv(root / "root_factors.csv")
np.testing.assert_allclose(
    factors.p, chi2.sf(factors.q / factors.c_hat, factors.df), rtol=1e-10, atol=1e-300
)
assert len(factors) == 8 * 3 * 5
for _, row in summary.iterrows():
    stem = f"{row['case']}-{row.seed}"
    edges = pd.read_csv(source / f"{stem}-edges.csv")
    siblings = pd.read_csv(source / f"{stem}-siblings.csv")
    assert len(edges) == 990 and len(siblings) == 495
    assert int(edges.Child_Parent_Divergence_Significant.sum()) == row.edge_rejections
    assert int(siblings.reject.sum()) == row.sibling_rejections
    assert not (edges.true_mean_null & edges.Child_Parent_Divergence_Significant).any()
    assert not (siblings.true_mean_null & siblings.reject).any()
    supported = siblings[siblings.supported]
    np.testing.assert_allclose(
        supported.loc[~supported.true_mean_null, "weighted_stat"].sum()
        / supported.weighted_stat.sum(),
        row.nonnull_stat_share,
    )
    assert row.support_enforced_status == "undefined_sparse_context"

siblings = pd.read_csv(source / "cancer_pc_hub_disconnected-20261106-siblings.csv")
support = siblings[siblings.supported]
family_c = (
    support.weighted_stat.sum() / (support.sibling_null_weight * support.degrees_of_freedom).sum()
)
assert len(support) == 494 and (~support.true_mean_null).sum() == 11
assert support.calibration_dependency_group.nunique() == 2
np.testing.assert_allclose(family_c, 14.84385583768377)
pure = support[support.true_mean_null]
pure_c = pure.weighted_stat.sum() / (pure.sibling_null_weight * pure.degrees_of_freedom).sum()
calibration = pd.read_csv(root / "calibration_factors.csv")
saved_pure = calibration[
    (calibration.case == "cancer_pc_hub_disconnected")
    & (calibration.seed == 20261106)
    & (calibration.variant == "known_mean_null_support_family_only")
].iloc[0]
np.testing.assert_allclose(pure_c, saved_pure.c_hat, rtol=1e-12)

traversal = pd.read_csv(source / "cancer_pc_hub_connected-20261104-traversal.csv")
boundary = traversal[traversal.node_id == "N985"].iloc[0]
members = [
    int(label.removeprefix("S")) for label in ast.literal_eval(boundary.descendant_leaf_signature)
]
truth = np.repeat(np.arange(6), [256, 64, 32, 16, 64, 64])
counts = np.bincount(truth[members], minlength=6)
np.testing.assert_array_equal(counts, [256, 0, 0, 0, 64, 64])
assert boundary.decision == "boundary" and boundary.true_mean_norm == 0

pairs = pd.read_csv(source / "cancer_pc_hub_disconnected-20261106-pairs.csv")
assert len(pairs) == 183
np.testing.assert_allclose(pairs.native_q_times_branch_factor, 30, atol=4e-5, rtol=0)
report = {
    "source_hashes_verified": len(manifest["source_sha256"]),
    "root_p_values_recomputed": len(factors),
    "edge_records_checked": 8 * 990,
    "sibling_records_checked": 8 * 495,
    "native_mean_null_rejections": 0,
    "contaminated_support_records": 11,
    "support_dependency_groups": 2,
    "restricted_family_c": pure_c,
    "connected_boundary_population_counts": counts.tolist(),
    "pair_q_times_branch_min": float(pairs.native_q_times_branch_factor.min()),
    "pair_q_times_branch_max": float(pairs.native_q_times_branch_factor.max()),
}
(root / "independent_verification.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
