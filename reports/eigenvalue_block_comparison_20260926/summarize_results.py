"""Validate and summarize this saved experiment, including rank sensitivity."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from sklearn.metrics import adjusted_rand_score
from tree_break_selection.space_separation import separate_adaptive_cosine_space

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "manifest.json").read_text())
assert manifest["complete"]
tables = {name: pd.read_csv(ROOT / f"{name}.csv") for name in (
    "cases", "spectra", "blocks", "block_pairs", "geometries", "runs", "assignments", "run_pairs",
)}
cases, spectra, blocks, pairs, geometries, runs, assignments, run_pairs = tables.values()
assert len(cases) == manifest["completed_cases"] == 11
assert len(runs) == len(blocks) * 16
assert not runs.duplicated(["case", "block", "geometry", "method"]).any()
for key, a in assignments.groupby(["case", "block", "geometry", "method"]):
    row = runs.set_index(["case", "block", "geometry", "method"]).loc[key]
    assert row.status == "ok"
    assert not a['sample'].duplicated().any()
    assert len(a) == int(cases.set_index("case").loc[key[0], "n_samples"])
    assert a.cluster.ge(0).all()
    assert a.cluster.nunique() == row.n_clusters
    assert np.isclose(adjusted_rand_score(a.truth, a.cluster), row.ari_truth)
for file, digest in manifest["source_sha256"].items():
    assert hashlib.sha256(Path(file).read_bytes()).hexdigest() == digest, file

rank_rows = []
for case in get_test_cases_by_suite("method_proof"):
    data = prepare_case_inputs(case, selected_methods=[]).data
    old = blocks[blocks.case.eq(case["name"])]
    record = cases.set_index("case").loc[case["name"]]
    assert hashlib.sha256(data.to_csv().encode()).hexdigest() == record.data_sha256
    new = separate_adaptive_cosine_space(data, weighting="binary", max_rank=80,
                                        min_segment_length=4, max_segments=8)
    old_eigenvalues = spectra[spectra.case.eq(case["name"])].eigenvalue.to_numpy()
    np.testing.assert_allclose(new.eigenvalues[:len(old_eigenvalues)], old_eigenvalues, rtol=1e-12)
    old_groups = np.empty(len(old_eigenvalues), dtype=int)
    new_groups = np.empty(len(new.eigenvalues), dtype=int)
    for row in old.itertuples():
        old_groups[row.block_start - 1:row.block_end] = row.block_id
    for block in new.blocks:
        new_groups[block.block_start - 1:block.block_end] = block.block_id
    rank_rows.append({
        "case": case["name"], "rank_40": len(old_groups), "rank_80": len(new_groups),
        "blocks_40": len(old), "blocks_80": len(new.blocks),
        "mode_grouping_ari_on_shared_modes": adjusted_rand_score(old_groups, new_groups[:len(old_groups)]),
        "block_ranges_40": "; ".join(f"{x.block_start}-{x.block_end}" for x in old.itertuples()),
        "block_ranges_80": "; ".join(f"{x.block_start}-{x.block_end}" for x in new.blocks),
        "energy_fraction_40": old.retained_energy_fraction.iloc[0],
        "energy_fraction_80": new.eigenvalues.sum() / len(data),
    })
rank = pd.DataFrame(rank_rows)
rank.to_csv(ROOT / "rank_sensitivity.csv", index=False)

# Adjacent boundaries only; nonadjacent pairs would inflate the apparent separation.
boundaries = []
for case, b in blocks.groupby("case", sort=False):
    b = b.sort_values("block_start")
    for i in range(len(b) - 1):
        left, right = b.iloc[i], b.iloc[i + 1]
        boundaries.append({"case": case, "left_block": left.block_name,
                           "right_block": right.block_name,
                           "boundary_ratio": left.eigenvalue_min / right.eigenvalue_max,
                           "left_width_ratio": left.eigenvalue_max / left.eigenvalue_min,
                           "right_width_ratio": right.eigenvalue_max / right.eigenvalue_min,
                           "common_boundary": left.block_type == "common_mode"})
boundaries = pd.DataFrame(boundaries)
boundaries.to_csv(ROOT / "adjacent_boundaries.csv", index=False)

lookup = runs.set_index(["case", "block", "geometry", "method"])
run_pairs["both_single_cluster"] = [
    lookup.loc[(r.case, r.left_block, r.left_geometry, r.left_method), "n_clusters"] == 1
    and lookup.loc[(r.case, r.right_block, r.right_geometry, r.right_method), "n_clusters"] == 1
    for r in run_pairs.itertuples()
]
run_pairs.to_csv(ROOT / "run_pairs_with_collapse.csv", index=False)
summary_rows = []
for case in cases.case:
    r = runs[runs.case.eq(case)]
    p = run_pairs[run_pairs.case.eq(case)]
    bb = boundaries[boundaries.case.eq(case) & ~boundaries.common_boundary]
    for geometry in ["raw", "adaptive_diffusion"]:
        rr = r[r.geometry.eq(geometry)]
        pp = p[p.comparison.eq("blocks") & p.left_geometry.eq(geometry)]
        active = pp[~pp.both_single_cluster & pp.ari_between_partitions.notna()]
        summary_rows.append({
            "case": case, "geometry": geometry,
            "blocks": len(blocks[blocks.case.eq(case)]),
            "median_noncommon_boundary_ratio": bb.boundary_ratio.median(),
            "ok_runs": rr.status.eq("ok").sum(), "total_runs": len(rr),
            "one_cluster_runs": rr.n_clusters.eq(1).sum(),
            "tree_pairs": pp.normalized_rooted_rf.count(),
            "median_between_block_rooted_rf": pp.normalized_rooted_rf.median(),
            "partition_pairs": pp.ari_between_partitions.count(),
            "median_between_block_partition_ari": pp.ari_between_partitions.median(),
            "pairs_with_at_least_one_multicluster_output": len(active),
            "median_ari_excluding_both_single_cluster": active.ari_between_partitions.median(),
        })
summary = pd.DataFrame(summary_rows)
summary.to_csv(ROOT / "summary_by_case.csv", index=False)
runs.groupby(["geometry", "method", "status"]).size().rename("runs").reset_index().to_csv(ROOT / "status_by_method.csv", index=False)

def table(frame):
    def cell(value):
        if isinstance(value, (float, np.floating)):
            return "—" if np.isnan(value) else f"{value:.3f}"
        return str(value)
    header = "| " + " | ".join(frame.columns) + " |"
    rule = "| " + " | ".join("---" for _ in frame.columns) + " |"
    return "\n".join([header, rule] + ["| " + " | ".join(map(cell, row)) + " |" for row in frame.itertuples(index=False, name=None)])

ordinary = boundaries[~boundaries.common_boundary]
bp = run_pairs[run_pairs.comparison.eq("blocks")]
tp = run_pairs[run_pairs.comparison.eq("tree_methods")]
ok = runs[runs.status.eq("ok")]
null_runs = runs[runs.case.isin(cases.loc[cases.true_clusters.eq(1), "case"])]
null_summary = null_runs.groupby(["case", "geometry"]).agg(
    attempts=("status", "size"), completed=("status", lambda v: v.eq("ok").sum()),
    multi_cluster_outputs=("n_clusters", lambda v: v.gt(1).sum()),
    maximum_clusters=("n_clusters", "max"),
).reset_index()
null_summary.to_csv(ROOT / "single_population_results.csv", index=False)
changed_rank = int((rank.mode_grouping_ari_on_shared_modes < 1 - 1e-12).sum())
body = f"""# Eigenvalue block and tree inference comparison

All 11 `method_proof` examples were evaluated on the same generated data and seeds across
seven linkage methods (average, complete, weighted, single, centroid, median, Ward) and
neighbor joining with MAD rooting. Each eigenvalue block was evaluated using raw Euclidean
coordinates and adaptive diffusion. There are {len(blocks)} blocks and {len(runs)} attempted runs.

## Findings

- Non-common-mode adjacent boundaries have a median eigenvalue ratio of
  {ordinary.boundary_ratio.median():.3f}, ranging from {ordinary.boundary_ratio.min():.3f} to
  {ordinary.boundary_ratio.max():.3f} across {len(ordinary)} boundaries. These are mostly
  changes in decay regime, not isolated groups separated by large empty eigenvalue intervals.
- Between blocks, the median normalized rooted tree distance is {bp.normalized_rooted_rf.median():.3f}
  over {bp.normalized_rooted_rf.count()} available pairs. Identical final labels need separate
  interpretation: {int((bp.both_single_cluster & bp.ari_between_partitions.notna()).sum())} of
  {bp.ari_between_partitions.count()} available partition pairs are both single-cluster outputs.
- Changing the tree inference method within the same block and geometry has median rooted
  tree distance {tp.normalized_rooted_rf.median():.3f} over {tp.normalized_rooted_rf.count()} pairs.
- {len(ok)} runs returned clusters; {int(ok.n_clusters.eq(1).sum())} returned one cluster.
  Status counts: `{runs.status.value_counts().to_dict()}`. Failure and unsupported rows are
  excluded from partition comparisons, never converted to one-cluster assignments.
- Single-population examples produced more than one cluster in
  {int(null_runs.n_clusters.gt(1).sum())}/{int(null_runs.status.eq('ok').sum())} completed
  runs. These are dependent method/block evaluations, not an estimated false-positive rate.
- Raising the target rank from 40 to 80 changes grouping among the already-retained modes
  in {changed_rank}/11 cases. Rank-80 tree inference was not run. See `rank_sensitivity.csv`.

## Figures and tables

- [Eigenvalue blocks](eigenvalue_blocks.png): every spectrum with existing block boundaries.
- [Benchmark cluster results](benchmark_clusters.pdf): one page per case, with spectra,
  TBS cluster counts and adjusted Rand index against generator labels for every method.
- [First binary example](binary_example.png).
- `summary_by_case.csv`: per-case results and denominators.
- `block_pairs.csv`: all pairwise eigenvalue separations; `adjacent_boundaries.csv`: adjacent gaps only.
- `geometry_pairs.csv`: matched sample-pair distance comparisons.
- `run_pairs_with_collapse.csv`: tree and partition comparisons with single-cluster matches flagged.
- `assignments.csv`, `trees.csv`, `runs.csv`: sample assignments, rooted clades and explicit outcomes.

## Per-case comparison

The tree distance is zero for identical rooted clades and one for no shared nontrivial clades.
Partition agreement is adjusted Rand index: one means identical partitions. All medians below
are descriptive; multiple pairs from one dataset are not independent observations.

{table(summary[['case', 'geometry', 'ok_runs', 'total_runs', 'median_between_block_rooted_rf', 'pairs_with_at_least_one_multicluster_output', 'median_ari_excluding_both_single_cluster']])}

### Single-population examples

{table(null_summary)}

## Measurement contract

The existing cosine operator is formed from row-normalized raw feature values. The name
`binary` in the weighting API means raw values, including signed continuous features.
Blocks are the existing contiguous, tie-preserving BIC segments of the log eigenspectrum.
The target rank is 40, minimum block length 4, maximum segments 8. Full-operator energy
fractions use trace = number of nonzero sample rows; retained-energy fractions are separate.

For an earlier block A and later block B, boundary ratio = min(A)/max(B).
Log-Wasserstein distance gives equal mass to each eigenvalue and compares log eigenvalues;
for disjoint ordered blocks it equals the difference of their mean log eigenvalues. It is
unchanged by a common rescaling. Within-block widths are reported separately.
Euclidean coordinates are Q_B sqrt(Lambda_B); distance Spearman correlation and RMS-normalized
distance error compare the same sample pairs. Unequal block dimensions remain visible.
Disjoint eigenspaces are already orthogonal, so principal angles would not measure useful
differences in their sample clustering.

Tree distance = size of the symmetric difference of nontrivial rooted clade sets divided by
the sum of their sizes. Neighbor joining includes its MAD-rooting choice. Cluster comparisons
use all original samples in their original order. Ground-truth labels are used only for
evaluation; no block, method, alpha or geometry was selected using them.

TBS gates use the original feature matrix and its FeatureSpace, edge alpha 0.001 and sibling
alpha 0.01. This follows the existing runner defaults with support-threshold enforcement off
and branch-length variance policy `none`. It does not establish selected-tail calibration,
false-positive control, or a winning method. `calibration_counts` remain in `runs.csv`.
Adaptive diffusion uses K=15, time=3, 30 components, Euclidean distance, bandwidth
`-1/(d+2)` and median epsilon. Failed geometry is recorded explicitly.

IQ-TREE takes state alignments rather than these block distances and is excluded. NNLS changes
branch lengths on a fixed topology; branch-time variants and other gate profiles were not run.
This is the complete method_proof panel, not the repository's full benchmark suite, and there
is one seed per example. The rank sensitivity check does not establish resampling stability.

## Reproduction and verification

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg uv run --no-sync python -m benchmarks.diagnostics.spectral.adaptive_cosine.compare_spectral_blocks --suite method_proof --output-dir <new-output-directory>
```

`manifest.json` records full case definitions, package versions, Git state and source hashes.
Verification regenerated all 11 datasets, matched their hashes and retained eigenvalues,
checked every successful assignment's sample coverage, cluster count and adjusted Rand index,
and confirmed source hashes did not change during inference. Analytic checks covered
log-distance values, common-scale invariance, proportional and constant distance vectors,
and rooted-clade distance. The 18 space-separation and 16 construction/rooting tests passed.
No production method, inference default, commit or push was changed by this analysis.
"""
(ROOT / "README.md").write_text(body)
print(body[:2400])
print("All saved-result checks and 11 rank-sensitivity comparisons passed.")
