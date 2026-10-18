"""Reproduce matched before/after comparisons; never modify production sources."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OLD = ROOT / "reports/eigenvalue_block_comparison_20260926"
sys.path.insert(0, str(ROOT))

from benchmarks.shared.util.case_inputs import prepare_case_inputs
from tree_break_selection.eigenvalue_groups import (
    extend_rank_through_ties,
    tied_eigenvalue_neighbors,
)
from tree_break_selection.space_separation import adaptive_cosine, diffusion

KEY = ["case", "seed", "block", "geometry", "method"]
GEOMETRY_KEY = KEY[:-1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def load_snapshot(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def historical_sources(manifest: dict) -> tuple[object, object]:
    """Recover old modules only when their bytes match the historical manifest.

    Diffusion is exactly the recorded Git HEAD. Cosine included an uncommitted
    tie fix; reverse only its later helper extraction, then verify the exact
    historical hash. An incorrect reconstruction fails before it is executed.
    """
    folder = OUT / "historical_sources"
    folder.mkdir(exist_ok=True)
    cosine_path = "tree_break_selection/space_separation/adaptive_cosine.py"
    diffusion_path = "tree_break_selection/space_separation/diffusion.py"
    old_diffusion = subprocess.check_output(
        ["git", "show", manifest["git_head"] + ":" + diffusion_path],
        cwd=ROOT,
    )
    old_cosine = (
        (ROOT / cosine_path)
        .read_text()
        .replace(
            "from tree_break_selection.eigenvalue_groups import (\n"
            "    extend_rank_through_ties,\n    tied_eigenvalue_neighbors,\n)\n\n",
            "",
        )
    )
    old_cosine = old_cosine.replace(
        "tied_eigenvalue_neighbors(eigenvalues)",
        "_tied_eigenvalue_neighbors(eigenvalues)",
    )
    helper = '''def _tied_eigenvalue_neighbors(eigenvalues: np.ndarray) -> np.ndarray:
    """Identify adjacent numerical ties, joined transitively into whole groups.

    Use relative tolerance 1e-10 and absolute tolerance 1e-12 times the
    spectral scale so eigensolver roundoff does not split repeated eigenvalues.
    """

    scale = float(np.max(np.abs(eigenvalues))) if len(eigenvalues) else 0.0
    return np.isclose(
        eigenvalues[:-1], eigenvalues[1:], rtol=1e-10, atol=1e-12 * scale
    )


'''
    old_cosine = (
        old_cosine.replace(
            "def cosine_eigendecomposition(",
            helper + "def cosine_eigendecomposition(",
        )
        .replace(
            "    rank = extend_rank_through_ties(eigenvalues, max_rank)\n",
            "    rank = min(max_rank, len(eigenvalues))\n"
            "    tied = _tied_eigenvalue_neighbors(eigenvalues)\n"
            "    while 0 < rank < len(eigenvalues) and tied[rank - 1]:\n"
            "        rank += 1\n",
        )
        .encode()
    )
    for name, content in [(cosine_path, old_cosine), (diffusion_path, old_diffusion)]:
        assert hashlib.sha256(content).hexdigest() == manifest["source_sha256"][name], name
        target = folder / Path(name).name
        if target.exists():
            assert target.read_bytes() == content
        else:
            target.write_bytes(content)
    return (
        load_snapshot("historical_cosine_for_reassessment", folder / "adaptive_cosine.py"),
        load_snapshot("historical_diffusion_for_reassessment", folder / "diffusion.py"),
    )


def geometry_replay() -> None:
    manifest = read_json(OLD / "manifest.json")
    source_paths = [
        Path(adaptive_cosine.__file__),
        Path(diffusion.__file__),
        ROOT / "tree_break_selection/eigenvalue_groups.py",
    ]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in source_paths}
    if (OUT / "geometry_replay.json").exists():
        assert read_json(OUT / "geometry_replay.json")["current_source_sha256"] == hashes, (
            "Source differs from the working tree recorded in geometry_replay.json; "
            "no commit contains that tree, so this replay cannot run on later source."
        )
        return
    before_cosine, before_diffusion = historical_sources(manifest)
    saved_cases = pd.read_csv(OLD / "cases.csv").set_index("case")
    saved_spectra = pd.read_csv(OLD / "spectra.csv")
    rows = []
    start = perf_counter()
    for case in manifest["cases"]:
        data = prepare_case_inputs(case, selected_methods=[]).data
        data_hash = hashlib.sha256(data.to_csv().encode()).hexdigest()
        assert data_hash == saved_cases.loc[case["name"], "data_sha256"]
        # Both recorded modules name raw weighting "binary"; the current API calls it "raw".
        old_space, new_space = [
            module.separate_adaptive_cosine_space(
                data,
                weighting="binary",
                max_rank=40,
                min_segment_length=4,
                max_segments=8,
            )
            for module in (before_cosine, adaptive_cosine)
        ]
        np.testing.assert_allclose(
            old_space.eigenvalues,
            saved_spectra.loc[saved_spectra.case == case["name"], "eigenvalue"],
            rtol=1e-12,
            atol=1e-12,
        )
        assert [vars(b) for b in old_space.blocks] == [vars(b) for b in new_space.blocks]
        np.testing.assert_array_equal(old_space.eigenvalues, new_space.eigenvalues)
        for old_block, new_block in zip(old_space.blocks, new_space.blocks, strict=True):
            old_coords, new_coords = (
                old_space.coordinates(old_block),
                new_space.coordinates(new_block),
            )
            for geometry in ("raw", "adaptive_diffusion"):
                row = {
                    "case": case["name"],
                    "seed": case["seed"],
                    "block": old_block.block_name,
                    "geometry": geometry,
                    "data_sha256": data_hash,
                    "cosine_coordinates_exact": np.array_equal(old_coords, new_coords),
                }
                distances = []
                for label, module, coordinates in (
                    ("before", before_diffusion, old_coords),
                    ("after", diffusion, new_coords),
                ):
                    try:
                        if geometry == "raw":
                            distance = pdist(coordinates)
                            dimension = coordinates.shape[1]
                        else:
                            result = module.block_adaptive_diffusion_geometry(
                                coordinates,
                                k_neighbors=15,
                                diffusion_time=3,
                                n_components=30,
                                metric="euclidean",
                                bandwidth_type="-1/(d+2)",
                                epsilon="median",
                            )
                            distance, dimension = (
                                result.distance_condensed,
                                result.coordinates.shape[1],
                            )
                        row[label + "_status"] = "ok"
                        row[label + "_dimension"] = dimension
                        row[label + "_distance_sha256"] = hashlib.sha256(
                            distance.tobytes()
                        ).hexdigest()
                        distances.append(distance)
                    except Exception as exc:
                        row[label + "_status"] = "failed_geometry"
                        row[label + "_error"] = repr(exc)
                        distances.append(None)
                a, b = distances
                if a is not None and b is not None:
                    row["distances_exact"] = np.array_equal(a, b)
                    row["distances_close"] = np.allclose(a, b, rtol=1e-10, atol=1e-12)
                    row["max_absolute_distance_change"] = float(np.max(np.abs(a - b)))
                    scale = float(np.sqrt(np.mean(a * a)))
                    row["distance_rmse_over_before_rms"] = float(
                        np.sqrt(np.mean((a - b) ** 2)) / scale
                    )
                    row["distance_spearman"] = float(spearmanr(a, b).statistic)
                rows.append(row)
        print("Geometry replay:", case["name"], flush=True)
    for path in source_paths:
        assert digest(path) == hashes[str(path.relative_to(ROOT))]
    pd.DataFrame(rows).to_csv(OUT / "geometry_comparison.csv", index=False)
    write_json(
        OUT / "geometry_replay.json",
        {
            "historical_manifest": str(OLD / "manifest.json"),
            "historical_manifest_sha256": digest(OLD / "manifest.json"),
            "historical_snapshots_exact_sha256_match": True,
            "current_source_sha256": hashes,
            "cases": len(manifest["cases"]),
            "geometry_rows": len(rows),
            "wall_seconds": perf_counter() - start,
            "closeness": {"relative_tolerance": 1e-10, "absolute_tolerance": 1e-12},
        },
    )


def leading_mode_details() -> None:
    """Record why the changed geometries differ without altering solver results."""
    cases = {c["name"]: c for c in read_json(OLD / "manifest.json")["cases"]}
    geometry = pd.read_csv(OUT / "geometry_comparison.csv")
    records = []
    for row in geometry.loc[~geometry.distances_close].itertuples():
        data = prepare_case_inputs(cases[row.case], selected_methods=[]).data
        # Recorded name for raw weighting; the current API calls it "raw".
        space = adaptive_cosine.separate_adaptive_cosine_space(
            data,
            weighting="binary",
            max_rank=40,
            min_segment_length=4,
            max_segments=8,
        )
        block = next(b for b in space.blocks if b.block_name == row.block)
        captured = {}
        original = diffusion.eigh

        def capture(matrix):
            eigenvalues, eigenvectors = original(matrix)
            captured["spectrum"] = eigenvalues[::-1]
            captured["components"] = int(
                connected_components(
                    matrix > 0,
                    directed=False,
                    return_labels=False,
                )
            )
            return eigenvalues, eigenvectors

        diffusion.eigh = capture
        try:
            diffusion.block_adaptive_diffusion_geometry(
                space.coordinates(block),
                k_neighbors=15,
                diffusion_time=3,
                n_components=30,
                metric="euclidean",
                bandwidth_type="-1/(d+2)",
                epsilon="median",
            )
        finally:
            diffusion.eigh = original
        values = captured["spectrum"]
        leading = extend_rank_through_ties(values, 1)
        start = 0 if tied_eigenvalue_neighbors(values)[0] else 1
        end = extend_rank_through_ties(values, start + 30)
        records.append(
            {
                "case": row.case,
                "block": row.block,
                "leading_numerical_group_size": leading,
                "leading_eigenvalues": values[: max(leading, 5)].tolist(),
                "positive_weight_connected_components": captured["components"],
                "old_one_based_modes": [2, 31],
                "new_one_based_modes": [start + 1, end],
            }
        )
    (OUT / "changed_geometry_leading_modes.json").write_text(json.dumps(records, indent=2) + "\n")


def compare_runs() -> None:
    old_m, new_m = read_json(OLD / "manifest.json"), read_json(OUT / "manifest.json")
    assert old_m["complete"] and new_m["complete"]
    assert old_m["cases"] == new_m["cases"]
    for k in ["edge_alpha", "sibling_alpha", "branch_length_variance_policy", "packages"]:
        assert old_m[k] == new_m[k], k
    for k in old_m["arguments"].keys() - {"output_dir"}:
        assert old_m["arguments"][k] == new_m["arguments"][k], k
    old, new = {}, {}
    for name in ["cases", "spectra", "blocks", "geometries", "runs", "trees", "assignments"]:
        old[name], new[name] = [pd.read_csv(folder / (name + ".csv")) for folder in (OLD, OUT)]
    for name in ["cases", "spectra", "blocks"]:
        pd.testing.assert_frame_equal(old[name], new[name], check_exact=True)
    for side in (old, new):
        assert not side["runs"].duplicated(KEY).any()
        assert not side["trees"].duplicated(KEY).any()
        assert not side["assignments"].duplicated(KEY + ["sample"]).any()
    merged = old["runs"].merge(
        new["runs"],
        on=KEY,
        suffixes=("_before", "_after"),
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    assert merged._merge.eq("both").all()
    merged = merged.drop(columns="_merge")
    merged["status_changed"] = merged.status_before.ne(merged.status_after)
    merged["jointly_ok"] = merged.status_before.eq("ok") & merged.status_after.eq("ok")
    merged["both_single_cluster"] = merged.n_clusters_before.eq(1) & merged.n_clusters_after.eq(1)
    merged["cluster_count_change"] = merged.n_clusters_after - merged.n_clusters_before
    merged["ari_truth_change"] = merged.ari_truth_after - merged.ari_truth_before
    geometry = pd.read_csv(OUT / "geometry_comparison.csv")
    geometry["distance_geometry_changed"] = ~geometry.distances_close
    merged = merged.merge(
        geometry[GEOMETRY_KEY + ["distance_geometry_changed", "max_absolute_distance_change"]],
        on=GEOMETRY_KEY,
        how="left",
        validate="many_to_one",
    )
    assert merged.distance_geometry_changed.notna().all()
    tree_maps, label_maps = [], []
    for side in (old, new):
        tree_maps.append(
            {
                tuple(row[k] for k in KEY): {frozenset(c) for c in json.loads(row.clades)}
                for _, row in side["trees"].iterrows()
            }
        )
        label_maps.append(
            {
                key: group.sort_values("sample").set_index("sample")
                for key, group in side["assignments"].groupby(KEY, sort=False)
            }
        )
        run_map = side["runs"].set_index(KEY)
        case_map = side["cases"].set_index("case")
        for key, labels in label_maps[-1].items():
            row = run_map.loc[key]
            assert row.status == "ok"
            assert len(labels) == case_map.loc[key[0], "n_samples"]
            assert labels.cluster.nunique() == row.n_clusters
            np.testing.assert_allclose(
                adjusted_rand_score(labels.truth, labels.cluster),
                row.ari_truth,
                rtol=1e-12,
                atol=1e-12,
            )
    for i, row in merged.iterrows():
        key = tuple(row[k] for k in KEY)
        trees = [mapping.get(key) for mapping in tree_maps]
        merged.loc[i, "tree_available_before"] = trees[0] is not None
        merged.loc[i, "tree_available_after"] = trees[1] is not None
        if all(t is not None for t in trees):
            a, b = trees
            merged.loc[i, "tree_identical"] = a == b
            merged.loc[i, "normalized_rooted_rf_before_after"] = (
                len(a ^ b) / (len(a) + len(b)) if a or b else 0.0
            )
        if row.jointly_ok:
            a, b = [mapping[key] for mapping in label_maps]
            pd.testing.assert_index_equal(a.index, b.index)
            np.testing.assert_array_equal(a.truth, b.truth)
            agreement = adjusted_rand_score(a.cluster, b.cluster)
            merged.loc[i, "ari_between_partitions"] = agreement
            merged.loc[i, "partition_identical"] = abs(agreement - 1.0) < 1e-12
    merged.to_csv(OUT / "run_comparison.csv", index=False)
    merged.loc[
        merged.status_changed | (merged.jointly_ok & merged.partition_identical.eq(False))
    ].to_csv(
        OUT / "changed_runs.csv",
        index=False,
    )
    merged.groupby(["status_before", "status_after"]).size().rename("count").reset_index().to_csv(
        OUT / "status_transitions.csv",
        index=False,
    )
    summary_rows = []
    for (case, geometry), group in merged.groupby(["case", "geometry"], sort=False):
        joint = group[group.jointly_ok]
        informative = joint[~joint.both_single_cluster]
        summary_rows.append(
            {
                "case": case,
                "geometry": geometry,
                "attempts": len(group),
                "ok_before": int(group.status_before.eq("ok").sum()),
                "ok_after": int(group.status_after.eq("ok").sum()),
                "jointly_ok": len(joint),
                "status_changes": int(group.status_changed.sum()),
                "changed_trees": int(group.tree_identical.eq(False).sum()),
                "changed_partitions": int(joint.partition_identical.eq(False).sum()),
                "single_cluster_before": int(group.n_clusters_before.eq(1).sum()),
                "single_cluster_after": int(group.n_clusters_after.eq(1).sum()),
                "jointly_ok_both_single_cluster": int(joint.both_single_cluster.sum()),
                "jointly_ok_at_least_one_multicluster": len(informative),
                "ari_truth_improved": int(joint.ari_truth_change.gt(1e-12).sum()),
                "ari_truth_worsened": int(joint.ari_truth_change.lt(-1e-12).sum()),
                "cluster_count_changed": int(joint.cluster_count_change.ne(0).sum()),
            }
        )
    pd.DataFrame(summary_rows).to_csv(OUT / "summary_by_case_geometry.csv", index=False)
    sources = []
    for path in sorted(old_m["source_sha256"].keys() | new_m["source_sha256"].keys()):
        before, after = old_m["source_sha256"].get(path), new_m["source_sha256"].get(path)
        current = digest(ROOT / path) if (ROOT / path).is_file() else None
        sources.append(
            {
                "path": path,
                "before_sha256": before,
                "after_sha256": after,
                "source_changed": before != after,
                "unchanged_during_run": after == current,
            }
        )
    source_frame = pd.DataFrame(sources)
    assert source_frame.unchanged_during_run.all()
    source_frame.to_csv(OUT / "source_comparison.csv", index=False)
    geometries = pd.read_csv(OUT / "geometry_comparison.csv")
    jointly_ok = merged[merged.jointly_ok]
    result = {
        "comparison_script_sha256": digest(Path(__file__)),
        "baseline_manifest_sha256": digest(OLD / "manifest.json"),
        "rerun_manifest_sha256": digest(OUT / "manifest.json"),
        "attempts": len(merged),
        "cases": len(new["cases"]),
        "blocks": len(new["blocks"]),
        "historical_status_counts": old["runs"].status.value_counts().to_dict(),
        "new_status_counts": new["runs"].status.value_counts().to_dict(),
        "status_changes": int(merged.status_changed.sum()),
        "jointly_ok": len(jointly_ok),
        "jointly_available_trees": int(merged.tree_identical.notna().sum()),
        "changed_trees": int(merged.tree_identical.eq(False).sum()),
        "changed_partitions_jointly_ok": int(jointly_ok.partition_identical.eq(False).sum()),
        "changed_partitions_with_unchanged_tree": int(
            (jointly_ok.partition_identical.eq(False) & jointly_ok.tree_identical.eq(True)).sum()
        ),
        "changed_partitions_with_changed_tree": int(
            (jointly_ok.partition_identical.eq(False) & jointly_ok.tree_identical.eq(False)).sum()
        ),
        "single_cluster_before": int(old["runs"].n_clusters.eq(1).sum()),
        "single_cluster_after": int(new["runs"].n_clusters.eq(1).sum()),
        "both_single_cluster_jointly_ok": int(jointly_ok.both_single_cluster.sum()),
        "at_least_one_multicluster_jointly_ok": int((~jointly_ok.both_single_cluster).sum()),
        "ari_truth_improved_jointly_ok": int(jointly_ok.ari_truth_change.gt(1e-12).sum()),
        "ari_truth_worsened_jointly_ok": int(jointly_ok.ari_truth_change.lt(-1e-12).sum()),
        "cluster_count_changed_jointly_ok": int(jointly_ok.cluster_count_change.ne(0).sum()),
        "geometries_compared": len(geometries),
        "geometries_exact": int(geometries.distances_exact.eq(True).sum()),
        "geometries_close": int(geometries.distances_close.eq(True).sum()),
        "diffusion_dimension_changes": int(
            geometries.before_dimension.ne(geometries.after_dimension).sum()
        ),
        "input_spectra_blocks_exact_match": True,
        "cases_seeds_flags_packages_match": True,
        "all_recorded_sources_stable_during_run": True,
        "source_changes": source_frame.loc[source_frame.source_changed, "path"].tolist(),
        "invocation": read_json(OUT / "invocation.json"),
    }
    write_json(OUT / "comparison_summary.json", result)
    render_report(result, merged, pd.DataFrame(summary_rows), source_frame)
    print(json.dumps(result, indent=2), flush=True)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in frame[columns].iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in columns) + " |")
    return "\n".join(lines)


def render_report(
    result: dict, runs: pd.DataFrame, cases: pd.DataFrame, sources: pd.DataFrame
) -> None:
    single_population = runs[runs.case.eq("mp_spike_below_bbp_continuous")]
    before_null_ok = single_population.status_before.eq("ok")
    after_null_ok = single_population.status_after.eq("ok")
    null_before = int((before_null_ok & single_population.n_clusters_before.gt(1)).sum())
    null_after = int((after_null_ok & single_population.n_clusters_after.gt(1)).sum())
    transitions = pd.read_csv(OUT / "status_transitions.csv")
    changed_sources = sources.loc[sources.source_changed, "path"].tolist()
    report = f"""# Matched spectral-block reassessment — 2026-09-28

The complete same 11-case `method_proof` panel was rerun after the numerical
eigenvalue-group rank repair: 47 blocks, eight tree inference methods, both raw
and adaptive-diffusion geometries, and 752 attempts. The historical results in
`reports/eigenvalue_block_comparison_20260926/` were preserved.

## Matched outcomes

- Status counts before: `{result["historical_status_counts"]}`; after:
  `{result["new_status_counts"]}`. There are {result["status_changes"]} status changes.
- {result["jointly_ok"]} runs returned partitions in both versions;
  {result["changed_partitions_jointly_ok"]} of these partitions changed, including
  {result["cluster_count_changed_jointly_ok"]} changes in cluster count.
- {result["jointly_available_trees"]} trees are available in both versions;
  {result["changed_trees"]} rooted clade sets changed.
- {result["changed_partitions_with_unchanged_tree"]} changed partitions occur on
  unchanged rooted trees; {result["changed_partitions_with_changed_tree"]} occur
  on changed trees. The former changes arise downstream of tree inference.
- Single-cluster outputs: {result["single_cluster_before"]} before and
  {result["single_cluster_after"]} after. Among jointly completed comparisons,
  {result["both_single_cluster_jointly_ok"]} are single-cluster outputs in both
  versions; {result["at_least_one_multicluster_jointly_ok"]} have at least one
  multi-cluster output. Partition agreement in the former group is not evidence
  of useful cluster recovery.
- Against generator labels, ARI increased in
  {result["ari_truth_improved_jointly_ok"]} and decreased in
  {result["ari_truth_worsened_jointly_ok"]} jointly completed runs. These dependent
  block/method counts are descriptive; no pooled quality score is reported.

The one-population `mp_spike_below_bbp_continuous` example had
{null_before}/{int(before_null_ok.sum())} completed outputs with multiple clusters
before and {null_after}/{int(after_null_ok.sum())} after. This is one generated
dataset evaluated repeatedly, not an estimated null error rate or proof of
calibration. The phylogenetic Brownian example's generator labels have 16 groups
and are not recoded as a one-population partition truth.

### Explicit outcome transitions

{markdown_table(transitions, ["status_before", "status_after", "count"])}

### Per-case and geometry counts

{markdown_table(cases, ["case", "geometry", "ok_before", "ok_after", "changed_trees", "changed_partitions", "single_cluster_before", "single_cluster_after"])}

## Geometry and provenance

Input case definitions, seeds, data hashes, retained spectra and all block records
match the historical output exactly. Package versions and all benchmark flags
also match. Direct replay compared all 94 geometries, using historical cosine and
diffusion source snapshots whose SHA-256 hashes exactly match the historical
manifest. All 47 raw geometries are exact; 44 of 47 adaptive-diffusion geometries
are exact. The three changed adaptive geometries are recorded in
`geometry_comparison.csv`. No component-count changes occurred in this panel.

The three changed geometries retain complete repeated leading groups: four
connected components for `support_selected_only_fail_closed` modes 2–6, two for
`phylo_brownian_null_16taxa` modes 2–5, and six for
`traversal_deep_branch_recovery_stress` modes 2–6. Their old diffusion selection
used modes 2–31; the repaired selection uses modes 1–30. The complete repeated
leading group is retained. See `changed_geometry_leading_modes.json`.

The source comparison records {len(sources)} paths; all remain unchanged during
the inference run. Differences from the historical source manifest are:

{chr(10).join("- `" + path + "`" for path in changed_sources)}

The historical diffusion snapshot was recovered from its recorded Git commit.
The historical cosine module contained an uncommitted tie fix: the comparison
helper reverses the later shared-helper extraction, then requires an exact
historical SHA-256 match before executing that reconstructed snapshot. It never
replaces working-tree source files.

## Reproduction and verification

The actual invocation, environment overrides, start/end times and exit code are
in `invocation.json`; all runner output is retained in `run.log`. Inference wall
time was {result["invocation"]["wall_seconds"]:.1f} seconds. Other concurrent work
means this is not a controlled speed comparison. The existing benchmark command
creates its output directory and refuses to overwrite an earlier run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg uv run --no-sync python -m benchmarks.diagnostics.spectral.adaptive_cosine.compare_spectral_blocks --suite method_proof --max-rank 40 --methods average complete weighted single centroid median ward neighbor_joining --geometries raw adaptive_diffusion --output-dir <new-output-directory>
```

Recompute these comparisons against the preserved baseline with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg .venv/bin/python reports/eigenspace_reassessment_20260928/spectral_blocks/compare_before_after.py
```

The helper first requires the three source hashes recorded in
`geometry_replay.json`. That working tree (commit `f34f59ca` plus uncommitted
edits) matches no commit, so on later source the helper stops at this check
instead of recomputing. It keeps the recorded weighting name `binary`, which
the current API calls `raw`.

The helper checks identical settings, input/spectrum/block records, unique
matched run keys, successful assignment coverage, truth labels, cluster counts,
recomputed ARI, paired tree availability and source stability. Failed and
unsupported runs retain their statuses and are excluded from partition metrics;
they are never turned into single-cluster assignments. The shared CLI calls
diffusion failures `failed_geometry`; none occurred in either version of this
panel. The local helper was linted successfully.

## Scope

The target rank remains 40, cosine weighting is raw values, minimum segment
length four and maximum segments eight. Diffusion uses K=15, time three,
30 target components, Euclidean metric, bandwidth `-1/(d+2)` and median epsilon.
The tree methods are average, complete, weighted, single, centroid, median,
Ward and neighbor joining with MAD rooting. Edge alpha is 0.001 and sibling
alpha is 0.01, with support-threshold enforcement false, branch-time policy
`none`, and contrast-energy adaptive Wald rank disabled in both runs.

This is a matched behavioral comparison, not held-out calibration evidence.
Failures, unsupported cases, one-cluster agreement and generator-label ARI must
remain separate. No inference about a method's general superiority or error
control follows from these dependent one-seed-per-case comparisons. IQ-TREE,
NNLS/branch-time variants and the full repository benchmark suite are outside
this panel.

Detailed changes: `run_comparison.csv`, `changed_runs.csv`,
`summary_by_case_geometry.csv`, `source_comparison.csv` and
`comparison_summary.json`. Original-format outputs (`runs.csv`, `trees.csv`,
`assignments.csv`, spectra, blocks and pair tables) remain alongside them.
"""
    (OUT / "README.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-only", action="store_true")
    args = parser.parse_args()
    geometry_replay()
    leading_mode_details()
    if not args.geometry_only:
        compare_runs()


if __name__ == "__main__":
    main()
