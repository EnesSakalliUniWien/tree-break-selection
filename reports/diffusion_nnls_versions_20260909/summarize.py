"""Preserve compact results and matched historical comparisons for this study."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score


def parse_params(value: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in re.split(r", (?=[A-Za-z_]+[A-Za-z_0-9]*=)", value))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    out = Path(__file__).resolve().parent
    status = json.loads((run / "continuation_status.json").read_text())
    assert status["all_configurations_attempted"] and status["exit_code"] == 0
    assert json.loads((run / "support_review/status.json").read_text())["exit_code"] == 0
    df = pd.read_csv(run / "full_benchmark_comparison.csv")
    summary = pd.read_csv(run / "configuration_summary.csv")
    summary.insert(0, "configuration", [f"C{i:02}" for i in range(1, len(summary) + 1)])
    ids = dict(zip(summary.run_id, summary.configuration))
    configs = dict(zip(summary.run_id, summary.method))
    errors = [
        json.loads(line) for line in (run / "execution_errors.jsonl").read_text().splitlines()
    ]
    columns = [
        "case_id",
        "method",
        "run_id",
        "source_family",
        "feature_representation",
        "samples",
        "features",
        "true_clusters",
        "found_clusters",
        "status",
        "ari",
        "nmi",
        "unsupported_reason_code",
        "unsupported_admissible_support_count",
        "skip_reason",
    ]
    outcomes = df[columns].copy()
    outcomes.loc[outcomes.status.ne("ok"), "found_clusters"] = float("nan")
    error_rows = pd.DataFrame(
        {
            "case_id": row["case_id"],
            "method": row["method"],
            "run_id": row["run_id"],
            "status": "execution_error",
            "error": row["error"],
        }
        for row in errors
    )
    outcomes = pd.concat([outcomes, error_rows], ignore_index=True)
    case_metadata = df.drop_duplicates("case_id").set_index("case_id")
    for column in [
        "source_family",
        "feature_representation",
        "samples",
        "features",
        "true_clusters",
    ]:
        outcomes[column] = outcomes[column].fillna(outcomes.case_id.map(case_metadata[column]))
    outcomes.insert(0, "configuration", outcomes.run_id.map(ids))
    assert not outcomes.duplicated(["configuration", "case_id"]).any()
    assert len(outcomes) == 1220
    assert outcomes.groupby("configuration").size().eq(122).all()
    manifest = json.loads((run / "run_manifest.json").read_text())
    expected_case_ids = {case["name"] for case in manifest["cases"]}
    for _, group in outcomes.groupby("configuration"):
        assert set(group.case_id) == expected_case_ids
    assert df.loc[df.status.ne("ok"), "ari"].isna().all()
    outcomes.drop(columns=["run_id", "method"]).sort_values(["configuration", "case_id"]).to_csv(
        out / "outcomes.csv", index=False
    )
    summary.to_csv(out / "configuration_summary.csv", index=False)
    representation_summary = []
    for (configuration, representation), group in outcomes.groupby(
        ["configuration", "feature_representation"]
    ):
        ok = group[group.status.eq("ok")]
        representation_summary.append(
            {
                "configuration": configuration,
                "feature_representation": representation,
                "requested": len(group),
                "ok": len(ok),
                "unsupported": int(group.status.eq("unsupported").sum()),
                "skip": int(group.status.eq("skip").sum()),
                "execution_error": int(group.status.eq("execution_error").sum()),
                "mean_ari_on_ok": ok.ari.mean(),
            }
        )
    pd.DataFrame(representation_summary).to_csv(out / "representation_summary.csv", index=False)
    reason_counts = (
        df[df.status.eq("unsupported")]
        .groupby(["run_id", "unsupported_reason_code"])
        .size()
        .rename("count")
        .reset_index()
    )
    reason_counts.insert(0, "configuration", reason_counts.run_id.map(ids))
    reason_counts.drop(columns="run_id").to_csv(out / "unsupported_counts.csv", index=False)
    error_counts = (
        pd.DataFrame(errors).groupby(["run_id", "error"]).size().rename("count").reset_index()
    )
    error_counts.insert(0, "configuration", error_counts.run_id.map(ids))
    error_counts.drop(columns="run_id").to_csv(out / "execution_error_counts.csv", index=False)

    historical = [
        ("run_20260801_130648Z_full_post_crossfit", "tbs_diffusion_adaptive_nnls"),
        (
            "run_20260706_170815Z_full_graphtools_tree_consensus",
            "tbs_diffusion_graphtools_adaptive_nnls",
        ),
        ("run_20260630_130921Z_full", "tbs_diffusion_graphtools_nnls"),
        ("run_20260811_selected_law_variants_full", "tbs_diffusion_adaptive_nnls"),
    ]
    comparisons, history_summary, source_hashes = [], [], {}
    for directory, method in historical:
        source = run.parent / directory / "full_benchmark_comparison.csv"
        source_hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        old = pd.read_csv(source)
        old = old[old.method.eq(method)].copy()
        old["parsed"] = old.params.map(parse_params)
        for run_id, current_method in configs.items():
            if current_method != method:
                continue
            current = df[df.run_id.eq(run_id)]
            config = summary[summary.run_id.eq(run_id)].iloc[0]
            matching = old[
                old.parsed.map(
                    lambda params: (
                        params.get("tree_builder", "linkage") == config.tree_builder
                        and params.get("tree_linkage_method", "average")
                        == config.tree_linkage_method
                    )
                )
            ]
            assert not matching.case_id.duplicated().any()
            ok = matching[matching.status.eq("ok")]
            history_summary.append(
                {
                    "source": directory,
                    "configuration": ids[run_id],
                    "method": method,
                    "tree_builder": config.tree_builder,
                    "tree_linkage_method": config.tree_linkage_method,
                    "rows": len(matching),
                    "ok": len(ok),
                    "unsupported": int(matching.status.eq("unsupported").sum()),
                    "skip": int(matching.status.eq("skip").sum()),
                    "mean_ari_on_ok": ok.ari.mean(),
                }
            )
            joined = matching.merge(
                current, on="case_id", suffixes=("_old", "_current"), validate="one_to_one"
            )
            for row in joined.to_dict("records"):
                old_params, new_params = row["parsed"], parse_params(row["params_current"])
                differences = {
                    key: [old_params.get(key), new_params.get(key)]
                    for key in sorted(old_params.keys() | new_params.keys())
                    if old_params.get(key) != new_params.get(key)
                }
                comparisons.append(
                    {
                        "source": directory,
                        "configuration": ids[run_id],
                        "case_id": row["case_id"],
                        "old_status": row["status_old"],
                        "current_status": row["status_current"],
                        "old_ari": row["ari_old"],
                        "current_ari": row["ari_current"],
                        "old_found_clusters": row["found_clusters_old"],
                        "current_found_clusters": row["found_clusters_current"],
                        "identical_recorded_params": not differences,
                        "parameter_differences": json.dumps(differences, sort_keys=True),
                        "matching_samples_features_truth": all(
                            row[f"{col}_old"] == row[f"{col}_current"]
                            for col in ["samples", "features", "true_clusters"]
                        ),
                        "matching_source_representation": all(
                            row[f"{col}_old"] == row[f"{col}_current"]
                            for col in ["source_family", "feature_representation"]
                        ),
                    }
                )
    pd.DataFrame(comparisons).to_csv(out / "historical_matched_cases.csv", index=False)
    pd.DataFrame(history_summary).to_csv(out / "historical_configuration_summary.csv", index=False)
    old_panel = pd.read_csv(
        run.parent / "run_20260801_130648Z_full_post_crossfit/full_benchmark_comparison.csv"
    )
    legacy_pair = old_panel[old_panel.method.eq("tbs_diffusion_adaptive_nnls")].merge(
        old_panel[old_panel.method.eq("tbs")],
        on="case_id",
        suffixes=("_diffusion_nnls", "_tbs"),
        validate="one_to_one",
    )
    legacy_pair = legacy_pair[
        legacy_pair.status_diffusion_nnls.eq("ok") & legacy_pair.status_tbs.eq("ok")
    ].copy()
    for column in [
        "samples",
        "features",
        "true_clusters",
        "source_family",
        "feature_representation",
    ]:
        assert legacy_pair[f"{column}_diffusion_nnls"].eq(legacy_pair[f"{column}_tbs"]).all()
    for name in ["benchmark_seed", "edge_alpha", "sibling_alpha"]:
        assert (
            legacy_pair.params_diffusion_nnls.map(parse_params)
            .map(lambda p: p[name])
            .eq(legacy_pair.params_tbs.map(parse_params).map(lambda p: p[name]))
            .all()
        )
    legacy_pair["ari_change"] = legacy_pair.ari_diffusion_nnls - legacy_pair.ari_tbs
    legacy_pair[
        [
            "case_id",
            "ari_diffusion_nnls",
            "ari_tbs",
            "ari_change",
            "found_clusters_diffusion_nnls",
            "found_clusters_tbs",
        ]
    ].to_csv(out / "legacy_diffusion_vs_tbs.csv", index=False)
    for source, name in [
        ("run_manifest.json", "run_manifest.json"),
        ("execution_status.json", "initial_execution_status.json"),
        ("continuation_manifest.json", "continuation_manifest.json"),
        ("continuation_status.json", "continuation_status.json"),
        ("parallel_manifest.json", "parallel_manifest.json"),
        ("parallel_exits.json", "parallel_exits.json"),
        ("support_review/results.csv", "support_review_results.csv"),
        ("support_review/manifest.json", "support_review_manifest.json"),
        ("support_review/status.json", "support_review_status.json"),
        ("support_review/initial_attempt/status.json", "support_review_initial_status.json"),
        ("support_review/initial_attempt/manifest.json", "support_review_initial_manifest.json"),
    ]:
        shutil.copyfile(run / source, out / name)
        source_hashes[str(run / source)] = hashlib.sha256((run / source).read_bytes()).hexdigest()
    for source in [run / "full_benchmark_comparison.csv", run / "execution_errors.jsonl"]:
        source_hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    (out / "evidence_sha256.json").write_text(json.dumps(source_hashes, indent=2) + "\n")
    for name, expected in manifest["source_sha256"].items():
        assert (
            hashlib.sha256((Path(manifest["cwd"]) / name).read_bytes()).hexdigest() == expected
        ), name
    scored = 0
    for row in df[df.status.eq("ok")].to_dict("records"):
        directory = (
            run
            / "cell_audits"
            / f"{row['test_case']:03d}_{hashlib.sha256(row['run_id'].encode()).hexdigest()[:12]}"
        )
        identity = json.loads((directory / "identity.json").read_text())
        assert identity == {"case_id": row["case_id"], "run_id": row["run_id"]}
        labels = pd.read_csv(directory / "labels.csv")
        assert abs(adjusted_rand_score(labels.truth, labels.cluster) - row["ari"]) < 1e-14
        assert labels.cluster.nunique() == row["found_clusters"]
        scored += 1
    signatures = {}
    signature_files = [
        run / "input_signatures.jsonl",
        *run.glob("parallel_shards/*/input_signatures.jsonl"),
    ]
    for path in signature_files:
        for line in path.read_text().splitlines():
            value = json.loads(line)
            case_id = value["case_id"]
            if case_id in signatures:
                assert value == signatures[case_id]
            signatures[case_id] = value
    verification = {
        "status": "passed",
        "unique_outcomes": len(outcomes),
        "configurations": len(summary),
        "cases_per_configuration": 122,
        "all_expected_case_ids_present": True,
        "independently_rescored_saved_label_files": scored,
        "unchanged_source_files": len(manifest["source_sha256"]),
        "distinct_input_signatures": len(signatures),
        "repeated_input_signatures_match": True,
        "input_signature_limit": "The first case's ten rows predate the continuation and have no input signature. They are all unsupported and are covered by original invocation/source/seed records.",
        "ok_cluster_count_distribution": df[df.status.eq("ok")]
        .found_clusters.value_counts()
        .to_dict(),
    }
    (out / "result_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(summary.drop(columns=["run_id"]).to_string(index=False))
    print("Unsupported reasons:", df.unsupported_reason_code.value_counts().to_dict())
    print("Execution errors:", pd.DataFrame(errors).error.value_counts().to_dict())


if __name__ == "__main__":
    main()
