"""Verify and compact the final-source restoration benchmark and matched history."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from reports.diffusion_nnls_versions_20260909.summarize import parse_params
from sklearn.metrics import adjusted_rand_score


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    out = Path(__file__).resolve().parent
    run = root / "benchmarks/results/run_20260909_calibration_restoration_verified"
    prior = root / "benchmarks/results/run_20260909_diffusion_nnls_versions_full"
    historical = root / "benchmarks/results/run_20260801_130648Z_full_post_crossfit"
    method = "tbs_diffusion_adaptive_nnls"
    manifest = json.loads((run / "run_manifest.json").read_text())
    for name, expected in manifest["source_sha256"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    status = json.loads((run / "continuation_status.json").read_text())
    assert status["exit_code"] == 0 and status["execution_errors"] == 0
    restored = pd.read_csv(run / "full_benchmark_comparison.csv")
    assert len(restored) == 122 and restored.case_id.is_unique
    assert set(restored.method) == {method}

    def read_method(directory: Path) -> pd.DataFrame:
        frame = pd.read_csv(directory / "full_benchmark_comparison.csv")
        return frame.loc[frame.method.eq(method)].copy()

    before, old = read_method(prior), read_method(historical)
    assert before.case_id.is_unique and old.case_id.is_unique
    verified_labels = 0
    for identity in sorted((run / "cell_audits").glob("*/identity.json")):
        key = json.loads(identity.read_text())
        labels = pd.read_csv(identity.parent / "labels.csv")
        row = restored.loc[restored.case_id.eq(key["case_id"])].iloc[0]
        assert row.run_id == key["run_id"] and row.status == "ok"
        assert len(labels) == row.samples
        assert labels.cluster.nunique() == row.found_clusters
        assert np.isclose(adjusted_rand_score(labels.truth, labels.cluster), row.ari)
        verified_labels += 1
    assert verified_labels == int(restored.status.eq("ok").sum())

    def signatures(directory: Path) -> dict[str, set[tuple[str, str]]]:
        result: dict[str, set[tuple[str, str]]] = {}
        for line in (directory / "input_signatures.jsonl").read_text().splitlines():
            item = json.loads(line)
            result.setdefault(item["case_id"], set()).add(
                (item["data_sha256"], item["truth_sha256"])
            )
        return result

    final_signatures, prior_signatures = signatures(run), signatures(prior)
    assert len(final_signatures) == 122
    matched_signatures = 0
    for case_id, hashes in final_signatures.items():
        assert len(hashes) == 1, case_id
        if case_id in prior_signatures:
            assert hashes == prior_signatures[case_id], case_id
            matched_signatures += 1

    summaries = {}
    comparisons = []
    columns = [
        "case_id",
        "status",
        "ari",
        "found_clusters",
        "params",
        "samples",
        "features",
        "true_clusters",
    ]
    for name, frame in [("august_1", old), ("before_restoration", before), ("restored", restored)]:
        ok = frame[frame.status.eq("ok")]
        summaries[name] = {
            "rows": len(frame),
            "statuses": frame.status.value_counts().to_dict(),
            "mean_ari_on_ok": float(ok.ari.mean()),
            "exact_k_on_ok": int(ok.found_clusters.eq(ok.true_clusters).sum()),
        }
        if name == "restored":
            continue
        matched = frame[columns].merge(
            restored[columns], on="case_id", suffixes=("_old", "_restored"), validate="one_to_one"
        )
        matched["source"] = name
        matched["identical_recorded_params"] = [
            parse_params(a) == parse_params(b)
            for a, b in zip(matched.params_old, matched.params_restored, strict=True)
        ]
        matched["matching_input_metadata"] = np.logical_and.reduce(
            [
                matched[f"{field}_old"].eq(matched[f"{field}_restored"])
                for field in ["samples", "features", "true_clusters"]
            ]
        )
        both_ok = matched.status_old.eq("ok") & matched.status_restored.eq("ok")
        matched["ari_delta"] = matched.ari_restored - matched.ari_old
        summaries[name]["matched"] = {
            "cases": len(matched),
            "both_ok": int(both_ok.sum()),
            "identical_params": int(matched.identical_recorded_params.sum()),
            "matching_input_metadata": int(matched.matching_input_metadata.sum()),
            "mean_old_ari_on_both_ok": float(matched.loc[both_ok, "ari_old"].mean()),
            "mean_restored_ari_on_both_ok": float(matched.loc[both_ok, "ari_restored"].mean()),
            "same_ari_and_k_on_both_ok": int(
                (
                    both_ok
                    & np.isclose(matched.ari_old, matched.ari_restored)
                    & matched.found_clusters_old.eq(matched.found_clusters_restored)
                ).sum()
            ),
            "status_transitions": matched.groupby(["status_old", "status_restored"])
            .size()
            .reset_index(name="count")
            .to_dict(orient="records"),
        }
        comparisons.append(matched)
    pd.concat(comparisons, ignore_index=True).to_csv(out / "matched_cases.csv", index=False)
    restored.loc[restored.case_id.str.contains("null")].to_csv(out / "null_cases.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    (out / "verification.json").write_text(
        json.dumps(
            {
                "completed_outcomes": len(restored),
                "independently_rescored_labels": verified_labels,
                "unchanged_source_hashes": len(manifest["source_sha256"]),
                "matched_prior_input_and_truth_hashes": matched_signatures,
                "missing_prior_input_hashes": sorted(set(final_signatures) - set(prior_signatures)),
                "historical_input_hashes_available": False,
            },
            indent=2,
        )
        + "\n"
    )
    for name in [
        "run_manifest.json",
        "configuration_summary.csv",
        "continuation_status.json",
        "full_benchmark_comparison.csv",
        "input_signatures.jsonl",
    ]:
        shutil.copyfile(run / name, out / name)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
