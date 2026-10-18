"""Finish the recorded benchmark through three isolated case shards."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    manifest = json.loads((run / "run_manifest.json").read_text())
    original = pd.read_csv(run / "full_benchmark_comparison.csv")
    errors = [
        json.loads(line) for line in (run / "execution_errors.jsonl").read_text().splitlines()
    ]
    shards, commands, processes = [], [], []
    started = time.monotonic()
    for index in range(3):
        shard = run / "parallel_shards" / str(index + 1)
        shard.mkdir(parents=True, exist_ok=False)
        indices = list(range(index + 1, manifest["case_count"] + 1, 3))
        cases = [manifest["cases"][number - 1] for number in indices]
        names = {case["name"] for case in cases}
        shard_manifest = dict(manifest)
        shard_manifest.update(
            {
                "case_count": len(cases),
                "cases": cases,
                "expected_rows": len(cases) * manifest["configuration_count"],
                "parent_run": str(run),
                "case_indices": indices,
            }
        )
        (shard / "run_manifest.json").write_text(json.dumps(shard_manifest, indent=2) + "\n")
        original[original.case_id.isin(names)].to_csv(
            shard / "full_benchmark_comparison.csv", index=False
        )
        (shard / "execution_errors.jsonl").write_text(
            "".join(json.dumps(error) + "\n" for error in errors if error["case_id"] in names)
        )
        command = [
            sys.executable,
            "-u",
            "-m",
            "reports.diffusion_nnls_versions_20260909.run_panel",
            "--run-dir",
            str(shard),
            "--case-indices",
            *map(str, indices),
        ]
        commands.append(command)
        shards.append(shard)
        log = (shard / "execution.log").open("w")
        processes.append(
            (
                subprocess.Popen(
                    command, cwd=manifest["cwd"], stdout=log, stderr=subprocess.STDOUT
                ),
                log,
            )
        )
    parallel_manifest = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "commands": commands,
        "case_partition": "Original one-based case index modulo three; all configurations for a case stay together.",
        "initial_completed_outcomes": len(original) + len(errors),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "worker_script_sha256": hashlib.sha256(
            Path(__file__).with_name("run_panel.py").read_bytes()
        ).hexdigest(),
        "worker_count": 3,
        "numerical_threads_per_worker": 1,
    }
    (run / "parallel_manifest.json").write_text(json.dumps(parallel_manifest, indent=2) + "\n")
    exits = []
    for process, log in processes:
        exits.append(process.wait())
        log.close()
    (run / "parallel_exits.json").write_text(json.dumps(exits) + "\n")
    assert exits == [0, 0, 0], exits
    combined = pd.concat(
        [pd.read_csv(shard / "full_benchmark_comparison.csv") for shard in shards],
        ignore_index=True,
    )
    all_errors = [
        json.loads(line)
        for shard in shards
        for line in (shard / "execution_errors.jsonl").read_text().splitlines()
    ]
    assert not combined.duplicated(["case_id", "run_id"]).any()
    assert len(combined) + len(all_errors) == manifest["expected_rows"]
    done = set(zip(original.case_id, original.run_id))
    added = combined[[key not in done for key in zip(combined.case_id, combined.run_id)]]
    assert set(zip(combined.case_id, combined.run_id)).issuperset(done)
    added.to_csv(run / "full_benchmark_comparison.csv", mode="a", header=False, index=False)
    error_done = {(row["case_id"], row["run_id"]) for row in errors}
    with (run / "execution_errors.jsonl").open("a") as stream:
        for row in all_errors:
            if (row["case_id"], row["run_id"]) not in error_done:
                stream.write(json.dumps(row) + "\n")
    for shard in shards:
        if (shard / "cell_audits").exists():
            for directory in (shard / "cell_audits").iterdir():
                shutil.copytree(directory, run / "cell_audits" / directory.name)
    # The existing worker's no-pending path verifies totals and writes the final summary.
    final_command = [
        sys.executable,
        "-u",
        "-m",
        "reports.diffusion_nnls_versions_20260909.run_panel",
        "--run-dir",
        str(run),
    ]
    subprocess.run(final_command, cwd=manifest["cwd"], check=True)
    status = json.loads((run / "continuation_status.json").read_text())
    status.update(
        {
            "parallel_elapsed_seconds": time.monotonic() - started,
            "parallel_worker_exit_codes": exits,
            "serial_stage_exit_code": 130,
            "summary_only_elapsed_seconds": status.pop("continuation_elapsed_seconds"),
        }
    )
    (run / "continuation_status.json").write_text(json.dumps(status, indent=2) + "\n")
    print(json.dumps(status), flush=True)


if __name__ == "__main__":
    main()
