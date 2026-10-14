#!/usr/bin/env python3
"""Decode saved per-contrast outcomes and check every aggregate table."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

directory = Path(__file__).resolve().parent
counts = pd.read_csv(directory / "counts.csv")
paired = pd.read_csv(directory / "paired_decisions.csv")
ranks = pd.read_csv(directory / "rank_counts.csv")
manifest = json.loads((directory / "manifest.json").read_text())
verified_counts = verified_pairs = verified_rank_counts = 0

for seed in sorted(counts.seed.unique()):
    archive = directory / f"replicate_outcomes_seed_{seed}.npz"
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest[
        "outcome_archive_sha256"
    ][archive.name]
    with np.load(archive) as arrays:
        for row in counts[counts.seed == seed].itertuples():
            prefix = f"{row.spectrum}__{row.selection}__{row.version}"
            rejected = np.unpackbits(
                arrays[prefix + "__alpha_" + str(row.alpha)], bitorder="little",
            )[:row.n_contrasts]
            assert int(rejected.sum()) == row.rejections
            assert arrays[prefix + "__rank"].shape == (row.n_contrasts,)
            verified_counts += 1
        for row in paired[paired.seed == seed].itertuples():
            prefix = f"{row.spectrum}__{row.selection}__"
            old, new = [
                np.unpackbits(
                    arrays[prefix + version + "__alpha_" + str(row.alpha)],
                    bitorder="little",
                )[:row.n_contrasts].astype(bool)
                for version in ("legacy_prefix", "current_complete_groups")
            ]
            assert (
                int((old & ~new).sum()), int((new & ~old).sum()),
                int((old & new).sum()), int((~old & ~new).sum()),
            ) == (
                row.old_only_rejected, row.new_only_rejected,
                row.both_rejected, row.neither_rejected,
            )
            assert row.rank_changed == int(np.sum(
                arrays[prefix + "legacy_prefix__rank"]
                != arrays[prefix + "current_complete_groups__rank"],
            ))
            verified_pairs += 1
        for row in ranks[ranks.seed == seed].itertuples():
            key = f"{row.spectrum}__{row.selection}__{row.version}__rank"
            assert int(np.sum(arrays[key] == row.rank)) == row.count
            verified_rank_counts += 1

kernel_checks = pd.read_csv(directory / "kernel_checks.csv")
assert bool(kernel_checks.matches.all())
rotation_checks = pd.read_csv(directory / "rotation_audit.csv")
current = rotation_checks[rotation_checks.version == "current_complete_groups"]
assert bool((current.rank_mismatches == 0).all())
for column in current.columns:
    if column.startswith("decision_mismatches_"):
        assert bool((current[column] == 0).all())
assert manifest["production_sources_unchanged_during_run"]

result = {
    "status": "passed",
    "verified_rejection_count_rows": verified_counts,
    "verified_paired_decision_rows": verified_pairs,
    "verified_rank_count_rows": verified_rank_counts,
    "kernel_checks_passed": len(kernel_checks),
    "current_rotation_rank_and_decision_mismatches": 0,
    "current_max_rotation_p_value_difference": float(current.max_abs_p_difference.max()),
    "outcome_archive_sha256_verified": True,
}
(directory / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
