"""Shard layout, S3 sync and checks shared by the AWS batch benchmark drivers."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd


def validate_shard_contract(*, shard_index: int, shard_count: int) -> None:
    """Validate shard dimensions before launching compute."""
    if shard_count <= 0:
        raise ValueError(f"shard_count must be positive; got {shard_count!r}.")
    if not 0 <= shard_index < shard_count:
        raise ValueError(
            f"shard_index must satisfy 0 <= index < shard_count; got "
            f"{shard_index!r} with shard_count={shard_count!r}."
        )


def validate_replicate_coverage(
    manifests: list[dict[str, object]],
    *,
    expected_replicates: int,
) -> None:
    """Require each replicate index exactly once across merged shard manifests."""
    observed: list[int] = []
    for manifest in manifests:
        observed.extend(int(index) for index in manifest["replicate_indices"])
    duplicated = sorted(index for index in set(observed) if observed.count(index) > 1)
    if duplicated:
        raise ValueError(f"Duplicate replicate index assignment after merge: {duplicated!r}")
    missing = sorted(set(range(expected_replicates)) - set(observed))
    if missing:
        raise ValueError(f"Missing replicate index assignment after merge: {missing!r}")


def iter_expected_shard_dirs(output_dir: Path, shard_count: int) -> tuple[Path, ...]:
    """Return the canonical shard directories for a completed array job."""
    if shard_count <= 0:
        raise ValueError(f"shard_count must be positive; got {shard_count!r}.")
    return tuple(output_dir / "shards" / f"shard_{index:04d}" for index in range(shard_count))


def sync_path_to_s3(local_path: Path, s3_uri: str) -> None:
    """Sync a local output directory to S3 with the AWS CLI."""
    subprocess.run(("aws", "s3", "sync", str(local_path), s3_uri), check=True)


def sync_s3_to_path(s3_uri: str, local_path: Path) -> None:
    """Sync S3 outputs into the local filesystem."""
    local_path.mkdir(parents=True, exist_ok=True)
    subprocess.run(("aws", "s3", "sync", s3_uri, str(local_path)), check=True)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_shard_table(shard_dir: Path, filename: str) -> pd.DataFrame:
    """Read one shard's CSV output, failing loudly when the shard did not write it."""
    path = shard_dir / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing shard table: {path}")
    return pd.read_csv(path)


def load_shard_manifest(shard_dir: Path, filename: str) -> dict[str, object]:
    """Read one shard's JSON manifest, failing loudly when it is absent."""
    path = shard_dir / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing shard manifest: {path}")
    return json.loads(path.read_text())


def resolve_synced_input(path: Path, input_work_dir: Path, *, use_synced_input: bool) -> Path:
    """Resolve a relative input path against the synced input directory when one is used."""
    if path.is_absolute() or not use_synced_input:
        return path
    return input_work_dir / path


__all__ = [
    "load_shard_manifest",
    "load_shard_table",
    "resolve_synced_input",
    "iter_expected_shard_dirs",
    "sync_path_to_s3",
    "sync_s3_to_path",
    "validate_replicate_coverage",
    "validate_shard_contract",
    "write_json",
]
