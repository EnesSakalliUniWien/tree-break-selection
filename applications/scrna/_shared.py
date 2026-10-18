"""Shared helpers for scRNA application reports and audits."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def add_generated_at(fig, generated_at: str) -> None:
    """Stamp a figure's lower-right corner with its generation time."""
    fig.text(
        0.99,
        0.01,
        f"Generated at: {generated_at}",
        ha="right",
        va="bottom",
        fontsize=7,
        color="#6b7280",
    )


def json_default(value: Any) -> Any:
    """Serialize standard scRNA report values without changing their schema."""
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def top_counts_text(values: pd.Series, n: int = 8) -> str:
    """Return the most frequent string values in stable report form."""
    counts = values.astype(str).value_counts()
    return "; ".join(f"{label}:{count}" for label, count in counts.head(n).items())


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest for a generated artifact."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_record(
    path: Path,
    *,
    relative_to: Path,
    role: str | None = None,
    generated_at: str | None = None,
) -> dict[str, object]:
    """Return a static-artifact provenance record."""

    record: dict[str, object] = {
        "path": str(path.relative_to(relative_to)),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if role is not None:
        record["role"] = role
    if path.suffix == ".csv":
        table = pd.read_csv(path, low_memory=False)
        record["rows"] = int(len(table))
        record["columns"] = int(len(table.columns))
        if generated_at is not None:
            record["generated_at"] = generated_at
        if "generated_at" in table.columns:
            record["generated_at_values"] = sorted(table["generated_at"].dropna().unique().tolist())
    return record


_BENCHMARK_ROOT = Path(__file__).resolve().parents[2] / "raw/assets/benchmark-results"


@dataclass(frozen=True)
class DatasetSpec:
    """One retained scRNA benchmark run audited by the analysis scripts."""

    key: str
    label: str
    output_dir: Path


@dataclass(frozen=True)
class AuditMethodSpec:
    """One TBS variant whose saved assignments and tree edges are audited."""

    geometry: str
    variant: str
    method_label: str
    assignment_key: str

    @property
    def edge_csv(self) -> str:
        return f"{self.assignment_key}_tree_edges.csv"


DATASETS = (
    DatasetSpec(
        key="adult_pancreas",
        label="Adult pancreas",
        output_dir=_BENCHMARK_ROOT / "pancreas_scrna_cluster_benchmark_20260623",
    ),
    DatasetSpec(
        key="goncalves_fetal",
        label="Goncalves fetal pancreas",
        output_dir=_BENCHMARK_ROOT / "goncalves_fetal_pancreas_progenitor_benchmark_20260624",
    ),
)

_SUFFIX = "projected_adaptive_k90_alpha0p01_edge0p001"
_LABEL_SUFFIX = "projected adaptive-k90 alpha=0.01 edge=0.001"
AUDIT_METHODS = (
    AuditMethodSpec(
        "pca_linkage", "topology_only", f"TBS topology {_LABEL_SUFFIX}", f"tbs_topology_{_SUFFIX}"
    ),
    AuditMethodSpec(
        "pca_linkage",
        "branch_time_nnls",
        f"TBS branch-time recomputed-NNLS {_LABEL_SUFFIX}",
        f"tbs_branch_time_recomputed_nnls_{_SUFFIX}",
    ),
    AuditMethodSpec(
        "pca_linkage",
        "branch_time_raw_linkage",
        f"TBS raw-linkage branch-time diagnostic {_LABEL_SUFFIX}",
        f"tbs_raw_linkage_branch_time_diagnostic_{_SUFFIX}",
    ),
    AuditMethodSpec(
        "adaptive_diffusion",
        "topology_only",
        f"TBS adaptive diffusion topology {_LABEL_SUFFIX}",
        f"tbs_adaptive_diffusion_topology_{_SUFFIX}",
    ),
    AuditMethodSpec(
        "adaptive_diffusion",
        "branch_time_nnls",
        f"TBS adaptive diffusion branch-time recomputed-NNLS {_LABEL_SUFFIX}",
        f"tbs_adaptive_diffusion_branch_time_recomputed_nnls_{_SUFFIX}",
    ),
    AuditMethodSpec(
        "adaptive_diffusion",
        "branch_time_raw_linkage",
        f"TBS adaptive diffusion raw-linkage branch-time diagnostic {_LABEL_SUFFIX}",
        f"tbs_adaptive_diffusion_raw_linkage_branch_time_diagnostic_{_SUFFIX}",
    ),
)
