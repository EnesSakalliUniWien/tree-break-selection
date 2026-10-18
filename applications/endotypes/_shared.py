"""Shared data and naming contracts for endotype application entry points."""

import csv
import json
import re
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from applications.output_dirs import timestamped_results_dir

_MYGENE_QUERY_ENDPOINT = "https://mygene.info/v3/query"
GO_ID_RE = re.compile(r"(GO:\d{7})")


def safe_name(value: object) -> str:
    """Return a filesystem-safe representation without changing case."""

    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(value))


def matrix_slug(input_path: Path, dataset_label: str | None = None) -> str:
    """Return the stable application slug for an input feature matrix."""

    raw = dataset_label or input_path.stem
    if raw.startswith("feature_matrix_"):
        raw = raw[len("feature_matrix_") :]
    return safe_name(raw).strip("_").lower() or "feature_matrix"


def default_output_dir(input_path: Path, dataset_label: str | None, suffix: str) -> Path:
    """Return a new ``results/analyses/<matrix slug>_<suffix>_<timestamp>`` path."""

    return timestamped_results_dir(f"{matrix_slug(input_path, dataset_label)}_{suffix}")


def parse_go_term(value: object) -> tuple[str, str]:
    """Split a ``term (GO:0000000)`` column label into ``(term, go_id)``; no id gives ``""``."""

    text = "" if pd.isna(value) else str(value)
    match = GO_ID_RE.search(text)
    go_id = match.group(1) if match else ""
    term = text.replace(f"({go_id})", "").strip() if go_id else text.strip()
    return term, go_id


def load_feature_matrix(path: Path) -> pd.DataFrame:
    """Load a numeric feature matrix and enforce cosine-space row support."""

    data = pd.read_csv(path, sep="\t", index_col=0)
    zero_columns = data.columns[(data.sum(axis=0) == 0).to_numpy()]
    if len(zero_columns):
        data = data.drop(columns=zero_columns)
    zero_rows = data.index[(data.sum(axis=1) == 0).to_numpy()]
    if len(zero_rows):
        raise ValueError(
            f"Rows with zero feature mass cannot enter cosine analysis: {list(zero_rows[:10])!r}"
        )
    return data.astype(float)


def load_binary_feature_matrix(
    path: Path,
    *,
    drop_zero_columns: bool = False,
    require_nonzero_rows: bool = False,
    non_binary_message: str = "contains non-binary values.",
    zero_rows_message: str = "Rows with no active GO terms",
) -> pd.DataFrame:
    """Load a tab-separated binary feature matrix."""

    data = pd.read_csv(path, sep="\t", index_col=0)
    data.index = data.index.astype(str)
    data.columns = data.columns.astype(str)
    data = data.apply(pd.to_numeric, errors="raise")
    values = data.to_numpy()
    if not np.isin(values, (0, 1)).all():
        raise ValueError(f"{path} {non_binary_message}")
    if drop_zero_columns:
        zero_columns = data.columns[(data.sum(axis=0) == 0).to_numpy()]
        if len(zero_columns):
            data = data.drop(columns=zero_columns)
    if require_nonzero_rows:
        zero_rows = data.index[(data.sum(axis=1) == 0).to_numpy()]
        if len(zero_rows):
            raise ValueError(f"{zero_rows_message}: {list(zero_rows[:10])!r}.")
    return data.astype(int)


def parse_reference_endotypes(path: Path) -> dict[str, dict[str, object]]:
    """Index a Julia reference-endotype table by Entrez gene ID."""

    entrez_to_endotype: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            first_cell = row[0].strip()
            if first_cell.startswith("#") or first_cell in {"SUM", "AVERAGE"}:
                continue
            if len(row) < 10:
                continue
            cluster_id_text = row[4].strip()
            if not cluster_id_text.isdigit():
                continue
            endotype = {
                "reference_cluster_id": int(cluster_id_text),
                "reference_cluster_color": row[5].strip() or "#808080",
                "reference_cluster_rank": row[6].strip(),
                "reference_cluster_name": row[7].strip() or f"Cluster {cluster_id_text}",
                "reference_cluster_type": row[8].strip() or "cluster",
            }
            for gene_id in (cell.strip() for cell in row[9:] if cell.strip()):
                entrez_to_endotype[gene_id] = endotype
    return entrez_to_endotype


def map_symbols_to_entrez(
    symbols: list[str],
    *,
    batch_size: int = 200,
) -> dict[str, str | None]:
    """Resolve human gene symbols to at most one Entrez ID per symbol."""

    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    mapped: dict[str, str | None] = {}
    for start in range(0, len(symbols), batch_size):
        batch = symbols[start : start + batch_size]
        body = urlencode(
            {
                "q": ",".join(batch),
                "scopes": "symbol",
                "fields": "symbol,entrezgene,taxid",
                "species": "human",
                "size": 1,
            }
        ).encode("utf-8")
        request = Request(
            _MYGENE_QUERY_ENDPOINT,
            data=body,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            method="POST",
        )
        try:
            with urlopen(request, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except URLError as exc:
            raise RuntimeError(
                "Failed to query mygene.info for gene symbol -> Entrez mapping."
            ) from exc

        records = [payload] if isinstance(payload, dict) else payload
        for record in records:
            query = record.get("query")
            if not query:
                continue
            # mygene's _id is the Entrez ID only when numeric; otherwise it is an Ensembl ID.
            record_id = str(record.get("_id", ""))
            entrez = record.get("entrezgene") or (record_id if record_id.isdigit() else None)
            mapped[str(query)] = None if entrez is None else str(entrez)
        for symbol in batch:
            mapped.setdefault(symbol, None)
    return mapped


def binary_entropy_bits(probabilities: np.ndarray) -> np.ndarray:
    """Return the Bernoulli entropy in bits of each probability in [0, 1]."""
    from scipy.special import entr

    p = np.asarray(probabilities, dtype=float)
    return (entr(p) + entr(1.0 - p)) / np.log(2.0)
