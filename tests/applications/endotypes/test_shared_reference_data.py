from __future__ import annotations

import json
from pathlib import Path
from urllib.error import URLError
from urllib.parse import parse_qs

import pytest

from applications.endotypes import _shared
from applications.output_dirs import PROJECT_ROOT


class _Response:
    def __init__(self, payload: object) -> None:
        self._payload = payload

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def read(self) -> bytes:
        return json.dumps(self._payload).encode("utf-8")


def test_parse_reference_endotypes_uses_one_canonical_schema(tmp_path) -> None:
    reference_path = tmp_path / "reference.tsv"
    reference_path.write_text(
        "\n".join(
            [
                "# ignored",
                "short\trow",
                "x\tx\tx\tx\t2\t\trank-2\t\t\t101\t102",
            ]
        ),
        encoding="utf-8",
    )

    expected = {
        "reference_cluster_id": 2,
        "reference_cluster_color": "#808080",
        "reference_cluster_rank": "rank-2",
        "reference_cluster_name": "Cluster 2",
        "reference_cluster_type": "cluster",
    }
    assert _shared.parse_reference_endotypes(reference_path) == {
        "101": expected,
        "102": expected,
    }


def test_map_symbols_to_entrez_batches_and_marks_unresolved(monkeypatch) -> None:
    requests = []
    payloads = iter(
        [
            [{"query": "A", "entrezgene": 101}, {"query": "B", "_id": "ENSG00000000001"}],
            {"query": "C", "_id": "303"},
        ]
    )

    def fake_urlopen(request, *, timeout):
        requests.append((request, timeout))
        return _Response(next(payloads))

    monkeypatch.setattr(_shared, "urlopen", fake_urlopen)

    assert _shared.map_symbols_to_entrez(["A", "B", "C"], batch_size=2) == {
        "A": "101",
        "B": None,
        "C": "303",
    }
    assert [parse_qs(request.data.decode())["q"] for request, _ in requests] == [
        ["A,B"],
        ["C"],
    ]
    assert all(parse_qs(request.data.decode())["size"] == ["1"] for request, _ in requests)
    assert [timeout for _, timeout in requests] == [60, 60]


def test_map_symbols_to_entrez_rejects_invalid_batch_size() -> None:
    with pytest.raises(ValueError, match="batch_size must be positive"):
        _shared.map_symbols_to_entrez(["A"], batch_size=0)


def test_map_symbols_to_entrez_wraps_transport_errors(monkeypatch) -> None:
    def fail_urlopen(request, *, timeout):
        del request, timeout
        raise URLError("offline")

    monkeypatch.setattr(_shared, "urlopen", fail_urlopen)

    with pytest.raises(RuntimeError, match="Failed to query mygene.info"):
        _shared.map_symbols_to_entrez(["A"])


def test_parse_go_term_returns_empty_id_when_label_has_none() -> None:
    assert _shared.parse_go_term("response to x (GO:0001234)") == ("response to x", "GO:0001234")
    assert _shared.parse_go_term(" plain label ") == ("plain label", "")
    assert _shared.parse_go_term(float("nan")) == ("", "")


def test_default_output_dir_is_a_timestamped_results_analyses_folder() -> None:
    path = _shared.default_output_dir(Path("feature_matrix_Demo.tsv"), None, "suffix")

    assert path.parent == PROJECT_ROOT / "results" / "analyses"
    assert path.name.startswith("demo_suffix_")
