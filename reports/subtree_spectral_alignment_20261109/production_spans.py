"""Compare surviving production-coordinate spans on the same saved child memberships."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_spectral_alignment as align,
)
from reports.subtree_spectral_alignment_20261109 import run

HERE = Path(__file__).resolve().parent


def coordinates(spec, features, operator):
    start, stop = run.spectral_window(spec, features)
    basis = spec.symmetric_vectors if operator == "symmetric" else spec.right_vectors
    weights = np.maximum(spec.eigenvalues[start:stop], 0) ** 3
    return basis[:, start:stop] * weights, int(np.count_nonzero(weights == 0))


def main():
    output = HERE / "parent_child/production_span_rows.jsonl.gz"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    before = run.source_hashes()
    with gzip.open(HERE / "parent_child/rows.jsonl.gz", "rt") as stream:
        original = [json.loads(line) for line in stream]
    cache, rows = {}, []
    for record in original:
        key = record["case"], record["replicate"]
        if key not in cache:
            x, _truth, _seed = run.dataset(*key)
            _geometry, parent = run.capture_diffusion(x)
            cache[key] = x, parent
        x, parent = cache[key]
        assert run.digest(x) == record["data_hash"]
        np.testing.assert_allclose(
            parent.degrees, record["result"]["diffusion"]["degrees"], rtol=0, atol=0
        )
        operator = record["operator"]
        parent_coordinates, parent_clipped = coordinates(parent, x.shape[1], operator)
        for index, child in enumerate(record["result"]["children"]):
            indices = np.asarray(child["rows"])
            induced = align.reversible_diffusion_spectrum(parent.weights[np.ix_(indices, indices)])
            geometry, rebuilt = run.capture_diffusion(x[indices], operator)
            result = {k: record[k] for k in ("case", "replicate", "operator", "seed", "data_hash")}
            result.update(child=index, parent_zero_weight_modes=parent_clipped)
            for name, spec in (("induced", induced), ("rebuilt", rebuilt)):
                child_coordinates, clipped = coordinates(spec, x.shape[1], operator)
                if name == "rebuilt":
                    np.testing.assert_allclose(
                        run.pdist(child_coordinates), geometry.distance_condensed, atol=1e-10
                    )
                inherited = align.restrict_parent_basis(parent_coordinates, indices)
                result[name] = {
                    "zero_weight_modes": clipped,
                    "inherited_parent_rank": inherited.parent_rank,
                    "rank_loss_on_restriction": inherited.rank_loss,
                    **asdict(align.compare_subspaces(inherited.basis, child_coordinates)),
                }
            rows.append(result)
    assert len(cache) == 80 and len(rows) == 320
    after = run.source_hashes()
    assert before == after
    with gzip.open(output, "xt") as stream:
        for record in rows:
            stream.write(json.dumps(record, default=diagnostic_json_default) + "\n")
    summary = {
        "status": "completed",
        "scope": "Same 80 datasets and saved 320 child memberships; no new native fits or independent observations",
        "child_rows": len(rows),
        "datasets": len(cache),
        "parent_windows_with_zero_weights": sum(
            row["parent_zero_weight_modes"] > 0 for row in rows
        ),
        "operators": {
            operator: {
                kind: {
                    "compared": len(group),
                    "child_windows_with_zero_weights": sum(
                        row[kind]["zero_weight_modes"] > 0 for row in group
                    ),
                    "same_span": sum(row[kind]["same_span"] for row in group),
                    "median_projector_gap": float(
                        np.median([row[kind]["projector_frobenius"] for row in group])
                    ),
                    "median_child_rank": float(np.median([row[kind]["rank_b"] for row in group])),
                }
                for kind in ("induced", "rebuilt")
            }
            for operator in ("symmetric", "markov")
            if (group := [row for row in rows if row["operator"] == operator])
        },
        "source_hashes": before,
    }
    (HERE / "parent_child/production_span_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "source_hashes"}, indent=2))


if __name__ == "__main__":
    main()
