"""Summarize saved native node-context comparisons without rerunning fits."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_spectral_alignment as align,
)
from reports.subtree_context_statistics_20261109.run import covariance_span_record
from reports.subtree_spectral_alignment_20261109.alternatives.run import dataset
from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NATIVE = HERE / "native"
PRIOR = HERE.parent / "subtree_spectral_alignment_20261109/parent_child/rows.jsonl.gz"
REFS = {"symmetric_average": "symmetric", "markov_ward": "markov"}
DF = "collected_sibling_degrees_of_freedom"
T = "collected_sibling_stat"
KEYS = ["case", "replicate", "reference_method", "alternative_method", "reference_node_id"]


def read(path, lines=False):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return [json.loads(line) for line in stream] if lines else json.load(stream)


def describe(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    return (
        dict(n=len(x), min=float(x.min()), median=float(np.median(x)), max=float(x.max()))
        if len(x)
        else dict(n=0, min=None, median=None, max=None)
    )


def paired(frame, field):
    a, b = (
        pd.to_numeric(frame[f"{arm}_{field}"], errors="coerce")
        for arm in ("reference", "alternative")
    )
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    delta = b - a
    changed = ~np.isclose(a, b, atol=1e-10, rtol=1e-10)
    denominator = np.maximum(np.maximum(np.abs(a), np.abs(b)), np.finfo(float).tiny)
    return describe(delta) | {
        "changed": int(changed.sum()),
        "median_abs_delta": float(np.median(np.abs(delta))) if len(delta) else None,
        "changed_abs_delta": describe(np.abs(delta[changed])),
        "changed_relative_delta": describe(np.abs(delta[changed]) / denominator[changed]),
    }


def transitions(frame, field):
    a = frame["reference_" + field].astype("string").fillna("not_recorded")
    b = (
        frame["alternative_" + field]
        .astype("string")
        .fillna("not_recorded")
        .where(frame.node_exists, "absent")
    )
    return [
        {"reference": left, "alternative": right, "n": int(count)}
        for (left, right), count in pd.DataFrame({"a": a, "b": b}).value_counts().items()
    ]


def context(frame):
    return {
        "comparisons": len(frame),
        "datasets": len(frame[["case", "replicate"]].drop_duplicates()),
        "T": paired(frame, T),
        "df": paired(frame, DF),
        "T_over_df": paired(frame, "T_over_df"),
        "branch_factor": paired(frame, "sibling_branch_factor"),
        "branch_unscaled_T": paired(frame, "branch_unscaled_T"),
        "inflation": paired(frame, "Sibling_Inflation_Factor"),
        **{
            name: int(frame[name].eq(True).sum())
            for name in (
                "same_incoming_test_rows",
                "same_incoming_basis_rows",
                "same_outgoing_basis_rows",
            )
        },
        **{
            field: transitions(frame, field)
            for field in ("sibling_gate_state", "visited", "traversal_decision", "final_boundary")
        },
    }


def structure(frame):
    assert not frame.duplicated(KEYS).any()
    frame["root"] = frame.original_parent.isna()
    assert not frame.loc[frame.root, ["case_1", "case_3", "case_4"]].any().any()
    assert frame.case_4.eq(frame.case_2 & frame.case_3).all()
    assert frame.case_1.eq(~frame.node_exists).all()
    assert frame.loc[~frame.node_exists, "alternative_" + T].isna().all()
    for arm in ("reference", "alternative"):
        stat = pd.to_numeric(frame[f"{arm}_{T}"], errors="coerce")
        df = pd.to_numeric(frame[f"{arm}_{DF}"], errors="coerce")
        frame[f"{arm}_T_over_df"] = stat / df.where(df > 0)
        frame[f"{arm}_branch_unscaled_T"] = stat * pd.to_numeric(
            frame[f"{arm}_sibling_branch_factor"], errors="coerce"
        )
    cases = {
        name: int(frame[name].sum())
        for name in ("case_1", "case_2", "case_3", "case_4", "at_least_one_child_absent")
    }
    cases["case_2_present_node"] = int((frame.case_2 & frame.node_exists).sum())
    cases["root_case_2"] = int((frame.root & frame.case_2).sum())
    present = frame.loc[frame.node_exists & ~frame.root]
    contexts = {}
    for parent in (True, False):
        for children in (True, False):
            name = f"parent_{'same' if parent else 'changed'}_children_{'same' if children else 'changed'}"
            mask = present.direct_parent_preserved.eq(
                parent
            ) & present.direct_children_preserved.eq(children)
            contexts[name] = context(present.loc[mask])
    assert sum(stats["comparisons"] for stats in contexts.values()) == len(present)
    same_children = present.loc[present.direct_children_preserved.eq(True)]
    size_presence = []
    for lower, upper in ((2, 3), (4, 7), (8, 15), (16, 31), (32, 63), (64, 64)):
        selected = frame.loc[frame.node_size.between(lower, upper)]
        size_presence.append(
            {
                "range": [lower, upper],
                "comparisons": len(selected),
                "present": int(selected.node_exists.sum()),
                "absent": int((~selected.node_exists).sum()),
            }
        )
    roots = {
        key: context(frame.loc[frame.root & frame.direct_children_preserved.eq(keep)])
        for key, keep in (("root_children_same", True), ("root_children_changed", False))
    }
    examples = {}
    cols = KEYS + [
        "node",
        "original_parent",
        "alternative_parent",
        "original_children",
        "alternative_children",
    ]
    cols += [
        f"{arm}_{field}"
        for arm in ("reference", "alternative")
        for field in (
            T,
            DF,
            "sibling_branch_factor",
            "Sibling_Inflation_Factor",
            "sibling_gate_state",
            "visited",
        )
    ]
    numeric = frame.loc[
        np.isfinite(frame.reference_T_over_df) & np.isfinite(frame.alternative_T_over_df)
    ]
    for name, mask in (
        ("case4", numeric.case_4),
        (
            "same_children_changed_parent",
            numeric.node_exists
            & ~numeric.root
            & numeric.direct_parent_preserved.eq(False)
            & numeric.direct_children_preserved,
        ),
    ):
        rows = numeric.loc[mask].sort_values(KEYS)
        if len(rows):
            examples[name] = json.loads(rows[cols].iloc[:1].to_json(orient="records"))[0]
    return {
        "comparisons": len(frame),
        "root_comparisons": int(frame.root.sum()),
        "present_nonroot": len(present),
        "cases": cases,
        "nonroot_2x2": contexts,
        "same_children_nonroot_branch_unscaled_T": paired(same_children, "branch_unscaled_T"),
        "presence_by_node_size": size_presence,
        "root": roots,
        "all_gate_transitions": transitions(frame, "sibling_gate_state"),
        "all_visit_transitions": transitions(frame, "visited"),
        "finite_raw_p_without_active_alternative_BH": int(
            (
                frame.node_exists
                & np.isfinite(frame.alternative_collected_sibling_p_value)
                & ~frame.alternative_sibling_gate_state.isin(["active_reject", "active_nonreject"])
            ).sum()
        ),
        "examples": examples,
    }


def span_group_summary(table):
    return [
        {
            "method": method,
            "root": bool(root),
            "n": len(rows),
            "same_span": int(rows.same.sum()),
            "rank_mismatch": int(rows.rank_mismatch.sum()),
            "projector_gap": describe(rows.gap),
            "equal_rank_projector_gap": describe(rows.loc[~rows.rank_mismatch, "gap"]),
            "equal_rank_different_span": int((~rows.rank_mismatch & ~rows.same).sum()),
            "local_share": describe(rows.local_share),
        }
        for (method, root), rows in table.groupby(["method", "root"])
    ]


def trace_summary(attempts):
    previous = {(r["case"], r["replicate"], r["operator"]): r for r in read(PRIOR, lines=True)}
    regenerated, covariance = {}, {}
    parity, nodes, paths = [], [], []
    corrected_nodes, sidecar = [], []
    worst = 0.0
    traced_nodes = 0
    for item in attempts:
        case, rep, method = item["case"], item["replicate"], item["method"]
        key = case, rep
        if key not in regenerated:
            x, _, seed = dataset(case, rep)
            regenerated[key] = (x, seed, hashlib.sha256(x.tobytes()).hexdigest())
        x, seed, data_hash = regenerated[key]
        assert item["data_hash"] == data_hash and item["seed"] == seed
        path = NATIVE / f"fit_{case}_{rep}_{method}.json.gz"
        if not path.is_file():
            assert item["status"] != "ok"
            continue
        paths.append(path)
        saved = read(path)
        assert (
            saved["data_hash"] == data_hash
            and saved["seed"] == seed
            and saved["status"] == item["status"]
        )
        basis = None
        if method == "cov_full":
            basis = np.eye(x.shape[1])
        elif "covariance_blocks" in saved:
            columns = [i for start, stop in saved["covariance_blocks"] for i in range(start, stop)]
            basis = align.covariance_spectrum(x).eigenvectors[:, columns]
        for node in saved["nodes"].values():
            traced_nodes += 1
            leaves = tuple(node["leaves"])
            cache_key = case, rep, leaves
            if cache_key not in covariance:
                cov = np.cov(x[list(leaves)], rowvar=False, ddof=1)
                covariance[cache_key] = np.maximum(np.linalg.eigvalsh(cov)[::-1], 0)
            observed = np.asarray(node["raw_covariance_eigenvalues"])
            expected = covariance[cache_key]
            worst = max(worst, float(np.max(np.abs(expected - observed))))
            np.testing.assert_allclose(observed, expected, atol=1e-11, rtol=1e-11)
            sub = node.get("inherited_vs_rebuilt_covariance_span")
            if sub is not None:
                positive = int(
                    np.count_nonzero(
                        observed
                        > max(observed[0], 0) * np.finfo(float).eps * max(len(leaves), x.shape[1])
                    )
                )
                ties = tied_eigenvalue_neighbors(observed)
                rebuilt = sub["rank_b"]
                nodes.append(
                    {
                        "method": method,
                        "case": case,
                        "replicate": rep,
                        "size": len(leaves),
                        "root": node["parent"] is None,
                        "gap": sub["projector_frobenius"],
                        "same": sub["same_span"],
                        "rank_mismatch": sub["rank_a"] != sub["rank_b"],
                        "zero_local_trace": bool(observed.sum() <= 0),
                        "positive_cutoff_slices_tie": bool(
                            0 < positive < len(observed) and ties[positive - 1]
                        ),
                        "rebuilt_cutoff_slices_tie": bool(
                            0 < rebuilt < len(observed) and ties[rebuilt - 1]
                        ),
                        "local_share": node["inherited_local_variance_share"],
                    }
                )
                assert basis is not None
                local = x[list(leaves)]
                rebuilt_cov = align.covariance_spectrum(local)
                corrected = covariance_span_record(rebuilt_cov, basis, local @ basis)
                new = corrected["inherited_vs_rebuilt_covariance_span"]
                share = corrected["inherited_local_variance_share"]
                np.testing.assert_allclose(
                    share, node["inherited_local_variance_share"], rtol=1e-10, atol=1e-12
                )
                new_rank = new["rank_b"]
                assert not (0 < new_rank < len(observed) and ties[new_rank - 1])
                corrected_nodes.append(
                    nodes[-1]
                    | {
                        "gap": new["projector_frobenius"],
                        "same": new["same_span"],
                        "rank_mismatch": new["rank_a"] != new_rank,
                        "local_share": share,
                    }
                )
                sidecar.append(
                    {
                        "case": case,
                        "replicate": rep,
                        "method": method,
                        "data_hash": data_hash,
                        "native_node_id": node["native_node_id"],
                        "leaves": leaves,
                        "original_rebuilt_rank": rebuilt,
                        "original_projector_gap": sub["projector_frobenius"],
                        **corrected,
                    }
                )
        if method not in REFS:
            continue
        old = previous[(case, rep, REFS[method])]
        assert old["data_hash"] == data_hash and old["seed"] == seed
        prior, current = old["result"]["native"], saved["metrics"]
        root = next((node for node in saved["nodes"].values() if node["parent"] is None), None)
        a, b = current.get("root_sibling_p"), prior.get("root_sibling_p")
        parity.append(
            {
                "case": case,
                "replicate": rep,
                "method": method,
                "status_equal": current["status"] == prior["status"],
                "root_p_equal": a is not None
                and b is not None
                and bool(np.isclose(a, b, rtol=1e-8, atol=0, equal_nan=True)),
                "root_split_equal": root is not None
                and sorted(root["children"])
                == sorted(child["rows"] for child in old["result"]["children"]),
                "root_reject_equal": current.get("root_sibling_reject")
                == prior.get("root_sibling_reject"),
                "partition_ari": adjusted_rand_score(current["labels"], prior["labels"])
                if "labels" in current and "labels" in prior
                else None,
                "ami_difference": current["ami"] - prior["ami"]
                if "ami" in current and "ami" in prior
                else None,
            }
        )
    table = pd.DataFrame(nodes)
    gaps = span_group_summary(table)
    pairs = pd.DataFrame(parity)
    return (
        {
            "datasets": len(regenerated),
            "trace_files": len(paths),
            "traced_nodes": traced_nodes,
            "unique_dataset_clades": len(covariance),
            "max_abs_raw_covariance_eigenvalue_error": worst,
            "native_prior": {
                "pairs": len(parity),
                **{
                    field: int(pairs[field].sum())
                    for field in (
                        "status_equal",
                        "root_p_equal",
                        "root_split_equal",
                        "root_reject_equal",
                    )
                },
                "partitions_ari_one": int(pairs.partition_ari.eq(1).sum()),
                "ami_difference": describe(pairs.ami_difference),
                "detail": parity,
            },
            "inherited_vs_rebuilt": span_group_summary(pd.DataFrame(corrected_nodes)),
            "initial_inherited_vs_rebuilt": gaps,
            "corrected_span_sidecar": {
                "rows": len(sidecar),
                "rebuilt_rank_changed": sum(
                    row["original_rebuilt_rank"]
                    != row["inherited_vs_rebuilt_covariance_span"]["rank_b"]
                    for row in sidecar
                ),
                "rebuilt_cutoff_slices_tie": 0,
            },
            "inherited_boundary_audit": {
                "comparisons": len(table),
                **{
                    name: int(table[name].sum())
                    for name in (
                        "zero_local_trace",
                        "positive_cutoff_slices_tie",
                        "rebuilt_cutoff_slices_tie",
                    )
                },
            },
        },
        paths,
        sidecar,
    )


def main():
    verification_file = NATIVE / "verification.json"
    verification = read(verification_file)
    assert verification["attempts"] == 160 and not verification["source_files_changed"]
    attempts_file, comparisons_file = NATIVE / "attempts.json", NATIVE / "comparisons.jsonl.gz"
    attempts = read(attempts_file)
    assert len({(a["case"], a["replicate"], a["method"]) for a in attempts}) == len(attempts) == 160
    comparisons = read(comparisons_file, lines=True)
    assert len(comparisons) == verification["node_comparisons"] == 16128
    comparison_summary = structure(pd.DataFrame(comparisons))
    traces, trace_paths, sidecar = trace_summary(attempts)
    assert len(sidecar) == 8064
    sidecar_path = HERE / "native_covariance_spans_corrected.jsonl.gz"
    with gzip.open(sidecar_path, "wt") as stream:
        for row in sidecar:
            stream.write(json.dumps(row, default=diagnostic_json_default, allow_nan=False) + "\n")
    archive = HERE / "source_during_native_run.py.txt"
    initial_runner_hash = hashlib.sha256(archive.read_bytes()).hexdigest()
    manifest = read(NATIVE / "manifest.json")
    runner = HERE / "run.py"
    assert initial_runner_hash == manifest["source_hashes"][str(runner.relative_to(ROOT))]
    paths = [
        Path(__file__),
        archive,
        runner,
        sidecar_path,
        verification_file,
        attempts_file,
        comparisons_file,
        NATIVE / "manifest.json",
        PRIOR,
        HERE.parent / "subtree_spectral_alignment_20261109/alternatives/run.py",
        *trace_paths,
    ]
    summary = {
        "verification": verification,
        "structure": comparison_summary,
        "traces": traces,
        "source_history": {
            "initial_runner": archive.name,
            "initial_runner_sha256": initial_runner_hash,
            "corrected_helper": "run.py:covariance_span_record",
            "corrected_helper_source_sha256": hashlib.sha256(runner.read_bytes()).hexdigest(),
            "clustering_refits": 0,
        },
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths
        },
    }
    (HERE / "native_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
