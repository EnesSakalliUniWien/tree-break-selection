"""Reanalyse all surviving-node paths in the saved native context panel."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict
from math import comb
from pathlib import Path

import numpy as np

from benchmarks.diagnostics.calibration.reporting import diagnostic_json_default
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_path_context import (
    compare_surviving_paths,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    topology_from_linkage,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_context import (
    _index_tree,
)
from reports.subtree_spectral_alignment_20261109.alternatives.run import REGIMES, dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NATIVE = ROOT / "reports/subtree_context_statistics_20261109/native"
REFS = ("symmetric_average", "markov_ward")
TOL = 1e-10


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with gzip.open(path, "rt") as stream:
        return json.load(stream)


class SavedFit:
    """Index recorded clades, actual half-row masks, branch fits and traversal."""

    def __init__(self, values, saved):
        self.tree = topology_from_linkage(np.asarray(saved["linkage"]))
        self.root, self.parents, self.children = _index_tree(self.tree)
        self.nodes = {tuple(n["leaves"]): n for n in saved["nodes"].values()}
        self.ids = {clade: node["native_node_id"] for clade, node in self.nodes.items()}
        masks = {}
        for node in self.nodes.values():
            for edge in node["outgoing_edges"]:
                clade = tuple(edge["leaves"])
                self.ids[clade] = edge["annotations"]["node_id"]
                masks[clade] = tuple(edge["test_rows"])
                assert masks[clade] == tuple(sorted(clade)[::2])
        assert self.ids.keys() == self.parents.keys()
        self.clades = {node: clade for clade, node in self.ids.items()}
        self.annotations = {self.clades[a["node_id"]]: a for a in saved["trace"]["annotations"]}
        self.visits = {
            self.clades[t["node_id"]]: t for t in saved["trace"]["traversal"]["traversal_trace"]
        }
        lengths = {
            (self.clades[parent], self.clades[child]): length
            for parent, child, length in saved["trace"]["branches"]
        }
        self.full_means = {c: values[list(c)].mean(axis=0) for c in self.parents}
        self.child_means = {c: values[list(rows)].mean(axis=0) for c, rows in masks.items()}
        self.parent_means = {
            p: values[[i for child in children for i in masks[child]]].mean(axis=0)
            for p, children in self.children.items()
        }
        self.potential = {self.root: np.zeros(values.shape[1])}
        self.prefix = {self.root: np.zeros(5)}
        self.context_gap = {}
        for node in sorted(self.parents, key=len, reverse=True):
            parent = self.parents[node]
            if parent is None:
                continue
            delta = self.child_means[node] - self.parent_means[parent]
            self.potential[node] = self.potential[parent] + delta
            a = self.annotations[node]
            length = lengths[parent, node]
            self.prefix[node] = self.prefix[parent] + np.array(
                [
                    length,
                    a["Child_Parent_Divergence_Tested"],
                    a["Child_Parent_Divergence_Ancestor_Blocked"],
                    a["Child_Parent_Divergence_Significant"],
                    length == 0,
                ]
            )
            if node in self.children:
                self.context_gap[node] = float(
                    np.linalg.norm(self.parent_means[node] - self.child_means[node])
                )

    def input_mean(self, node, neighbor):
        return self.child_means[node] if self.parents[node] == neighbor else self.parent_means[node]

    def path_metrics(self, path, lca):
        start, end = path[0], path[-1]
        aggregate = self.prefix[start] + self.prefix[end] - 2 * self.prefix[lca]
        half_delta = self.potential[end] - self.potential[start]
        endpoint_half_delta = self.input_mean(end, path[-2]) - self.input_mean(start, path[1])
        return {
            "edge_count": len(path) - 1,
            "branch_sum": float(aggregate[0]),
            "bh_tested_edges": int(aggregate[1]),
            "ancestor_blocked_edges": int(aggregate[2]),
            "significant_edges": int(aggregate[3]),
            "zero_length_edges": int(aggregate[4]),
            "visited_vertices": sum(node in self.visits for node in path),
            "start_visited": start in self.visits,
            "end_visited": end in self.visits,
            "half_context_switches": sum(
                self.context_gap.get(node, 0) > TOL for node in path[1:-1] if node != lca
            ),
            "half_continuity_error": float(np.linalg.norm(half_delta - endpoint_half_delta)),
        }

    def node_statistics(self, node):
        record = self.nodes.get(node)
        annotation = self.annotations[node]
        trace = self.visits.get(node)
        return {
            "native_node_id": self.ids[node],
            "leaves": node,
            "parent": self.parents[node],
            "children": self.children.get(node),
            "incoming_edge_T": None
            if node == self.root
            else annotation["Child_Parent_Divergence_Test_Statistic"],
            "incoming_edge_df": None
            if node == self.root
            else annotation["Child_Parent_Divergence_df"],
            "incoming_edge_tested": None
            if node == self.root
            else annotation["Child_Parent_Divergence_Tested"],
            "incoming_ancestor_blocked": None
            if node == self.root
            else annotation["Child_Parent_Divergence_Ancestor_Blocked"],
            "sibling_raw": None if record is None else record["raw_sibling"],
            "sibling_gate_state": "not_applicable_leaf"
            if record is None
            else record["sibling_gate_state"],
            "sibling_branch_factor": None if record is None else record["sibling_branch_factor"],
            "visited": trace is not None,
            "decision": None if trace is None else trace["decision"],
            "final_boundary": None if trace is None else trace["final_boundary"],
        }


def row_for(context, reference, alternative):
    c, r, a = context, reference, alternative
    rp, ap = c.reference_path, c.alternative_path
    rm, am = r.path_metrics(rp, c.reference_lca), a.path_metrics(ap, c.alternative_lca)
    internal = sum(len(node) > 1 for node in (c.start, c.end))
    old_neighbors = set().union(*c.reference_neighborhood.values())
    new_neighbors = set().union(*c.alternative_neighborhood.values())
    np.testing.assert_array_equal(
        r.full_means[c.end] - r.full_means[c.start],
        a.full_means[c.end] - a.full_means[c.start],
    )
    row = {
        "start_reference_id": r.ids[c.start],
        "end_reference_id": r.ids[c.end],
        "start_alternative_id": a.ids[c.start],
        "end_alternative_id": a.ids[c.end],
        "start_size": len(c.start),
        "end_size": len(c.end),
        "relation": c.relation,
        "endpoint_kind": ("leaf_leaf", "internal_leaf", "internal_internal")[internal],
        "reference_path": ";".join(r.ids[n] for n in rp),
        "alternative_path": ";".join(a.ids[n] for n in ap),
        "missing_intermediates": ";".join(r.ids[n] for n in c.missing_intermediates),
        "surviving_intermediates_off_path": ";".join(
            r.ids[n] for n in c.surviving_intermediates_off_path
        ),
        "missing_neighborhood": ";".join(r.ids[n] for n in c.missing_neighborhood),
        "surviving_neighborhood_displaced": ";".join(
            r.ids[n] for n in c.surviving_neighborhood_displaced
        ),
        "missing_intermediate_count": len(c.missing_intermediates),
        "inserted_intermediate_count": len(c.inserted_intermediates),
        "missing_neighborhood_count": len(c.missing_neighborhood),
        "same_path": rp == ap,
        "same_neighborhood": old_neighbors == new_neighbors,
        "lca_changed": c.reference_lca != c.alternative_lca,
        "reference_lca_absent": c.reference_lca not in a.parents,
        "half_path_increment_change": float(
            np.linalg.norm(
                (r.potential[c.end] - r.potential[c.start])
                - (a.potential[c.end] - a.potential[c.start])
            )
        ),
    }
    row.update({"reference_" + k: v for k, v in rm.items()})
    row.update({"alternative_" + k: v for k, v in am.items()})
    return row


def update_counts(counts, row):
    flags = {
        "paths": True,
        "missing_intermediate": row["missing_intermediate_count"] > 0,
        "all_reference_intermediates_missing": row["reference_edge_count"] > 1
        and row["missing_intermediate_count"] == row["reference_edge_count"] - 1,
        "survivor_moved_off_path": bool(row["surviving_intermediates_off_path"]),
        "missing_neighborhood": row["missing_neighborhood_count"] > 0,
        "surviving_neighborhood_displaced": bool(row["surviving_neighborhood_displaced"]),
        "same_path": row["same_path"],
        "same_path_and_neighborhood": row["same_path"] and row["same_neighborhood"],
        "lca_changed": row["lca_changed"],
        "reference_lca_absent": row["reference_lca_absent"],
        "path_edge_count_changed": row["reference_edge_count"] != row["alternative_edge_count"],
        "path_branch_sum_changed": not np.isclose(
            row["reference_branch_sum"], row["alternative_branch_sum"], atol=TOL, rtol=TOL
        ),
        "endpoint_visit_changed": any(
            row["reference_" + k] != row["alternative_" + k]
            for k in ("start_visited", "end_visited")
        ),
        "bh_tested_edge_count_changed": row["reference_bh_tested_edges"]
        != row["alternative_bh_tested_edges"],
        "half_path_increment_changed": row["half_path_increment_change"] > TOL,
        "reference_half_non_telescope": row["reference_half_continuity_error"] > TOL,
        "alternative_half_non_telescope": row["alternative_half_continuity_error"] > TOL,
        "reference_half_context_switch": row["reference_half_context_switches"] > 0,
        "alternative_half_context_switch": row["alternative_half_context_switches"] > 0,
    }
    counts.update({k: int(v) for k, v in flags.items()})


def run(output=HERE / "panel"):
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    input_paths = sorted(NATIVE.glob("fit_*.json.gz"))
    assert len(input_paths) == 160
    sources = [
        Path(__file__),
        ROOT / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_path_context.py",
        ROOT / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_context.py",
        ROOT
        / "benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_candidate.py",
        ROOT / "reports/subtree_spectral_alignment_20261109/alternatives/run.py",
    ]
    before = {str(p.relative_to(ROOT)): digest(p) for p in input_paths + sources}
    (output / "manifest.json").write_text(json.dumps({"sha256": before}, indent=2))
    counts = defaultdict(Counter)
    examples, pair_counts, edge_rows = {}, [], []
    writer = None
    with gzip.open(output / "paths.csv.gz", "wt", compresslevel=1, newline="") as stream:
        for case in REGIMES:
            for replicate in range(4):
                values, _, seed = dataset(case, replicate)
                fits = {}
                for path in NATIVE.glob(f"fit_{case}_{replicate}_*.json.gz"):
                    saved = read(path)
                    assert saved["status"] == "ok" and saved["seed"] == seed
                    assert saved["data_hash"] == hashlib.sha256(values.tobytes()).hexdigest()
                    fits[saved["method"]] = SavedFit(values, saved)
                assert len(fits) == 10
                for reference_name in REFS:
                    r = fits[reference_name]
                    for alternative_name, a in sorted(fits.items()):
                        if not alternative_name.startswith("cov_"):
                            continue
                        metadata = {
                            "case": case,
                            "replicate": replicate,
                            "reference_method": reference_name,
                            "alternative_method": alternative_name,
                        }
                        n = 0
                        for c in compare_surviving_paths(r.tree, a.tree):
                            row = {**metadata, **row_for(c, r, a)}
                            if writer is None:
                                writer = csv.DictWriter(stream, fieldnames=list(row))
                                writer.writeheader()
                            writer.writerow(row)
                            group = c.relation + ":" + row["endpoint_kind"]
                            update_counts(counts["all"], row)
                            update_counts(counts[group], row)
                            if c.missing_intermediates:
                                update_counts(counts[group + ":missing_intermediate"], row)
                            if row["same_path"]:
                                update_counts(counts[group + ":same_path"], row)
                            conditions = {
                                "ancestor_missing": c.relation == "ancestor_descendant"
                                and bool(c.missing_intermediates),
                                "cross_lca_absent": row["reference_lca_absent"],
                                "cross_lca_survives_off_path": c.relation == "cross_branch"
                                and row["lca_changed"]
                                and not row["reference_lca_absent"],
                                "same_path_changed_neighborhood": row["same_path"]
                                and not row["same_neighborhood"],
                                "same_path_changed_visit": row["same_path"]
                                and any(
                                    row["reference_" + k] != row["alternative_" + k]
                                    for k in ("start_visited", "end_visited")
                                ),
                            }
                            if row["endpoint_kind"] == "internal_internal":
                                for label, selected in conditions.items():
                                    if selected and label not in examples:
                                        examples[label] = {
                                            **row,
                                            "context": asdict(c),
                                            "reference_endpoints": [
                                                r.node_statistics(v) for v in (c.start, c.end)
                                            ],
                                            "alternative_endpoints": [
                                                a.node_statistics(v) for v in (c.start, c.end)
                                            ],
                                        }
                            n += 1
                        common = r.parents.keys() & a.parents.keys()
                        assert n == comb(len(common), 2)
                        pair_counts.append({**metadata, "shared_clades": len(common), "paths": n})
                        for child in sorted(common):
                            parent = r.parents[child]
                            if parent is None or a.parents[child] != parent:
                                continue
                            left, right = r.annotations[child], a.annotations[child]
                            edge = {**metadata, "parent": parent, "child": child}
                            for suffix in (
                                "Test_Statistic",
                                "df",
                                "Tested",
                                "Ancestor_Blocked",
                                "Significant",
                            ):
                                field = "Child_Parent_Divergence_" + suffix
                                edge["reference_" + suffix] = left[field]
                                edge["alternative_" + suffix] = right[field]
                            edge_rows.append(edge)
                print(f"{case} {replicate}: paths recorded", flush=True)
    changed = [p for p, h in before.items() if digest(ROOT / p) != h]
    assert not changed, changed
    summary = {
        "datasets": 16,
        "saved_fits": 160,
        "new_fits": 0,
        "tree_pairs": len(pair_counts),
        "tolerance": TOL,
        "counts": dict(counts),
        "pairs": pair_counts,
        "matched_edges": len(edge_rows),
        "matched_edge_raw_T_changed": sum(
            not np.isclose(
                e["reference_Test_Statistic"], e["alternative_Test_Statistic"], atol=TOL, rtol=TOL
            )
            for e in edge_rows
        ),
        "matched_edge_df_changed": sum(e["reference_df"] != e["alternative_df"] for e in edge_rows),
        "matched_edge_bh_tested_changed": sum(
            e["reference_Tested"] != e["alternative_Tested"] for e in edge_rows
        ),
        "examples": examples,
        "source_and_inputs_changed": changed,
    }
    with gzip.open(output / "matched_edges.json.gz", "wt") as stream:
        json.dump(edge_rows, stream, default=diagnostic_json_default)
    (output / "summary.json").write_text(
        json.dumps(summary, default=diagnostic_json_default, indent=2)
    )
    print(
        json.dumps(
            {
                "paths": counts["all"]["paths"],
                "tree_pairs": len(pair_counts),
                "examples": list(examples),
            }
        )
    )


if __name__ == "__main__":
    run()
