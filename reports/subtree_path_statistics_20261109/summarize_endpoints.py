"""Compare endpoint decisions and missing-neighborhood roles in saved paths."""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict

from reports.subtree_path_statistics_20261109.run import HERE, NATIVE, digest, read


def summarize():
    panel = HERE / "panel"
    fits, hashes = {}, {}
    for path in sorted(NATIVE.glob("fit_*.json.gz")):
        r = read(path)
        parents = {child: parent for parent, child, _ in r["trace"]["branches"]}
        children = defaultdict(set)
        for child, parent in parents.items():
            children[parent].add(child)
        fits[r["case"], r["replicate"], r["method"]] = {
            "parents": parents,
            "children": children,
            "trace": {t["node_id"]: t for t in r["trace"]["traversal"]["traversal_trace"]},
            "nodes": {n["native_node_id"]: n for n in r["nodes"].values()},
        }
        hashes[str(path)] = digest(path)
    path_file = panel / "paths.csv.gz"
    hashes[str(path_file)] = digest(path_file)
    counts, examples = defaultdict(Counter), {}
    n = 0
    with gzip.open(path_file, "rt", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["endpoint_kind"] != "internal_internal":
                continue
            key = row["case"], int(row["replicate"])
            r, a = (fits[*key, row[side + "_method"]] for side in ("reference", "alternative"))
            endpoints = {}
            for side, fit in (("reference", r), ("alternative", a)):
                endpoints[side] = []
                for end in ("start", "end"):
                    node_id = row[end + "_" + side + "_id"]
                    visit = fit["trace"].get(node_id)
                    endpoints[side].append(
                        {
                            "node_id": node_id,
                            "visited": visit is not None,
                            "decision": None if visit is None else visit["decision"],
                            "final_boundary": None if visit is None else visit["final_boundary"],
                            "sibling_gate_state": fit["nodes"][node_id]["sibling_gate_state"],
                        }
                    )
            flags = {"paths": True}
            for field in ("visited", "decision", "final_boundary", "sibling_gate_state"):
                flags[field + "_changed"] = any(
                    x[field] != y[field]
                    for x, y in zip(endpoints["reference"], endpoints["alternative"], strict=True)
                )
            old_path = set(row["reference_path"].split(";"))
            missing = set(row["missing_neighborhood"].split(";")) - {""}
            neighbor_roles = {
                "parent": {r["parents"].get(node) for node in old_path} - {None},
                "child": {c for node in old_path for c in r["children"].get(node, ())},
                "sibling": {
                    sibling
                    for node in old_path
                    for sibling in r["children"].get(r["parents"].get(node), ())
                    if sibling != node
                },
            }
            for role, members in neighbor_roles.items():
                flags["missing_neighbor_" + role] = bool((members - old_path) & missing)
            groups = ["all", row["relation"]]
            if row["same_path"] == "True":
                groups.append("same_path")
                if row["same_neighborhood"] == "True":
                    groups.append("same_path_and_neighborhood")
            if int(row["missing_intermediate_count"]):
                groups.append("missing_intermediate")
            if missing:
                groups.append("missing_neighborhood")
            for group in groups:
                counts[group].update({k: int(v) for k, v in flags.items()})
            if flags["decision_changed"]:
                label = "same_path" if row["same_path"] == "True" else "changed_path"
                if label not in examples:
                    examples[label] = {"path": row, "endpoints": endpoints}
            n += 1
    main = json.loads((panel / "summary.json").read_text())
    expected = sum(
        v["paths"] for k, v in main["counts"].items() if k.endswith(":internal_internal")
    )
    assert n == expected
    assert counts["all"]["visited_changed"] == sum(
        v["endpoint_visit_changed"]
        for k, v in main["counts"].items()
        if k.endswith(":internal_internal")
    )
    assert counts["all"]["visited_changed"] <= counts["all"]["decision_changed"]
    assert all(digest(path) == h for path, h in hashes.items())
    out = {
        "scope": "internal/internal endpoints",
        "counts": dict(counts),
        "examples": examples,
        "sha256": hashes,
    }
    target = panel / "endpoint_summary.json"
    with target.open("x") as stream:
        json.dump(out, stream, indent=2)
    print(json.dumps(dict(counts), indent=2))


if __name__ == "__main__":
    summarize()
