"""Small development panel, separate from any final validation gate."""

import hashlib
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.stats import binomtest
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_family_reference as family,
    subtree_gaussian_reference as gaussian,
    subtree_statistic_attribution as attribution,
    subtree_topology_tests as topology,
)

FACTOR = np.array([[1, 0.55, 0.1], [0, 0.8, 0.25], [0, 0, 1.3]])
CELLS = (
    ("attribution_null", 16, 3, 0, ()),
    ("attribution_signal", 16, 3, 6, ()),
    ("fixed_beta", 6, 2, 0, ()),
    ("family_null", 4, 2, 0, ()),
    ("family_selected_child", 8, 2, 14, (0,)),
    ("family_mixed_root", 8, 2, 4, ()),
)


def normalized_score(matrix, rows, left, right):
    if len(rows) == 2:
        return 0.0
    local = matrix[rows, 0]
    denom = ((local - local.mean()) ** 2).sum()
    delta = matrix[left, 0].mean() - matrix[right, 0].mean()
    return float(delta**2 / ((1 / len(left) + 1 / len(right)) * denom))


def beta_score(matrix):
    center = matrix - matrix.mean(axis=0)
    u, _, _ = np.linalg.svd(center, full_matrices=False)
    e = np.repeat([1.0, -1.0], 3) / np.sqrt(6)
    return float((u[:, 0] @ e) ** 2)


def trial(i, replicate):
    name, n, p, shift, path = CELLS[i]
    data_seed, chain_seed = np.random.SeedSequence([7811260911, i, replicate]).spawn(2)
    seed = int(chain_seed.generate_state(1, dtype=np.uint64)[0])
    factor = FACTOR[:p, :p]
    x = np.random.default_rng(data_seed).normal(size=(n, p)) @ factor
    means = np.zeros(n)
    if shift:
        x[n // 2 :, 0] += shift
        means[n // 2 :] = shift
    out = dict(
        cell=name,
        replicate=replicate,
        seed=seed,
        data_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
    )
    if i < 2:
        result = attribution.native_sibling_attribution(
            x, range(n), oracle_covariance=factor.T @ factor, gram_rank=1
        )
        np.testing.assert_allclose(
            result["unit_branch_raw"], result["branch_factor"] * result["native_raw"], rtol=1e-12
        )
        assert result["unit_branch_df"] == result["native_df"]
        if result["native_display"] is not None:
            np.testing.assert_allclose(
                result["native_display"], result["annotated_display"], rtol=1e-12
            )
        return out | {"attribution": result}
    if i == 2:
        exact = gaussian.fixed_subtree_contrast_test(
            x, subtree_indices=range(n), left_indices=range(3), rank=1
        )
        np.testing.assert_allclose(exact["statistic"], beta_score(x), atol=1e-13)
        out["beta_p"] = exact["p_value"]
    comparisons = []
    for steps in (99, 999):
        start = perf_counter()
        if i == 2:
            r = gaussian.subtree_gaussian_chain_test(
                x, subtree_indices=range(n), score=beta_score, steps=steps, seed=seed
            )
            record = dict(p_value=r["p_value"], decision=r["decision"], moves=r["accepted_moves"])
        else:
            r = family.ward_subtree_family_test(
                x, path=path, score=normalized_score, steps=steps, seed=seed
            )
            nodes = []
            for node in r["nodes"]:
                indices = node["subtree_indices"]
                is_null = np.ptp(means[indices]) == 0
                equal_contrast = means[node["left"]].mean() == means[node["right"]].mean()
                nodes.append(
                    dict(
                        rows=indices,
                        null=bool(is_null),
                        equal_contrast=bool(equal_contrast),
                        decision=node["family_decision"],
                        p_value=node["p_value"],
                        holm=node["holm_p_value"],
                        moves=node["accepted_moves"],
                    )
                )
            record = dict(
                decision=r["decision"],
                nodes=nodes,
                hypotheses=r.get("hypotheses"),
                minimum_steps=r.get("minimum_steps_for_holm"),
                resolution=r.get("resolution_sufficient"),
            )
        comparisons.append(record | dict(steps=steps, seconds=perf_counter() - start))
    return out | {"comparisons": comparisons}


def interval(events):
    n = len(events)
    hits = sum(events)
    ci = binomtest(hits, n).proportion_ci()
    return dict(hits=hits, n=n, ci=[ci.low, ci.high])


def summarize(rows):
    summary = {}
    for name, *_ in CELLS:
        panel = [r for r in rows if r["cell"] == name]
        assert [r["replicate"] for r in panel] == list(range(50))
        if name.startswith("attribution"):
            a = [r["attribution"] for r in panel]
            summary[name] = dict(
                n=50,
                display_supported=sum(r["native_display"] is not None for r in a),
                medians={
                    k: float(np.median([r[k] for r in a if r[k] is not None]))
                    for k in a[0]
                    if k not in ("child_sizes", "annotated_display")
                    and any(r[k] is not None for r in a)
                },
            )
            continue
        if name == "fixed_beta":
            stats = dict(beta=interval([r["beta_p"] <= 0.05 for r in panel]))
        else:
            stats = {}
        for j, steps in enumerate((99, 999)):
            runs = [r["comparisons"][j] for r in panel]
            st = dict(
                rejected=interval([r["decision"] == "reject" for r in runs]),
                unresolved=sum(r["decision"] == "unresolved" for r in runs),
                seconds_median=float(np.median([r["seconds"] for r in runs])),
            )
            if name == "fixed_beta":
                st["median_absolute_p_difference"] = float(
                    np.median([abs(a["beta_p"] - b["p_value"]) for a, b in zip(panel, runs)])
                )
            else:
                for null, endpoint in ((True, "false_family"), (False, "true_family")):
                    st[endpoint] = interval(
                        [
                            any(n["null"] == null and n["decision"] == "reject" for n in r["nodes"])
                            for r in runs
                        ]
                    )
                st.update(
                    resolution_failures=sum(not r["resolution"] for r in runs),
                    median_moves=float(np.median([n["moves"] for r in runs for n in r["nodes"]])),
                    unresolved_nodes=sum(
                        n["decision"] == "unresolved" for r in runs for n in r["nodes"]
                    ),
                    heterogeneous_zero_contrast=sum(
                        not n["null"] and n["equal_contrast"] for r in runs for n in r["nodes"]
                    ),
                )
            stats[str(steps)] = st
        summary[name] = stats
    return summary


def main():
    files = [Path(mod.__file__) for mod in (family, gaussian, attribution, topology)] + [
        Path(__file__),
        HERE / "protocol.md",
    ]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    rows = []
    with threadpool_limits(limits=1), (HERE / "results.jsonl").open("x") as stream:
        for i, (name, *_) in enumerate(CELLS):
            for j in range(50):
                result = trial(i, j)
                rows.append(result)
                stream.write(json.dumps(result, allow_nan=False) + "\n")
                stream.flush()
            print(name, flush=True)
    assert all(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    (HERE / "summary.json").write_text(json.dumps(summarize(rows), indent=2) + "\n")
    (HERE / "manifest.json").write_text(
        json.dumps(dict(sha256=hashes, attempts=len(rows)), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
