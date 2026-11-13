"""Use the same conditional trajectory for each declared score control."""

import json
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from run import FACTOR, HERE, attribution, gaussian, topology

KEYS = ("native_raw", "unit_branch_raw", "oracle_covariance_raw", "gram_rank_raw", "native_display")


def run(i, steps):
    data, seed = np.random.SeedSequence([7811260911, i, 0]).spawn(2)
    seed = int(seed.generate_state(1, dtype=np.uint64)[0])
    x = np.random.default_rng(data).normal(size=(16, 3)) @ FACTOR
    if i:
        x[8:, 0] += 6
    rows, event = topology.ward_subtree_selection(x, keep_split=True)
    scored = []

    def score(matrix):
        evidence = attribution.native_sibling_attribution(
            matrix, rows, oracle_covariance=FACTOR.T @ FACTOR, gram_rank=1
        )
        scored.append({k: evidence[k] for k in KEYS})
        return evidence["native_raw"]

    start = perf_counter()
    r = gaussian.subtree_gaussian_chain_test(
        x, subtree_indices=rows, score=score, selection_event=event, seed=seed, steps=steps
    )
    out = dict(cell=i, steps=steps, seed=seed, seconds=perf_counter() - start, result=r)
    if r["status"] != "complete":
        return out
    anchor = r["anchor_index"]
    cursor = iter(scored)
    observed = next(cursor)
    states = [None] * (steps + 1)
    states[anchor] = observed
    for positions in (range(anchor - 1, -1, -1), range(anchor + 1, steps + 1)):
        current = observed
        for j in positions:
            if r["chain_moves"][j] is not None:
                current = next(cursor)
            states[j] = current
    assert next(cursor, None) is None
    assert [s["native_raw"] for s in states] == r["chain_scores"]
    ranks = {}
    for k in KEYS:
        if any(s[k] is None for s in states):
            ranks[k] = None
        else:
            ranks[k] = sum(s[k] >= observed[k] for s in states) / len(states)
    assert ranks["native_raw"] == r["p_value"]
    return out | dict(score_states=states, rank_p_values=ranks)


if __name__ == "__main__":
    with threadpool_limits(limits=1), (HERE / "variant.jsonl").open("x") as stream:
        for i in range(2):
            for steps in (99, 999):
                row = run(i, steps)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
                print(i, steps, row.get("rank_p_values"), row["seconds"], flush=True)
