"""Do different cluster pairs separate in different subspaces? Held-out check.

For every pair of true clusters (g, h): split their rows in half at random,
select the pair's separating features on half A (Welch t, BH q < 0.01), then on
half B compare how well a 2-cluster average-linkage cut separates g from h using
(1) all features, (2) the union of all pairs' features selected on half A, and
(3) only this pair's features. Using held-out rows avoids rewarding features
that were picked because they separate the very rows they are scored on.
"""

import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.metrics import adjusted_rand_score

from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.util.case_inputs import prepare_case_inputs

warnings.filterwarnings("ignore")


def select(x, a, b):
    p = ttest_ind(x[a], x[b], axis=0, equal_var=False).pvalue
    p = np.where(np.isfinite(p), p, 1.0)
    return false_discovery_control(p, method="bh") < 0.01


def two_cut(x, rows, truth, binary):
    if x.shape[1] == 0:
        return np.nan
    v = x[rows]
    if binary:
        d = pdist(v, "hamming")
    else:
        var = v.var(0, ddof=1)
        var[var <= 0] = 1.0
        d = pdist(v, "seuclidean", V=var)
    return adjusted_rand_score(truth, fcluster(linkage(d, "average"), 2, "maxclust"))


HERE = Path(__file__).parent
summary = pd.read_csv(HERE / "failing_case_subspace_summary.csv")
names = summary[summary.diagnosis.str.contains("partly|misses|no tree")].case_id
cases = {c["name"]: c for c in get_test_cases_by_suite("full")}
rng = np.random.default_rng(0)
out = []
for name in names:
    prep = prepare_case_inputs(cases[name], ["tbs"])
    x = prep.data.to_numpy(float)
    y = np.asarray(prep.labels)
    binary = summary.set_index("case_id").loc[name, "representation"] != "continuous"
    groups = {g: rng.permutation(np.flatnonzero(y == g)) for g in np.unique(y)}
    half_a = {g: r[: len(r) // 2] for g, r in groups.items()}
    half_b = {g: r[len(r) // 2 :] for g, r in groups.items()}
    pairs = list(combinations(groups, 2))
    if len(pairs) > 120:
        pairs = [pairs[i] for i in rng.choice(len(pairs), 120, replace=False)]
    masks = {pr: select(x, half_a[pr[0]], half_a[pr[1]]) for pr in pairs}
    union = np.any(list(masks.values()), axis=0)
    for (g, h), m in masks.items():
        rows = np.concatenate([half_b[g], half_b[h]])
        truth = np.r_[np.zeros(len(half_b[g])), np.ones(len(half_b[h]))]
        out.append(
            {
                "case_id": name,
                "pair": f"{g}-{h}",
                "pair_features": int(m.sum()),
                "full": two_cut(x, rows, truth, binary),
                "union": two_cut(x[:, union], rows, truth, binary),
                "pair_only": two_cut(x[:, m], rows, truth, binary),
            }
        )
    print(name, flush=True)
res = pd.DataFrame(out)
res.to_csv(HERE / "pair_subspace_heldout.csv", index=False)
agg = res.groupby("case_id").agg(
    pairs=("pair", "size"),
    pair_features=("pair_features", "median"),
    full=("full", "mean"),
    union=("union", "mean"),
    pair_only=("pair_only", "mean"),
)
print(agg.round(3).sort_values("pair_only", ascending=False).to_string())
