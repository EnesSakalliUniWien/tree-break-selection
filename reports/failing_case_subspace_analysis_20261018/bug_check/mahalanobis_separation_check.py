"""Does mahalanobis_time keep a large mean separation between two Gaussian clusters?"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

from tree_break_selection.tree.continuous_distance import continuous_time_distance_condensed
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

rng = np.random.default_rng(3)
n = 120
rows = []
for p in (40, 80, 110, 119, 150, 300):
    for shift in (1.0, 3.0):
        y = np.repeat([0, 1], n // 2)
        x = rng.normal(size=(n, p))
        x[y == 1] += shift  # every feature shifted by `shift` standard deviations
        space = continuous_feature_space_from_columns([f"F{i}" for i in range(p)])
        d = continuous_time_distance_condensed(x, space)
        rows.append(
            {
                "p": p,
                "shift_sd": shift,
                "mahalanobis_time_cv": d.std() / d.mean(),
                "mahalanobis_time_ari": adjusted_rand_score(
                    y, fcluster(linkage(d, "average"), 2, "maxclust")
                ),
                "euclidean_ari": adjusted_rand_score(
                    y, fcluster(linkage(pdist(x), "average"), 2, "maxclust")
                ),
            }
        )
table = pd.DataFrame(rows)
table.to_csv(Path(__file__).parent / "mahalanobis_separation_check.csv", index=False)
print(table.round(4).to_string(index=False))
