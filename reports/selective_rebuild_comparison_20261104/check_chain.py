"""Fixed-law sanity checks of the radial chain, independent of diffusion."""

import json
from pathlib import Path

import numpy as np
from run_comparison import chain_rank
from scipy.stats import chi

rng = np.random.default_rng(2026110442)
results = {}
for name in ["all_positive", "two_disconnected_intervals"]:
    if name == "all_positive":
        observations = np.sqrt(rng.chisquare(3, 4000))

        def event(radius):
            return radius > 0

    else:
        intervals = np.array([[3.0, 3.5], [5.0, 6.0]])
        mass = chi.cdf(intervals[:, 1], 3) - chi.cdf(intervals[:, 0], 3)
        selected = rng.choice(2, size=4000, p=mass / mass.sum())
        lower = chi.cdf(intervals[selected, 0], 3)
        upper = chi.cdf(intervals[selected, 1], 3)
        observations = chi.ppf(lower + rng.random(4000) * (upper - lower), 3)

        def event(radius):
            return (3 <= radius <= 3.5) or (5 <= radius <= 6)

    p_values, moves = [], []
    for index, observed in enumerate(observations):
        result, _ = chain_rank(observed, 3, event, steps=99, seed=2026110500 + index)
        p_values.append(result["p_value"])
        moves.append(result["moves"])
    p_values = np.asarray(p_values)
    count = int(np.count_nonzero(p_values <= 0.05))
    assert count < 200 + 6 * np.sqrt(4000 * 0.05 * 0.95)
    results[name] = dict(
        draws=4000,
        steps=99,
        rejections_005=count,
        rejections_001=int(np.count_nonzero(p_values <= 0.01)),
        median_moves=float(np.median(moves)),
    )
(Path(__file__).parent / "chain_null_control.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
