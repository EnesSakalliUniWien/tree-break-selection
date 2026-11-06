"""Reproduce complete diffusion-merge-sequence reentry; no calibrated integration.

uv run --no-sync python raw/assets/diffusion-selection-event-20261102/reentry_probe.py
The finite scan and transition brackets do not exhaust the radial selection set.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from pydiffmap import kernel

from tree_break_selection.space_separation.diffusion import feature_family_diffusion_geometry
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

OUTPUT = Path(__file__).resolve().parent
ORIGINAL_KERNEL = kernel.Kernel


def replay(values):
    captured = []

    def factory(*args, **kwargs):
        fitted = ORIGINAL_KERNEL(*args, **kwargs)
        captured.append(fitted)
        return fitted

    data = pd.DataFrame(values, index=[f'S{i}' for i in range(len(values))])
    space = continuous_feature_space_from_columns(data.columns)
    with patch('pydiffmap.kernel.Kernel', factory):
        geometry = feature_family_diffusion_geometry(data, space)
    fitted = captured[0]
    neighbors = fitted.neigh.kneighbors(values, return_distance=False).tolist()
    built = build_tree(data, geometry.distance_condensed, linkage_method='average')
    signature = [sorted(map(int, row[:2])) for row in built.linkage_matrix]
    return dict(signature=signature, epsilon=geometry.metadata['epsilon'],
                components=geometry.metadata['diffusion_components'],
                kde_dimension=float(fitted.kde.d), kde_epsilon=float(fitted.kde.epsilon_fitted),
                neighbors=neighbors, metadata=geometry.metadata)


def main():
    seed, n, p = 20261105, 9, 1
    x = np.random.default_rng(seed).normal(size=(n, p))
    observed = replay(x)
    groups = {i: [i] for i in range(n)}
    for step, (a, b) in enumerate(observed['signature']):
        groups[n + step] = groups[a] + groups[b]
    left, right = observed['signature'][-1]
    eta = np.zeros(n)
    eta[groups[left]], eta[groups[right]] = 1 / len(groups[left]), -1 / len(groups[right])
    e = eta / np.linalg.norm(eta)
    z = e @ x
    radius = np.linalg.norm(z)
    u = z / radius
    offset = x - np.outer(e, z)
    np.testing.assert_allclose(offset + np.outer(e, radius * u), x, atol=1e-14)
    np.testing.assert_allclose(e @ offset, np.zeros(p), atol=1e-14)

    def evaluate(fraction):
        state = replay(offset + np.outer(e, fraction * radius * u))
        return dict(fraction=float(fraction), same=state['signature'] == observed['signature'], **state)

    grid = [evaluate(f) for f in np.unique(np.r_[np.linspace(0, 4, 121), np.geomspace(4.1, 100, 20)])]
    assert evaluate(1.)['same'] and not evaluate(1.1)['same'] and evaluate(1.2)['same']
    brackets = []
    for first, second in zip(grid[:-1], grid[1:]):
        if first['same'] == second['same']:
            continue
        lo, hi = first['fraction'], second['fraction']
        left_state = first['same']
        for _ in range(30):
            mid = (lo + hi) / 2
            if evaluate(mid)['same'] == left_state:
                lo = mid
            else:
                hi = mid
        brackets.append(dict(lower=evaluate(lo), upper=evaluate(hi)))
    source_paths = [Path('tree_break_selection/space_separation/diffusion.py'), Path(kernel.__file__)]
    result = dict(seed=seed, n=n, p=p, sigma=1, radius=float(radius),
                  root_groups=[groups[left], groups[right]], x=x.tolist(), e=e.tolist(), u=u.tolist(),
                  offset=offset.tolist(), observed=observed, grid=grid, transition_brackets=brackets,
                  versions={name: importlib.metadata.version(name) for name in ['numpy','scipy','scikit-learn','pydiffmap']},
                  source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
                  limitation='Finite replay and local transition brackets, not exhaustive intervals or a selected p-value.')
    (OUTPUT / 'reentry.json').write_text(json.dumps(result, indent=2) + '\n')
    print('radius', radius, 'root_groups', result['root_groups'])
    for bracket in brackets:
        print('transition', bracket['lower']['fraction'], bracket['upper']['fraction'],
              bracket['lower']['same'], bracket['upper']['same'])
    print('verified matching r/r_obs=1 and1.2 separated by nonmatching1.1; output', OUTPUT / 'reentry.json')


if __name__ == '__main__':
    main()
