"""Matched internal-distribution and branch-variance diagnostics; no production edits."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import multiprocessing
import platform
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.shared.runners.method_registry import FIXED_TOPOLOGY_NNLS_BRANCH_TIME_PARAMS
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from reports.independent_diffusion_nnls_20260928.run_benchmark import (
    ALPHAS,
    PROFILES,
    digest,
    frozen_hashes,
    json_hash,
)
from scipy.linalg import solve_triangular
from sklearn.metrics import adjusted_rand_score
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.statistics.independent_gaussian import (
    IndependentGaussianDiffusion,
    _annotate_tests,
    fit_independent_gaussian_tree,
    test_independent_gaussian_tree,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.energy_selected_reference import (
    energy_selected_gaussian_sf,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.tree_estimator import (
    _precompute_branch_length_internal_states,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.tree_helpers import (
    precompute_descendant_internal_nodes,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    decide_empirical_null_calibration,
    fit_empirical_null_inflation_model,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.record_collection import (
    collect_sibling_pair_records,
)
from tree_break_selection.space_separation import adaptive_diffusion_geometry
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

ROOT = Path(__file__).resolve().parents[2]
CONTEXTS = ("leaf_only", "empirical_barycenter", "branch_length_state")
SOURCES = ("diffusion", "diffusion_nnls")
POLICIES = ("none", "normalized_branch_length")


def make_data(profile, signal, replicate, namespace):
    """Exactly the previous evaluation's seed and Gaussian generator contract."""
    dimension = {"identity2": 2, "ar1_8": 8, "identity40": 40}[profile]
    sigma = (0.6 ** np.abs(np.subtract.outer(np.arange(dimension), np.arange(dimension)))
             if profile == "ar1_8" else np.eye(dimension))
    seeds = np.random.SeedSequence([namespace, PROFILES.index(profile), int(signal), replicate]).spawn(2)
    truth = np.zeros(32, dtype=int)
    means = np.zeros((32, dimension))
    if signal:
        truth[16:] = 1
        means[:, :2] = np.where(truth[:, None] == 0, -1.0, 1.0)
    matrices = [pd.DataFrame(
        means + np.random.default_rng(s).multivariate_normal(np.zeros(dimension), sigma, 32),
        index=[f"S{i}" for i in range(32)], columns=[f"X{i}" for i in range(dimension)],
    ) for s in seeds]
    return *matrices, sigma, truth


def augment_discovery(model, discovery, context):
    """Report-only adapter: orthonormal discovery basis, original candidate-rank cap."""
    if context == "leaf_only":
        return model
    values = solve_triangular(model.covariance_cholesky, discovery.to_numpy().T, lower=True).T
    tree = model.tree.copy()
    positions = {label: i for i, label in enumerate(discovery.index)}
    descendants = tree.compute_descendant_sets(use_labels=True)
    for n in tree:
        tree.nodes[n]["distribution"] = values[[positions[x] for x in descendants[n]]].mean(0)
    states = (_precompute_branch_length_internal_states(tree, feature_count=values.shape[1])
              if context == "branch_length_state" else
              {n: tree.nodes[n]["distribution"] for n in tree if tree.out_degree(n)})
    internal = precompute_descendant_internal_nodes(tree)
    nodes = []
    for node in model.nodes:
        rows = values[list(node.left_rows + node.right_rows)]
        extra = [states[n] for n in internal[node.parent]]
        if not extra:
            nodes.append(node)
            continue
        eigen = eigendecompose_covariance(np.vstack([rows, *extra]), compute_eigenvectors=True)
        projection, eigenvalues = build_pca_projection(
            eigen, projection_dimension=len(node.projection), n_features_total=values.shape[1],
        )
        assert len(projection) <= len(rows) - 1
        np.testing.assert_allclose(projection @ projection.T, np.eye(len(projection)), atol=1e-8)
        projection.setflags(write=False)
        eigenvalues.setflags(write=False)
        nodes.append(replace(node, projection=projection, eigenvalues=eigenvalues))
    return replace(model, nodes=tuple(nodes))


def partition_metrics(tree, decomposition, truth, labels):
    predicted, _ = labels_and_report_from_decomposition(decomposition, labels)
    descendants = tree.compute_descendant_sets(use_labels=True)
    positions = {x: i for i, x in enumerate(labels)}
    false = 0
    split_nodes = [r['node_id'] for r in decomposition['traversal_trace'] if r['decision'] == 'split']
    for parent in split_nodes:
        a, b = tree.successors(parent)
        li, ri = [positions[x] for x in descendants[a]], [positions[x] for x in descendants[b]]
        false += int(truth[li].sum()) * len(ri) == int(truth[ri].sum()) * len(li)
    return dict(clusters=decomposition['num_clusters'], ari=adjusted_rand_score(truth, predicted),
                labels=json.dumps(np.asarray(predicted).tolist()), n_splits=len(split_nodes),
                n_false_splits=int(false), any_false_split=bool(false))


def tree_hashes(tree):
    return dict(topology_sha256=json_hash(sorted((str(u), str(v)) for u, v in tree.edges)),
                branches_sha256=json_hash([(str(u), str(v), float(a['branch_length']))
                                           for u, v, a in tree.edges(data=True)]))


def run_dataset(job):
    profile, signal, replicate, namespace = job
    discovery, inference, sigma, truth = make_data(*job)
    key = dict(profile=profile, signal=signal, replicate=replicate,
               discovery_sha256=digest(discovery), inference_sha256=digest(inference))
    runs, nodes, metadata = [], [], []
    feature_space = continuous_feature_space_from_columns(discovery.columns)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        geometry = adaptive_diffusion_geometry(discovery, k_neighbors=10, diffusion_time=3,
                    n_components=30, metric='euclidean', bandwidth_type='-1/(d+2)', epsilon='median')
    metadata.append({**key, 'kind': 'geometry', 'metadata': geometry.metadata,
                     'warnings': [str(w.message) for w in caught]})
    for source in SOURCES:
        nnls = dict(FIXED_TOPOLOGY_NNLS_BRANCH_TIME_PARAMS) if source == 'diffusion_nnls' else {}
        nnls.pop('edge_branch_length_variance_policy', None)
        for context in CONTEXTS:
            for policy in POLICIES:
                base = {**key, 'pipeline': 'current', 'source': source, 'context': context, 'policy': policy}
                for alpha in ALPHAS:
                    try:
                        with warnings.catch_warnings(record=True) as caught:
                            warnings.simplefilter('always')
                            result = run_tbs_on_distance(
                                discovery, geometry.distance_condensed, alpha, edge_alpha=.001,
                                feature_space=feature_space, tree_linkage_method='average',
                                branch_length_data_df=discovery, **nnls,
                                spectral_include_internal_barycenters=context != 'leaf_only',
                                spectral_internal_distribution_mode=(context if context != 'leaf_only' else 'empirical_barycenter'),
                                edge_branch_length_variance_policy=policy,
                                allow_linkage_ultrametric_branch_time=source == 'diffusion',
                                passthrough=True, trace_level='full',
                            )
                        bundle = result.extra['gate_bundle']
                        cfg = asdict(bundle.metadata.config)
                        assert cfg['spectral_include_internal_barycenters'] == (context != 'leaf_only')
                        assert cfg['spectral_internal_distribution_mode'] == (context if context != 'leaf_only' else 'empirical_barycenter')
                        assert cfg['edge_branch_length_variance_policy'] == policy
                        tree, annotations = result.extra['tree'], result.extra['annotations']
                        root = tree.root()
                        row = {**base, 'alpha': alpha, 'status': result.status,
                               'reason': str(result.unsupported_reason or result.skip_reason or ''),
                               'root_p': annotations.loc[root, 'Sibling_Divergence_P_Value'],
                               **tree_hashes(tree)}
                        if result.status == 'ok':
                            row.update(partition_metrics(tree, result.extra['decomposition'], truth, list(discovery.index)))
                        if alpha == .01:
                            spectral = bundle.edge_gate_result.spectral_context
                            records, _ = collect_sibling_pair_records(tree, bundle.edge_gate_result.annotated_df,
                                sibling_projection_dimensions_from_edge_comparisons=spectral.test_projection_dimensions_by_node,
                                parent_principal_component_projections=spectral.principal_component_projections_by_node,
                                parent_principal_component_eigenvalues=spectral.principal_component_eigenvalues_by_node,
                                feature_space=feature_space)
                            inflation = fit_empirical_null_inflation_model(records)
                            by_parent = {record.parent: record for record in records}
                            metadata.append({**base, 'kind': 'current_configuration', 'actual_config': cfg,
                                'warnings': [str(w.message) for w in caught],
                                'nnls': {k: v for k, v in result.extra.items() if k.startswith('branch_length_optimization_')}})
                            for n in spectral.test_projection_dimensions_by_node:
                                if tree.out_degree(n) == 0:
                                    continue
                                attrs = annotations.loc[n].to_dict()
                                record = by_parent[n]
                                decision = decide_empirical_null_calibration(inflation, record, enforce_support_thresholds=False)
                                projection = spectral.principal_component_projections_by_node[n]
                                eigenvalues = spectral.principal_component_eigenvalues_by_node[n]
                                nodes.append({**base, 'node': n, 'is_root': n == root,
                                    'p_value': attrs.get('Sibling_Divergence_P_Value'),
                                    'statistic': attrs.get('Sibling_Test_Statistic'),
                                    'raw_statistic': record.stat, 'raw_p_value': record.p_value,
                                    'c_hat': decision.c_hat, 'calibration_status': decision.status,
                                    'null_weight': record.sibling_null_weight,
                                    'edge_blocked': record.is_edge_blocked,
                                    'rank': attrs.get('Sibling_Projection_Dimension'),
                                    'candidate_rank': len(projection),
                                    'spectrum_sha256': digest(eigenvalues),
                                    'projection_sha256': digest(projection),
                                    'independent_rows': spectral.effective_independent_rows_by_node[n],
                                    'mp_rows': spectral.mp_threshold_rows_by_node[n],
                                    'internal_rows': spectral.internal_distribution_row_counts_by_node[n],
                                    'matrix_rows': spectral.spectral_matrix_row_counts_by_node[n],
                                    'annotations_json': json.dumps(attrs, default=str)})
                        runs.append(row)
                    except Exception as exc:
                        runs.append({**base, 'alpha': alpha, 'status': 'error', 'reason': repr(exc)})
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                fitted = fit_independent_gaussian_tree(discovery, known_covariance=sigma,
                    observation_ids=[f'd{i}' for i in range(32)],
                    diffusion=IndependentGaussianDiffusion(fit_nnls=source == 'diffusion_nnls'))
            metadata.append({**key, 'source': source, 'kind': 'independent_fit',
                             'warnings': [str(w.message) for w in caught],
                             'nnls': json.loads(fitted.branch_length_metadata_json)})
        except Exception as exc:
            for context in CONTEXTS:
                for alpha in ALPHAS:
                    runs.append({**key, 'pipeline': 'independent', 'source': source, 'context': context,
                        'policy': 'known_covariance', 'alpha': alpha, 'status': 'error', 'reason': repr(exc)})
            continue
        for context in CONTEXTS:
            base = {**key, 'pipeline': 'independent', 'source': source, 'context': context, 'policy': 'known_covariance'}
            try:
                model = augment_discovery(fitted, discovery, context)
                before = frozen_hashes(model)
                result = test_independent_gaussian_tree(model, inference, observation_ids=[f't{i}' for i in range(32)])
                assert before == frozen_hashes(model)
                mean_length = compute_mean_branch_length(model.tree)
                assert mean_length is not None and mean_length > 0
                for item, node in zip(result.tests.to_dict('records'), model.nodes, strict=True):
                    assert item['parent'] == node.parent
                    multiplier = 1 + sum(model.tree.edges[node.parent, c]['branch_length']
                                         for c in (node.left, node.right)) / (2 * mean_length)
                    kwargs = dict(candidate_dimension=item['candidate_dimension'],
                        projection_dimension=item['projection_dimension'],
                        previous_group_dimension=item['previous_group_dimension'], energy_fraction=.9)
                    scaled_stat = item['statistic'] / multiplier
                    naive_p = energy_selected_gaussian_sf(scaled_stat, **kwargs)
                    corrected_p = energy_selected_gaussian_sf(scaled_stat * multiplier, **kwargs)
                    assert abs(corrected_p - item['p_value']) <= 1e-10
                    assert naive_p >= item['p_value'] - 1e-12
                    nodes.append({**base, 'node': node.parent, 'is_root': node.parent == model.tree.root(),
                        'p_value': item['p_value'], 'statistic': item['statistic'],
                        'rank': item['projection_dimension'], 'candidate_rank': len(node.projection),
                        'spectrum_sha256': digest(node.eigenvalues), 'projection_sha256': digest(node.projection),
                        'independent_rows': len(node.left_rows) + len(node.right_rows),
                        'variance_multiplier': multiplier, 'incompatible_scaled_p': naive_p,
                        'corrected_scaled_p': corrected_p})
                metadata.append({**base, 'kind': 'independent_context', **before, 'frozen_unchanged': True})
                for alpha in ALPHAS:
                    _, decomposition = _annotate_tests(model, result.tests, alpha=alpha, edge_alpha=.001, passthrough=True)
                    runs.append({**base, 'alpha': alpha, 'status': 'ok', 'reason': '',
                        **tree_hashes(model.tree), 'root_p': result.tests.loc[model.tree.root(), 'p_value'],
                        **partition_metrics(model.tree, decomposition, truth, list(discovery.index))})
            except Exception as exc:
                for alpha in ALPHAS:
                    runs.append({**base, 'alpha': alpha, 'status': 'error', 'reason': repr(exc)})
    return runs, nodes, metadata


def source_hashes():
    paths = sorted((ROOT / 'tree_break_selection').rglob('*.py')) + sorted((ROOT / 'benchmarks/shared').rglob('*.py'))
    paths += [Path(__file__), Path(__file__).with_name('protocol.md'), ROOT / 'uv.lock',
              ROOT / 'reports/independent_diffusion_nnls_20260928/run_benchmark.py']
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=['pilot', 'evaluation'], required=True)
    parser.add_argument('--null-replicates', type=int, default=250)
    parser.add_argument('--signal-replicates', type=int, default=50)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    hashes = source_hashes()
    namespace = 2026092820 if args.phase == 'pilot' else 2026092813
    manifest = dict(arguments={**vars(args), 'output_dir': str(args.output_dir)}, seed_namespace=namespace,
                    source_sha256=hashes, python=sys.version, platform=platform.platform(), complete=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    jobs = [(p, s, r, namespace) for p in PROFILES for s in (False, True)
            for r in range(args.signal_replicates if s else args.null_replicates)]
    rows, root_nodes = [], []
    start = perf_counter()
    with (ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('spawn')) as pool,
          gzip.open(args.output_dir / 'nodes.jsonl.gz', 'wt') as nf,
          (args.output_dir / 'contexts.jsonl').open('w') as mf):
        for count, (runs, nodes, metadata) in enumerate(pool.map(run_dataset, jobs), 1):
            rows.extend(runs)
            root_nodes.extend(n for n in nodes if n['is_root'])
            for n in nodes:
                nf.write(json.dumps(n, default=str) + '\n')
            for m in metadata:
                mf.write(json.dumps(m, default=str) + '\n')
            if count % 25 == 0 or count == len(jobs):
                pd.DataFrame(rows).to_csv(args.output_dir / 'runs.csv', index=False)
                print(f'{count}/{len(jobs)} datasets; {len(rows)} runs; {perf_counter()-start:.1f}s', flush=True)
    assert hashes == source_hashes(), 'Sources changed during evaluation'
    runs = pd.DataFrame(rows)
    write_diagnostic_bundle(output_dir=args.output_dir,
        tables={'runs': runs, 'root_nodes': pd.DataFrame(root_nodes)},
        filenames={'runs': 'runs.csv', 'root_nodes': 'root_nodes.csv'},
        manifest={**manifest, 'complete': True, 'paired_datasets': len(jobs),
                  'elapsed_seconds': perf_counter()-start, 'source_hashes_unchanged': True,
                  'status_counts': runs.status.value_counts().to_dict()})
    print(runs.groupby(['pipeline', 'status']).size().to_string(), flush=True)


if __name__ == '__main__':
    main()
