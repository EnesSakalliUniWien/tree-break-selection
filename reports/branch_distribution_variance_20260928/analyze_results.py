"""Verify matched inputs/settings and summarize the branch-distribution benchmark."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from benchmarks.validation.contracts.report_contract import wilson_interval
from scipy.stats import binomtest
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

KEY = ['profile', 'signal', 'replicate']
VARIANT = ['pipeline', 'source', 'context', 'policy']


def node_effects(items):
    frame = pd.DataFrame(items)
    groups = {k: g.set_index('node') for k, g in frame.groupby(VARIANT)}
    out = []
    for key, variant in groups.items():
        pipe, source, context, policy = key
        baselines = []
        if context != 'leaf_only':
            baselines.append(('internal_rows', (pipe, source, 'leaf_only', policy)))
        if context == 'branch_length_state':
            baselines.append(('branch_weighting', (pipe, source, 'empirical_barycenter', policy)))
        if policy == 'normalized_branch_length':
            baselines.append(('edge_variance', (pipe, source, context, 'none')))
        if source == 'diffusion_nnls':
            baselines.append(('nnls', (pipe, 'diffusion', context, policy)))
        for effect, baseline_key in baselines:
            if baseline_key not in groups:
                continue
            baseline = groups[baseline_key]
            assert set(variant.index) == set(baseline.index)
            baseline = baseline.loc[variant.index]
            row = {**dict(zip(VARIANT, key)), **variant.iloc[0][KEY].to_dict(),
                   'effect': effect, 'nodes': len(variant)}
            for col in ['p_value', 'raw_p_value', 'statistic', 'raw_statistic', 'rank', 'c_hat', 'null_weight']:
                a = pd.to_numeric(baseline.get(col), errors='coerce')
                b = pd.to_numeric(variant.get(col), errors='coerce')
                finite = np.isfinite(a) & np.isfinite(b)
                row[col + '_compared'] = int(np.sum(finite))
                row[col + '_changed'] = int(np.sum(finite & ~np.isclose(a, b, rtol=1e-10, atol=1e-12)))
                row[col + '_availability_changed'] = int(np.sum(np.isfinite(a) != np.isfinite(b)))
            for col in ['projection_sha256', 'spectrum_sha256']:
                row[col + '_changed'] = int((variant[col] != baseline[col]).sum())
            assert np.array_equal(variant.independent_rows, baseline.independent_rows)
            out.append(row)
    return out


def analyze(directory):
    runs = pd.read_csv(directory / 'runs.csv')
    manifest = json.loads((directory / 'manifest.json').read_text())
    assert manifest['complete'] and manifest['source_hashes_unchanged']
    assert not runs.duplicated(KEY + VARIANT + ['alpha']).any()
    assert len(runs) == manifest['paired_datasets'] * 54
    assert runs.groupby(KEY).size().eq(54).all()
    ok = runs[runs.status.eq('ok')]
    for _, g in ok.groupby(KEY):
        assert g.discovery_sha256.nunique() == g.inference_sha256.nunique() == 1
        assert g.topology_sha256.nunique() == 1
        assert g.groupby('source').branches_sha256.nunique().eq(1).all()
    summary = []
    for key, group in runs.groupby(['profile', 'signal'] + VARIANT + ['alpha']):
        good = group[group.status.eq('ok')]
        row = dict(zip(['profile', 'signal'] + VARIANT + ['alpha'], key))
        row.update(attempts=len(group), completed=len(good), unavailable=len(group)-len(good),
                   mean_ari=good.ari.mean(), mean_clusters=good.clusters.mean(),
                   any_false_split_count=int(good.any_false_split.sum()))
        if not key[1] and len(good):
            count = int((good.clusters > 1).sum())
            low, high = wilson_interval(count, len(good))
            row.update(false_partitions=count, false_partition_rate=count / len(good),
                       wilson95_low=low, wilson95_high=high)
        summary.append(row)
    pd.DataFrame(summary).to_csv(directory / 'summary.csv', index=False)
    run_effect_rows = []
    for key, group in runs.groupby(VARIANT):
        pipe, source, context, policy = key
        baselines = []
        if context != 'leaf_only':
            baselines.append(('internal_rows', (pipe, source, 'leaf_only', policy)))
        if context == 'branch_length_state':
            baselines.append(('branch_weighting', (pipe, source, 'empirical_barycenter', policy)))
        if policy == 'normalized_branch_length':
            baselines.append(('edge_variance', (pipe, source, context, 'none')))
        if source == 'diffusion_nnls':
            baselines.append(('nnls', (pipe, 'diffusion', context, policy)))
        for effect, baseline_key in baselines:
            baseline = runs[(runs[VARIANT] == pd.Series(dict(zip(VARIANT, baseline_key)))).all(axis=1)]
            paired = group.merge(baseline, on=KEY + ['alpha'], suffixes=('', '_baseline'), validate='one_to_one')
            assert len(paired) == len(group)
            for record in paired.to_dict('records'):
                row = {k: record[k] for k in KEY + VARIANT + ['alpha']}
                row.update(effect=effect, status=record['status'], baseline_status=record['status_baseline'])
                if record['status'] == record['status_baseline'] == 'ok':
                    row.update(delta_ari=record['ari']-record['ari_baseline'],
                        delta_clusters=record['clusters']-record['clusters_baseline'],
                        partition_changed=adjusted_rand_score(json.loads(record['labels']), json.loads(record['labels_baseline'])) < 1-1e-12)
                run_effect_rows.append(row)
    pd.DataFrame(run_effect_rows).to_csv(directory / 'run_effects.csv', index=False)
    roots = pd.read_csv(directory / 'root_nodes.csv')
    tails = []
    for key, group in roots[~roots.signal].groupby(['profile'] + VARIANT):
        for field in (['p_value', 'incompatible_scaled_p'] if key[1] == 'independent' else ['p_value', 'raw_p_value']):
            p = group[field].dropna()
            for cutoff in [.001, .005, .01, .05]:
                count = int((p <= cutoff).sum())
                low, high = wilson_interval(count, len(p)) if len(p) else (np.nan, np.nan)
                tails.append({**dict(zip(['profile'] + VARIANT, key)), 'field': field,
                    'cutoff': cutoff, 'n': len(p), 'unavailable': len(group)-len(p),
                    'rejections': count, 'rate': count/len(p) if len(p) else np.nan,
                    'wilson95_low': low, 'wilson95_high': high,
                    'binomial_excess_p': binomtest(count, len(p), cutoff, alternative='greater').pvalue if len(p) else np.nan})
    tails = pd.DataFrame(tails)
    mask = tails.pipeline.eq('independent') & tails.field.eq('p_value')
    assert mask.sum() == 72
    tails.loc[mask, 'holm_p_72'] = multipletests(tails.loc[mask, 'binomial_excess_p'], method='holm')[1]
    tails.to_csv(directory / 'root_tails.csv', index=False)
    # Stream dataset groups: do not keep the large annotation payload in memory.
    effects, items, previous, count = [], [], None, 0
    max_scale_error = 0.0
    with gzip.open(directory / 'nodes.jsonl.gz', 'rt') as f:
        for line in f:
            item = json.loads(line)
            dataset = tuple(item[k] for k in KEY)
            if previous is not None and dataset != previous:
                effects.extend(node_effects(items))
                items = []
            if item['pipeline'] == 'independent':
                max_scale_error = max(max_scale_error, abs(item['p_value'] - item['corrected_scaled_p']))
                assert item['incompatible_scaled_p'] >= item['p_value'] - 1e-12
            items.append(item)
            previous = dataset
            count += 1
    effects.extend(node_effects(items))
    effects = pd.DataFrame(effects)
    edge_effect = effects[effects.effect.eq('edge_variance')]
    assert edge_effect[['raw_statistic_changed', 'raw_p_value_changed', 'rank_changed',
                        'projection_sha256_changed', 'spectrum_sha256_changed']].eq(0).all().all()
    effects.to_csv(directory / 'node_effects_by_dataset.csv', index=False)
    effects.groupby(VARIANT + ['effect']).sum(numeric_only=True).drop(columns=['signal', 'replicate']).to_csv(directory / 'node_effects.csv')
    contexts, warning_counts = Counter(), Counter()
    for line in (directory / 'contexts.jsonl').open():
        item = json.loads(line)
        contexts[item['kind']] += 1
        warning_counts.update(item.get('warnings', []))
        if item['kind'] == 'independent_context':
            assert item['frozen_unchanged']
    verify = dict(run_count=len(runs), dataset_count=manifest['paired_datasets'],
        status_counts=runs.status.value_counts().to_dict(), node_count=count,
        contexts=dict(contexts), warning_counts=dict(warning_counts),
        corrected_scale_max_absolute_p_error=max_scale_error,
        shared_data_topology_and_per_source_lengths=True, independent_row_counts_preserved=True,
        edge_variance_preserves_raw_sibling_statistics_and_spectra=True,
        analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if manifest['arguments']['phase'] == 'evaluation':
        old = pd.read_csv(Path(__file__).resolve().parents[1] / 'independent_diffusion_nnls_20260928/evaluation/runs.csv')
        checks = []
        for source in ['diffusion', 'diffusion_nnls']:
            for pipeline in ['current', 'independent']:
                chosen = runs[runs.source.eq(source) & runs.pipeline.eq(pipeline) & runs.context.eq('leaf_only')]
                if pipeline == 'current':
                    chosen = chosen[chosen.alpha.eq(.01) & chosen.policy.eq('none' if source == 'diffusion' else 'normalized_branch_length')]
                    prior = old[old.geometry.eq('tbs_diffusion_adaptive' + ('_nnls' if source.endswith('nnls') else ''))]
                else:
                    prior = old[old.geometry.eq(source) & old.reference.eq('conditional')]
                pair = chosen.merge(prior, on=KEY + ['alpha'], suffixes=('_new', '_old'), validate='one_to_one')
                assert len(pair) == len(chosen)
                for col in ['status', 'discovery_sha256', 'inference_sha256']:
                    assert pair[col+'_new'].equals(pair[col+'_old'])
                for col in ['clusters', 'ari'] + (['root_p'] if pipeline == 'independent' else []):
                    np.testing.assert_allclose(pair[col+'_new'], pair[col+'_old'], atol=1e-12, rtol=1e-10)
                checks.append(dict(source=source, pipeline=pipeline, matched_rows=len(pair)))
        verify['previous_benchmark_baselines_reproduced'] = checks
    (directory / 'verification.json').write_text(json.dumps(verify, indent=2) + '\n')
    print(json.dumps(verify, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    analyze(parser.parse_args().directory)
