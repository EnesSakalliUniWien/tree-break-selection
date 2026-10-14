"""Matched independent Gaussian tree/rank selection replacement benchmark."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import multiprocessing
import platform
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.decomposition import labels_and_report_from_decomposition
from benchmarks.validation.contracts.report_contract import wilson_interval
from sklearn.metrics import adjusted_rand_score
from tree_break_selection.hierarchy_analysis.statistics.independent_gaussian import (
    _annotate_tests,
    fit_independent_gaussian_tree,
    test_independent_gaussian_tree,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

ROOT = Path(__file__).resolve().parents[2]
METHODS = ('average','complete','weighted','single','centroid','median','ward','neighbor_joining')
PROFILES = ('identity2','ar1_8','identity40')
SEED = 2026092811


def digest(x):
    return hashlib.sha256(np.ascontiguousarray(x, dtype=float).tobytes()).hexdigest()


def source_hashes():
    paths = sorted((ROOT/'tree_break_selection').rglob('*.py'))
    paths += sorted((ROOT/'benchmarks/shared').rglob('*.py'))
    paths += [Path(__file__), Path(__file__).with_name('protocol.md'), ROOT/'uv.lock']
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run_dataset(job):
    profile, signal, replicate = job
    dimension = {'identity2':2,'ar1_8':8,'identity40':40}[profile]
    sigma = .6**np.abs(np.subtract.outer(np.arange(dimension),np.arange(dimension))) if profile=='ar1_8' else np.eye(dimension)
    seed = np.random.SeedSequence([SEED, PROFILES.index(profile), int(signal), replicate])
    discovery_seed, inference_seed = seed.spawn(2)
    truth = np.zeros(32, dtype=int)
    if signal:
        truth[16:]=1
    means = np.zeros((32,dimension))
    if signal:
        means[:,:2] = np.where(truth[:,None]==0,-1.,1.)
    labels = [f'S{i}' for i in range(32)]
    columns = [f'X{i}' for i in range(dimension)]
    discovery = pd.DataFrame(means + np.random.default_rng(discovery_seed).multivariate_normal(np.zeros(dimension),sigma,32),index=labels,columns=columns)
    inference = pd.DataFrame(means + np.random.default_rng(inference_seed).multivariate_normal(np.zeros(dimension),sigma,32),index=labels,columns=columns)
    key = dict(profile=profile,signal=signal,replicate=replicate,discovery_sha256=digest(discovery),inference_sha256=digest(inference))
    runs, nodes = [], []
    for method in METHODS:
        start = perf_counter()
        alphas = (.005,.01,.05) if method in ('average','neighbor_joining') else (.01,)
        try:
            model = fit_independent_gaussian_tree(discovery, known_covariance=sigma,
                observation_ids=[f'd{i}' for i in range(32)], energy_fraction=.9,
                tree_builder='neighbor_joining' if method=='neighbor_joining' else 'linkage',
                linkage_method='average' if method=='neighbor_joining' else method)
            result = test_independent_gaussian_tree(model,inference,observation_ids=[f't{i}' for i in range(32)])
            root = model.tree.root()
            true_null = {x.parent: int(truth[list(x.left_rows)].sum())*len(x.right_rows)==int(truth[list(x.right_rows)].sum())*len(x.left_rows) for x in model.nodes}
            for row in result.tests.to_dict('records'):
                nodes.append({**key,'method':method,**row,'is_root':row['parent']==root,'true_null':true_null[row['parent']]})
            for variant in ('independent_conditional','independent_nominal'):
                tests = result.tests.copy()
                if variant=='independent_nominal':
                    tests['p_value']=tests.nominal_p_value
                    tests['reference_law']='nominal_chi_square_after_energy_selection'
                for alpha in alphas:
                    if variant=='independent_conditional' and alpha==.01:
                        decomposition=result.decomposition
                    else:
                        _, decomposition = _annotate_tests(model,tests,alpha=alpha,edge_alpha=.001,passthrough=True)
                    predicted,_=labels_and_report_from_decomposition(decomposition,labels)
                    splits=[r['node_id'] for r in decomposition['traversal_trace'] if r['decision']=='split']
                    false_splits=sum(true_null[node] for node in splits)
                    runs.append({**key,'method':method,'variant':variant,'alpha':alpha,'status':'ok','reason':'',
                        'clusters':decomposition['num_clusters'],'ari':adjusted_rand_score(truth,predicted),
                        'any_false_split':bool(false_splits),'n_false_splits':false_splits,'n_splits':len(splits),
                        'root_p':float(tests.loc[root,'p_value']),'root_true_null':true_null[root],
                        'root_rank':int(tests.loc[root,'projection_dimension']),
                        'root_candidate_dimension':int(tests.loc[root,'candidate_dimension']),
                        'elapsed_seconds':perf_counter()-start})
        except Exception as exc:
            status = 'unsupported_tree' if isinstance(exc, ValueError) and 'merge distances must be nondecreasing' in str(exc) else 'error'
            for variant in ('independent_conditional','independent_nominal'):
                for alpha in alphas:
                    runs.append({**key,'method':method,'variant':variant,'alpha':alpha,'status':status,'reason':repr(exc),'elapsed_seconds':perf_counter()-start})
        old_name={'average':'tbs','complete':'tbs_complete','neighbor_joining':'tbs_neighbor_joining'}.get(method)
        if old_name:
            params=dict(METHOD_SPECS[old_name].param_grid[0])
            params['tree_distance_metric']='standardized_euclidean'
            try:
                old=run_clustering_result(discovery,old_name,params,seed=int(seed.generate_state(1)[0]),
                    significance_level=.01,edge_alpha=.001,feature_space=continuous_feature_space_from_columns(columns))
                row={**key,'method':method,'variant':'current_single_matrix','alpha':.01,'status':old.status,
                     'reason':str(old.unsupported_reason or old.skip_reason or '')}
                if old.status=='ok':
                    row.update(clusters=old.found_clusters,ari=adjusted_rand_score(truth,old.labels))
                runs.append(row)
            except Exception as exc:
                runs.append({**key,'method':method,'variant':'current_single_matrix','alpha':.01,'status':'error','reason':repr(exc)})
    return runs,nodes


def summarize(runs):
    rows=[]
    for key,g in runs.groupby(['profile','signal','method','variant','alpha']):
        ok=g[g.status=='ok']
        row=dict(zip(['profile','signal','method','variant','alpha'],key))
        row.update(attempts=len(g),completed=len(ok),unavailable=len(g)-len(ok),mean_ari=ok.ari.mean(),mean_clusters=ok.clusters.mean())
        if not key[1] and len(ok):
            count=int((ok.clusters>1).sum())
            lo,hi=wilson_interval(count,len(ok))
            row.update(false_partition_count=count,false_partition_rate=count/len(ok),wilson95_low=lo,wilson95_high=hi)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--null-replicates',type=int,default=250)
    parser.add_argument('--signal-replicates',type=int,default=50)
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    hashes=source_hashes()
    manifest=dict(arguments=vars(args),seed=SEED,source_sha256=hashes,python=sys.version,platform=platform.platform(),complete=False)
    (args.output_dir/'manifest.json').write_text(json.dumps(manifest,default=str,indent=2)+'\n')
    jobs=[(profile,signal,rep) for profile in PROFILES for signal in (False,True)
          for rep in range(args.signal_replicates if signal else args.null_replicates)]
    rows=[]
    start=perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool, gzip.open(args.output_dir/'nodes.csv.gz','wt') as node_file:
        for completed,(runs,nodes) in enumerate(pool.map(run_dataset,jobs),1):
            rows.extend(runs)
            if nodes:
                pd.DataFrame(nodes).to_csv(node_file,index=False,header=node_file.tell()==0)
            if completed%50==0 or completed==len(jobs):
                pd.DataFrame(rows).to_csv(args.output_dir/'runs.csv',index=False)
                print(f'{completed}/{len(jobs)} paired datasets; {len(rows)} runs; {perf_counter()-start:.1f}s',flush=True)
    runs=pd.DataFrame(rows)
    summary=summarize(runs)
    roots=runs[runs.alpha.eq(.01)&runs.root_true_null.eq(True)&runs.status.eq('ok')]
    tails=[]
    for key,g in roots.groupby(['profile','signal','method','variant']):
        for alpha in (.001,.005,.01,.05):
            count=int((g.root_p<=alpha).sum())
            lo,hi=wilson_interval(count,len(g))
            tails.append({**dict(zip(['profile','signal','method','variant'],key)),'cutoff':alpha,'n':len(g),'rejections':count,'rate':count/len(g),'wilson95_low':lo,'wilson95_high':hi})
    assert hashes==source_hashes(), 'Sources changed during benchmark.'
    write_diagnostic_bundle(output_dir=args.output_dir,tables={'runs':runs,'summary':summary,'root_tails':pd.DataFrame(tails)},
        filenames={'runs':'runs.csv','summary':'summary.csv','root_tails':'root_tails.csv'},
        manifest={**manifest,'complete':True,'paired_datasets':len(jobs),'elapsed_seconds':perf_counter()-start,'source_hashes_unchanged':True,'status_counts':runs.status.value_counts().to_dict()})
    print(runs.groupby(['variant','status']).size().to_string(),flush=True)


if __name__=='__main__':
    main()
