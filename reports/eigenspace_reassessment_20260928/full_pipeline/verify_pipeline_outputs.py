#!/usr/bin/env python3
"""Independent decoding audit of a completed frozen pipeline evaluation.

Does not modify the driver, protocol, raw evaluation tables or production code.
Writes only new verification.json and deterministic_rerun_checks.csv artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import runpy
import sys
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "TBS_N_JOBS"):
    os.environ[variable] = "1"

import numpy as np
import pandas as pd
from scipy.stats import binom, chi2, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
KEYS = ["family", "signal", "replicate", "method", "sibling_alpha"]
CELL = ["family", "signal", "method", "sibling_alpha"]
DATA_KEYS = ["family", "signal", "replicate"]


def near(actual, expected, label, *, atol=2e-11, rtol=2e-10):
    if not np.allclose(np.asarray(actual, float), np.asarray(expected, float),
                       atol=atol, rtol=rtol, equal_nan=True):
        raise AssertionError(f"{label}: actual={actual!r}, expected={expected!r}")


def wilson(events, count):
    if count == 0:
        return np.nan, np.nan
    z = norm.ppf(.975)
    p = events/count
    denominator = 1+z*z/count
    center = (p+z*z/(2*count))/denominator
    margin = z*np.sqrt(p*(1-p)/count+z*z/(4*count*count))/denominator
    return max(0.,center-margin), min(1.,center+margin)


def key(row):
    return tuple(row[name] for name in KEYS)


def reconstruct_summaries(runs, summary, tails):
    checks = 0
    for row in summary.to_dict("records"):
        frame = runs.loc[np.logical_and.reduce([runs[c].eq(row[c]) for c in CELL])]
        completed = frame.loc[frame.status.eq("ok")]
        audited = completed.loc[completed.audit_error.fillna("").eq("")]
        event = int(audited.n_false_splits.gt(0).sum())
        lo, hi = wilson(event,len(audited))
        expected = {
            "attempts":len(frame), "completed":len(completed), "audited":len(audited),
            "unsupported":int(frame.status.eq("unsupported").sum()),
            "skipped":int(frame.status.eq("skip").sum()),
            "runtime_errors":int(frame.status.eq("runtime_error").sum()),
            "audit_errors":int(frame.audit_error.fillna("").ne("").sum()),
            "any_false_split_count":event,
            "any_false_split_rate":event/len(audited) if len(audited) else np.nan,
            "any_false_split_wilson_low":lo,"any_false_split_wilson_high":hi,
            "mean_split_fdp":audited.split_fdp.mean(),
            "mean_ari":completed.ari.mean() if row["signal"] else np.nan,
            "mean_clusters":completed.found_clusters.mean(),
            "mean_pass_throughs":audited.pass_throughs.mean(),
            "mean_root_test_rank":audited.root_test_rank.mean(),
            "mean_max_test_rank":audited.max_test_rank.mean(),
            "datasets_with_nonnull_calibration_records":int(audited.nonnull_calibration_records.gt(0).sum()),
            "calibration_records":int(audited.calibration_records.sum()),
            "nonnull_calibration_records":int(audited.nonnull_calibration_records.sum()),
            "sum_calibration_base_weight":audited.calibration_base_weight.sum(),
            "sum_nonnull_calibration_base_weight":audited.nonnull_calibration_base_weight.sum(),
        }
        if not row["signal"] and len(completed):
            count = int(completed.false_partition.sum())
            lo,hi = wilson(count,len(completed))
            expected.update({
                "false_partition_count":count,"false_partition_rate":count/len(completed),
                "false_partition_wilson_low":lo,"false_partition_wilson_high":hi,
                "observed_false_partition_per_attempt":count/len(frame),
                "unavailable_outcome_bound_upper":(count+len(frame)-len(completed))/len(frame),
                "binomial_inflation_p":binom.sf(count-1,len(completed),row["sibling_alpha"]),
            })
        for name,value in expected.items():
            near(row[name],value,f"summary {tuple(row[c] for c in CELL)} {name}")
            checks += 1
    for row in tails.to_dict("records"):
        frame = runs.loc[np.logical_and.reduce([runs[c].eq(row[c]) for c in CELL])]
        completed = frame.loc[frame.status.eq("ok")]
        audited = completed.loc[completed.audit_error.fillna("").eq("")]
        eligible = audited.loc[audited.root_p.notna() & audited.root_true_null.eq(True)]
        count = int(eligible.root_p.le(row["cutoff"]).sum())
        lo,hi = wilson(count,len(eligible))
        expected = {
            "attempts":len(frame),"completed":len(completed),"audited":len(audited),
            "eligible_root_nulls":len(eligible),"rejections":count,
            "rate":count/len(eligible) if len(eligible) else np.nan,
            "wilson_low":lo,"wilson_high":hi,
            "binomial_inflation_p":binom.sf(count-1,len(eligible),row["cutoff"]) if len(eligible) else np.nan,
        }
        for name,value in expected.items():
            near(row[name],value,f"root_tail {tuple(row[c] for c in CELL)} {row['cutoff']} {name}")
            checks += 1
    for table,mask,size in ((summary,~summary.signal,34),(tails,np.ones(len(tails),bool),272)):
        assert int(np.sum(mask)) == size
        assert table.holm_family_size.eq(size).all()
        p = table.loc[mask,"binomial_inflation_p"].fillna(1.).to_numpy()
        order = np.argsort(p)
        q_ordered = np.minimum(1.,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
        q = np.empty_like(q_ordered)
        q[order] = q_ordered
        valid = table.loc[mask,"binomial_inflation_p"].notna().to_numpy()
        near(table.loc[mask,"holm_p"].to_numpy()[valid],q[valid],f"Holm family {size}")
        assert table.loc[mask,"holm_p"].to_numpy()[~valid].size == int(np.sum(~valid))
        assert np.isnan(table.loc[mask,"holm_p"].to_numpy()[~valid]).all()
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=HERE/"evaluation")
    args = parser.parse_args()
    path = args.input_dir.resolve()
    manifest = json.loads((path/"manifest.json").read_text())
    if manifest.get("complete") is not True:
        raise RuntimeError("Evaluation is incomplete; no result verification is written.")
    assert manifest["arguments"]["phase"] == "evaluation"
    assert manifest["planned_datasets"] == manifest["completed_datasets"] == 4200
    assert manifest["method_alpha_runs"] == 23800
    for relative,expected in manifest["source_sha256"].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
    frozen = runpy.run_path(str(HERE/"run_reassessment.py"))
    runs = pd.read_csv(path/"runs.csv",dtype={"seed":"uint64"})
    assert len(runs)==23800 and not runs.duplicated(KEYS).any()
    assert set(runs.status) <= {"ok","skip","unsupported","runtime_error"}
    completed_mask = runs.status.eq("ok")
    unavailable = runs.loc[~completed_mask]
    assert unavailable[["found_clusters","ari","false_partition"]].isna().all().all()
    completed = runs.loc[completed_mask]
    null_complete = completed.loc[~completed.signal]
    near(null_complete.false_partition,null_complete.found_clusters.gt(1).astype(int),"global null K endpoint")
    assert completed.loc[completed.signal,"false_partition"].isna().all()
    datasets = runs.groupby(DATA_KEYS,sort=False)
    assert len(datasets)==4200
    for name in ("seed","data_sha256","n_samples","n_features"):
        assert datasets[name].nunique().eq(1).all(),name
    dataset_rows = runs.drop_duplicates(DATA_KEYS).to_dict("records")
    seeds = []
    dataset_cache = {}
    for row in dataset_rows:
        expected_seed = int(np.random.SeedSequence([
            20260928,17,frozen["FAMILIES"].index(row["family"]),int(row["signal"]),row["replicate"],
        ]).generate_state(1,dtype=np.uint64)[0])
        assert row["seed"]==expected_seed
        seeds.append(expected_seed)
        data,truth,space = frozen["dataset"](row["family"],row["signal"],expected_seed)
        assert data.shape==(row["n_samples"],row["n_features"])
        assert hashlib.sha256(np.ascontiguousarray(data,dtype=np.float64).tobytes()).hexdigest()==row["data_sha256"]
        if row["replicate"]==0:
            dataset_cache[tuple(row[c] for c in DATA_KEYS)] = (data,truth,space)
    assert len(set(seeds))==4200
    for group,frame in datasets:
        expected = {(method,alpha) for method in frozen["METHODS"] for alpha in (
            (.005,.01,.05) if group[0]=="bernoulli_iid" and method in {"tbs","tbs_diffusion_adaptive_nnls"}
            else (.01,)
        )}
        assert set(zip(frame.method,frame.sibling_alpha,strict=True))==expected
    for row in runs.loc[runs.actual_config.notna()].to_dict("records"):
        config = json.loads(row["actual_config"])
        assert config["edge"]["alpha"]==.001
        assert config["sibling"]["alpha"]==row["sibling_alpha"]
        assert config["config"]["adaptive_projection_dimension_energy_fraction"] is None
        assert config["config"]["spectral_minimum_dimension"]==2

    # Select deterministic examples before decoding nodes: one per family/method/status,
    # plus a signal tbs run in each family to cover mixed-population truth.
    selected = runs.sort_values(KEYS).drop_duplicates(["family","method","status"])
    signal_examples = runs.loc[runs.signal & runs.method.eq("tbs")].sort_values(KEYS).drop_duplicates("family")
    selected = pd.concat([selected,signal_examples]).drop_duplicates(KEYS)
    selected_keys = {key(row) for row in selected.to_dict("records")}
    selected_nodes = {}
    group_parts,roots,identities = [],[],[]
    node_count = raw_p_checks = adjusted_p_checks = 0
    for chunk in pd.read_csv(path/"nodes.csv.gz",dtype={"seed":"uint64"},chunksize=100000):
        node_count += len(chunk)
        identities.append(pd.util.hash_pandas_object(chunk[KEYS+["node"]],index=False).to_numpy())
        assert chunk.left_n.gt(0).all() and chunk.right_n.gt(0).all()
        assert chunk.left_population1.between(0,chunk.left_n).all()
        assert chunk.right_population1.between(0,chunk.right_n).all()
        null = chunk.left_population1*chunk.right_n == chunk.right_population1*chunk.left_n
        assert chunk.true_null.eq(null).all()
        assert chunk.loc[~chunk.signal,["left_population1","right_population1"]].eq(0).all().all()
        population_sum = chunk.left_population1+chunk.right_population1
        assert chunk.pure_parent.eq(population_sum.eq(0)|population_sum.eq(chunk.left_n+chunk.right_n)).all()
        near(chunk.df,chunk.projection_dimension,"node df versus projection rank")
        assert chunk.df.ge(0).all() and chunk.df.eq(np.floor(chunk.df)).all()
        assert chunk.df.le(chunk.left_n+chunk.right_n-1).all()
        assert np.isfinite(chunk.raw_statistic).all() and chunk.raw_statistic.ge(0).all()
        assert np.isfinite(chunk.reference_scale).all() and chunk.reference_scale.gt(0).all()
        expected_raw = np.where(chunk.df.eq(0),1.,chi2.sf(chunk.raw_statistic/chunk.reference_scale,chunk.df))
        near(chunk.raw_p,expected_raw,"raw node chi-square law")
        raw_p_checks += len(chunk)
        finite = chunk.adjusted_p.notna()
        expected_adjusted = np.where(chunk.loc[finite,"df"].eq(0),1.,chi2.sf(chunk.loc[finite,"adjusted_statistic"],chunk.loc[finite,"df"]))
        near(chunk.loc[finite,"adjusted_p"],expected_adjusted,"executed node chi-square law")
        adjusted_p_checks += int(finite.sum())
        assert (chunk.loc[finite,"adjusted_statistic"] <= (chunk.loc[finite,"raw_statistic"]/chunk.loc[finite,"reference_scale"])+1e-8).all()
        has_bh = chunk.bh_p.notna()
        assert (chunk.loc[has_bh,"bh_p"]+1e-12>=chunk.loc[has_bh,"adjusted_p"]).all()
        chunk["n_splits"] = chunk.decision.eq("split").astype(int)
        chunk["n_false_splits"] = (chunk.decision.eq("split") & chunk.true_null).astype(int)
        chunk["root_false_count"] = (chunk.n_false_splits.eq(1)&chunk.is_root).astype(int)
        chunk["descendant_false_count"] = (chunk.n_false_splits.eq(1)&~chunk.is_root).astype(int)
        chunk["pass_throughs"] = chunk.decision.eq("pass_through").astype(int)
        split = chunk.decision.eq("split")
        assert chunk.loc[split,"status"].eq("ok").all()
        assert chunk.loc[split,"calibration_role"].eq("active_traversal_sibling_gate").all()
        assert chunk.loc[split,"bh_p"].le(chunk.loc[split,"sibling_alpha"]).all()
        assert not chunk.loc[split,"calibration_record"].any()
        chunk["calibration_records"] = chunk.calibration_record.astype(int)
        chunk["nonnull_calibration_records"] = (chunk.calibration_record&~chunk.true_null).astype(int)
        chunk["calibration_weight"] = np.where(chunk.calibration_record,chunk.calibration_base_weight,0.)
        chunk["nonnull_calibration_weight"] = np.where(chunk.calibration_record&~chunk.true_null,chunk.calibration_base_weight,0.)
        chunk["node_count"] = 1
        sums = ["node_count","n_splits","n_false_splits","root_false_count","descendant_false_count","pass_throughs","calibration_records","nonnull_calibration_records","calibration_weight","nonnull_calibration_weight"]
        group_parts.append(chunk.groupby(KEYS,sort=False)[sums].sum())
        roots.append(chunk.loc[chunk.is_root].copy())
        for values,frame in chunk.groupby(KEYS,sort=False):
            if values in selected_keys:
                selected_nodes.setdefault(values,[]).append(frame.copy())
    all_ids = np.concatenate(identities)
    assert len(np.unique(all_ids))==len(all_ids)
    aggregate = pd.concat(group_parts).groupby(level=list(range(len(KEYS)))).sum()
    root_nodes = pd.concat(roots)
    assert not root_nodes.duplicated(KEYS).any()
    root_lookup = {key(row):row for row in root_nodes.to_dict("records")}
    audited = runs.loc[runs.status.eq("ok")&runs.audit_error.fillna("").eq("")]
    for row in audited.to_dict("records"):
        k = key(row)
        computed = aggregate.loc[k]
        assert int(computed.node_count)==row["n_samples"]-1
        for name in ("n_splits","n_false_splits","pass_throughs","calibration_records","nonnull_calibration_records"):
            near(row[name],computed[name],f"run aggregate {k} {name}")
        near(row["calibration_base_weight"],computed.calibration_weight,f"run calibration weights {k}")
        near(row["nonnull_calibration_base_weight"],computed.nonnull_calibration_weight,f"run nonnull weights {k}")
        near(row["split_fdp"],computed.n_false_splits/max(1,computed.n_splits),f"run split FDP {k}")
        near(row["root_false_split"],int(computed.root_false_count>0),f"root split {k}")
        near(row["descendant_false_split"],int(computed.descendant_false_count>0),f"descendant split {k}")
        near(row["found_clusters"],1+computed.n_splits+computed.pass_throughs,f"K from traversal actions {k}")
        root = root_lookup[k]
        assert root["left_n"]+root["right_n"]==row["n_samples"]
        expected_ones = 0 if not row["signal"] else row["n_samples"]//(4 if row["family"] in {"bernoulli_dependent","gaussian_ar1"} else 2)
        assert root["left_population1"]+root["right_population1"]==expected_ones
        near(row["root_p"],root["adjusted_p"],f"root p {k}")
        assert row["root_true_null"]==root["true_null"]
    summary = pd.read_csv(path/"summary.csv")
    tails = pd.read_csv(path/"root_tails.csv")
    assert len(summary)==68 and len(tails)==272
    summary_checks = reconstruct_summaries(runs,summary,tails)

    rerun_checks = []
    numeric_fields = ["found_clusters","ari","false_partition","n_splits","n_false_splits","split_fdp","root_false_split","descendant_false_split","pass_throughs","root_p","root_test_rank","max_test_rank","n_rank_above_mp_floor","calibration_records","nonnull_calibration_records","calibration_base_weight","nonnull_calibration_base_weight"]
    node_numeric = ["left_n","right_n","left_population1","right_population1","raw_statistic","raw_p","reference_scale","adjusted_statistic","adjusted_p","bh_p","df","projection_dimension","calibration_base_weight"]
    for row in selected.to_dict("records"):
        k = key(row)
        dkey = tuple(row[c] for c in DATA_KEYS)
        data,truth,space = dataset_cache.get(dkey,(None,None,None))
        if data is None:
            data,truth,space = frozen["dataset"](row["family"],row["signal"],int(row["seed"]))
        run_key = {name:row[name] for name in (*DATA_KEYS,"seed","n_samples","n_features","data_sha256")}
        fresh,node_rows = frozen["method_run"](data,truth,space,run_key,row["method"],row["sibling_alpha"])
        assert fresh["status"]==row["status"]
        assert (fresh["audit_error"] or "")==("" if pd.isna(row["audit_error"]) else row["audit_error"])
        for field in numeric_fields:
            near(fresh[field],row[field],f"deterministic rerun {k} {field}",atol=2e-8,rtol=2e-7)
        old_nodes = pd.concat(selected_nodes[k]).set_index("node") if k in selected_nodes else pd.DataFrame()
        assert len(node_rows)==len(old_nodes)
        for node in node_rows:
            old = old_nodes.loc[node["node"]]
            for field in ("decision","is_root","true_null","pure_parent","calibration_record"):
                assert node[field]==old[field],(k,node["node"],field)
            for field in node_numeric:
                near(node[field],old[field],f"node rerun {k} {node['node']} {field}",atol=2e-8,rtol=2e-7)
        rerun_checks.append({**{name:row[name] for name in KEYS},"status":row["status"],"nodes_checked":len(node_rows),"matches":True})
        print("rerun verified",k,row["status"],flush=True)
    for relative,expected in manifest["source_sha256"].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
    pd.DataFrame(rerun_checks).to_csv(path/"deterministic_rerun_checks.csv",index=False)
    result = {
        "status":"passed", "completed_manifest_required":True,
        "datasets_and_input_hashes_regenerated":len(dataset_rows),"unique_dataset_seeds":len(set(seeds)),
        "method_alpha_runs_checked":len(runs),"node_rows_checked":node_count,
        "raw_p_value_law_checks":raw_p_checks,"executed_p_value_law_checks":adjusted_p_checks,
        "audited_run_aggregates_reconstructed":len(audited),
        "summary_and_tail_numeric_checks":summary_checks,
        "holm_partition_family_size":34,"holm_root_tail_family_size":272,
        "unavailable_results_kept_separate":True,"actual_split_and_pass_through_distinguished":True,
        "integer_population_mixture_truth_verified":True,"rank_equals_df_verified_all_nodes":True,
        "deterministic_reruns":len(rerun_checks),
        "source_hashes_unchanged":True,
        "limitations":[
            "Node tables contain sibling ranks, not all edge spectral ranks or raw MP counts; edge-rank aggregate fields are verified on deterministic reruns only.",
            "This verifies extraction and calculations, not a mathematical error-control theorem.",
        ],
        "raw_artifact_sha256":{name:hashlib.sha256((path/name).read_bytes()).hexdigest() for name in ("runs.csv","nodes.csv.gz","summary.csv","root_tails.csv")},
    }
    (path/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
