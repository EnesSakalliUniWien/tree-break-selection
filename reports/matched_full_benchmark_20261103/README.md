# Matched full benchmark rerun — 2026-11-03

Baseline: `run_20261019_064234Z_full`, commit `8e0093d3`. All 1,342 case/run pairs and recorded parameters, seeds, sample/feature counts, true K and benchmark metadata match. Case/generator source and lockfile match the baseline commit; original-run matrix hashes were not recorded. New reconstructed input hashes are saved.

| method                      |   ok_before |   ok_after |   mean_ari_before |   mean_ari_after |   mean_ari_delta |   exact_k_before |   exact_k_after |   ari_improved |   ari_degraded |   null_split_before |   null_split_after |
|:----------------------------|------------:|-----------:|------------------:|-----------------:|-----------------:|-----------------:|----------------:|---------------:|---------------:|--------------------:|-------------------:|
| dbscan                      |         122 |        122 |          0.576000 |         0.576000 |         0.000000 |               71 |              71 |              0 |              0 |                   0 |                  0 |
| hdbscan                     |         122 |        122 |          0.603754 |         0.603754 |         0.000000 |               56 |              56 |              0 |              0 |                   4 |                  4 |
| kmeans                      |         122 |        122 |          0.857050 |         0.857050 |         0.000000 |              122 |             122 |              0 |              0 |                   0 |                  0 |
| leiden                      |         122 |        122 |          0.807762 |         0.807762 |         0.000000 |               92 |              92 |              0 |              0 |                   4 |                  4 |
| louvain                     |         122 |        122 |          0.804345 |         0.804345 |         0.000000 |               93 |              93 |              0 |              0 |                   4 |                  4 |
| optics                      |         122 |        122 |          0.528018 |         0.528018 |         0.000000 |               66 |              66 |              0 |              0 |                   0 |                  0 |
| spectral                    |         122 |        122 |          0.826922 |         0.826922 |         0.000000 |              122 |             122 |              0 |              0 |                   0 |                  0 |
| tbs                         |         122 |        122 |          0.794135 |         0.805953 |         0.011818 |               91 |              93 |              3 |              0 |                   0 |                  0 |
| tbs_diffusion               |         109 |        109 |          0.530765 |         0.537087 |         0.006322 |               39 |              41 |              3 |              0 |                   0 |                  0 |
| tbs_diffusion_adaptive_nnls |         101 |        101 |          0.821914 |         0.823085 |         0.001171 |               85 |              84 |              1 |              0 |                   0 |                  0 |
| tbs_family_diffusion        |         122 |        122 |          0.824874 |         0.825843 |         0.000969 |               99 |              98 |              1 |              0 |                   0 |                  0 |

## Changed results

| case_id                           | method_before               |   ari_before |   ari_after |   found_clusters_before |   found_clusters_after |
|:----------------------------------|:----------------------------|-------------:|------------:|------------------------:|-----------------------:|
| cat_highcard_20cat_4c             | tbs                         |     0.000000 |    0.716061 |                       1 |                      4 |
| phylo_conserved_4taxa             | tbs                         |     0.000000 |    0.707965 |                       1 |                      3 |
| phylo_large_64taxa                | tbs_diffusion               |     0.180729 |    0.500613 |                      43 |                     54 |
| cat_highcard_20cat_4c             | tbs_diffusion               |     0.654333 |    0.896517 |                       3 |                      4 |
| cat_simplex_face_rare_20cat       | tbs_diffusion               |     0.208377 |    0.335422 |                       3 |                      4 |
| cat_dirichlet_overdispersed_20cat | tbs_diffusion_adaptive_nnls |     0.214239 |    0.332460 |                       4 |                      7 |
| cat_dirichlet_overdispersed_20cat | tbs_family_diffusion        |     0.214239 |    0.332460 |                       4 |                      7 |
| phylo_large_64taxa                | tbs                         |     0.982208 |    1.000000 |                      63 |                     64 |

## Interpretation limits

The current working tree includes both the self-distance correction and the prior empirical-inflation repair. This is an integrated implementation comparison, not an isolated causal estimate of the self-distance change. K-means and spectral clustering receive true K. Skipped cases are excluded from quality means and counted separately. Six suite null cases are not a calibrated Type-I-error study. Stage timings omit parts of diffusion/NNLS and do not provide comparable end-to-end method runtimes.

## Execution and reproduction

The full run completed all 122 cases and 1,342 rows with exit code 0 in 654.13 seconds (10 minutes 54 seconds), including plots and report generation. No case failed; all 210 recorded source hashes and the baseline CSV hash remained unchanged. The runner baseline tests passed 8 tests before execution.

The complete recorded environment is in `execution.json`; data recipes, method grids, dependency versions and source hashes are in `matched-run-manifest.json`. `input-fingerprints.csv` records reconstructed matrices and labels. The working-tree source changes are preserved in `source-changes.patch`.

Run the existing full entry point from the repository root, using a fresh output directory, with the environment in `execution.json`:

```bash
TBS_RUN_DIR=/absolute/path/to/fresh-output uv run --no-sync python -u -m benchmarks.full.run
```

The original defaults are full cases, the 11 default methods, plots and UMAP case isolation enabled, manifold and tree consensus disabled, relationship analysis/plots enabled, 1 spectral worker, and a 1,800-second case timeout. Numeric thread variables are all 1. No threshold was tuned.

To repeat the matched checks on the recorded run:

```bash
.venv/bin/python reports/matched_full_benchmark_20261103/compare_matched.py benchmarks/results/run_20261103_matched_full_ixx1jf0m
```

The checks require the recorded source state and local original result paths. Full CSV snapshots are retained as local-only generated outputs under `generated/`, with hashes in `reports/artifacts.json`; compact paired results and summaries are retained alongside this report. The original full report and 122 per-case PDFs remain under `benchmarks/results/run_20261103_matched_full_ixx1jf0m/`. Nothing was published or committed.

## Main tradeoff

Family-diffusion TBS mean ARI rises from 0.824874 to 0.825843, while exact K falls from 99/122 to 98/122. Its only changed result is `cat_dirichlet_overdispersed_20cat`: ARI 0.214239 to 0.332460, but K 4 to 7 when true K is 4. The adaptive-NNLS variant has the same change. No ARI declines occur at absolute tolerance 1e-12; this does not mean every performance measure improved. Both full-coverage TBS variants still return one cluster on all six suite nulls. The outlier-shell weakness remains (family-diffusion ARI 0). Non-TBS methods match the baseline.
