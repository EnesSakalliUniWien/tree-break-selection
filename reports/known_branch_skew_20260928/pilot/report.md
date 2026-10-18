# Known-branch skew diagnostic

Generated for diagnostic attribution only. It does not establish that the selected production hierarchy is calibrated and does not change production methods or defaults.

## Protocol

Each replicate is one independent dataset with four leaf groups, each containing 32 observations in 8 dimensions. Population covariance profiles are identity and AR(1) with correlation 0.6. Each replicate is reused across matched local variants. For each fixed sibling node, the saved left and right population means and their vector contrast are computed directly from the named leaf-group means; `population_null` means exact equality of those population vectors. The JSON vector columns preserve the actual contrast, including nodes whose descendant groups differ while their parent means coincide.

| Stage | Fixed component | Changed component | Interpretation |
|---|---|---|---|
| Covariance | fixed node membership; full whitened coordinates | known population covariance vs sample pooled covariance | The sample pooled covariance includes between-sibling mean variation; it is labeled as such and is not silently called a within-child estimator. |
| Projection/MP | fixed node membership; covariance held within matched rows | fixed full basis vs PCA and production MP rank from all raw descendant observations | Production Cholesky whitening puts PCA rows and the tested contrast in the same coordinates. Eigenvalue ties are retained by production rank logic; this controlled diagnostic uses MP floor 1. |
| Selection | same dataset; oracle covariance and full basis | fixed root membership vs same-data average-linkage root membership | The selected root is a separate node estimand. Its exact composition-weighted population means are saved; selected roots are not pooled with fixed-node null calibration. |
| FDR/traversal | known tree; exact Gaussian local p-values | production traversal-aligned sibling BH and `TreeDecomposition`, with passthrough disabled/enabled | Every edge prerequisite is deliberately open. Outcomes are empirical for this controlled regime and count only reached `split` decisions. |

## Planted branch patterns

| Pattern | True sibling nodes | Leaf population means |
|---|---|---|
| `global_null` | R, A, B null | all four leaf means are zero |
| `root_only` | R signal; A and B null | A0=A1=+mu; B0=B1=-mu |
| `A_deep_only` | R and B null; A signal | A0=+mu; A1=-mu; B0=B1=0 |
| `separate_deep_branches` | R null; A and B signal | A0=+mu; A1=-mu; B0=+nu; B1=-nu |
| `strong_root_weak_A` | R strong; A weak; B null | A0/A1 average +mu with contrast 0.5 nu; B0=B1=-mu |

`mu = 1.25 e1` and `nu = e2`; all remaining coordinates are zero. Every internal node has a stored population contrast vector, not only a purity or label summary.

## Outcome

- `ar1_0.6`, passthrough=False: observed any reached false split 0.0000 (n=5); global-null ideal reference 0.0100.
- `ar1_0.6`, passthrough=True: observed any reached false split 0.0000 (n=5); global-null ideal reference 0.0199.
- `identity`, passthrough=False: observed any reached false split 0.0000 (n=5); global-null ideal reference 0.0100.
- `identity`, passthrough=True: observed any reached false split 0.0000 (n=5); global-null ideal reference 0.0199.

For the ideal global-null traversal reference, known-covariance full-basis sibling contrasts at R, A, and B are independent in this balanced Gaussian design. The root family is a singleton and the deeper family is BH over A/B. With passthrough disabled, only a root split can be reached; with passthrough enabled, a deeper BH rejection can be reached through a root pass-through. This analytic reference relies on exact Gaussian known-covariance p-values and deliberately open edges; it is not evidence about empirical production p-values or selected hierarchies.

## Run record

- Replicates per pattern/profile: 5
- Seed: 20260928
- Alpha: 0.01
- Source commit: `f34f59cacd09b19c34fafa3a4a6e783dd3973f72`
- Rows include `ok`, `unsupported`, and `runtime_error` status plus failure reason.
- Selected-tree root p-values are descriptive after same-data hierarchy selection.
- Raw-statistic skew is reported beside the expected chi-square skew at each observed degree of freedom; p-value rejection rates and quantiles are reported separately.
