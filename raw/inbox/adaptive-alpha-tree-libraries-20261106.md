# Libraries and rules for adaptive alphas on the tree

Checked on 2026-11-06 by three web-research agents, then re-checked directly where
marked "verified". Nothing was installed.

## Rules taken from the primary papers (verified, arXiv full text)

- Alpha-investing (Foster and Stine 2008, JRSS-B 70(2):429-444,
  https://doi.org/10.1111/j.1467-9868.2007.00643.x), as tabulated in Ramdas,
  Yang, Wainwright and Jordan, "Online control of the false discovery rate with
  decaying memory" (https://arxiv.org/abs/1710.00499), Table 2: level
  `alpha_t = phi_t / (1 + phi_t)`, penalty `phi_t <= W(t-1)`, reward
  `psi_t = phi_t + B0`. Their Table 1 lists its dependence condition as
  SuperCoAD (null p-values super-uniform conditional on all earlier decisions)
  and its error rate as mFDR. They state that Aharoni and Rosset showed
  `W0 = B0 = alpha` suffices for mFDR_1.
- AdaPT (Lei and Fithian 2018, JRSS-B 80(4):649-679,
  https://doi.org/10.1111/rssb.12274; https://arxiv.org/abs/1609.06035):
  `FDP_hat = (1 + A_t) / (R_t v 1)`. Theorem 2.3 assumes null p-values
  independent of each other and of the non-nulls, uniform or mirror-conservative.
- STAR (Lei, Ramdas and Fithian, Biometrika 108(2):253-267,
  https://doi.org/10.1093/biomet/asaa064; https://arxiv.org/abs/1710.02776):
  `FDP_hat = (h(1) + sum_{i in R} h(p_i)) / (1 + |R|)`, `h` non-decreasing with
  integral one, so `h(1) >= 1`. Theorem 2.3 assumes null p-values mutually
  independent and independent of the non-nulls conditional on the covariates,
  with non-decreasing null densities. The paper names hierarchical testing as a
  case where a rejected node's parent must also be rejected.

Consequence: at alpha 0.05, AdaPT needs at least 20 rejections and STAR at least
19 before either can stop (39 with the paper's example `h(p) = 2 * 1{p >= 0.5}`).

## Packages

| Method | Package | Language | Status | Tree constraint | Source |
| --- | --- | --- | --- | --- | --- |
| Alpha-investing, LORD, SAFFRON, ADDIS | `online-fdr` | Python | 1.0.0, uploaded 2026-10-19, BSD-3-Clause (verified on PyPI) | no | https://pypi.org/project/online-fdr/ |
| Same family | `onlineFDR` | R / Bioconductor | 2.20.0 (agent-reported); `Alpha_investing` uses SAFFRON's update with `lambda_i = alpha_i`, initial wealth alpha/2 by default (agent-reported) | no | https://bioconductor.org/packages/onlineFDR |
| AdaPT | `adaptMT` | R | archived from CRAN on 2026-02-24 (verified); GitHub lihualei71/adaptMT | no | https://CRAN.R-project.org/package=adaptMT |
| STAR | lihualei71/STAR (paper code, not a package) | R | last push 2019-05-13, no license (verified via GitHub API); `R/STAR_tree.R`, `R/STAR_DAG.R`, `R/generic_STAR.R` | yes | https://github.com/lihualei71/STAR |
| IHW | `IHW` | R / Bioconductor | agent-reported 1.38 to 1.40 | no | https://bioconductor.org/packages/IHW |
| AdaFDR | `adafdr` | Python | PyPI 0.1.8, 2020 (agent-reported) | no | https://github.com/martinjzhang/adafdr |
| Selective test, hierarchical clustering | `clusterpval` | R, GitHub only | DESCRIPTION 1.0.1, MIT (agent-reported) | n/a | https://github.com/lucylgao/clusterpval |
| Count splitting | `countsplit` | R / CRAN | 4.0.0 (agent-reported) | n/a | https://cran.r-project.org/package=countsplit |
| Weighted Bonferroni | none found in Python; `gMCP` in R (agent-reported) | | | | |
| Hierarchical FWER / FDR | `hierinf` (Bioconductor), `structSSI` (CRAN) (agent-reported) | R | | tree | |

Agents found no Python port of AdaPT, STAR or clusterpval. Neither online-FDR
package documents choosing the next hypothesis from earlier decisions, which tree
descent does; the tree use rests on the SuperCoAD condition instead.

## Choice for the prototype

- Alpha-investing: implemented directly from the verified Table 2 rule (a few
  lines, with a hand-computed self-check); `online-fdr` is the reference package if
  a cross-check is wanted later.
- AdaPT and STAR: not run, because a family has at most 7 nodes and both need at
  least about 20 rejections.
- Weighted Bonferroni with alphas fixed by each node's event: implemented
  directly; no Python library exists.
