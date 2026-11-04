# Root conditional sampler development validation

Fixed before validation, 2026-11-04.

The maintained diagnostic API is `selected_diffusion_root_chain_test` in
`benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py`.
It uses the root-split random-position reversible chi chain from the preceding
comparison. No production TBS gate is replaced. Known spherical Gaussian scale
is mandatory; default steps=999, seed=0, alpha=0.05; proposal scale remains one.

Validation has three parts:

1. Replay all 24 previously evaluated nine-by-three datasets with the same
   sampler seeds. Require exact agreement of stored radii, ranks and moves for
   the root-split chain, not only aggregate counts.
2. A fresh fixed panel: 40 null, ten moderate, ten strong nine-by-three datasets.
   Data seeds 2026110701 + replicate, chain seeds 2026111701 + replicate.
   Unit noise; first four versus last five means -shift/+shift in feature one,
   with shifts 0, 1, 2. Use 999 steps and alpha 0.05. Report 0.01 as a secondary
   fixed threshold without choosing a different sampler.
3. The registered `gauss_clear_medium_continuous` and
   `gauss_moderate_3c_continuous` cases at their existing sizes and seeds, with
   their declared spherical noise scales and unchanged input matrices. Use
   999 steps, seeds 2026112701 and 2026112702. This checks the supported raw
   Gaussian benchmark path; it does not validate fitted-PC cancer inputs.

Archive all inputs and complete results, including unresolved cases. Do not
resample failed cases or drop them from denominators. Recompute every complete
rank from the saved chain. Check reconstruction, observed membership, and a
small fixed set of chain memberships through native replay. Preserve all prior
reports. Tests and the existing ordered suite are run separately; their status
does not establish the statistical null law. No deployment, installation or
full cancer rerun is authorized by this development step.
