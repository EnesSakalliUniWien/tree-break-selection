# Frozen score and family development panel

Lock before development draws. Use 50 datasets per cell from independent
NumPy SeedSequences [7811260911, cell index, replicate] and separate spawned
integer chain seeds. Use paired 99 and 999 transitions on the exact same
inputs and seeds. The path rule and score stay fixed. Nominal alpha=.05.
These are development diagnostics. Final calibration and recovery gates
are separate.

- Cells 0 and 1: n=16,p=3 correlated Gaussian rows, common mean or shift
  6 in the last eight rows. Fit native Ward/NNLS and record raw, displayed
  calibrated statistic where eligible, unit branch factor, oracle
  covariance control, and prespecified leading Gram rank 1. Require
  unit-branch score to equal native raw times g and rank to agree,
  conditional on the same fitted components. Use true simulation covariance
  only for the oracle. No claim of analytical validity for these variants.
- Cell 2: n=6,p=2 homogeneous. Fixed three-versus-three contrast, leading
  Gram rank 1. Compare the exact Beta helper to both orbit-chain budgets
  using the same normalized projection score. Record score identity and
  rejection rates. Exact Beta and finite-path ranks need not coincide.
- Cells 3 through 5: homogeneous n=4 root, homogeneous first child of an n=8
  Ward root with the final four rows shifted by 14, and mixed n=8 root with
  the final four rows shifted by 4. At each internal node, use the
  branch-free squared first-feature normalized contrast, or the constant
  zero at two-row nodes. Replay the common Ward selection event and apply
  Holm. Classify a true null by common generating mean, not by equality of
  selected child average means. Record false and true family rejections,
  unresolved outcomes, attainable p grid, movement, and runtime. No
  seed substitution or optional stopping.

Data are Gaussian with feature factor [[1,.55,.1],[0,.8,.25],[0,0,1.3]]
truncated to p dimensions, and covariance factor-transpose times factor.
Report exact binomial intervals over all attempted datasets.

Validation correction after the fourth attribution attempt. Production adds a
fixed ridge after variance scaling, so the simple g multiplication is an
unregularized identity. Preserve the three written attempts and the original
runner. Record the exact g=1 re-evaluation, its difference from the
unregularized identity, and any rank changes instead of requiring the
identity. Retain unsupported native fits in the fixed denominator. Restart
on the original data and seeds, with the same score, cells and budgets.
