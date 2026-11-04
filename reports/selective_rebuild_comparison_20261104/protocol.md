# Conditional selection versus whole-procedure rebuilding

Frozen before running the panel, 2026-11-04.

## Question and scope

Evaluate two candidate replacements on the same native symmetric diffusion /
average-linkage trees. Keep production code unchanged. This is a root-only,
known-identity Gaussian mechanism study, not cancer-data validation or a port
of the SHC package. The observed tree remains fixed for reporting; synthetic
trees are rebuilt to define reference distributions.

## Fixed panel

- Nine rows, three features; independent unit Gaussian noise.
- Null: 12 datasets. Moderate and strong: six each, adding means -1/+1 or
  -2/+2 in feature one to the first four / last five rows.
- Data seeds 2026110401 + replicate, reusing the previous prototype's fixed
  cases. No cases are selected according to the new results.
- Known Gaussian scale one; all three input features; no fitted preprocessing,
  covariance, projected rank, branch multiplier, empirical inflation or traversal.
- Report thresholds 0.01 and 0.05; no error-rate claim from this small panel.

## Conditional comparisons

Let e be the unit observed root contrast, r=||e'X|| and u=e'X/r.
Rebuild native diffusion at X(t)=(I-ee')X+t e u'. Under the fixed known-I
Gaussian contrast null, t has a chi(3) distribution before selection.

Compare two events: the complete observed merge history and only the unordered
pair of observed root children. Every history match must also match the root.
Coarser root conditioning is a root test; it does not justify the entire tree's
descendant decisions.

1. Existing rejection sampling: 499 chi proposals and the observation-inclusive
   rank (1 + accepted values >= r)/(1 + accepted count).
2. Importance sampling: 999 draws from 0.5 chi(3) + 0.5 N(r,1), including
   negative proposals with zero target density. Weight by chi density/proposal
   density. Record the ratio estimate, effective sample size and approximate
   simulation standard error. This is an approximation, not an exact finite
   Monte Carlo p-value.
3. Reversible-chain rank: 999 random-walk steps, standard deviation one, with
   chi-density Metropolis acceptance and event membership. Place the observation
   at an independent uniformly random position among the 1,000 path positions,
   and simulate both directions with the same reversible kernel. Count all ties
   in the upper tail. Under an exact event oracle this rank is conservative under
   the conditional null without a mixing assumption. Poor movement can still
   destroy power; report moves and unique states. Run separately for both events.

No event is approximated by a single threshold or interval. Native deterministic
merge signatures define the numerical event, including native tie breaking;
these are not exact-real interval certificates. Geometry failures stop a case
rather than becoming nonselection. The original strict replay is checked on a
fixed fixture separately; differing numerical membership policies are explicit.

## Whole-procedure rebuilding

At each rebuilt root, use BSS/TSS, the fraction of total squared variation
explained by its two children. This is one minus the usual cluster index.
Generate 999 independent known-I null datasets once and rebuild diffusion and
average linkage for each. Compare each observation by (1 + null scores at least
as large)/(1 + 999). The shared reference bank is independent of observations;
counts across observations share simulation noise and are not independent.

Also generate 499 null datasets per observation from its fitted full sample
covariance, rebuilding the same steps and score. This plug-in variant has no
finite-sample guarantee. Neither arm recreates production TBS inference: they
replace its score and root decision. The tested null is one Gaussian population,
which is stronger than equal averages in the two selected children.

## Checks and outputs

Before the panel: root extraction under reordered commuting merges; BSS identity;
importance estimates against a known truncated-chi law; random-position chain
on an exactly enumerable two-state model; source/environment hashes. Baseline
native replay and conditional-rank tests passed 16 tests. Save per-case records,
proposal and chain arrays, inputs, compact summaries and a comparison figure.
Report numerical failure, low simulation support and hypothesis differences.
