# Wolfram query provenance

These are the exact `code` strings sent to
`mcp__codex_apps__wolfram_wolframlanguageevaluator`, recovered from the
tool history. Both calls used `timeConstraint: 20`. The original returned
payloads are retained in [wolfram_evidence.json](wolfram_evidence.json).

## Branch scaling and the cancelling-contrast example

```wolfram
With[{mu={-a,a,-a,a}, e={1,1,-1,-1}/2, h={{1,-1,0,0}/Sqrt[2],{1,1,-2,0}/Sqrt[6],{1,1,1,-3}/Sqrt[12]}}, FullSimplify[{(1+(c l+c r)/(2 c b))-(1+(l+r)/(2 b)), (Norm[{u,v}/Sqrt[g]]^2-Norm[{u,v}]^2/g), e.mu, h.mu, (h.mu).(h.mu)}, Assumptions->c>0&&b>0&&g>0&&Element[{a,l,r,u,v},Reals]]]
```

Saved output:

```text
Out[1]= {0, 0, 0, {-(Sqrt[2]*a), Sqrt[2/3]*a, (-2*a)/Sqrt[3]}, 4*a^2}
```

The first zero verifies invariance of the normalized branch multiplier
under common positive length scaling. The second verifies the ridge-free
quadratic scaling under a fixed projection. The third shows that the
first-half versus second-half contrast cancels in this heterogeneous
four-row mean vector. Its Helmert mean still has squared norm 4a², so
the cancelling contrast does not imply the common-mean null.

## Fixed ridge and Holm grid resolution

```wolfram
{FullSimplify[1/(v+lam)-g/(g v+lam),Assumptions->v>0&&lam>0&&g>0], Ceiling[127/(1/100)]-1, Ceiling[7/(1/20)]-1}
```

Saved output:

```text
Out[1]= {(lam + v)^(-1) - g/(lam + g*v), 12699, 139}
```

The evaluator returned the ridge difference in its original two-term form.
Combining denominators in the report gives
`lam (1-g)/((v+lam)(g v+lam))`. This explains why the ridge-free
multiplicative identity need not hold with a fixed additive ridge.
The integers are the minimum step counts for the p-value grid to reach
Holm's first cutoff in the two specified families. They do not guarantee
power or successful numerical replay.
