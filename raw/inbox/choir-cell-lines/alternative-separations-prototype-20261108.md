# Prototype: alternative separations where a subtree's top eigenvalues tie

Prototype written and run on 2026-11-08. Local folder: `results/analyses/alternative_separations_20261108`
(git-ignored); code, runner, summary and output in `alternative-separations-prototype-scripts-20261108.py.txt`.
No protocol file: settings were fixed in the code before the first run; the real-data pass on smaller subtrees
was added after it and is exploratory.

## What it does

For each subtree (simulated: >= 60 cells; real: >= 500 cells with a smaller child of >= 50):

1. Spectrum in the 30 test PCs as in `subtree-spectrum-topology-20261108.md`: eigenvectors above the noise band,
   multiplicity groups by the block equality test. Flag the subtree when its top group has multiplicity m >= 2.
2. Cut the subtree along the tree into m + 1 blocks (split the block with the largest merge height first; peels
   smaller than max(20, 2% of the subtree) are walked through and their cells left out).
3. Report every two-sided split of the blocks: strength n_S n_R / n * d' W^-1 d (W pooled within blocks), whether
   the tree holds it (one side is the block set under a tree node), and its strength relative to the tree's own
   split of the subtree. Splits the tree does not hold are the alternative separations.

With m = 2 there are 3 blocks and a tree holds all 3 splits, so alternatives need m >= 3.

Self-check (assert-based, in the module): 4 groups of 250 cells at equal distances in 30 dimensions, Ward tree:
the root is flagged with m = 3 and 4 blocks, each block is one group, 5 of 7 splits are held, the two dropped 2+2
splits are 0.998 and 0.990 as strong as the tree's split; groups at unequal spacing are not flagged at the root.

## Simulated data (185 datasets, both trees)

| Data | Flagged subtrees per dataset (blocks from >= 2 groups / within one group) | Alternatives at flagged subtrees with blocks from >= 2 groups |
| --- | --- | --- |
| Balanced 4 groups | 2.0-2.3 / 6.9-12.2 | root flagged in 10/10 per tree, m = 3, blocks pure; median 1-2 alternatives at 1.03-1.05 of the tree's split |
| Imbalanced, singleton, density variation | 0 / 6.9-27.4 | none |
| No groups (sd 0.25, 0.42, 1.0) | 0 / 0.2-3.9 | none |

Flags inside one population: m = 2 in 554, m = 3 in 315, m = 4-10 in 498 subtrees; for m >= 3 they report
alternatives too (best at 0.88-0.98 of the tree's split, medians by m).

## Real CHOIR data (48,879 cells, 190 lines)

- Subtrees >= 500 cells with a smaller child >= 50: 11 flagged (default) and 10 (Markov/Ward), all with m = 2, so
  no alternative separation can be reported.
- Exploratory, subtrees >= 100 cells with a smaller child >= 20: 206 flagged (default; m = 2: 151, 3: 43, 4: 12)
  and 308 (Markov/Ward; m = 2: 187, 3: 85, 4: 23, 5: 9, 6: 4). Of the m >= 3 flags, 8 (default) and 2
  (Markov/Ward) have blocks led by different lines, and none of those has all blocks >= 90% one line. The rest
  are inside one line (for example four blocks all NCIH2110_LUNG).

## Reading

The prototype finds the designed case: when separate groups sit in one eigenspace, it names the separations
the tree had to drop, and they are about as strong as the one it kept. Inside subtrees the tree cut out, it
also flags ties within a single population, because a subtree chosen by the tree has a spectrum shaped by that
choice and the noise band read from its own bulk does not account for it. In the real data it finds no clean
tie between cell lines. Before its flags can be trusted below the root, the noise band of a selected subtree
needs a reference that includes the tree's choice.
