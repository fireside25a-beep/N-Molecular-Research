# Methods and mathematical definitions

## Computational experimental unit

One paired comparison consists of two independently initialized 48-tick native SNN runs under matched declared seeds and stimulation policies, with one computational graph input differing by a specified intervention. Deterministic seeds are not biological replicates.

## Held-out classification

\[\mathrm{Accuracy} = N_{\mathrm{correct}}/N_{\mathrm{evaluated}}.\]

The main cycle-breaking case reports 196/256 correct labels; the ring-preserving control reports 116/256. These percentages are classification scores on computational graph labels.

## Mask disagreement

\[D=\frac{\sum_{t=0}^{47}\operatorname{popcount}(A_t\oplus B_t)}{48\times22}.\]

A and B are target-neuron firing masks from paired original native executions. Recorded changed positions sum over 48 steps; each step has 22 target positions. The package provides every per-comparison row and condition aggregation.

## Permutation probabilities

The historical six-case experiment records a 4,096-permutation sign-swap procedure with a plus-one empirical probability:

\[p=(1+\#\{a_{\mathrm{null}}\geq a_{\mathrm{observed}}\})/(4096+1).\]

Source results include Holm step-down adjustments for documented hypothesis families. These are computational conditional analyses, not inferences from repeated physical biological samples.

## Model boundaries

The weighted RNA pairing routine uses computational base-pair scores. The later RNA model uses selected published energy parameters and simplified loop terms; pseudoknot search results are heuristic crossing-stem candidates. The atomistic path uses MMFF94 force-field minimization and finite-difference gradient checks on small molecules and short RNA fragments. Source documents describe the exact implementation.

## Audited example arithmetic

- Graph classification: 196 / 256 = 0.765625 = 76.5625%.
- Ring-preserving contrast: 116 / 256 = 0.453125 = 45.3125%.
- The paired disagreement denominator for one run pair is 48 × 22 = 1,056 temporal positions.
- The later-stage paired records aggregate 2,048 + 3,328 + 7,936 = 13,312 comparisons, but their distinct protocols are not pooled into a single performance metric.

For complete network update equations and the actual C/C++/Python algorithm, use the native source archive and historical methods rather than inferring them from this summary.
