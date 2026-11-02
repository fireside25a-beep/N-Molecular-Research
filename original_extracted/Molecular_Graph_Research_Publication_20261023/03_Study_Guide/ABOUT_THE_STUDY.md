# Understanding the study

## What question did the researchers investigate?

The study asks whether a compiled, deterministically driven spiking neural network retains measurable information about selected changes in computational molecular graph connections, and whether recurrent feedback and learning alter those measurements.

A molecular graph represents sites and relationships by nodes and edges. Some experimental inputs are based on real deposited DNA/RNA sequences or chemically valid SMILES records. **The original graph rewires were constructed computationally**, not measured atomic refolding experiments.

## Why investigate it?

A spiking network evolves over time, so a structural change in its inputs can influence later activity. The study tests that possibility with matched resets and independent computational seeds, using readouts that can be reconstructed from recorded outputs. This is a study of *graph representations and neural dynamics*, not a clinical or pharmaceutical prediction test.

## How was the experiment performed?

The preserved native scaffold contains 118 neurons and 324 directed synapses. Study runs advanced it through 48 deterministic ticks. Graph interventions changed selected connection patterns while matching declared seed and stimulation conditions. Four operating settings crossed native-output-conditioned or matched-fixed feedback with enabled or frozen synaptic plasticity.

Two different measurements appear in the results:

1. **Graph-label classification:** the number of correct held-out intact-versus-rewired labels divided by all evaluated labels.
2. **Neural disagreement:** the fraction of target-neuron firing positions that differ between independently reset paired native runs. Its denominator is 48 ticks × 22 recorded target slots.

These measurements must not be averaged together or compared as if they measured the same scientific quantity.

## What happened in the main structural experiment?

In the constructed graph associated with the deposited PDB identifier **1KF1**, the cycle-breaking contrast yielded **196/256 correct held-out labels (76.56%)**. A separate cycle-preserving contrast yielded **116/256 (45.31%)**. These are different interventions, and the cycle-preserving experiment was introduced as a later control. The results support selectivity in this particular computational setting, not a proof that the network recognizes all topological equivalences.

Other archived experiments included non-significant classifications and conditions with null or reversed response differences. They are retained alongside positive results.

## How were chemistry and RNA inputs added?

Chemical bond connectivity from SMILES, calculated nested RNA base pairs, and RDKit-derived spatial contacts were converted to native input connections. A later method used selected RNA energy parameters and searched for crossing-stem candidates. Atomic-force checks used RDKit MMFF94 on eight small molecules and five four-nucleotide RNA fragments, not complete RNA molecules. These inputs affected the SNN's *computational output*; measured binding energies and experimental RNA folding accuracy were not established.

## How much recorded computational evidence is present?

The independently audited later-stage data contain **13,312 paired comparison rows**, **208 condition summaries**, and **100 historical classifier-statistic rows**. **81 named final test cases** appear as passing in preserved historical logs. Original native execution and source build evidence is supplied separately from this textual summary.

## What the results cannot establish

The available records do not establish complete physical RNA folding, a molecule's drug efficacy, ligand affinity, clinical benefit, or superiority to published structure predictors on identical experimentally validated test sets. They establish the narrower computational question described in the article.

## Where do the numbers come from?

See `../04_Data_And_Results/` for machine-readable source-derived records, `../06_Experimental_Evidence/` for historic logs and original native outputs, and `../10_Provenance_And_Integrity/` for the original filename-to-current-location and checksum ledger. All quantitative claims in this guide are drawn from the preserved manuscript or its underlying recorded tables.
