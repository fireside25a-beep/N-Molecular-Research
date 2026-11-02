# Methods and numerical contract

## Research question

The native experiment measures how a recurrent spiking neural network responds when the ligand chemical bond graph is enabled versus disabled, with all ligand atoms retained. It asks about computational sensitivity, not physical drug binding. The source of truth is `src/bridge.cpp` and the immutable historical NeuroMantle core in `core/`.

## Input admission

The entry point `tools/prepare_smiles.py` uses RDKit 2025.09.4 in the tested environment. Sanitization rejects impossible valences, malformed SMILES, disconnected multi-component inputs, duplicate identifiers and unsupported bonds. An isomeric canonical SMILES is reparsed before atom records are serialized. Atomic number, formal charge, aromaticity, ring membership, isotope, chirality, bond multiplicity and bond stereo survive the native graph exchange. The exchange format includes the SHA-256 of its exact serialized graph. The native engine independently checks its boundaries, indices, duplicate bonds, connectivity, and digest. This is an integrity check, not cryptographic proof of a molecular assay.

## SNN configuration

The unchanged 118 original scaffold neurons and 324 directed scaffold edges are installed first. Another 22 neurons represent the experimentally deposited 1KF1 parent sequence as categorical labels, and additional neurons represent heavy atoms from an RDKit parsed ligand graph. The three four-membered contact rings are manually declared in the source as a computational hypothesis. They are **not** atomically extracted 1KF1 G-tetrad interactions. Computational relay connections from ligand atoms to target and scaffold are likewise explicit hypotheses and do not represent measured chemical bonds, forces or ligand contacts.

The original native `nm_add_neuron`, `nm_add_edge`, `nm_inject` and `nm_step` functions are called for each independent 48 tick trial. The original integer leaky integrate and fire transitions and bounded STDP are preserved. The built in capacity remains 1,024 neurons and 16,384 edges. The user receives an explicit refusal on capacity overflow, not silent truncation.

## Factorial controls

For each ligand and deterministic entropy seed, eight fresh native network states are computed. The three binary factors are presence versus ablation of the internal ligand bond edges, closed versus matched-current fixed feedback, and native STDP versus frozen synaptic weights. The same Surreal entropy tape drives the same exogenous ligand and target injections within every matched four-condition block. Both feedback types inject current 12 for ticks 1 through 47; only the selection policy differs.

Closed-loop feedback at tick `t` depends on the preceding `nm_step` output and the native synaptic weights on the original 324-edge scaffold, which prevents a trivial change in the extra ligand edge count from automatically changing feedback. The target is `(spike_signature + scaffold_weight_signature) mod 118`. Both signature inputs and every feedback target are preserved in the receipts. Fixed schedule targets are `(17 + 31*t) mod 118` and do not inspect SNN state.

## Observable and interpretation

For the 22 categorical target nodes, let `A_t` and `B_t` be 22 bit firing masks in a matched intact/ablated pair. The reported sensitivity is

\[
D = \frac{1}{48\cdot22}\sum_{t=0}^{47}\operatorname{popcount}(A_t\oplus B_t).
\]

This fraction ranges from zero to one. The independently implemented Python audit reconstructs it from both 48-element mask arrays. It also checks target spike totals, synaptic updates, same-seed exogenous tape hashes, controller arithmetic, scaffold capacity, the four control conditions and input commitments. This measures computational graph-intervention sensitivity only. A larger number is not a prediction of higher ligand potency, affinity, selectivity, stabilization, accessibility or toxicity.

The dispersion over deterministic seeds is a description of the computational stimulation schedule. These seeds are neither replicate chemical assays nor biological replicates. No confidence claim about physical binding follows from them.

## Physical structure evaluation

PDB 1KF1 is a deposited 2.10 Å diffraction structure of a parallel human telomeric G-quadruplex that includes potassium ions. 3D coordinates are physically meaningful only when read from a real deposited PDB/mmCIF file. `tools/validate_coordinates.py` checks the actual 22-residue target sequence, finiteness and bounds of atom coordinates, per-chain sequences and counts, deposited potassium sites, guanine O6 coordinates and simple Cartesian distances. It calculates observed coordinate summaries from the supplied file, separately from the SNN.

Neither RDKit 2D SMILES coordinates nor graphical SVG layout is used as 1KF1 atomic ground truth. The coordinate validator is deliberately excluded from the SNN learning path. The current release has no measured ligand complex or binding labels and therefore cannot validate docking, screen sensitivity/specificity, lead enrichment or the rejection of inactive compounds.

## Known statistical boundary

Earlier v6 decoding of the hand-authored 1KF1 contact intervention was sensitive to four-membered versus two-membered graph cycles. Under a ring-preserving intervention it did not retain the original apparent signal. That negative control remains preserved in `research/v6/`. In v7, none of four additional RNA targets passed its corrected primary held-out classification test. The current ligand screening pipeline does not recycle those old decoder percentages as a measure of chemical efficacy.

## Additional information observability

`tools/entropy_audit.py` uses the native 48-tick target bitmasks to form an empirical distribution over the 22 target neuron identities for each paired graph state. It reports base-2 Shannon entropies and Jensen Shannon divergence, with the usual unit bound of 1 bit for two distributions. Both measurements summarize observed native spike patterns and remain distinct from molecular thermal entropy or measured ligand binding free energy.
