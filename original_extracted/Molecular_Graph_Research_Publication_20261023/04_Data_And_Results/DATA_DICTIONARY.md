# Machine-readable dataset and data dictionary

**The authoritative files are the archived CSVs and original JSON/JSONL/logs.** This dictionary describes recorded fields, not invented biological measurement units.

Each CSV below was carried forward byte-for-byte from the authoritative publication release.

## `ALL_CONDITIONS_FROM_RAW.csv`
- **Rows:** 208
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `branch` — Computational branch recorded by the original analysis
  - `structure_id` — Molecular or PDB-labelled computational input identifier
  - `feedback` — Recorded feedback policy, including state-dependent and fixed modes
  - `plasticity` — Recorded learning status, including STDP and frozen weights
  - `paired_rows` — Number of paired computational records in the aggregate
  - `total_changed_spikes` — Sum of recorded changed target-neuron firing positions within an aggregate
  - `ticks` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `target_slots` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `fraction_from_hits` — Changed positions / (paired rows × tick count × target slots)
  - `mean_of_recorded_fractions` — Mean of archived per-pair disagreement fractions
  - `mask_errors` — Number of flagged reconstruction mismatches in a group
  - `contacts` — Number of graph contacts under the declared computed-contact rule
  - `source_jsonl` — Original relative path of native result records within the historical release


## `ALL_V9_V10_V11_PAIRED_COMPARISONS.csv`
- **Rows:** 13312
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `branch` — Computational branch recorded by the original analysis
  - `structure_id` — Molecular or PDB-labelled computational input identifier
  - `seed` — Deterministic stimulation seed identifier for this computational pair
  - `feedback` — Recorded feedback policy, including state-dependent and fixed modes
  - `plasticity` — Recorded learning status, including STDP and frozen weights
  - `tick_count` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `target_slots` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `changed_spikes` — Recorded changed target-neuron firing positions
  - `fraction_recorded` — Original stored firing-mask disagreement fraction
  - `fraction_recomputed` — Independently recomputed disagreement fraction
  - `mask_xor_count` — Popcount of recorded paired firing-mask differences
  - `mask_mismatch` — Audit mismatch marker for recorded vs recomputed masks
  - `contacts` — Number of graph contacts under the declared computed-contact rule
  - `source_jsonl` — Original relative path of native result records within the historical release


## `ATOMIC_FORCE_RESULTS.csv`
- **Rows:** 13
- **Columns:**

  - `structure_id` — Molecular or PDB-labelled computational input identifier
  - `atoms_h_included` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `energy_before` — MMFF94 computed energy before minimization, in original force-field units
  - `energy_after` — MMFF94 computed energy after minimization, in original force-field units
  - `final_gradient_norm` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `max_abs_gradient_delta` — Largest absolute analytical-versus-numerical gradient difference among tested components
  - `components_checked` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `converged` — Original minimizer convergence flag
  - `coordinates_sha256` — Cryptographic checksum of archived generated coordinates
  - `source_jsonl` — Original relative path of native result records within the historical release


## `CLAIM_EVIDENCE_MATRIX.csv`
- **Rows:** 9
- **Columns:**

  - `claim` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `status` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `source` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `interpretation` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `FIGURE_PROVENANCE.csv`
- **Rows:** 13
- **Columns:**

  - `figure` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `title` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `csv_source` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `transformation` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `png_sha256` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `FIGURE_PROVENANCE_EXTENDED.csv`
- **Rows:** 9
- **Columns:**

  - `figure` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `title` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `sources` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `FIGURE_SHA256_EXTENDED.csv`
- **Rows:** 22
- **Columns:**

  - `image_file` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `sha256` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `source_type` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `GRAPH_CLASSIFICATION_CONDITIONS.csv`
- **Rows:** 100
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `family` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `pdb` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `readout` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `intervention` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `original_accuracy` — Original held-out label accuracy for the archived readout
  - `correct` — Number of correct computational classification labels for a specified readout
  - `total` — Number of evaluated computational classification labels
  - `fresh_accuracy` — Fresh-seed held-out label accuracy for the archived readout
  - `null_exceedances` — Recorded permutation null draws at least as extreme as observed
  - `permutations` — Number of permutations used in the archived test
  - `permutation_p` — Empirical p value as recorded by original statistical code
  - `holm_p` — Archived multiplicity-adjusted probability
  - `feedback` — Recorded feedback policy, including state-dependent and fixed modes
  - `plasticity` — Recorded learning status, including STDP and frozen weights
  - `source_json` — Original relative path of source results within the historical release


## `LITERATURE_BENCHMARKS_VERIFIED.csv`
- **Rows:** 2
- **Columns:**

  - `method` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `dataset` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `f1` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `source` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `PUBLICATION_FIGURE_PROVENANCE.csv`
- **Rows:** 14
- **Columns:**

  - `figure` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `title` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `derived_from` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `caption` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `sha256` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `RECORDED_TIMING_SECONDS.csv`
- **Rows:** 14
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `stage` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `seconds` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `source_json` — Original relative path of source results within the historical release


## `RNA_FOLDING_STRUCTURES.csv`
- **Rows:** 10
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `pdb` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `length` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `sequence` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `nested_pairs` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `model_energy_kcal_mol` — Computational model energy in kcal/mol; not measured physical free energy
  - `weighted_score` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `pseudoknot_candidate_selected` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `dot_bracket` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `source_jsonl` — Original relative path of native result records within the historical release


## `SCIENTIFIC_REFERENCES.csv`
- **Rows:** 8
- **Columns:**

  - `publication_or_provider` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `role` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `url` — Archived field name; meaning must be verified from the corresponding source method before reuse.


## `TEST_CASE_INDEX.csv`
- **Rows:** 140
- **Columns:**

  - `version` — Original internal experiment identifier, retained solely for provenance
  - `suite` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `test_id` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `status` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `log_path` — Original path of historic test log used by the index


## `VALIDATION_REFERENCE_STANDARDS.csv`
- **Rows:** 7
- **Columns:**

  - `resource` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `what_it_documents` — Archived field name; meaning must be verified from the corresponding source method before reuse.
  - `url` — Archived field name; meaning must be verified from the corresponding source method before reuse.

