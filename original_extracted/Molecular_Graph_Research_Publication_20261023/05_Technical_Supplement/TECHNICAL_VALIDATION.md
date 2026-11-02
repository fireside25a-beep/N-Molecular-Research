# Technical validation and traceability

## Source-controlled evidence and negative results

The manuscript combines distinct computational tasks but never equates their output measures: the original graph-label classifier, paired firing-mask ablations, RNA algorithm outputs, and MMFF94 minimization. The original 48-tick mask calculations are reported as independently checked for 13,312 rows; the condition-level summaries count 208 and the historical classifier-statistic table contains 100 rows.

**Negative findings remain scientifically relevant:** four later RNA-labelled graph cases did not pass their recorded multiple-testing correction. The ring-preserving 1KF1-labelled control did not demonstrate reliable graph-label classification. Ethanol's selected spatial-contact definition produced zero qualifying contacts and no resulting firing-mask change. Some frozen-plasticity comparisons exceeded the analogous learning-enabled response.

## Source-dependent claims that must remain distinct

- Passing tests and deterministic replay establish software consistency within the tested environment, not chemical or biological accuracy.
- The presence of experimentally solved PDB structures does not imply that deposited atomic coordinates were used in all original graph interventions.
- RNA pairing model energies are not experimental ΔG measurements.
- MMFF94 minimization of the first four nucleotides is not optimization of the full RNA structure; identical initial four-base sequences are not independent confirmations.
- A sampled finite-difference force-gradient check verifies implementation consistency for sampled Cartesian components, not external validation of physical force fields.
- Original software logs may include successful test subsections in a terminated wider test command; inspect actual exit status before claiming complete-suite success.

## Supporting records

- `../04_Data_And_Results/ALL_V9_V10_V11_PAIRED_COMPARISONS.csv` — per-pair computations.
- `../04_Data_And_Results/ALL_CONDITIONS_FROM_RAW.csv` — 208 condition aggregates.
- `../04_Data_And_Results/GRAPH_CLASSIFICATION_CONDITIONS.csv` — historical classifier summaries.
- `../04_Data_And_Results/ATOMIC_FORCE_RESULTS.csv` — 13 minimization records.
- `../06_Experimental_Evidence/original_logs/` — original build, native, sanitizer, regression and statistical logs.
- `../08_Executable_Source_Archive/` — original native source release, binary results, tests and their metadata.
- `../10_Provenance_And_Integrity/ORIGINAL_PATHS_AND_SHA256.csv` — all exact original paths, publication paths, sizes and digests.

## Data availability and reproducibility status

The complete computational data accompany this local release. A publicly accessible archival DOI or journal record has **not** been established by this package and must not be invented. The upstream published journal articles cited in the manuscript are *methodological context*, not experiments executed here.

## Historical manifest discrepancy discovered during this release

When the original folder structure was reconstructed from the unmodified, authoritative input ZIP, `reproduce/VERIFY_EXTENDED_PORTFOLIO.py` passed: 22 figure digests, 13,312 paired rows, 208 condition summaries, and 100 classifier rows. The next historical script, `reproduce/VERIFY_PORTFOLIO.py`, stopped with a missing-file check. Its *older* `PACKAGE_SHA256.csv` and `PUBLICATION_SHA256.csv` reference `publication_methods/__pycache__/build_research_manuscript.cpython-313.pyc`, which is **not included** in that source ZIP. This is one missing legacy generated Python bytecode entry, not a newly measured experimental failure.

A separate later manifest, `PACKAGE_SHA256_PUBLICATION.csv`, matches all 374 entries it covers in the original release. The current organized release verifies all **376 actual original files** by its own exact byte hashes and independently reconstructs the archived comparisons. The original old manifest files remain unchanged under `09_Legacy_Edition/Previous_Release_Manifests_And_Instructions/`. See `10_Provenance_And_Integrity/HISTORICAL_VERIFIER_REPLAY.log` for the actual failure output. No missing bytecode was invented or regenerated as historical evidence.
