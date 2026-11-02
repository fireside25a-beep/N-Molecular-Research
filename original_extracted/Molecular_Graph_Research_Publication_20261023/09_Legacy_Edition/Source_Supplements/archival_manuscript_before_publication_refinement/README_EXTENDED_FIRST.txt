N MOLECULAR RESEARCH PORTFOLIO V3–V11 — EXTENDED SOURCE EDITION

START HERE
- MOLECULAR_RESEARCH_PORTFOLIO_V3-V11_EXTENDED_LATEX.pdf : new compiled scientific LaTeX report with full chart appendix
- MOLECULAR_RESEARCH_PORTFOLIO_V3-V11.pdf : original 18-page professional portfolio, preserved
- INDEX.html : offline portfolio with added figures
- COMPLETE_AUDIT_WORKBOOK.xlsx : 17-tab evidence workbook
- figures/ : original 13 charts drawn from computation records
- figures_extra/ : additional charts from archived data, with external literature chart separately labelled
- data/ : all source-derived CSVs, chart provenance, literature reference table
- latex/portfolio_extended.tex : LaTeX source of new report
- reproduce/VERIFY_EXTENDED_PORTFOLIO.py : quick file/hash and count audit
- reproduce/VERIFY_PORTFOLIO.py : original portfolio validation routine
- evidence_and_logs/ : historical run logs, independent verifier and evidence bundle
- prior_portfolios/ : byte-preserved earlier portfolio ZIPs and evidence archive
- proof_of_build/ : LaTeX compile logs and reproducibility audit

SCIENTIFIC SCOPE
The plotted results are recorded computational experiments, not biological binding or physical folding accuracy.
External literature scores are from their own published datasets, shown separately.

FULL NATIVE FRAMEWORK
- native_source_release/: original V11 source, binaries, scripts and historic logs in unchanged verified ZIP.

VERIFY FIGURES
Run: python3 reproduce/VERIFY_EXTENDED_PORTFOLIO.py

LATEX REBUILD
From latex/: pdflatex -halt-on-error portfolio_extended.tex (twice).
