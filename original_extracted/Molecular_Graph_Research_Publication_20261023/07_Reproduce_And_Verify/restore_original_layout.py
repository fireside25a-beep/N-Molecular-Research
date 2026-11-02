#!/usr/bin/env python3
"""Copy the exact archived bytes into their previous source-relative paths for legacy scripts."""
import csv,shutil,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
if len(sys.argv)!=2: raise SystemExit('Usage: python3 restore_original_layout.py UNUSED_OUTPUT_DIRECTORY')
T=Path(sys.argv[1]).expanduser().resolve()
if T.exists(): raise SystemExit('Refusing existing directory: '+str(T))
T.mkdir(parents=True)
with (R/'10_Provenance_And_Integrity/ORIGINAL_PATHS_AND_SHA256.csv').open(newline='') as f:rows=list(csv.DictReader(f))
for r in rows:
 dst=T/r['original_path'];dst.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(R/r['publication_path'],dst)
print('Reconstructed',len(rows),'historical files in',T)
