#!/usr/bin/env python3
"""Verify recorded molecular research and public manuscript after extraction. No simulation reruns."""
from pathlib import Path
import csv,hashlib,subprocess,json
R=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
for name in ['reproduce/VERIFY_EXTENDED_PORTFOLIO.py','reproduce/VERIFY_PORTFOLIO.py']:
 subprocess.run(['python3',str(R/name)],check=True,cwd=R)
with (R/'data/PUBLICATION_FIGURE_PROVENANCE.csv').open() as f: figs=list(csv.DictReader(f))
assert len(figs)==14
for rec in figs:
 p=R/rec['figure']; assert p.exists() and sha(p)==rec['sha256'],str(p)
with (R/'PACKAGE_SHA256_PUBLICATION.csv').open() as f:rows=list(csv.DictReader(f))
for rec in rows:
 p=R/rec['path']
 assert p.exists() and int(rec['bytes'])==p.stat().st_size and sha(p)==rec['sha256'],str(p)
print('PASS publication integrity:',len(figs),'publication figures;',len(rows),'source-record files')
