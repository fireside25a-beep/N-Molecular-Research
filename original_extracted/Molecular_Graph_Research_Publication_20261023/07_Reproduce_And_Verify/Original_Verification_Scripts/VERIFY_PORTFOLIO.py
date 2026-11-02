#!/usr/bin/env python3
"""Read-only integrity test for the extracted professional portfolio."""
from pathlib import Path
import csv, hashlib, json, collections
root=Path(__file__).resolve().parent.parent
with (root/'PACKAGE_SHA256.csv').open(newline='') as f: master=list(csv.DictReader(f))
errors=[]
for row in master:
 p=root/row['relative_path']
 if not p.is_file() or p.stat().st_size!=int(row['size_bytes']): errors.append('missing or size mismatch '+row['relative_path']);continue
 h=hashlib.sha256(p.read_bytes()).hexdigest()
 if h!=row['sha256']:errors.append('hash mismatch '+row['relative_path'])
with (root/'data/FIGURE_PROVENANCE.csv').open(newline='') as f:figs=list(csv.DictReader(f))
for x in figs:
 h=hashlib.sha256((root/'figures'/x['figure']).read_bytes()).hexdigest()
 if h!=x['png_sha256']:errors.append('figure hash mismatch '+x['figure'])
with (root/'data/ALL_V9_V10_V11_PAIRED_COMPARISONS.csv').open(newline='') as f: rr=list(csv.DictReader(f))
if len(rr)!=13312:errors.append('incorrect number of paired records')
if any(int(x['mask_mismatch'])!=0 or int(x['changed_spikes'])!=int(x['mask_xor_count']) for x in rr):errors.append('row-level popcount mismatch')
for x in rr:
 if abs(float(x['fraction_recomputed'])-int(x['mask_xor_count'])/(int(x['tick_count'])*int(x['target_slots'])))>1e-10:
  errors.append('rate mismatch '+x['version']);break
print(('FAIL' if errors else 'PASS'), 'verified',len(master),'file SHA256/size values;',len(figs),'figure hashes;',len(rr),'paired comparisons')
if errors:
 print('\n'.join(errors[:20]));raise SystemExit(1)
