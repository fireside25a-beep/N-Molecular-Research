#!/usr/bin/env python3
"""Read-only publication verification; numerical audits only for recorded outputs."""
import csv,json,hashlib,sys,collections
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()
with (R/'10_Provenance_And_Integrity/ORIGINAL_PATHS_AND_SHA256.csv').open(newline='') as f:mapping=list(csv.DictReader(f))
for r in mapping:
 p=R/r['publication_path']
 if not p.is_file() or p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:
  raise SystemExit('FILE INTEGRITY FAILURE '+r['publication_path'])
D=R/'04_Data_And_Results'
with (D/'ALL_V9_V10_V11_PAIRED_COMPARISONS.csv').open(newline='') as f:pairs=list(csv.DictReader(f))
assert len(pairs)==13312
assert all(int(r['mask_mismatch'])==0 for r in pairs),'mask mismatch in archived pair'
for r in pairs:
 assert int(r['changed_spikes'])==int(r['mask_xor_count'])
 assert abs(float(r['fraction_recomputed'])-int(r['changed_spikes'])/(int(r['tick_count'])*int(r['target_slots'])))<1e-12
 assert abs(float(r['fraction_recorded'])-float(r['fraction_recomputed']))<1e-12
with (D/'ALL_CONDITIONS_FROM_RAW.csv').open(newline='') as f:conditions=list(csv.DictReader(f))
assert len(conditions)==208
buckets=collections.defaultdict(list)
for r in pairs:buckets[(r['version'],r['branch'],r['structure_id'],r['feedback'],r['plasticity'])].append(r)
for r in conditions:
 k=(r['version'],r['branch'],r['structure_id'],r['feedback'],r['plasticity'])
 group=buckets[k]
 assert len(group)==int(r['paired_rows'])
 assert sum(int(q['changed_spikes']) for q in group)==int(r['total_changed_spikes'])
 denom=sum(int(q['tick_count'])*int(q['target_slots']) for q in group)
 assert abs(int(r['total_changed_spikes'])/denom-float(r['fraction_from_hits']))<1e-12
with (D/'GRAPH_CLASSIFICATION_CONDITIONS.csv').open(newline='') as f:cla=list(csv.DictReader(f))
assert len(cla)==100
for r in cla:
 if r['correct'] and r['total'] and r['fresh_accuracy']:
  assert abs(int(r['correct'])/int(r['total'])-float(r['fresh_accuracy']))<1e-12
 if r['null_exceedances'] and r['permutations'] and r['permutation_p']:
  assert abs((int(r['null_exceedances'])+1)/(int(r['permutations'])+1)-float(r['permutation_p']))<1e-12
with (D/'PUBLICATION_FIGURE_PROVENANCE.csv').open(newline='') as f:fig=list(csv.DictReader(f))
assert len(fig)==14
original_to_current={r['original_path']:r['publication_path'] for r in mapping}
for r in fig:
 p=R/original_to_current[r['figure']]
 assert sha(p)==r['sha256'], r['figure']
with (D/'TEST_CASE_INDEX.csv').open(newline='') as f:tests=list(csv.DictReader(f))
assert sum(t['version']=='v11' and t['status']=='ok' for t in tests)==81
print('PASS: original bytes',len(mapping),'paired masks',len(pairs),'condition summaries',len(conditions),'classification summaries',len(cla),'publication figures',len(fig),'historical recorded tests 81')
