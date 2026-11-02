from pathlib import Path
import csv,json,hashlib
root=Path(__file__).resolve().parents[1]
data=root/'data'
with (data/'FIGURE_SHA256_EXTENDED.csv').open() as f: rows=list(csv.DictReader(f))
for r in rows:
    p=root/r['image_file']
    assert p.is_file(),p
    assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
with (data/'ALL_V9_V10_V11_PAIRED_COMPARISONS.csv').open() as f:
    n=sum(1 for _ in f)-1
assert n==13312,n
with (data/'ALL_CONDITIONS_FROM_RAW.csv').open() as f:
    n=sum(1 for _ in f)-1
assert n==208,n
with (data/'GRAPH_CLASSIFICATION_CONDITIONS.csv').open() as f:
    n=sum(1 for _ in f)-1
assert n==100,n
print(f'PASS: {len(rows)} source-derived figure hashes, 13,312 paired rows, 208 conditions, 100 classifier rows')
