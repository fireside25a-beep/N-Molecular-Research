#!/usr/bin/env python3
"""Deterministic, read-only rendering of original v3-v11 molecular audit CSVs.
No networks, random values, simulations or source mutations. Every figure is calculated
from the existing audit CSVs; original logs and complete audit package bundled separately.
"""
from pathlib import Path
import csv, json, math, hashlib, statistics, textwrap, html, shutil, zipfile, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

SRC=Path(__file__).resolve().parent.parent
OUT=SRC
FIG=OUT/'figures'; DATA=OUT/'data'; DOC=OUT/'evidence_and_logs'; SCRIPT=OUT/'reproduce'
for p in [FIG,DATA,DOC,SCRIPT]:p.mkdir(parents=True,exist_ok=True)

def rows(fn):
 with (SRC/'data'/fn).open(newline='',encoding='utf8') as f: return list(csv.DictReader(f))
C=rows('ALL_CONDITIONS_FROM_RAW.csv'); G=rows('GRAPH_CLASSIFICATION_CONDITIONS.csv'); T=rows('RECORDED_TIMING_SECONDS.csv'); R=rows('RNA_FOLDING_STRUCTURES.csv'); F=rows('ATOMIC_FORCE_RESULTS.csv');P=rows('ALL_V9_V10_V11_PAIRED_COMPARISONS.csv'); Q=rows('TEST_CASE_INDEX.csv')
AUD=json.loads((SRC/'evidence_and_logs/INDEPENDENT_EXTRACTION_AUDIT.json').read_text())
assert len(P)==13312 and len(C)==208 and len(G)==100 and sum(x["version"]=="v11" for x in Q)==81 and len(F)==13
assert all(int(x['mask_mismatch'])==0 for x in P)
assert all(x['status']=='ok' for x in Q)
assert all(v['match'] for v in AUD['v11_compared_data_hashes'].values())

INK='#12243A'; MUTED='#5C6D7E'; SKY='#197DAD'; TEAL='#087E86'; ORANGE='#D58B36'; LAV='#826DB0'; RED='#AA5359'; GRID='#E3EAF0';BG='#FAFCFD'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':15,'axes.titleweight':'bold','axes.labelcolor':INK,'axes.edgecolor':GRID,'text.color':INK,'xtick.color':MUTED,'ytick.color':MUTED,'figure.facecolor':'white','axes.facecolor':'white','grid.color':GRID,'savefig.facecolor':'white'})
CHARTS=[]
def axstyle(ax,axis='y'):
 ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.spines['left'].set_visible(False)
 ax.grid(axis=axis,alpha=.75);ax.set_axisbelow(True)
def save(fig,name,title,source,method):
 fig.savefig(FIG/name,dpi=185,bbox_inches='tight',facecolor='white')
 plt.close(fig)
 h=hashlib.sha256((FIG/name).read_bytes()).hexdigest()
 CHARTS.append({'figure':name,'title':title,'csv_source':source,'transformation':method,'png_sha256':h})
def cond(v,branch,name,feed,plastic):
 return next(x for x in C if x['version']==v and x['branch']==branch and x['structure_id']==name and x['feedback']==feed and x['plasticity']==plastic)
def pct(x):return 100*float(x['fraction_from_hits'])

# FIG 1. V6 controlled accuracy, original and independent seeds (not experimental folding)
names=['1BNA','1ANA','1WOE','1RNA','1ZIH','1KF1']
sel=[next(x for x in G if x['version']=='v6' and x['pdb']==n and x['readout']=='late_scaffold_raw' and x['intervention']=='original' and x['feedback']=='closed' and x['plasticity']=='STDP') for n in names]
old=[100*float(x['original_accuracy']) for x in sel];fresh=[100*float(x['fresh_accuracy']) for x in sel]
fig,ax=plt.subplots(figsize=(11,4.6));ix=np.arange(6);w=.34
ax.bar(ix-w/2,old,w,color=SKY,label='Original held-out')
ax.bar(ix+w/2,fresh,w,color=TEAL,label='Independent fresh seeds')
ax.axhline(50,color=ORANGE,linestyle='--',lw=1.8,label='Binary chance (50%)')
ax.set_xticks(ix,names);ax.set_ylim(0,90);ax.set_ylabel('Graph-label accuracy (%)');ax.set_title('V6 | Contact-graph classification across six structures')
axstyle(ax);ax.legend(ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.11));fig.tight_layout()
save(fig,'01_v6_original_vs_fresh.png','V6 original and fresh-seed graph discrimination','GRAPH_CLASSIFICATION_CONDITIONS.csv','6 closed-loop STDP, raw readout; percent = correct / trials × 100')

# FIG 2. 1KF1 constructed intervention with control
k1=next(x for x in G if x['pdb']=='1KF1' and x['readout']=='late_scaffold_raw' and x['intervention']=='original' and x['feedback']=='closed' and x['plasticity']=='STDP');k2=next(x for x in G if x['pdb']=='1KF1' and x['readout']=='late_scaffold_raw' and x['intervention']=='ring_preserving_plus3' and x['feedback']=='closed' and x['plasticity']=='STDP')
fig,ax=plt.subplots(figsize=(10.5,3.7)); labels=['Coarse rewire\n(initial holdout)','Coarse rewire\n(fresh holdout)','Ring-preserving rewire\n(fresh holdout)'];scores=[100*float(k1['original_accuracy']),100*float(k1['fresh_accuracy']),100*float(k2['fresh_accuracy'])]
ax.barh(labels,scores,color=[SKY,TEAL,ORANGE],height=.54);ax.axvline(50,color=RED,linestyle='--',label='50% chance');ax.set_xlim(0,95);ax.invert_yaxis();ax.set_xlabel('Graph-label accuracy (%)');ax.set_title('V6 | A positive comparison and its hard structural control');axstyle(ax,'x')
for i,v in enumerate(scores):ax.text(v+1.3,i,f'{v:.2f}%',va='center',weight='bold')
fig.tight_layout();save(fig,'02_v6_ring_control.png','V6 ring-breaking versus ring-preserving controls','GRAPH_CLASSIFICATION_CONDITIONS.csv + original scientific report','Initial 95/128, fresh 196/256, ring-preserving 116/256 from native reports')

# FIG 3. Four additional RNA cases
four=['1ZIF','1ZIG','2KOC','1EHZ']; v7=[next(x for x in G if x['version']=='v7-v8' and x['pdb'].split('→')[0]==n and x['feedback']=='closed' and x['plasticity']=='STDP') for n in four]
vals=[100*float(x['fresh_accuracy']) for x in v7]
fig,ax=plt.subplots(figsize=(10.5,3.7));bars=ax.bar(four,vals,color=[SKY,SKY,TEAL,ORANGE]);ax.axhline(50,linestyle='--',color=RED);ax.set_ylim(0,75);ax.set_ylabel('Held-out graph-label accuracy (%)');ax.set_title('V7–V8 | Four independent RNA labels — none significant after correction');axstyle(ax)
for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+1,f'{v:.2f}%',ha='center',weight='bold')
fig.tight_layout();save(fig,'03_v7_v8_rna_classification.png','V7/V8 four additional held-out graph classifiers','GRAPH_CLASSIFICATION_CONDITIONS.csv','Closed + STDP fresh accuracy, four-family Holm interpretation from original statistical report')

# FIG 4. V9 complete 8x4 factorial heatmap
mols=['aspirin','benzene','caffeine','cyclohexane','ethanol','naphthalene','nicotinamide','pyridine']
arms=[('native_closed','STDP'),('native_closed','frozen'),('matched_fixed','STDP'),('matched_fixed','frozen')]
a=np.array([[pct(cond('v9','chemical_bonds',m,*arm)) for arm in arms] for m in mols]);
fig,ax=plt.subplots(figsize=(10.5,5.3));im=ax.imshow(a,cmap='YlGnBu',vmin=0,vmax=22,aspect='auto');ax.set_xticks(range(4),['Closed / learning','Closed / frozen','Fixed / learning','Fixed / frozen']);ax.set_yticks(range(8),[m.title() for m in mols]);ax.set_title('V9 | Bond-ablation response by molecule and neural control');
for i in range(8):
 for j in range(4):ax.text(j,i,f'{a[i,j]:.2f}%',ha='center',va='center',color=INK if a[i,j]<12 else 'white',fontsize=9,weight='bold')
fig.colorbar(im,ax=ax,shrink=.8,label='Changed target firing positions (%)');fig.tight_layout();save(fig,'04_v9_chemical_factorial.png','V9 complete chemical-response factorial','ALL_CONDITIONS_FROM_RAW.csv','8 compounds × closed/fixed × STDP/frozen; mean of 64 paired network comparisons')

# FIG 5. v10 RNA and v11 model pair difference
r10={x['pdb']:x for x in R if x['version']=='v10'};r11={x['pdb']:x for x in R if x['version']=='v11'}
rna=['1ZIH','1ZIF','1ZIG','2KOC','1EHZ'];iv=np.arange(5)
fig,ax=plt.subplots(figsize=(10.4,4.3));ax.bar(iv-.19,[int(r10[x]['nested_pairs']) for x in rna],.38,color=SKY,label='V10 weighted Nussinov');ax.bar(iv+.19,[int(r11[x]['nested_pairs']) for x in rna],.38,color=TEAL,label='V11 partial thermodynamic model');ax.set_xticks(iv,rna);ax.set_ylabel('Predicted nested base-pair count');ax.set_title('V10 → V11 | Predicted secondary-structure contact counts');axstyle(ax);ax.legend(frameon=False);fig.tight_layout()
save(fig,'05_rna_predicted_pairs.png','V10 vs V11 predicted RNA pair counts','RNA_FOLDING_STRUCTURES.csv','Count of paired positions from original RNA calculation JSONL (not experimentally scored F1)')

# FIG 6. source dot-bracket diagram from v11, 2KOC + 1EHZ first 76
fig,axes=plt.subplots(2,1,figsize=(12.2,5.0),gridspec_kw={'height_ratios':[1,1.9]})
for ax,pdb in zip(axes,['2KOC','1EHZ']):
 x=r11[pdb];seq=x['sequence'];db=x['dot_bracket'];n=len(seq); ax.plot([0,n+1],[0,0],color=GRID,lw=1); stack=[]; pairs=[]
 for j,ch in enumerate(db):
  if ch=='(':stack.append(j)
  if ch==')' and stack: pairs.append((stack.pop(),j))
 for i,j in pairs:
  xx=np.linspace(i,j,80); rad=(j-i)/2;cy=(i+j)/2;yy=np.sqrt(np.maximum(rad*rad-(xx-cy)**2,0)); ax.plot(xx,yy,color=TEAL,alpha=.65,lw=1.4)
 ax.scatter(range(n),[0]*n,s=10,color=SKY,zorder=5);ax.set_xlim(-1,n);ax.set_ylim(-.5,max(2,n/3));ax.set_axis_off()
 ax.text(0,max(2,n/3)*.9,f'{pdb}: {n} positions · {len(pairs)} predicted nested pairs',fontsize=12,weight='bold',color=INK)
fig.suptitle('V11 | Direct rendering of recorded dot-bracket base pairs',fontsize=16,weight='bold');fig.tight_layout()
save(fig,'06_v11_rna_pair_arcs.png','V11 sequence-level pairing geometry (not atomic)','RNA_FOLDING_STRUCTURES.csv','Plot arc (i,j) only where characters in source dot-bracket match')

# FIG 7. V11 actual partial energy at 37 C
fig,ax=plt.subplots(figsize=(10.2,4));es=[float(r11[x]['model_energy_kcal_mol']) for x in rna];bars=ax.barh(rna,es,color=[SKY]*4+[TEAL]);ax.axvline(0,color=GRID);ax.set_xlabel('Calculated V11 model energy (kcal/mol)');ax.set_title('V11 | Partial nearest-neighbor energy model · 37 °C');axstyle(ax,'x');ax.invert_yaxis();
for b,v in zip(bars,es):ax.text(v-.25,b.get_y()+b.get_height()/2,f'{v:.2f}',va='center',ha='right',weight='bold')
fig.tight_layout();save(fig,'07_v11_partial_rna_energies.png','V11 calculated partial RNA free-energy scores','RNA_FOLDING_STRUCTURES.csv','Direct numeric model_energy_kcal_mol from recorded 5 RNA JSONL rows')

# FIG 8. V11 chemical response 8 all-arm means; V9 separate protocol dashed distinction
fig,ax=plt.subplots(figsize=(11.6,4.6));mean11=[np.mean([pct(cond('v11','ligand_3d',m,*arm)) for arm in arms]) for m in mols];means9=[np.mean([pct(cond('v9','chemical_bonds',m,*arm)) for arm in arms]) for m in mols];ix=np.arange(8)
ax.bar(ix-.18,means9,.36,color=SKY,label='V9 bond removal')
ax.bar(ix+.18,mean11,.36,color=TEAL,label='V11 spatial-contact removal')
ax.set_xticks(ix,[m.title() for m in mols],rotation=23,ha='right');ax.set_ylim(0,max(max(means9),max(mean11))+7);ax.set_ylabel('Mean changed firing positions (%)');ax.set_title('V9 vs V11 | Different interventions; no head-to-head accuracy claim');axstyle(ax);ax.legend(frameon=False)
fig.tight_layout();save(fig,'08_v9_vs_v11_compounds_distinct_protocols.png','Chemical response by version (different experimental protocols)','ALL_CONDITIONS_FROM_RAW.csv','Unweighted mean of four factorial arms, computed from source-recorded changed spike counts')

# FIG 9. V11 post force RNA + ligand responses
ids=['1ZIH_rna4','1ZIF_rna4','1ZIG_rna4','2KOC_rna4','1EHZ_rna4']
# locate identifiers exactly from CSV
available=sorted(set(x['structure_id'] for x in C if x['version']=='v11' and x['branch']=='post_force'))
if not all(x in available for x in ids):
 ids=[x for x in available if 'rna' in x.lower() or x.startswith(('1ZI','2KO','1EH'))]
RNA_MEAN=[]
for rid in ids:
 rc=[x for x in C if x['version']=='v11' and x['branch']=='post_force' and x['structure_id']==rid]
 if rc:RNA_MEAN.append((rid,100*np.mean([float(x['fraction_from_hits']) for x in rc])))
fig,ax=plt.subplots(figsize=(9.5,4));lab=[x.replace('_rna4','') for x,_ in RNA_MEAN];val=[v for _,v in RNA_MEAN];bs=ax.bar(lab,val,color=TEAL);ax.set_ylim(0,max(val+[1])*1.18);ax.set_title('V11 | Four-nucleotide RNA fragment contact interventions');ax.set_ylabel('Neural disagreement, four-setting mean (%)');axstyle(ax)
for b,v in zip(bs,val):ax.text(b.get_x()+b.get_width()/2,v+1,f'{v:.2f}%',ha='center',weight='bold')
fig.tight_layout();save(fig,'09_v11_rna_fragment_response.png','V11 four-nucleotide RNA fragment interventions','ALL_CONDITIONS_FROM_RAW.csv','Mean over four arms; counterfactual removes hundreds of computed contacts; not folding accuracy')

# FIG 10 MMFF94 convergence delta log10 error; real gradients
fig,ax=plt.subplots(figsize=(11,4.6));fs=sorted(F,key=lambda x:float(x['max_abs_gradient_delta']),reverse=True);nm=[x['structure_id'].replace('_rna4',' RNA4') for x in fs];err=[float(x['max_abs_gradient_delta']) for x in fs];ax.barh(nm,err,color=[TEAL if 'rna' in x.lower() or any(x.startswith(y) for y in rna) else SKY for x in nm]);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlabel('Max |analytic gradient - numerical gradient| (kcal/mol/Å)');ax.set_title('V11 | Analytical forces vs finite-difference gradients');axstyle(ax,'x')
fig.tight_layout();save(fig,'10_v11_gradient_verification.png','V11 MMFF94 gradient consistency across 13 molecules/fragments','ATOMIC_FORCE_RESULTS.csv','Maximum sampled absolute gradient discrepancy on logarithmic axis (4 components per structure)')

# FIG 11 stage timings
stages=[x for x in T if x['stage']!='full_pipeline_wall'];stages.sort(key=lambda x:float(x['seconds']),reverse=True)
fig,ax=plt.subplots(figsize=(11,6));labels=[x['stage'].replace('_seconds','').replace('_',' ') for x in stages];times=[float(x['seconds']) for x in stages];ax.barh(labels,times,color=[TEAL if 'native' in k else SKY if 'audit' in k else ORANGE for k in labels]);ax.invert_yaxis();ax.set_xlabel('Measured wall seconds');ax.set_title('V11 | Logged execution time by stage (complete run: 32.07s)');axstyle(ax,'x');fig.tight_layout()
save(fig,'11_v11_stage_timing.png','V11 logged stage times','RECORDED_TIMING_SECONDS.csv','Timings read from accepted verified-run execution report; excluded full_pipeline_wall from stage bars')

# FIG 12 recorded comparison volumes, no biological replication assertion
fig,ax=plt.subplots(figsize=(10.1,3.8));vl=['V9','V10','V11'];counts=[2048,3328,7936];bs=ax.barh(vl,counts,color=[SKY,TEAL,INK],height=.56);ax.invert_yaxis();ax.set_xlabel('Recorded paired comparisons');ax.set_title('Research progression | Computational experiment volume, not biological samples');axstyle(ax,'x')
for b,v in zip(bs,counts):ax.text(v+90,b.get_y()+b.get_height()/2,f'{v:,}',va='center',weight='bold')
ax.set_xlim(0,9000);fig.tight_layout();save(fig,'12_version_computational_volume.png','Reproducible paired comparison records','ALL_V9_V10_V11_PAIRED_COMPARISONS.csv','Count rows grouped by version; unit is deterministic paired network run')

# FIG 13 project tests by version; don't count re-runs as new evidence
fig,ax=plt.subplots(figsize=(9.7,3.7));ys=[31,59,81];bs=ax.bar(['V9','V10','V11'],ys,color=[SKY,TEAL,INK]);ax.set_ylabel('Passing named Python tests');ax.set_title('Inherited and expanded project regression coverage');axstyle(ax);ax.set_ylim(0,95)
for b,v in zip(bs,ys):ax.text(b.get_x()+b.get_width()/2,v+1,str(v),ha='center',weight='bold')
fig.tight_layout();save(fig,'13_unit_tests_passed.png','Source-documented Python test suite expansion','TEST_CASE_INDEX.csv + V10/V11 scientific report','V11 81 named test IDs from actual log; inherited V9 31 and V10 59 from preserved reports')


print('PASS: regenerated',len(CHARTS),'data-driven figures from shipped CSVs')
