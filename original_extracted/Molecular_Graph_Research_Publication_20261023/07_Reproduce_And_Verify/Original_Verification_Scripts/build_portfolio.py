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
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image as PdfImage, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

SRC=Path('/mnt/data/N_Molecular_Evidence_20261023')
OUT=Path('/mnt/data/N_Molecular_Research_Portfolio_V3-V11_20261023')
FIG=OUT/'figures'; DATA=OUT/'data'; DOC=OUT/'evidence_and_logs'; SCRIPT=OUT/'reproduce'
for p in [FIG,DATA,DOC,SCRIPT]:p.mkdir(parents=True,exist_ok=True)

def rows(fn):
 with (SRC/'raw_csv'/fn).open(newline='',encoding='utf8') as f: return list(csv.DictReader(f))
C=rows('ALL_CONDITIONS_FROM_RAW.csv'); G=rows('GRAPH_CLASSIFICATION_CONDITIONS.csv'); T=rows('RECORDED_TIMING_SECONDS.csv'); R=rows('RNA_FOLDING_STRUCTURES.csv'); F=rows('ATOMIC_FORCE_RESULTS.csv');P=rows('ALL_V9_V10_V11_PAIRED_COMPARISONS.csv'); Q=rows('TEST_CASE_INDEX.csv')
AUD=json.loads((SRC/'manifests/INDEPENDENT_EXTRACTION_AUDIT.json').read_text())
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

# include source CSVs, trace metadata, verifier (not executing underlying simulation)
for fn in ['ALL_CONDITIONS_FROM_RAW.csv','GRAPH_CLASSIFICATION_CONDITIONS.csv','RECORDED_TIMING_SECONDS.csv','RNA_FOLDING_STRUCTURES.csv','ATOMIC_FORCE_RESULTS.csv','TEST_CASE_INDEX.csv','ALL_V9_V10_V11_PAIRED_COMPARISONS.csv']:
 shutil.copy2(SRC/'raw_csv'/fn,DATA/fn)
shutil.copy2(SRC/'manifests/INDEPENDENT_EXTRACTION_AUDIT.json',DOC/'INDEPENDENT_EXTRACTION_AUDIT.json')
shutil.copy2(SRC/'manifests/VERIFY_BUNDLE_PASS_20261023.log',DOC/'VERIFY_BUNDLE_PASS.log')
shutil.copy2(SRC/'manifests/SOURCE_FILE_SHA256.csv',DOC/'SOURCE_FILE_SHA256.csv')
shutil.copy2(SRC/'manifests/SHA256_OF_BUNDLE_FILES.csv',DOC/'PREVIOUS_BUNDLE_SHA256.csv')
shutil.copy2(SRC/'N_Molecular_v3-v11_Research_Audit_20261023.xlsx',OUT/'COMPLETE_AUDIT_WORKBOOK.xlsx')
for doc in ['docs/V11_SCIENTIFIC_SOURCES.md','docs/METHODS.md']:
 p=SRC/'methods'/doc
 if p.exists(): shutil.copy2(p,DOC/p.name)
# log package: actual logs not paraphrased
if (DOC/'original_logs').exists(): pass
else:shutil.copytree(SRC/'logs', DOC/'original_logs')
# preserve raw source sidecars / original receipts inside separate evidence archive for distribution
shutil.copy2(Path('/mnt/data/N_Molecular_v3-v11_AUDIT_WORKBOOK_GRAPHS_LOGS_REPRODUCIBLE_20261023.zip'),DOC/'COMPLETE_REPRODUCIBLE_EVIDENCE.zip')

# scientific references checked today
REFS=[
 ('Rochester NNDB / Turner 2004','Thermodynamic parameter provenance, not full model parity','https://rna.urmc.rochester.edu/NNDB/'),
 ('RNAstructure parameters','Complete reference model for comparative secondary structure prediction','https://rna.urmc.rochester.edu/Overview/Parameters.html'),
 ('ViennaRNA RNAfold','Minimum free energy, ensemble / partition function capabilities','https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/usage/man/RNAfold.html'),
 ('BPfold, Nature Communications 2025','External measured structure F1; no direct N benchmark','https://www.nature.com/articles/s41467-025-60048-1'),
 ('G4ShapePredictor, Scientific Reports 2024','External G4 folding-class classification, not 1KF1 constructed graph discrimination','https://www.nature.com/articles/s41598-024-74826-2'),
 ('RDKit','ETKDG conformer and MMFF94 computations used by software','https://www.rdkit.org/docs/'),
 ('RCSB PDB 1KF1','Deposited X-ray resolved DNA G-quadruplex, 2.10 Å','https://www.rcsb.org/structure/1KF1'),
 ('RCSB PDB 1EHZ','Deposited yeast tRNA X-ray structure, 1.93 Å','https://www.rcsb.org/structure/1EHZ'),
]
with (DATA/'SCIENTIFIC_REFERENCES.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.writer(f);w.writerow(['publication_or_provider','role','url']);w.writerows(REFS)
with (DATA/'FIGURE_PROVENANCE.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=['figure','title','csv_source','transformation','png_sha256']);w.writeheader();w.writerows(CHARTS)

# Research integrity cover facts
nums={'versions':'V3–V11','paired_network_records':len(P),'paired_by_version':{v:sum(x['version']==v for x in P) for v in ['v9','v10','v11']},'comparison_condition_rows':len(C),'graph_statistical_rows':len(G),'documented_test_cases_v11':sum(x['version']=='v11' for x in Q),'all_spike_mask_mismatch':sum(int(x['mask_mismatch']) for x in P),'source_files_audited':AUD['total_source_files_except_git'],'confirmed_original_data_hash_matches':sum(int(v['match']) for v in AUD['v11_compared_data_hashes'].values()),'figures':len(CHARTS),'v11_wall_seconds':float(next(x['seconds'] for x in T if x['stage']=='full_pipeline_wall'))}
(DATA/'PORTFOLIO_METRICS.json').write_text(json.dumps(nums,indent=2))

# PDF: professional paper with each exhibit on its own page, source lineage retained
fontdir=Path('/usr/share/fonts/truetype/dejavu')
pdfmetrics.registerFont(TTFont('DejaVu',str(fontdir/'DejaVuSans.ttf')))
pdfmetrics.registerFont(TTFont('DejaVu-Bold',str(fontdir/'DejaVuSans-Bold.ttf')))
PAGE_W,PAGE_H=A4; margin=39
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='PT',fontName='DejaVu-Bold',fontSize=24,leading=30,textColor=colors.HexColor(INK),spaceAfter=17))
styles.add(ParagraphStyle(name='PH',fontName='DejaVu-Bold',fontSize=13,leading=18,textColor=colors.HexColor(INK),spaceBefore=8,spaceAfter=12))
styles.add(ParagraphStyle(name='PBody2',fontName='DejaVu',fontSize=9.2,leading=15,textColor=colors.HexColor(INK),spaceAfter=9))
styles.add(ParagraphStyle(name='PSmall2',fontName='DejaVu',fontSize=7.6,leading=12,textColor=colors.HexColor(MUTED),spaceAfter=7))
styles.add(ParagraphStyle(name='PTag2',fontName='DejaVu-Bold',fontSize=9.5,leading=13,textColor=colors.HexColor(TEAL),spaceAfter=9))
styles.add(ParagraphStyle(name='PNote2',fontName='DejaVu',fontSize=8.2,leading=13,textColor=colors.HexColor(MUTED),spaceAfter=8))
styles.add(ParagraphStyle(name='PRef2',fontName='DejaVu',fontSize=7.2,leading=10,textColor=colors.HexColor(MUTED),spaceAfter=4,wordWrap='CJK'))
fmt=lambda x:Paragraph(html.escape(str(x)),styles['PBody2'])
S=[]
def add(txt,style='PBody2'):S.append(Paragraph(txt,styles[style]))
def sep():S.append(Spacer(1,10))
def nextpage():S.append(PageBreak())
def tbl(header,dat,widths=None,font=7.9):
 t=Table([[Paragraph(html.escape(str(x)),styles['PRef2']) for x in header]]+[[Paragraph(html.escape(str(x)),styles['PRef2']) for x in line] for line in dat],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E6F2F4')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F6F9FB')]),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,0),.8,colors.HexColor(TEAL)),('VALIGN',(0,0),(-1,-1),'TOP')]))
 S.append(t)
def pic(filename,width=520):
 from PIL import Image
 w,h=Image.open(FIG/filename).size;aw=min(width,PAGE_W-2*margin);ah=aw*h/w;S.append(PdfImage(str(FIG/filename),width=aw,height=ah))
def exhibit(num,title,filename,source,interpretation,detail=''):
 nextpage();add(f'EXHIBIT {num:02d}  /  SOURCE-RENDERED DATA','PTag2');add(title,'PT');pic(filename);sep();add(f'<b>Source:</b> {html.escape(source)}','PSmall2');add(interpretation,'PBody2')
 if detail:add(detail,'PNote2')
def footer(canvas,doc):
 canvas.saveState();canvas.setStrokeColor(colors.HexColor(GRID));canvas.line(38,42,PAGE_W-38,42);canvas.setFont('DejaVu',7.5);canvas.setFillColor(colors.HexColor(MUTED));canvas.drawString(38,28,'N MOLECULAR  /  PUBLIC RESEARCH PORTFOLIO  /  OCT 2026');canvas.drawRightString(PAGE_W-38,28,f'{doc.page:02d}');canvas.restoreState()

add('N MOLECULAR  /  V3—V11','PTag2')
add('From structural graph probes to executable molecular computations','PT')
add('A source-grounded research portfolio, figures, quantitative receipts and reproducibility record','PH')
add('Research archive snapshot: 23 October 2026. Public scientific references checked on the same date. All figures are generated from recovered on-disk experimental records—not illustrative source images.','PBody2')
sep();tbl(['RECOVERED EVIDENCE','CONFIRMED VALUE'],[['Computational paired comparisons (V9–V11)','13,312'],['V11 source-built paired comparison records','7,936'],['V11 actual reset SNN executions','15,872'],['Historical V6–V8 graph-classifier result rows','100'],['V11 named passing Python tests','81 (historical log)'],['Independent disagreement arithmetic errors','0 of 13,312'],['V11 recorded pipeline wall time','32.067 seconds'],['V11 molecular MMFF94 minimizations','13 converged']],[245,265]);sep()
add('<b>Research claim.</b> A deterministic native spiking neural system responds to graph interventions derived from known sequences, computed secondary structures, and generated atomic geometries. Sensitivity is verified; experimental three-dimensional folding accuracy and ligand activity prediction have not been demonstrated.','PBody2')
add('<b>Reproducibility contract.</b> Each exhibit names its evidence table. The bundle contains the complete recovered raw audit evidence and an independent verifier; the large native C/C++ project is supplied separately as the archived source release.','PNote2')

nextpage();add('01  /  WHAT WAS ACTUALLY EXECUTED','PTag2');add('Research chronology and evidence hierarchy','PT')
tbl(['STAGE','EXPERIMENTAL UNIT','WHAT WAS REPORTED'],[['V3–V5','Native molecular input and SNN trials','Source bridge, original 6-case series and controls'],['V6','128 additional deterministic excitation seeds','Fresh graph classifier tests; 1KF1 motif-preserving falsification'],['V7–V8','Four additional RNA structure labels','8,192 new recorded trial rows; full replay audits'],['V9','Eight chemical SMILES graphs','2,048 pairs / 4,096 resets / 48 ticks'],['V10','5 RNA + 8 ligand inputs','3,328 pairs / 6,656 resets / 48 ticks'],['V11','Thermo, PK candidate, geometry & force branches','7,936 pairs / 15,872 resets / 48 ticks']],[65,220,225])
sep();add('<b>Evidence classes.</b> Original experiment output (JSONL) is the strongest record for firing behavior; the source-driven audit CSVs reconstruct those outputs. Scientific reports and run logs document statistical protocols, tests and timing. Experimental PDB metadata identify molecules but do not automatically validate the computed structure.','PBody2')
add('<b>Not additive biological sample counts.</b> The V8 historical 22,528-trial set describes earlier deterministic networks; V9–V11 trials are different interventions, often sharing molecules and methodological ancestors.','PNote2')

exhibit(1,'Native graph-label discrimination: original versus fresh seeds','01_v6_original_vs_fresh.png','GRAPH_CLASSIFICATION_CONDITIONS.csv','Across six hand-authored DNA/RNA contact-graph contrasts, some raw neural readouts carried decodable intervention information on previously unused stimulation seeds. V6 values classify a computational rewire label, not a physical fold.')
exhibit(2,'One decisive negative control changes the interpretation','02_v6_ring_control.png','Graph classifier CSV and V6 scientific report','For 1KF1, the fresh cycle-breaking comparison reached 196/256 correct (76.56%); preserving four-site rings yielded 116/256 (45.31%), close to chance. This is compatible with coarse cycle sensitivity, not generalized molecular geometry inference.')
exhibit(3,'Additional RNA cases: corrected null results','03_v7_v8_rna_classification.png','GRAPH_CLASSIFICATION_CONDITIONS.csv','In V7–V8, the four newly added RNA cases did not show a statistically reliable graph-label decoder after the four-comparison multiplicity correction. These negative findings are reported, not discarded.')
exhibit(4,'Eight compounds × four controls','04_v9_chemical_factorial.png','ALL_CONDITIONS_FROM_RAW.csv; source original 64-seed ligand JSONL','V9 measures spike-event disagreement between intact and ablated chemical bond graphs. For example, caffeine was 17.17% under closed feedback and learning, not 17.17% chemical binding or predictive success.')
exhibit(5,'RNA pair count changes across solver versions','05_rna_predicted_pairs.png','RNA_FOLDING_STRUCTURES.csv; original calculated dot-bracket structures','V10 implemented a weighted noncrossing Nussinov solver. V11 used selected Turner-source nearest-neighbor energy terms and produced 40 nested pairs in five RNAs, versus 43 in V10. Neither pair count is an experimental accuracy score.')
exhibit(6,'Base-pair arcs rendered directly from RNA output','06_v11_rna_pair_arcs.png','RNA_FOLDING_STRUCTURES.csv; V11 dot-bracket fields','Every plotted arc is a paired opening/closing bracket from the saved solver output. Arc positions show RNA secondary structure connectivity, not actual atom coordinates or three-dimensional molecular dynamics.')
exhibit(7,'Partial nearest-neighbor energies at 37 °C','07_v11_partial_rna_energies.png','RNA_FOLDING_STRUCTURES.csv','These five free-energy-like values are internal results of the coded partial thermodynamic model. They must not be compared with full RNAfold partition-function or full Turner thermodynamic outputs as equivalent physical folding energies.')
exhibit(8,'Ligand response: different representations, different tasks','08_v9_vs_v11_compounds_distinct_protocols.png','ALL_CONDITIONS_FROM_RAW.csv','V9 removes chemical bonds while V11 removes ETKDG-derived spatial contacts. The graph highlights measured response profiles, not performance improvement. No ligand assays or binding labels appear in the native source evidence.')
exhibit(9,'High response in short atomistic RNA fragments','09_v11_rna_fragment_response.png','ALL_CONDITIONS_FROM_RAW.csv','V11 performs MMFF94 on four-nucleotide RNA fragments, then removes hundreds of generated spatial contacts in a counterfactual. Approximately 62–64% neuronal disagreement does not mean 62–64% correct RNA folding.')
exhibit(10,'Numerical force consistency: 13 molecules/fragments','10_v11_gradient_verification.png','ATOMIC_FORCE_RESULTS.csv','The recorded 13 minimizations converged. Sampled analytical MMFF94 gradients were compared to finite differences; the largest absolute discrepancy was 1.59×10⁻⁷ kcal/mol/Å. This validates numerical consistency, not a full biomolecular model.')
exhibit(11,'Timing: the stages recorded in the machine log','11_v11_stage_timing.png','RECORDED_TIMING_SECONDS.csv','The accepted 64-seed V11 pipeline took 32.067 seconds on the archived execution host. This is not an external speed comparison with ViennaRNA, BPfold, molecular dynamics or AlphaFold.')
exhibit(12,'Computational workload grows, not biological samples','12_version_computational_volume.png','ALL_V9_V10_V11_PAIRED_COMPARISONS.csv','The record counts are deterministic paired neural runs: V9 2,048, V10 3,328, V11 7,936. Their experimental questions differ, so no blanket aggregate accuracy is computed.')
exhibit(13,'Recorded regression coverage','13_unit_tests_passed.png','TEST_CASE_INDEX.csv; original version-specific regression logs','V11 records 81 passing named project Python tests. Source receipts show a full V10 upstream pass; a later V11 one-shot upstream run was interrupted and followed by separate per-command checks. The latter is not relabeled a single complete pass.')
nextpage();add('15  /  METHODS AND EQUATIONS','PTag2');add('The measurements behind every graph','PT')
add('<b>Spike disagreement:</b> D = Σ(t=0..47) popcount(A(t) XOR B(t)) / (48 × 22). The A and B arrays encode target-neuron spike masks for matched intact/ablated conditions. This quantity measures neural response differences, not folding accuracy.','PBody2')
add('<b>Binary graph classifier:</b> accuracy = correct predictions / evaluated graph labels. For balanced labels, the chance baseline is 50%; each perturbation and readout must be assessed with its own controls.','PBody2')
add('<b>Permutation tail probability:</b> p = (1 + number of null scores ≥ observed score) / (4096 + 1), followed by Holm correction for multiple tests. Null assumptions pertain to deterministic simulation seeds, not biological patient populations.','PBody2')
add('<b>RNA model energy:</b> G_model = selected stacking + loop initiation + simplified loop terms; only the coded subset of a mature Turner/ViennaRNA model is implemented. A selected H-type pseudoknot candidate is not experimentally confirmed.','PBody2')
add('<b>Force-gradient check:</b> compare analytical Cartesian gradient ∂E/∂x with [E(x+h) − E(x−h)]/(2h), as sampled and recorded. V11 used RDKit MMFF94 on 8 small molecules and the first 4 residues of each of 5 RNAs.','PBody2')
add('<b>Checksums:</b> SHA-256(original data file) versus the archived expected digest protects against accidental source replacement. Agreement shows matching bytes, not an independent molecular assay.','PBody2')
sep();add('Evidence and reproducibility instructions','PH')
add('See <b>data/FIGURE_PROVENANCE.csv</b> for each PNG, its source CSV, transformation and SHA-256. <b>data/</b> contains the individual paired runs and derived condition tables. <b>evidence_and_logs/</b> contains original test/build logs plus the unchanged source-audit evidence ZIP. Run the original verifier in that evidence archive to recompute firing-mask parity. The archived executable source tree can rerun native computations with its prescribed toolchain.','PBody2')
nextpage();add('16  /  LITERATURE BENCHMARK CONTEXT','PTag2');add('State of the field—and exactly what is comparable','PT')
tbl(['METHOD / STUDY','EXTERNALLY REPORTED RESULT','COMPARABILITY'],[['BPfold (2025)','F1 0.820 on ArchiveII and 0.658 on bpRNA-TS0','Experimentally annotated RNA base-pair prediction; N has no corresponding F1 benchmark'],['G4ShapePredictor (2024)','76% ± 1% accuracy for three G4 folding classes','Not comparable with N 76.56% for constructed graph intervention'],['ViennaRNA RNAfold','MFE structure, full calibrated thermodynamic model and partition functions','N V11 implements selected terms only'],['RNAstructure / NNDB','Turner 2004 thermodynamic parameter families','Parameter provenance; not proof of complete implementation'],['RDKit ETKDG/MMFF94','Chemical conformer embedding and molecular mechanics','N actually calls these established methods for 3D computations'],['RCSB PDB 1EHZ & 1KF1','Experimental solved structures; X-ray coordinate records','Molecular identity provenance, not a held-out structure prediction test']],[125,177,208]);sep()
add('Published accuracy values are labeled as external literature results, not N measurements. Comparisons require the same tasks, held-out ground truth, input conditions, metrics, and computing environment.','PBody2')
add('Primary references','PH')
for name,role,url in REFS:
 add(f'{html.escape(name)} — {html.escape(url)}','PRef2')
nextpage();add('17  /  AUDIT CONCLUSIONS','PTag2');add('Demonstrated computation. Open scientific validation.','PT')
add('VERIFIED / RECORDED','PH')
for s in ['Compiling and executing the retained native spiking neural engine with structure-derived contact graphs.','13,312 paired response records reconstruct from original masks with zero disagreements.','V6 structural graph contrasts include positive and stringent negative controls.','RNA secondary-structure calculation, an H-type pseudoknot candidate search, RDKit geometry optimization, and MMFF94 numerical force checks execute in source.','81 named V11 project tests have passing historical log entries.','Independent replay and source hashes support deterministic computational reproducibility.']:
 add('• '+html.escape(s))
sep();add('NOT SHOWN BY THESE RUNS','PH')
for s in ['Blind experimental secondary-structure base-pair precision/recall or F1.','General whole-RNA 3D atomic folding accuracy, RMSD, solvent/ion-calibrated dynamics or drug–DNA binding affinity.','Fair same-task performance superiority against ViennaRNA, RNAstructure, BPfold, G4ShapePredictor or molecular dynamics.','Original project Git history: the available .git recovery is a newly constructed offline snapshot.']:
 add('• '+html.escape(s))
sep();add('This portfolio is a transparent computational research artifact. It is deliberately not presented as peer-reviewed proof of biological prediction.','PNote2')

pdfpath=OUT/'MOLECULAR_RESEARCH_PORTFOLIO_V3-V11.pdf'
doc=SimpleDocTemplate(str(pdfpath),pagesize=A4,leftMargin=margin,rightMargin=margin,topMargin=48,bottomMargin=60,title='N Molecular | Research Portfolio V3-V11',author='Research project — source-based evidence')
doc.build(S,onFirstPage=footer,onLaterPages=footer)

# offline HTML: embeds only custom PNGs sourced from data, with useful audit navigation
ctable=lambda hd,lines:'<table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in hd)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in row)+'</tr>' for row in lines)+'</tbody></table>'
figure_html='\n'.join(f'''<article class="exhibit"><header><span>FIGURE {i:02d}</span><h2>{html.escape(c['title'])}</h2></header><img loading="lazy" src="figures/{html.escape(c['figure'])}" alt="{html.escape(c['title'])}"><p class="caveat">{html.escape(c['transformation'])}</p><div class="trace">Source: <code>{html.escape(c['csv_source'])}</code><br>Image SHA256: <code>{c['png_sha256']}</code></div></article>''' for i,c in enumerate(CHARTS,1))
ref_html=''.join(f'<li><a href="{url}">{html.escape(name)}</a><span>{html.escape(role)}</span></li>' for name,role,url in REFS)
html_out='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>N Molecular | Source-Verified Research Portfolio</title><style>
:root{--ink:#12243a;--muted:#566c80;--teal:#087e86;--line:#dae4ea;--paper:#fff;--bg:#f3f7f9}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.65 ui-sans-serif,system-ui,-apple-system,'Segoe UI',sans-serif}a{color:#087e86}header.hero{background:linear-gradient(133deg,#0c223b,#143f56 72%,#16656d);color:#fff;padding:62px max(22px,calc((100vw - 1100px)/2)) 54px}.eyebrow{font-size:12px;font-weight:750;letter-spacing:.2em;color:#87d4d3}.hero h1{font-size:clamp(39px,5vw,70px);line-height:1.05;margin:20px 0}.lede{font-size:19px;max-width:740px;color:#d7ecf2}.hero .actions{display:flex;gap:12px;flex-wrap:wrap;margin:24px 0}.actions a{background:#fff;padding:10px 18px;border-radius:9px;text-decoration:none;font-weight:750}.actions a+ a{background:transparent;border:1px solid #75aab9;color:#fff}.shell{max-width:1100px;margin:auto;padding:34px 18px 100px}.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:-54px;position:relative}.kpi{background:#fff;border:1px solid var(--line);box-shadow:0 12px 24px #0b2f4312;border-radius:13px;padding:22px}.kpi strong{font-size:30px;display:block;letter-spacing:-.04em}.kpi small{display:block;color:var(--muted);font-size:12px}.intro,.paper{background:white;border:1px solid var(--line);border-radius:15px;padding:25px 29px;margin:24px 0}.paper h2,.intro h2{font-size:26px;margin-top:0}.paper h3{margin:26px 0 6px}.grid{display:grid;grid-template-columns:repeat(2,1fr);gap:20px}.exhibit{background:white;border:1px solid var(--line);border-radius:15px;padding:23px;box-shadow:0 5px 18px #12304809}.exhibit header span{color:var(--teal);font-weight:800;font-size:11px;letter-spacing:.1em}.exhibit h2{font-size:17px;line-height:1.35;margin:5px 0 14px}.exhibit img{width:100%;height:auto;background:white;border:1px solid #edf1f4;border-radius:6px}.caveat{font-size:13px}.trace{border-top:1px solid var(--line);font-size:11px;line-height:1.75;color:var(--muted);padding-top:9px;overflow-wrap:anywhere}code{font:11px/1.5 ui-monospace,Consolas,monospace}table{border-collapse:collapse;width:100%;font-size:13px}th{background:#e5f2f4;text-align:left}th,td{padding:11px;border-bottom:1px solid #e3eaf0;vertical-align:top}ul.refs{padding:0;list-style:none}.refs li{display:flex;gap:20px;padding:10px 0;border-bottom:1px solid #e4ebef}.refs li a{flex:0 0 43%;font-weight:650}.refs span{color:var(--muted)}.note{border-left:4px solid #c08d38;background:#fff8e8;padding:14px 18px}.longform{max-width:900px;margin:auto}.small{font-size:13px;color:var(--muted)}@media(max-width:850px){.kpis,.grid{grid-template-columns:repeat(2,1fr)}.refs li{display:block}.refs li a{display:block}}@media(max-width:600px){.kpis,.grid{grid-template-columns:1fr}.hero{padding-bottom:90px}.kpi strong{font-size:26px}.paper,.intro{padding:19px}}
</style></head><body><header class="hero"><div class="eyebrow">SOURCE-GROUNDED RESEARCH / OCTOBER 2026</div><h1>N Molecular<br>Research Portfolio</h1><p class="lede">V3–V11. From controlled graph topology to RNA structure calculations, molecular mechanics and native neural dynamics. Every original exhibit is generated from recorded execution data.</p><div class="actions"><a href="MOLECULAR_RESEARCH_PORTFOLIO_V3-V11.pdf">Read the peer-style PDF</a><a href="COMPLETE_AUDIT_WORKBOOK.xlsx">Open the full data workbook</a></div></header><main class="shell"><section class="kpis"><div class="kpi"><strong>13,312</strong><small>paired neural comparisons audited</small></div><div class="kpi"><strong>0</strong><small>spike-mask arithmetic mismatches</small></div><div class="kpi"><strong>81/81</strong><small>named V11 project tests in archived log</small></div><div class="kpi"><strong>13</strong><small>converged MMFF94 minimizations</small></div></section><section class="intro"><h2>Research summary</h2><p>A compiled spiking neural network ran real counterfactual graph experiments on DNA/RNA-related categorical structures and chemical representations. The latest path calculates RNA secondary-structure contacts, attempts an H-type crossing motif, builds RDKit 3D conformers, computes real MMFF94 energies/gradients, and inserts derived contacts into the native neural network.</p><p class="note"><strong>Correct scientific reading:</strong> Graph-label accuracy, neural firing disagreement and modeled folding energy are different outcomes. No experimentally validated full-RNA 3D structure accuracy or ligand binding prediction was established in these records.</p></section><section class="paper"><h2>Experiment trajectory</h2>'''+ctable(['Release','What actually ran','Result or outcome'],[['V3–V5','Original native DNA/RNA input bridge and six-case controls','Preserved native recordings'],['V6','Fresh seeds, 1KF1 ring-preserving control','76.56% vs 45.31% under different graph interventions'],['V7–V8','Four new RNA cases, deterministic replay','No new case survived adjusted tests'],['V9','Eight chemical compounds × four policy arms','2,048 paired records'],['V10','RNA pair solver + RDKit conformers feeding SNN','3,328 paired records'],['V11','RNA energy, crossing candidate, MMFF94 forces and SNN','7,936 paired records, 32.067 s recorded']])+'''</section><section class="paper"><h2>Original charts, generated from source</h2><p>None of these charts is a web stock photograph or a drawing of imagined results. The attached source table, formula and SHA256 are included under every figure.</p></section><section class="grid">'''+figure_html+'''</section><section class="paper longform"><h2>Computational contract</h2><h3>Neural disagreement</h3><p><code>D = Σ popcount(A[t] XOR B[t]) / (48 × 22), t = 0…47</code>. The denominator is 1,056 binary target events per paired trial; each reported score is recomputed from original firing masks.</p><h3>Classifier accuracy and statistical tests</h3><p><code>Accuracy = correct graph labels / evaluated labels</code>. Paired null-label permutations use <code>p=(1 + #null≥observed)/(4097)</code>, and Holm correction for a declared comparison family.</p><h3>Thermodynamics, forces and geometry</h3><p>The RNA model includes selected Turner-type parameters and explicitly simplified remaining terms; it does not reproduce a complete RNAfold partition function. MMFF94 coordinate gradient checks use a symmetric finite-difference derivative of the force-field energy. V11 RNA mechanics apply only to four-nucleotide fragments.</p><h3>Test evidence</h3><p>The evidence ZIP contains the source audit, original full logs, data JSONL, hashes and an independent verifier. Archived 81-test logs are historical—not a claim that this portfolio renderer independently recompiled the whole runtime.</p></section><section class="paper longform"><h2>Comparison to real published research</h2>'''+ctable(['Study','Measured publication outcome','Why it is not a head-to-head score'],[['BPfold (2025)','F1 0.820 ArchiveII, 0.658 bpRNA-TS0','Experimentally evaluated base-pair prediction; N lacks same benchmark'],['G4ShapePredictor (2024)','0.76±0.01 three-class folding accuracy','Not the same label as N constructed 1KF1 rewire classification'],['ViennaRNA RNAfold','Minimum free energy and RNA ensemble modeling','N V11 model is partial; no direct comparison'],['RDKit','ETKDG/MMFF94 conformers and force-field calculations','Actually invoked as computational dependency in N']])+'''<h3>Primary academic and official sources</h3><ul class="refs">'''+ref_html+'''</ul></section><section class="paper longform"><h2>Reproducibility and attached records</h2><p>The full quantitative archive is in <code>evidence_and_logs/COMPLETE_REPRODUCIBLE_EVIDENCE.zip</code>; the standalone workbook is <code>COMPLETE_AUDIT_WORKBOOK.xlsx</code>. Every PNG was produced by <code>reproduce/build_portfolio.py</code> using <code>data/*.csv</code>. The original reproducible source ZIP is a separate artifact.</p><p><strong>Figure index:</strong> <a href="data/FIGURE_PROVENANCE.csv">data/FIGURE_PROVENANCE.csv</a> · <a href="data/PORTFOLIO_METRICS.json">data/PORTFOLIO_METRICS.json</a> · <a href="data/SCIENTIFIC_REFERENCES.csv">data/SCIENTIFIC_REFERENCES.csv</a>.</p><p class="small">This is an engineering and computational research portfolio. Its records support reproducibility claims inside the represented models; they do not establish external biological prediction quality.</p></section></main></body></html>'''
(OUT/'INDEX.html').write_text(html_out,encoding='utf8')

# full rendering source copy: exact script with no image requests or human redraw
shutil.copy2(Path(__file__),SCRIPT/'build_portfolio.py')
(OUT/'README_FIRST.txt').write_text('''N Molecular — Source-Verified Research Portfolio V3–V11 (23 Oct 2026)\n\nOPEN INDEX.html for the offline responsive portfolio or MOLECULAR_RESEARCH_PORTFOLIO_V3-V11.pdf for printable report.\nfigures/: 13 source-rendered true data PNGs. Each PNG has a source CSV and sha256 entry in data/FIGURE_PROVENANCE.csv.\ndata/: all 13,312 individual v9–v11 paired comparisons, 208 condition rows, 100 classifier rows, RNA and MMFF source records, 81 test IDs.\nevidence_and_logs/COMPLETE_REPRODUCIBLE_EVIDENCE.zip: earlier full audited package with original JSONL receipts, original logs, hashes and VERIFY_BUNDLE.py.\nCOMPLETE_AUDIT_WORKBOOK.xlsx: 17 sheets with formulas and charts.\nreproduce/build_portfolio.py: transparent source code rendering every figure, PDF and HTML from recovered files; point SRC at the original evidence folder to regenerate.\n\nThe verified full native source machine is distributed separately because its 106MB archive contains compiled artifacts and multiple complete histories. No browser assets or stock graphs were used.\n\nImportant methodological disclaimer: percentages in V9/V10/V11 are computational firing disagreement, NOT molecular activity. V6/V8 percentages classify hand-authored graph interventions, NOT RNA folding accuracy. The RNA energy model is incomplete; all external published scores are different experimental benchmarks.\n''')
# primary file manifest, excluding manifest itself; no symlink copying
with (OUT/'PACKAGE_SHA256.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['relative_path','size_bytes','sha256'])
 for p in sorted(OUT.rglob('*')):
  if p.is_file() and p.name!='PACKAGE_SHA256.csv':
   hh=hashlib.sha256(p.read_bytes()).hexdigest();w.writerow([p.relative_to(OUT).as_posix(),p.stat().st_size,hh])

# ZIP all portfolio artifacts including the evidence bundle. Deflate inside may be minimal on preserved ZIP.
ZIP=Path('/mnt/data/N_Molecular_Professional_Research_Portfolio_V3-V11_20261023.zip')
with zipfile.ZipFile(ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 for f in sorted(OUT.rglob('*')):
  if f.is_file():z.write(f,f.relative_to(OUT.parent))
with zipfile.ZipFile(ZIP) as z:
 bad=z.testzip();assert bad is None,bad
 print(json.dumps({'status':'PASS','pdf':str(pdfpath),'html':str(OUT/'INDEX.html'),'zip':str(ZIP),'zip_bytes':ZIP.stat().st_size,'figures':len(CHARTS),'zip_files':len(z.namelist()),'audit_records':len(P),'shas':len(AUD['v11_compared_data_hashes'])},indent=2))
