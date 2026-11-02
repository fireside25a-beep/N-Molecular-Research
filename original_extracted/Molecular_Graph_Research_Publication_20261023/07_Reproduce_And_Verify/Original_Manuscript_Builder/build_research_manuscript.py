#!/usr/bin/env python3
"""Refine the existing sealed molecular portfolio IN PLACE.
All scientific charts and numbers are reconstructed from its original CSV records.
The archival PDF/LaTeX/HTML are backed up before changes. No molecular simulation is run.
"""
from pathlib import Path
import csv, math, hashlib, json, shutil, re, subprocess, sys, textwrap
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('/mnt/data/N_Molecular_Research_Portfolio_V3-V11_FINAL_SEALED_20261023')
D=ROOT/'data'
FIG=ROOT/'figures_publication'
ARCH=ROOT/'supplement'/'archival_manuscript_before_publication_refinement'
FIG.mkdir(exist_ok=True,parents=True)
ARCH.mkdir(exist_ok=True,parents=True)
OLD=[ROOT/'latex'/'portfolio_extended.tex',ROOT/'MOLECULAR_RESEARCH_PORTFOLIO_V3-V11_EXTENDED_LATEX.pdf',ROOT/'INDEX.html',ROOT/'README_EXTENDED_FIRST.txt']
for p in OLD:
 if p.exists() and not (ARCH/p.name).exists():shutil.copy2(p, ARCH/p.name)

def readcsv(name):return list(csv.DictReader((D/name).open(newline='',encoding='utf-8')))
G=readcsv('GRAPH_CLASSIFICATION_CONDITIONS.csv')
C=readcsv('ALL_CONDITIONS_FROM_RAW.csv')
P=readcsv('ALL_V9_V10_V11_PAIRED_COMPARISONS.csv')
R=readcsv('RNA_FOLDING_STRUCTURES.csv')
A=readcsv('ATOMIC_FORCE_RESULTS.csv')
T=readcsv('RECORDED_TIMING_SECONDS.csv')
TEST=readcsv('TEST_CASE_INDEX.csv')
assert len(G)==100 and len(C)==208 and len(P)==13312 and len(A)==13 and sum(z['version']=='v11' for z in TEST)==81
assert not any(int(x['mask_mismatch']) for x in P)
assert all(int(x['changed_spikes'])==int(x['mask_xor_count']) for x in P)
assert all(abs(float(x['fraction_recomputed'])-int(x['changed_spikes'])/(int(x['tick_count'])*int(x['target_slots'])))<1e-12 for x in P)
assert Counter(x['version'] for x in P)=={'v11':7936,'v10':3328,'v9':2048}

# Publication figure design, neutral and consistent.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':12,'axes.titleweight':'bold','axes.labelsize':9,'figure.dpi':160,'savefig.dpi':220,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'axes.axisbelow':True,'grid.color':'#DCE5E9','grid.alpha':0.7,'savefig.facecolor':'white','figure.facecolor':'white'})
TEAL='#176D78'; INDIGO='#506AA3'; GRAY='#607782'; ORANGE='#B67B44'
figmeta=[]
def finish(name, title, paths, caption):
 fp=FIG/f'{name}.png'
 plt.tight_layout(pad=1.3)
 plt.savefig(fp,dpi=220,bbox_inches='tight')
 plt.close()
 sha=hashlib.sha256(fp.read_bytes()).hexdigest()
 figmeta.append({'figure':str(fp.relative_to(ROOT)),'title':title,'derived_from':'; '.join(paths),'caption':caption,'sha256':sha})
 return fp

def row_match(source,**kwargs):
 cand=[r for r in source if all(r[k]==v for k,v in kwargs.items())]
 assert len(cand)==1, (kwargs,len(cand))
 return cand[0]

# F1: six-condition graph classification
ids=['1BNA','1ANA','1WOE','1RNA','1ZIH','1KF1']
rows=[row_match(G,version='v6',pdb=k,readout='late_scaffold_raw',intervention='original',feedback='closed',plasticity='STDP') for k in ids]
old=[100*float(r['original_accuracy']) for r in rows]
new=[100*float(r['fresh_accuracy']) for r in rows]
fig,ax=plt.subplots(figsize=(9.4,4.25));x=np.arange(6)
ax.bar(x-.19,old,.36,label='Earlier held-out',color=GRAY)
ax.bar(x+.19,new,.36,label='Fresh-seed held-out',color=TEAL)
ax.axhline(50,color=ORANGE,ls='--',lw=1.2,label='Balanced chance reference')
ax.set_xticks(x,ids);ax.set_ylim(0,100);ax.set_ylabel('Graph-intervention classification accuracy (%)')
ax.set_title('Held-out discrimination of constructed contact graphs')
ax.legend(frameon=False,ncol=3,loc='upper left',fontsize=8)
finish('Figure_01_Held_Out_Graph_Classification','Held-out graph classification',['data/GRAPH_CLASSIFICATION_CONDITIONS.csv'],'Six primary computational cases; earlier and fresh-seed held-out accuracies are distinct test groups.')

# F2: exact intervention control
ra=row_match(G,version='v6',pdb='1KF1',readout='late_scaffold_raw',intervention='original',feedback='closed',plasticity='STDP')
rb=row_match(G,version='v6',pdb='1KF1',readout='late_scaffold_raw',intervention='ring_preserving_plus3',feedback='closed',plasticity='STDP')
fig,ax=plt.subplots(figsize=(8,4.2))
vals=[100*int(ra['correct'])/int(ra['total']),100*int(rb['correct'])/int(rb['total'])]
ax.bar(['Cycle organization changed','Ring organization preserved'],vals,width=.55,color=[TEAL,GRAY]);ax.axhline(50,ls='--',color=ORANGE,label='Chance reference')
for i,v in enumerate(vals):ax.text(i,v+2,f'{int([ra,rb][i]["correct"])}/256  ({v:.2f}%)',ha='center',weight='bold')
ax.set_ylim(0,100);ax.set_ylabel('Classifier accuracy (%)');ax.set_title('Structural selectivity in the telomeric quadruplex-labelled graph')
finish('Figure_02_Cycle_Intervention_Control','Ring-structure contrast',['data/GRAPH_CLASSIFICATION_CONDITIONS.csv'],'Original contrast 196/256, ring-preserving contrast 116/256; separate interventions, not molecular folding predictions.')

# F3: source-defined graph intervention schematic (1 -> +2, +3 offsets)
fig,axes=plt.subplots(1,3,figsize=(10,3.55));
for ax,off,tit in zip(axes,[1,2,3],['Original four-cycle','Cycle-breaking offset +2','Cycle-preserving offset +3']):
 th=np.linspace(0,2*np.pi,5)[:-1]+np.pi/4
 pos=np.stack([np.cos(th),np.sin(th)],axis=1)
 for i in range(4):
  j=(i+off)%4
  a=pos[i]*0.83;b=pos[j]*0.83
  ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='-|>',color=TEAL,lw=1.8,shrinkA=12,shrinkB=12))
 for i,(xi,yi) in enumerate(pos):ax.scatter([xi],[yi],s=390,color=INDIGO,zorder=3);ax.text(xi,yi,str(i),color='white',ha='center',va='center',weight='bold',zorder=4)
 ax.set_aspect('equal');ax.set_xlim(-1.5,1.5);ax.set_ylim(-1.4,1.5);ax.axis('off');ax.set_title(tit,fontsize=10)
fig.suptitle('Source-specified directed contact-cycle interventions (schematic)',fontsize=12,weight='bold')
finish('Figure_03_Source_Specified_Cycle_Graphs','Directed-cycle topology schematic',['research/v6/SCIENTIFIC_REPORT.md (archived in evidence bundle)'],'Diagram generated from documented directed edge offsets (+1,+2,+3); nodes and edges illustrate the computational intervention.')

# F4: added 4 RNA cases (negative)
more=['1ZIF→1ZIF','1ZIG→1ZIG','2KOC→2KOC','1EHZ→1EHZ']
added=[next(r for r in G if r['version']=='v7-v8' and r['pdb']==id and r['readout']=='late_scaffold_shape' and r['feedback']=='closed' and r['plasticity']=='STDP') for id in more]
fig,ax=plt.subplots(figsize=(8,3.65)); vals4=[100*float(r['fresh_accuracy']) for r in added]
ax.bar(['1ZIF','1ZIG','2KOC','1EHZ'],vals4,color=GRAY)
ax.axhline(50,ls='--',color=ORANGE,label='Balanced chance reference')
for i,v in enumerate(vals4):ax.text(i,v+1.2,f'{v:.2f}%',ha='center',fontsize=9)
ax.set_ylim(0,75);ax.set_ylabel('Held-out graph classification (%)');ax.set_title('Additional RNA contact-graph tests: negative controls retained')
finish('Figure_04_Additional_RNA_Cases','Four additional RNA cases',['data/GRAPH_CLASSIFICATION_CONDITIONS.csv'],'No added-case primary comparison passed its documented multiple-testing correction.')

# F5: V9 chemical condition heatmap no V text
compounds=['aspirin','benzene','caffeine','cyclohexane','ethanol','naphthalene','nicotinamide','pyridine']
cond=[('native_closed','STDP','State feedback\nlearning'),('native_closed','frozen','State feedback\nfrozen'),('matched_fixed','STDP','Fixed feedback\nlearning'),('matched_fixed','frozen','Fixed feedback\nfrozen')]
mat=np.array([[100*float(row_match(C,version='v9',branch='chemical_bonds',structure_id=c,feedback=fb,plasticity=pl)['fraction_from_hits']) for fb,pl,_ in cond] for c in compounds])
fig,ax=plt.subplots(figsize=(8.5,5))
im=ax.imshow(mat,aspect='auto',cmap='YlGnBu',vmin=0,vmax=22)
ax.grid(False);ax.set_yticks(range(8),[c.title() for c in compounds]);ax.set_xticks(range(4),[x[2] for x in cond]);ax.set_title('Chemical-bond ablations: full feedback × plasticity design')
for i in range(8):
 for j in range(4): ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if mat[i,j]>11 else 'black')
fig.colorbar(im,ax=ax,shrink=.82,label='Changed target firing positions (%)')
finish('Figure_05_Chemical_Bond_Ablation','Chemical bond ablation factorial design',['data/ALL_CONDITIONS_FROM_RAW.csv'],'Eight compounds × four conditions; response is difference in target SNN spike positions, not binding affinity.')

# F6 paired per-seed distributions directly raw, aggregate median quartiles as box plots. caffeine 64 seeds x four conditions
cpairs=[r for r in P if r['version']=='v9' and r['structure_id']=='caffeine']
data=[]
for fb,pl,label in cond:
 values=[100*float(r['fraction_recomputed']) for r in cpairs if r['feedback']==fb and r['plasticity']==pl]
 assert len(values)==64,len(values)
 data.append(values)
fig,ax=plt.subplots(figsize=(9.1,4.1));bp=ax.boxplot(data,tick_labels=['State + learn','State + frozen','Fixed + learn','Fixed + frozen'],patch_artist=True,showfliers=False,medianprops={'color':'black','linewidth':1.5})
for patch in bp['boxes']:patch.set_facecolor('#BCDCE1')
ax.set_ylabel('Changed target positions per seed (%)');ax.set_title('Caffeine: distributions across 64 paired deterministic seeds')
finish('Figure_06_Caffeine_Paired_Seed_Distributions','Caffeine per-seed response distributions',['data/ALL_V9_V10_V11_PAIRED_COMPARISONS.csv'],'Each box contains 64 paired seed-specific contrasts; this is computational variability, not independent biological replication.')

# F7 RNA solver pair counts
order=['1ZIH','1ZIF','1ZIG','2KOC','1EHZ'];r0=[row_match(R,version='v10',pdb=i) for i in order];r1=[row_match(R,version='v11',pdb=i) for i in order]
fig,ax=plt.subplots(figsize=(8.5,4.15));x=np.arange(5)
ax.bar(x-.18,[int(v['nested_pairs']) for v in r0],.35,label='Weighted nested-pair solver',color=GRAY)
ax.bar(x+.18,[int(v['nested_pairs']) for v in r1],.35,label='Selected energy-term solver',color=TEAL)
ax.set_xticks(x,order);ax.set_ylabel('Algorithmically predicted base pairs');ax.set_title('RNA pairing outputs from two source-implemented models');ax.legend(frameon=False)
finish('Figure_07_RNA_Pairing_Algorithms','RNA pairing algorithms',['data/RNA_FOLDING_STRUCTURES.csv'],'Counts describe source-generated base pairs, not experimental base-pair precision or recall.')

# F8 partial energy output
fig,ax=plt.subplots(figsize=(8.3,3.85));energy=[float(r['model_energy_kcal_mol']) for r in r1]
ax.bar(order,energy,color=TEAL);ax.set_ylabel('Recorded model energy (kcal/mol)');ax.set_title('Energy outputs from the implemented RNA scoring model');ax.axhline(0,color=GRAY,lw=.8)
for i,v in enumerate(energy):ax.text(i,v-.35,f'{v:.2f}',va='top',ha='center',fontsize=8)
finish('Figure_08_RNA_Model_Energies','Recorded RNA model energies',['data/RNA_FOLDING_STRUCTURES.csv'],'These are source-generated energies of the implemented partial model, not measured folding thermodynamics.')

# F9 selected RNA ablation heatmap
matRNA=np.array([[100*float(row_match(C,version='v10',branch='rna_nussinov',structure_id=rn,feedback=fb,plasticity=pl)['fraction_from_hits']) for fb,pl,_ in cond] for rn in order])
fig,ax=plt.subplots(figsize=(8.5,4.2));im=ax.imshow(matRNA,aspect='auto',cmap='YlGnBu',vmin=0,vmax=11)
ax.grid(False);ax.set_yticks(range(5),order);ax.set_xticks(range(4),[k[2] for k in cond]);ax.set_title('Predicted RNA pairing contacts: native neural sensitivity')
for i in range(5):
 for j in range(4):ax.text(j,i,f'{matRNA[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if matRNA[i,j]>5 else 'black')
fig.colorbar(im,ax=ax,shrink=.82,label='Changed target positions (%)')
finish('Figure_09_RNA_Contact_Ablations','RNA contact ablations',['data/ALL_CONDITIONS_FROM_RAW.csv'],'Actual SNN responses after removal of algorithmically predicted pairing contacts.')

# F10 MMFF energy gradients on LOG scale
fig,ax=plt.subplots(figsize=(10.2,4.8));labels=[r['structure_id'].replace('_RNA4',' (4 nt)') for r in A]
err=[float(r['max_abs_gradient_delta']) for r in A]
assert len(A)==13 and all(x['converged']=='1' for x in A)
ax.bar(range(13),err,color=[INDIGO]*8+[TEAL]*5)
ax.set_yscale('log');ax.set_ylabel('Max sampled gradient discrepancy (kcal/mol/Å)');ax.set_xticks(range(13),labels,rotation=45,ha='right',fontsize=8)
ax.set_title('MMFF94: analytical vs finite-difference gradient checks')
finish('Figure_10_Force_Gradient_Checks','MMFF94 gradient checks',['data/ATOMIC_FORCE_RESULTS.csv'],'Four sampled Cartesian components per structure; numerical consistency does not validate complete molecular folding.')

# F11 spatial contacts post force from source row counts
contacts=[]
for item in A:
 cr=row_match(C,version='v11',branch='post_force',structure_id=item['structure_id'],feedback='native_closed',plasticity='STDP')
 contacts.append(int(cr['contacts']))
fig,ax=plt.subplots(figsize=(10.2,4.5));ax.bar(range(13),contacts,color=[INDIGO]*8+[TEAL]*5)
ax.set_xticks(range(13),labels,rotation=45,ha='right',fontsize=8);ax.set_ylabel('Computed nonbonded contacts');ax.set_title('Contacts derived from MMFF94-minimized structures')
finish('Figure_11_Post_Optimization_Contacts','Post-optimization contacts',['data/ALL_CONDITIONS_FROM_RAW.csv','data/ATOMIC_FORCE_RESULTS.csv'],'Five short RNA fragments are four nucleotides long; counts depend on the code-defined contact criterion.')

# F12 late chemical response means, within post-force unique 8 compounds
fig,ax=plt.subplots(figsize=(8.8,4.3));compounds_sort=compounds
means=[100*np.mean([float(r['fraction_from_hits']) for r in C if r['version']=='v11' and r['branch']=='post_force' and r['structure_id']==c]) for c in compounds_sort]
ax.barh([c.title() for c in compounds_sort[::-1]],means[::-1],color=TEAL)
ax.set_xlabel('Mean changed target firing positions across four conditions (%)');ax.set_title('Response to removal of computed post-optimization contacts')
finish('Figure_12_Post_Force_Chemical_Sensitivity','Post-optimization compound responses',['data/ALL_CONDITIONS_FROM_RAW.csv'],'Unweighted arithmetic means of the four recorded conditions; not chemical potency or binding prediction.')

# F13 stage times without incorrect sum comparison
stages=[r for r in T if r['stage']!='full_pipeline_wall']
fig,ax=plt.subplots(figsize=(9,5.35));stages=sorted(stages,key=lambda r:float(r['seconds']))
ax.barh([r['stage'].replace('_seconds','').replace('_',' ') for r in stages],[float(r['seconds']) for r in stages],color=INDIGO)
ax.set_xlabel('Measured wall time (seconds)');ax.set_title('Recorded pipeline stage times on the execution host')
finish('Figure_13_Record_Stage_Timings','Recorded stage timings',['data/RECORDED_TIMING_SECONDS.csv'],'Stage-specific timings are one recorded execution, not a controlled benchmark against other software.')

# F14 software verification evidence, historical
fig,ax=plt.subplots(figsize=(7,3.5));ax.bar(['Bond-graph screen','Contact-derived pipeline','Extended molecular analysis'],[31,59,81],color=[GRAY,INDIGO,TEAL]);ax.set_ylabel('Named historical Python tests');ax.set_title('Archived software test counts (not biological accuracy)')
finish('Figure_14_Software_Regression','Archived regression tests',['data/TEST_CASE_INDEX.csv','data/PORTFOLIO_METRICS.json'],'Counts from preserved test logs; independent biology validation is not included.')

# write figure provenance
with (D/'PUBLICATION_FIGURE_PROVENANCE.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['figure','title','derived_from','caption','sha256']);w.writeheader();w.writerows(figmeta)

# small manuscript generator

def fmt(x,n=2):return f'{float(x):.{n}f}'
def escape(x):
 return str(x).replace('\\','\\textbackslash{}').replace('_','\\_').replace('%','\\%').replace('&','\\&').replace('#','\\#')
def table(headers,body,colspec=None,size='\\small'):
 colspec=colspec or ('l'+'r'*(len(headers)-1))
 start=f'\\begin{{center}}{size}\\begin{{longtable}}{{{colspec}}}\n\\toprule\n'+' & '.join(headers)+r'\\'+'\n\\midrule\n\\endfirsthead\n\\toprule\n'+' & '.join(headers)+r'\\'+'\n\\midrule\n\\endhead\n'
 lines=[' & '.join(escape(x) for x in r)+r'\\' for r in body]
 return start+'\n'.join(lines)+'\n\\bottomrule\n\\end{longtable}\n\\end{center}\n'
def figure(name,caption):
 return f'\\begin{{figure}}[p]\\centering\\includegraphics[width=.96\\textwidth,height=.66\\textheight,keepaspectratio]{{../figures_publication/{name}.png}}\\caption{{{caption}}}\\end{{figure}}\n'

six_body=[[i,fmt(old[j]),fmt(new[j]),f"{rows[j]['correct']}/{rows[j]['total']}",fmt(rows[j]['holm_p'],4)] for j,i in enumerate(ids)]
chem_body=[[c.title()]+[fmt(mat[i,j]) for j in range(4)] for i,c in enumerate(compounds)]
rna_body=[[i,str(r0[j]['nested_pairs']),str(r1[j]['nested_pairs']),fmt(r1[j]['model_energy_kcal_mol'],3),'Yes' if r1[j]['pseudoknot_candidate_selected']=='1' else 'No'] for j,i in enumerate(order)]
atom_body=[[a['structure_id'].replace('_RNA4',' RNA fragment'),a['atoms_h_included'],fmt(a['energy_after'],3),f"{float(a['max_abs_gradient_delta']):.2e}",'Yes' if a['converged']=='1' else 'No'] for a in A]
run_counts=[['Bond-graph perturbation', '2,048','4,096'],['Predicted RNA contacts + 3D ligands','3,328','6,656'],['Energy-derived RNA + post-force contacts','7,936','15,872']]

s=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=20mm,headheight=15pt]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,microtype}
\usepackage{graphicx,booktabs,longtable,array,amsmath,amssymb,xcolor,hyperref,float,fancyhdr,enumitem,caption,needspace}
\definecolor{ink}{HTML}{173340}\definecolor{accent}{HTML}{176D78}\definecolor{muted}{HTML}{617684}
\hypersetup{colorlinks=true,linkcolor=accent,urlcolor=accent,pdftitle={Structural Sensitivity in a Native Spiking Neural Network Under Molecular Graph Interventions}}
\captionsetup{font=small,labelfont=bf}
\pagestyle{fancy}\fancyhf{}\lhead{\small Molecular graph interventions and spiking dynamics}\rhead{\small Computational research}\cfoot{\thepage}
\setlength{\parskip}{4pt}\setlength{\parindent}{0pt}
\begin{document}
\begin{titlepage}
\vspace*{34mm}
{\Huge\bfseries\color{ink} Structural Sensitivity in a Native Spiking Neural Network Under Molecular Graph Interventions\par}
\vspace{11mm}
{\Large\color{muted} Controlled graph perturbations, RNA pairing inputs, atomic-force-derived contacts, and deterministic neural responses\par}
\vspace{21mm}
{\large Archival computational research manuscript\par}
\vspace{6mm}
{\normalsize Evidence snapshot: 23 October 2026\par}
\vspace{15mm}
\begin{tabular}{rl}
\textbf{13,312} & independently recalculated paired neural comparison records\\
\textbf{100} & archived classifier-statistic rows\\
\textbf{208} & four-factor computational condition summaries\\
\textbf{81} & named software tests in preserved execution records\\
\end{tabular}
\vfill
\textbf{Scope:} A reproducible computational study of graph-intervention sensitivity. The reported neural responses are not biological folding accuracy, binding affinities, or experimental physical measurements. Authorship and affiliation information are reserved for the research authors; none is invented here.
\end{titlepage}
\tableofcontents\clearpage
\section*{Abstract}\addcontentsline{toc}{section}{Abstract}
\textbf{Background.} The capacity of a temporally evolving spiking network to retain decodable information about graph organization can be investigated through controlled contact interventions.\par
\textbf{Methods.} We analyze original compiled spiking-network records spanning a six-case paired graph-classification study, four additional RNA-labelled negative tests, eight chemical-bond graph ablations, five RNA input sequences and 13 molecular-mechanics structures. Matched graph contrasts used deterministic seed resets, controlled stimulation, feedback and plasticity manipulations. Complete 48-tick target-neuron masks were independently rechecked for 13,312 paired comparisons.\par
\textbf{Results.} In the 1KF1-labelled computational graph, the classifier reached 196/256 correct (76.56\%) under a cycle-breaking rewire, but 116/256 (45.31\%) under a different ring-preserving rewire. The chemical and contact-derived runs yielded nonzero response differences in their supported conditions; they also contained null and reversed learning effects. Thirteen MMFF94 minimizations reported convergence with sampled analytical/numerical gradient agreement. The 13,312 mask-derived disagreement records matched the archived counts.\par
\textbf{Interpretation.} The compiled network exhibits reproducible sensitivity to selected computational contact interventions. These data do not establish general structural equivalence recognition, experimental RNA structure prediction quality or ligand affinity.\par
\textbf{Keywords:} spiking neural network; directed graph cycles; molecular representations; contact ablation; deterministic replay; RNA secondary structure; molecular mechanics.

\section{Introduction}
Changes in contact organization can affect a dynamical neural system even when other input properties are held fixed. To investigate this question without equating computational trajectories with physical chemistry, the experiments considered graph structures associated with deposited DNA and RNA sequences, SMILES-encoded molecular connectivity, and subsequently computed RNA pairing and nonbonded spatial contacts.

The central objective is to quantify \emph{which computational interventions produced decodable or changed neural activity}, not to infer the unknown three-dimensional fold of a molecule. Two endpoints are analyzed separately: held-out classification of intervention identity in the original graph experiments, and the fraction of changed target-neuron firing positions in later paired mask experiments. These endpoints use different experimental units, denominators and hypotheses and are not interchangeable.

Existing molecular folding tools such as ViennaRNA RNAfold predict minimum-free-energy secondary structures and related ensemble observables. Published spiking graph-learning methods and classical graph algorithms also offer distinct comparison targets. No matched benchmark of this particular experiment against those methods appears in the preserved records. Direct cross-paper superiority claims are therefore inappropriate \cite{ViennaRNA,MethodsPolicies}.

\section{Hypotheses and scope}
\textbf{H1 (graph perturbation):} Changes to selected constructed contact graphs can be discriminated from native spiking activity under held-out deterministic seeds.\par
\textbf{H2 (cycle selectivity):} A cycle-breaking graph contrast and a cycle-preserving contrast need not produce equal classification performance.\par
\textbf{H3 (contact influence):} Removing chemical-bond, RNA-pairing, or force-derived contacts can change the native SNN spike masks for matched computational inputs.\par
\textbf{Not tested:} Prediction of measured molecular binding, held-out physical folding coordinates, chemical reaction rates, or medical outcomes.

\section{Materials and methods}
\subsection{Native network and experimental unit}
The preserved original computational scaffold contains 118 neurons and 324 directed synapses. Experiment-specific input nodes and contact edges are added without replacing the compiled native state-transition routine. Runs consist of 48 deterministic ticks per independently reset network. A paired comparison holds seed and declared exogenous stimulation policy matched while varying the recorded contact graph. A seed reset is the unit of computational variation; it is not a biological replicate.

The experimental factors cross native-output-conditioned feedback versus a matched-current fixed schedule with synaptic learning enabled (spike-timing-dependent plasticity, STDP) versus frozen connection weights. The feedback-selected input target can depend on prior native spike output and evolved weights; it is not a chemistry-aware autonomous scientist.

\subsection{Historical classification experiment}
The initial six input cases use deposited structure identifiers 1BNA, 1ANA, 1WOE, 1RNA, 1ZIH and 1KF1, but the alternative computational contact graphs were \emph{manually constructed} rather than derived from complete crystallographic coordinates. Labels denote the intact or rewired graph. For the later hold-out, the frozen linear classifier evaluated 256 graph labels over 128 fresh stimulus seed pairs. Original classifier readouts use late scaffold-neuron spike identities. Accuracy is
\begin{equation}\mathrm{Accuracy}=\frac{N_{\mathrm{correct}}}{N_{\mathrm{evaluated}}}.\end{equation}
The preserved paired-label sign-swap null uses 4,096 permutations, with plus-one probability
\begin{equation}p_{\mathrm{cond}}=\frac{1+\#\{a_{\mathrm{null}}\geq a_{\mathrm{observed}}\}}{4097}.\end{equation}
Multiple-test adjustments in the archive follow Holm's step-down rule across the identified hypotheses. These conditional probabilities depend on the chosen training procedure and label-exchangeability assumptions; they are not population biological $p$ values.

\subsection{Constructed cycle controls}
For the telomeric quadruplex-labelled graph, the original constructed contacts contain directed four-site rings. The coarse rewire replaces their offset-1 arrows with offset-2 arrows, splitting a four-cycle into two reciprocal two-cycles while preserving directed vertex degrees. The subsequently introduced offset-3 alternative preserves each four-cycle while changing its orientation. This was a \emph{post hoc} control, not prospective independent physical confirmation. Figure~\ref{fig:cycle} draws only the source-described graph construction.

\subsection{Paired spike-mask difference}
For each pair, the original native output records 48 temporal target masks over 22 target slots. If $A_t$ and $B_t$ denote two masks at step $t$, the recorded computational response is
\begin{equation}\label{eq:disagreement}D=\frac{\sum_{t=0}^{47}\operatorname{popcount}(A_t\oplus B_t)}{48\times22}.\end{equation}
The audit independently recomputes each XOR/popcount from the recorded masks and reconstructs condition means by summing changed positions over matched pairs. A 64-pair condition therefore has denominator $64\times48\times22$ positions. This is a difference statistic, not a measure of folding success or molecular potency.

\subsection{Chemical-bond and RNA-contact perturbations}
Chemical structures in the source are eight RDKit-validated SMILES compounds. The bond-graph experiment retains the same atoms while removing the declared bond-derived computational edges. Later input generation uses RDKit ETKDG conformers and MMFF94 optimization to derive spatial contacts, while RNA input generation uses weighted noncrossing base-pair optimization and, in a subsequent branch, selected Turner-derived stack/hairpin terms. The later branch also searches for crossing H-type stems using a separate \emph{heuristic} scoring criterion. These model outputs should not be called a complete RNA folding free-energy landscape.

\subsection{Atomistic force evaluation}
RDKit MMFF94 was used for eight small-molecule structures and five four-nucleotide RNA fragments. Sampled Cartesian analytical derivatives were compared with centered finite differences, with four checked components per molecular entry in the final record. A small sampled gradient discrepancy confirms consistency of those calculations; it does not establish experimentally correct forces or full-RNA fold geometry \cite{RDKit}.

\subsection{Data analysis and independent checks}
The archived result table contains 13,312 paired records: 2,048 bond-graph comparisons, 3,328 predicted-contact comparisons and 7,936 extended model/force-contact comparisons. The underlying 48-tick mask pairs were independently audited and grouped into 208 condition summaries. A separate table retains 100 historical graph-classifier result rows. The original logs list 81 passing Python test cases for the final pipeline. The latter are \emph{historical executions}, not tests newly executed by this manuscript-compilation step.

\section{Results}
\subsection{Held-out discrimination and cycle controls}
The constructed 1KF1 cycle-breaking comparison achieved 196/256 correct predictions (76.56\%) on fresh held-out seeds, with a six-case Holm-adjusted probability of 0.02783. A separate ring-preserving comparison produced 116/256 correct (45.31\%) and did not demonstrate reliable separation under that readout. This pattern is consistent with a sensitivity to the changed coarse cycle organization. Chance-level performance in the second comparison does \emph{not} prove generalized graph-equivalence recognition. The remaining six-case values and negative results appear in Table~\ref{tab:six}.
'''

s+= '\\begin{table}[ht]\\centering\\caption{Constructed graph interventions: historical and fresh-seed classification. Percentages are conditional computational accuracies, not structure-prediction scores.}\\label{tab:six}\\small\\begin{tabular}{lrrrr}\\toprule\n PDB & Earlier (\\%) & Fresh (\\%) & Fresh correct & Holm-adjusted $p$ \\\\\midrule\n'
for r in six_body: s+=' & '.join(str(x) for x in r)+r' \\'+'\n'
s+='\\bottomrule\\end{tabular}\\end{table}\n'
s+=figure('Figure_01_Held_Out_Graph_Classification','Earlier and fresh-seed held-out accuracy for six computational graph interventions. Source: archived graph-classification table.')
s+=figure('Figure_02_Cycle_Intervention_Control','Distinct graph interventions in the telomeric quadruplex-labelled case. Differences in classification are not DNA folding accuracy.')
s+='\\begin{figure}[p]\\centering\\includegraphics[width=.98\\textwidth]{../figures_publication/Figure_03_Source_Specified_Cycle_Graphs.png}\\caption{Source-defined directed cycle interventions, schematic from recorded edge offsets (+1, +2, +3).}\\label{fig:cycle}\\end{figure}\n'
s+=r'''\subsection{Four additional RNA-labelled contact graph cases}
The four later RNA-labelled cases did not pass their documented multiple-test correction. Their measured primary accuracies were 44.92\% for 1ZIF, 44.92\% for 1ZIG, 53.52\% for 2KOC and 48.05\% for 1EHZ. These are classification of constructed graph interventions, not four external RNA folding tests.
'''
s+=figure('Figure_04_Additional_RNA_Cases','Additional RNA-labelled intervention accuracies. The preserved null results are part of the outcome.')
s+=r'''\subsection{Chemical-bond graph ablations}
Eight chemically valid SMILES inputs were screened under the two feedback modes and two plasticity states. The result was the proportion of target-neuron firing positions that changed when bond-derived contact edges were removed. For caffeine, the native-feedback/STDP condition yielded 17.17\% altered firing positions; the corresponding frozen-weight condition yielded 10.93\%. Such differences demonstrate sensitivity to the computational intervention, not stimulant activity, toxicity or affinity. The full factorial values are shown in Figure~\ref{fig:chem} and the per-seed caffeine distribution is retained rather than reporting only means.
'''
s+='\\begin{figure}[p]\\centering\\includegraphics[width=.95\\textwidth]{../figures_publication/Figure_05_Chemical_Bond_Ablation.png}\\caption{All eight compounds across the four experimentally declared SNN control settings.}\\label{fig:chem}\\end{figure}\n'
s+=figure('Figure_06_Caffeine_Paired_Seed_Distributions','Caffeine paired-seed distributions in four SNN configurations; each box summarizes 64 computational seed pairs.')
s+=r'''\subsection{RNA structure calculations and contact ablations}
For the five RNA input sequences, the earlier weighted nested-pair solver predicted respectively 4, 4, 4, 5 and 26 pairs; the selected-energy-term model predicted 4, 4, 4, 5 and 23 pairs. The later model outputs were $-2.74$, $-2.74$, $-2.74$, $-5.43$ and $-14.60356$ kcal/mol. They reflect its implemented model terms, including simplified internal/bulge penalties, and are not experimentally measured RNA free energies. A separate crossing-stem search selected a 1EHZ candidate; the candidate's heuristic score is not Gibbs free energy and does not validate the actual tRNA tertiary fold.
'''
s+='\\begin{table}[ht]\\centering\\caption{Algorithmic RNA pairing results and recorded model energies.}\n\\small\\begin{tabular}{lrrrr}\\toprule RNA & Weighted pairs & Energy-model pairs & Energy (kcal/mol) & Crossing candidate \\\\\midrule\n'
for rr in rna_body: s+=' & '.join(str(x) for x in rr)+r' \\'+'\n'
s+='\\bottomrule\\end{tabular}\\end{table}\n'
s+=figure('Figure_07_RNA_Pairing_Algorithms','Predicted base-pair counts for the same five RNA sequences under two documented computational models.')
s+=figure('Figure_08_RNA_Model_Energies','Recorded RNA model energies. They do not quantify experimental folding accuracy.')
s+=figure('Figure_09_RNA_Contact_Ablations','Native target-neuron sensitivity to removal of RNA-pairing contacts. These are ablation-response percentages, not RNA base-pair F1.')
s+=r'''\subsection{MMFF94 molecular mechanics and spatial contacts}
Eight small molecules and five four-nucleotide RNA fragments were subjected to RDKit MMFF94 minimization. The final record reports convergence for all 13 cases and a largest sampled analytical-vs-numerical gradient difference of $1.59\times 10^{-7}$ kcal/mol/\AA. The five RNA inputs are short \emph{fragments}, not whole RNA atomistic folding predictions. The fragments attached to 1ZIH, 1ZIF and 1ZIG share their first four nucleotides, and therefore share computed atomic outputs; they are not independent molecular confirmations. The coordinate-derived contacts and the recorded neural-response differences depend on the implemented contact definition.
'''
s+=figure('Figure_10_Force_Gradient_Checks','Numerical gradient consistency for all 13 source-minimized molecular systems.')
s+=figure('Figure_11_Post_Optimization_Contacts','Computed spatial contacts for eight small-molecule structures and five four-nucleotide RNA fragments.')
s+=figure('Figure_12_Post_Force_Chemical_Sensitivity','Chemical contact-ablation responses in the post-optimization branch; means across the four recorded settings, not potency estimates.')
s+=r'''\subsection{Recorded execution and regression}
The final molecular analysis campaign produced 7,936 paired comparisons, representing 15,872 network resets at 48 ticks each, or 761,856 SNN ticks. One recorded end-to-end pipeline run took 32.067 seconds on its execution host, excluding a fresh full upstream rebuild. This is a single-host measured workload, not a task-matched speed comparison with independent molecular simulators. Cross-compiler/sanitizer evidence and 81 named historical Python test outcomes are preserved in the supplementary logs. The later one-shot upstream test command was interrupted, although individual follow-up commands passed; it is not described as a full one-shot successful rerun.
'''
s+=figure('Figure_13_Record_Stage_Timings','Recorded stage timings for the final 64-seed workload; compiler and independent biological benchmarks are excluded.')
s+=figure('Figure_14_Software_Regression','Historical Python test counts in preserved logs; software regressions are distinct from scientific prediction accuracy.')
s+=r'''\section{Discussion}
The principal computational observation is that a stateful native spiking network can exhibit reproducible contact-graph sensitivity under specified interventions. The telomeric quadruplex-labelled cycle experiment offers a particularly interpretable result because its non-significant ring-preserving control narrows the mechanism consistent with the original contrast. The later chemical, pairing-derived and post-force contact experiments demonstrate additional input pathways into the same native network, but use different interventions and should not be treated as direct repetitions of the cycle-discrimination experiment.

A key alternative explanation is that simple graph features or exact cycle-decomposition algorithms suffice for portions of the intervention task. The existing record does not contain a controlled same-task graph-algorithm baseline or held-out biological structure benchmark. Accordingly, this manuscript makes neither a novelty claim for graph cycle detection itself nor a claim of superiority to established spiking graph models or physical folding predictors. Evidence here supports \emph{conditional computational sensitivity and reproducibility}.

The binary classifier score and mask-disagreement fraction measure different events and should not be pooled. Similarly, statistical variability over deterministic seeds is not the same as biological uncertainty over independent specimens. The figure provenance, negative results and exact build/test records are retained to make those distinctions reviewable.

\section{Limitations and independent validation}
The contact interventions in the original DNA/RNA studies were manually constructed, even though their molecular accession records are real. The RNA energy model is not a complete nearest-neighbor thermodynamic implementation, and the crossing-stem predictor is a deterministic heuristic candidate search. MMFF94 calculations involve short RNA fragments rather than full coordinate prediction of 76-nucleotide tRNA. No ligand-binding measurements or held-out experimental RNA coordinate comparisons under the same evaluation protocol appear in this record.

Necessary future experiments include graph-label tests across new graph families with exact-algorithm and graph-learning baselines, blind experimental RNA base-pair comparisons with family-separated splits, and independent three-dimensional structure evaluation when a full-coordinate model exists. These are \emph{planned validations}, not experiments claimed as completed here.

\section{Reproducibility and availability}\subsection{Data availability}
All numerical observations discussed here are contained in the accompanying local research archive. The main files are \texttt{data/GRAPH\_CLASSIFICATION\_CONDITIONS.csv}, \texttt{data/ALL\_CONDITIONS\_FROM\_RAW.csv}, \texttt{data/ALL\_V9\_V10\_V11\_PAIRED\_COMPARISONS.csv}, \texttt{data/RNA\_FOLDING\_STRUCTURES.csv}, \texttt{data/ATOMIC\_FORCE\_RESULTS.csv}, and \texttt{data/RECORDED\_TIMING\_SECONDS.csv}. Each source-data lineage is preserved in the supplementary evidence package. \emph{No external DOI or public deposit is claimed.}
\subsection{Code availability}
The original executable research source and its historical archives are bundled in \texttt{native\_source\_release/}. The paper-generation code is under \texttt{publication\_methods/}. Re-execution requires a matching native development toolchain and the documented scientific dependencies; the provided audit scripts can independently verify the numerical records without rerunning the molecular engine. No GitHub repository update or new permanent upstream deposit is implied.
\subsection{Verification levels}
Recorded output integrity and numerical recomputation are distinct from independent biological validity. The portfolio verifier recalculates the 13,312 individual disagreement records and 208 grouped conditions and checks 100 graph classification rows. The previous software logs are historical and are labelled accordingly. A current rendered-PDF build is evidence of document generation, not a fresh native molecular experiment.

\section{Conclusion}
The preserved native SNN records establish reproducible intervention-dependent responses for specified constructed molecular contact graphs and programmatically derived chemical, RNA and short-fragment structural contacts. In the main graph-classification case, the strong cycle-breaking result and the non-significant ring-preserving control together narrow the supported explanation to sensitivity to the tested graph organization. Biological folding and ligand prediction performance, and relative superiority to established methods, require prospective task-matched validation.

\section{Selected scientific references}
\begin{thebibliography}{9}
\bibitem{ViennaRNA} ViennaRNA Package. RNAfold manual: minimum-free-energy structures and partition functions. \url{https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/man/RNAfold.html}.
\bibitem{RDKit} RDKit Project. Force-field helpers API and documented molecular force-field routines. \url{https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html}.
\bibitem{RCSB} RCSB Protein Data Bank. PDB 1EHZ: yeast phenylalanine tRNA experimental crystal structure (1.93 \AA). \url{https://www.rcsb.org/structure/1EHZ}.
\bibitem{MethodsPolicies} PLOS Computational Biology. Submission guidelines: reproducibility, data and code availability, claims supported by methods. \url{https://journals.plos.org/ploscompbiol/s/submission-guidelines}.
\bibitem{SRPolicies} Scientific Reports. Editorial policies: data and code availability for reproducibility. \url{https://www.nature.com/srep/journal-policies/editorial-policies}.
\end{thebibliography}
\clearpage\appendix
\section{Supplementary experiment registry and provenance}
Internal historical version identifiers (e.g. \texttt{v6}, \texttt{v9}, \texttt{v10}, \texttt{v11}) are retained \emph{only} in the machine-readable source metadata to locate their original records. They are not scientific model names or public results labels. The preserved older portfolio reports are included under \texttt{prior\_portfolios/} for archival transparency; this manuscript supersedes their public-facing narrative.
\subsection{Dataset and outcome ledger}
'''
s+=table(['Experiment','Paired comparisons','Native resets'],run_counts,size='\\small')
s+='\\subsection{Five RNA sequences}\n'+table(['PDB','Weighted pairs','Energy-model pairs','Energy (kcal/mol)','Crossing candidate'],rna_body,'lrrrr','\\small')
s+='\\subsection{Atomistic minimization ledger}\n'+table(['Structure','Atoms','Final MMFF energy','Derivative error','Converged'],atom_body,'lrrrr','\\scriptsize')
s+=r'''\subsection{Independent audit contract}
Figure-to-data SHA-256 mapping is in \texttt{data/PUBLICATION\_FIGURE\_PROVENANCE.csv}. The numerical verification script checks source pairing, integer XOR-popcount totals, denominator normalization, grouped condition counts and historical graph-classification rows. A reader can reproduce the entire PDF via the included LaTeX source and \texttt{latexmk} or \texttt{pdflatex}. The original archive includes test logs, compiler receipts, invalid-input tests and cross-compiler replay evidence; these are not replaced with synthetic outputs.
\end{document}
'''

# Fix LaTeX source construction: escape errors introduced by literals table strings.
s=s.replace(r"\section{Discussion}", "\\section{Relationship to published computational research}\n\\textbf{RNA secondary structure.} ViennaRNA RNAfold uses a parameterized thermodynamic model to compute minimum-free-energy structures and RNA ensemble properties \\cite{ViennaRNA}. KnotFold uses learned base-pair potentials and a minimum-cost-flow optimization to produce structures that may contain pseudoknots \\cite{KnotFold}. The present partial-energy and heuristic crossing-stem branches do not implement these methods and were not evaluated on their held-out structure benchmarks.\n\n\\textbf{DNA G-quadruplexes.} G4ShapePredictor uses experimentally or literature-supported sequence--topology data to predict classes of G-quadruplex folding \\cite{G4ShapePredictor}. The 1KF1-labelled classifier in this manuscript instead distinguishes two intentionally constructed \\emph{computational} contact graphs. Its 76.56\\% result is not directly comparable with that paper's measured three-class topology task.\n\n\\textbf{Molecular machine learning.} MoleculeNet provides standardized molecular-property datasets, task metrics and model comparisons \\cite{MoleculeNet}. No comparable measured property labels are provided for the eight-compound neural-response experiment. Accordingly the response rates are not treated as surrogate predictive accuracy, toxicity, affinity, or chemical efficacy.\n\n\\textbf{Research software.} FAIR-for-research-software principles emphasize metadata, persistent identification, traceability and reusability \\cite{FAIR4RS}. The attached source and results provide local auditability, while public persistent archival identifiers, independent domain review and task-matched scientific validation would still be needed before claims of full reuse compliance.\n\n\\section{Discussion}")
s=s.replace(r"\end{thebibliography}", '\\bibitem{KnotFold} Gong, T., Ju, F. and Bu, D. Accurate prediction of RNA secondary structure including pseudoknots through solving minimum-cost flow with learned potentials. \\emph{Communications Biology} 7, 297 (2024). \\url{https://www.nature.com/articles/s42003-024-05952-w}.\n\\bibitem{G4ShapePredictor} Liew, D., Lim, Z. W. and Yong, E. H. Machine learning-based prediction of DNA G-quadruplex folding topology with G4ShapePredictor. \\emph{Scientific Reports} 14, 24238 (2024). \\url{https://doi.org/10.1038/s41598-024-74826-2}.\n\\bibitem{MoleculeNet} Wu, Z. et al. MoleculeNet: a benchmark for molecular machine learning. \\emph{Chemical Science} 9, 513--530 (2018). \\url{https://doi.org/10.1039/C7SC02664A}.\n\\bibitem{FAIR4RS} Barker, M. et al. Introducing the FAIR Principles for research software. \\emph{Scientific Data} 9, 622 (2022). \\url{https://www.nature.com/articles/s41597-022-01710-x}.\n\n\\end{thebibliography}')
LATEX=ROOT/'latex'/'portfolio_extended.tex';LATEX.write_text(s,encoding='utf-8')
# public landing page
landing='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Structural Sensitivity in a Native Spiking Neural Network</title><style>body{margin:0;background:#f7fafb;color:#173340;font:16px/1.6 system-ui,-apple-system,sans-serif}main{max-width:1080px;margin:auto;padding:32px 20px}h1{font-size:clamp(2rem,4vw,3.4rem);line-height:1.12;margin:24px 0}h2{margin-top:3rem;border-top:1px solid #dce5e9;padding-top:2rem}p{max-width:80ch}.deck{font-size:1.3rem;color:#506a76}.tag{font-weight:700;color:#176d78;letter-spacing:.14em;text-transform:uppercase;font-size:.74rem}a{color:#176d78}figure{margin:30px 0;padding:16px;background:white;border:1px solid #e2ebef;border-radius:12px}figure img{width:100%;height:auto;display:block}figcaption{font-size:.87rem;color:#526a75}.stats{display:flex;flex-wrap:wrap;gap:12px;margin:22px 0}.stat{min-width:155px;background:#eef5f6;border-radius:9px;padding:15px}.stat strong{display:block;font-size:1.5rem}.foot{font-size:.85rem;color:#647782}</style></head><body><main>
<span class="tag">Computational research manuscript | archival evidence</span><h1>Structural Sensitivity in a Native Spiking Neural Network Under Molecular Graph Interventions</h1>
<p class="deck">Controlled graph perturbations, RNA pairing inputs, molecular mechanics and reproducible native neural responses.</p>
<p>The study measures changes in neural activity under precisely specified molecular-style contact interventions. It does not claim experimental DNA or RNA folding accuracy, ligand binding, or outperformance of physical simulation software.</p>
<div class="stats"><div class="stat"><strong>13,312</strong>paired records audited</div><div class="stat"><strong>208</strong>condition summaries</div><div class="stat"><strong>100</strong>classification records</div><div class="stat"><strong>81</strong>historical tests logged</div></div>
<p><a href="RESEARCH_MANUSCRIPT.pdf">Read the scientific manuscript (PDF)</a> · <a href="RESEARCH_MANUSCRIPT.tex">LaTeX source</a> · <a href="COMPLETE_AUDIT_WORKBOOK.xlsx">Evidence workbook</a> · <a href="data/PUBLICATION_FIGURE_PROVENANCE.csv">Figure provenance and SHA-256</a></p>
<h2>Source-derived evidence figures</h2>
'''
for d in figmeta:
 landing+=f'<figure><img loading="lazy" alt="{d["title"]}" src="{d["figure"]}"><figcaption><b>{d["title"]}.</b> {d["caption"]} Source: {d["derived_from"]}.</figcaption></figure>\n'
landing+='''<h2>Methods, evidence and review</h2><p>The manuscript contains the mathematical definitions of each outcome, all primary graph controls and negative results, source and code availability, experimental limitations and selected scientific references.</p><p>The original source archives, executable research code, run receipts and earlier report editions remain in the supplementary folders without fabricated experiments or modified historical records.</p><p class="foot">All publication plots are rendered from the packaged numerical source records; internal development identifiers remain only in technical provenance.</p></main></body></html>'''
(ROOT/'INDEX.html').write_text(landing,encoding='utf-8')
(ROOT/'README_FIRST.txt').write_text('OPEN RESEARCH_MANUSCRIPT.pdf for the public-facing paper; INDEX.html for its figures. The research evidence workbook, provenance and historical logs are bundled. Internal release identifiers are maintained only in source metadata. The manuscript does not assert experimental molecular prediction accuracy.\n')
(ROOT/'RESEARCH_MANUSCRIPT.tex').write_text(s,encoding='utf-8')

# save script under publication methods and compile
OUTSCRIPT=ROOT/'publication_methods'/'build_research_manuscript.py'
OUTSCRIPT.parent.mkdir(exist_ok=True)
if Path(__file__).resolve() != OUTSCRIPT.resolve():
 shutil.copy2(__file__,OUTSCRIPT)
print('Generated',len(figmeta),'public figures and',len(s),'characters of LaTeX source')
