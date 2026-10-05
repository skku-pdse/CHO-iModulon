"""Redraw stored discovery metadata associations; no statistical rerun."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from figure1_paths import DATA, RESULTS
ROOT = DATA.parents[1]
SOURCE = ROOT / 'discovery iModulon/conditional association/discovery_metadata_association_results.csv'
OUT = RESULTS / 'overview'
OUT.mkdir(parents=True, exist_ok=True)
FIELDS = ['Cell line', 'Temperature shift', 'Time point (day)', 'Phase', 'Culture', 'Media', 'Engineered']
COLORS = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#687887', '#8064A2']
LABELS = ['Cell line', 'Temperature shift', 'Timepoint', 'Phase', 'Culture', 'Media', 'Engineered']
d = pd.read_csv(SOURCE)
assert len(d) == 945 and not d.duplicated(['iModulon', 'metadata']).any()
assert np.isfinite(d[['omega2', 'q_value']]).all().all()
assert d.groupby('metadata').size().eq(105).all()
s = d[d.metadata.isin(FIELDS) & d.q_value.lt(.05)].copy()
plt.rcParams.update({'font.family': 'Arial', 'font.size': 10, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(7.1, 4.6))
fig.subplots_adjust(left=.24, right=.88, bottom=.16, top=.91)
coords = {}
for i, (field, color) in enumerate(zip(FIELDS, COLORS)):
    sub = s[s.metadata.eq(field)].sort_values(['omega2', 'iModulon'])
    placed = []
    for row in sub.itertuples():
        x = row.omega2
        # Deterministic compact swarm: x values are never changed.
        candidates = [0.] + [v for j in range(1, 12) for v in (j*.072, -j*.072)]
        yoff = next(y for y in candidates if all(((x-px)/.013)**2 + ((y-py)/.072)**2 >= 1 for px, py in placed))
        placed.append((x, yoff))
        coords[(field, row.iModulon)] = (x, i+yoff)
    ax.scatter([p[0] for p in placed], [i+p[1] for p in placed], s=17, color=color, edgecolors='white', linewidths=.32, alpha=.92, zorder=3)
    ax.text(1.04, i, str(len(sub)), va='center', ha='center', fontsize=10, clip_on=False)
ax.set_yticks(range(len(FIELDS)),LABELS)
ax.set_ylim(len(FIELDS)-.15,-.8)
ax.set_xlim(-.025,1.015)
ax.set_xticks(np.arange(0,1.01,.2))
ax.set_xlabel(r'Metadata association ($\omega^2$)',labelpad=8)
ax.grid(axis='x',color='#E5E8EB',linewidth=.65,zorder=0)
ax.tick_params(axis='y',length=0,pad=9)
ax.tick_params(axis='x',length=3,color='#666666')
for spine in ['top','right','left']: ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_color('#777777')
ax.text(1.04,-.65,'n',ha='center',va='center',fontstyle='italic')
fig.text(.025,.94,'F',fontsize=16,fontweight='bold')
fig.text(.24,.965,'Significant associations (BH-adjusted q < 0.05)',fontsize=9,color='#444444',va='top')
for ext in ['png','pdf','svg','tiff']:
    options = {'pil_kwargs': {'compression':'tiff_lzw'}} if ext == 'tiff' else {}
    fig.savefig(OUT/f'figure1F_metadata_association.{ext}',dpi=600 if ext=='tiff' else 300,facecolor='white',**options)
plt.close(fig)
s.to_csv(OUT/'plotted_associations.csv',index=False)
summary=pd.DataFrame({'metadata':FIELDS,'significant_iModulons':[sum(s.metadata.eq(f)) for f in FIELDS],'tested_iModulons':105,'color':COLORS})
summary.to_csv(OUT/'factor_summary.csv',index=False)
(OUT/'provenance.json').write_text(json.dumps({'source':str(SOURCE),'sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'source_tests':945,'displayed_factors':FIELDS,'omitted_factors':['Producer','Project'],'q_threshold':.05,'new_tests':False,'displayed_points':len(s)},indent=2))
(OUT/'caption_and_notes.md').write_text('## Suggested caption\n\n(F) Distribution of iModulon activity associations with seven metadata factors. Each point represents an iModulon with a significant association (BH-adjusted q < 0.05); horizontal position indicates one-way ANOVA omega-squared. Vertical offsets separate overlapping points and carry no quantitative meaning. Colors distinguish metadata factors, and n indicates the number of significant iModulons out of 105 tested per factor. An iModulon may appear in multiple rows.\n\n## Provenance and interpretation\n\nValues are taken unchanged from the existing version-2 discovery results (945 tests across nine factors). Global BH-adjusted q values are retained, including tests for Producer and Project, which are not displayed. The seven plotted factors follow the agreed panel scope. Timepoint uses quartile-binned culture time. Each factor was tested separately; associations are not adjusted for Project or other factors. Counts depend on sampling and power and do not measure factor importance. Original results and figures are preserved.\n',encoding='utf-8')
print(summary.to_string(index=False))
print(OUT)

