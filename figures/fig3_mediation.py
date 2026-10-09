# ============================================================================
# Figure 3 - Mediation NPS -> Brain EAA -> CDR (full sample, adjusted)
# Inputs: analysis dataset (Epi_cr3x.csv) for the a/b/c/c' path coefficients,
#         and the ACME Excel (full-sample mediation) for the standardized ACME
#         and its FDR. Paths a, b, c, c' are standardized betas from OLS.
#
# Set DATA_PATH and ACME_XLSX below.
# Output: Fig3_mediation.{pdf,png,svg}
# ============================================================================
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42; mpl.rcParams['ps.fonttype'] = 42; mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']})
import matplotlib.pyplot as plt, matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'
ACME_XLSX = '/content/drive/MyDrive/Doctorado/Mediation_NPS_ACME_FULL_country.xlsx'

df = pd.read_csv(DATA_PATH)
# --- variables (as in the analysis) ---
df['NPS']  = (df['npi_total'] >= 3).astype(int)     # clinically significant NPS
df['sexm'] = (df['demo_sex'] == 2).astype(int)      # 2 = male
EXPO, MED, OUT, DX = 'NPS', 'Brain_EAA', 'cdr_cdrtot', 'clinical_diagnosis'
COVS = ['demo_age', 'sexm', 'demo_smoke', 'cog_ed', 'meds_nps']
CTRY = 'C(demo_resid, Treatment(reference=4))'      # Colombia reference

req = [EXPO, MED, OUT, 'demo_age', 'sexm', 'demo_smoke', 'cog_ed', 'meds_nps', DX, 'demo_resid']
d = df.dropna(subset=req).copy()
rhs = ' + '.join(COVS + [f'C({DX})', CTRY])
ma = smf.ols(f'{MED} ~ {EXPO} + {rhs}', data=d).fit()               # a: NPS -> Brain
mb = smf.ols(f'Q("{OUT}") ~ {MED} + {EXPO} + {rhs}', data=d).fit()  # b and c'
mc = smf.ols(f'Q("{OUT}") ~ {EXPO} + {rhs}', data=d).fit()          # c (total)

def zb(m, term, xcol, ycol):
    return m.params[term] * d[xcol].std() / d[ycol].std(), m.pvalues[term]
ba, pa   = zb(ma, EXPO, EXPO, MED)
bb, pb   = zb(mb, MED,  MED,  OUT)
cpr, ppr = zb(mb, EXPO, EXPO, OUT)
ct, pt   = zb(mc, EXPO, EXPO, OUT)

raw = pd.read_excel(ACME_XLSX); raw.columns = [str(c).strip() for c in raw.columns]
sel = raw[(raw['Block'] == 'Full sample') &
          (raw['Mediator'].astype(str).str.contains('Brain')) &
          (raw['Outcome'].astype(str).str.contains('CDR'))].iloc[0]
def num(x): return float(str(x).replace('*', '').replace(',', '.').strip())
acme = num(sel['ACME (std)']); fdr = num(sel['ACME P-FDR'])

def st(p): return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else ''
print(f"βa={ba:.3f}{st(pa)}  βb={bb:.3f}{st(pb)}  c={ct:.3f}{st(pt)}  "
      f"c'={cpr:.3f}{st(ppr)}  ACME={acme:.3f}{st(fdr)}(FDR)  N={int(mb.nobs)}")

# ---------------- Figure (manuscript colors) ----------------
FILL_NPS, FILL_BRAIN, FILL_CDR = '#F0E442', '#009E73', '#A60000'
def tint(c, w=0.90):
    r, g, b = mcolors.to_rgb(c); return (r+(1-r)*w, g+(1-g)*w, b+(1-b)*w)
fig, ax = plt.subplots(figsize=(92/25.4, 56/25.4)); ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis('off')
def box(x, y, w, h, t, c):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle='round,pad=0.02,rounding_size=0.1',
                 linewidth=2.0, edgecolor=c, facecolor=tint(c), zorder=2))
    ax.text(x, y, t, ha='center', va='center', fontsize=7, color='black', zorder=3)
box(1.7, 1.25, 2.3, 0.95, 'NPS\n(significant)', FILL_NPS)
box(5.0, 4.7, 2.3, 0.95, 'Brain EAA', FILL_BRAIN)
box(8.3, 1.25, 2.3, 0.95, 'CDR\n(severity)', FILL_CDR)
def arr(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>,head_width=2.5,head_length=4',
                 mutation_scale=1, linewidth=0.9, color='black', zorder=1, shrinkA=5, shrinkB=5))
arr(2.2, 1.78, 4.2, 4.25); arr(5.8, 4.25, 7.8, 1.78); arr(2.95, 1.25, 7.05, 1.25)
ax.text(2.35, 3.15, f'$\\beta$a = {ba:.3f}{st(pa)}', fontsize=7, ha='center')
ax.text(7.65, 3.15, f'$\\beta$b = {bb:.3f}{st(pb)}', fontsize=7, ha='center')
ax.text(5.0, 1.5, f'c = {ct:.3f}{st(pt)}', fontsize=7, ha='center', bbox=dict(facecolor='white', edgecolor='none', pad=0.4))
ax.text(5.0, 0.72, f"c' = {cpr:.3f}{st(ppr)}", fontsize=7, ha='center')
ax.text(5.0, 0.25, f"Indirect effect (ACME-std) = {acme:.3f}{st(fdr)}", fontsize=7, ha='center')
plt.tight_layout()
for ext in ['pdf', 'png', 'svg']:
    fig.savefig(f'Fig3_mediation.{ext}', dpi=450, bbox_inches='tight', facecolor='white')
print('Figure 3 generated (PDF, PNG, SVG).')
