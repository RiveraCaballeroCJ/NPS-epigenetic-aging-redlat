# ============================================================================
# Figure - NPI-Q symptoms by group (two panels)
#   (a) Within the clinically significant NPS group: mean NPI-Q severity
#       CONTRIBUTION, Dementia (AD+FTD) vs Controls (CN)
#   (b) Clinically significant NPS vs Normal NPS - prevalence (%)
#   Horizontal grouped bars, symptoms ordered within the four Gonzalez domains.
#   Panel (a) metric = mean of each symptom's 0-3 severity over the subgroup
#   (absent = 0); summing across symptoms approximates the group's mean NPI-Q
#   total, so bar length shows how much each symptom contributes.
#   Width <= 18 cm. Colors: Dementia #a60000, Controls #1d5f22,
#   Significant #E69F00, Normal #0072B2.
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: Fig_symptom_prevalence.{pdf,png,svg}
# ============================================================================
import pandas as pd, numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42; mpl.rcParams['ps.fonttype'] = 42; mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 6})
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'

DOMAINS = [
 ('Psychosis',           [('npi_del', 'Delusions'), ('npi_hall', 'Hallucinations')]),
 ('Activation',          [('npi_agit', 'Agitation/aggression'), ('npi_disn', 'Disinhibition'), ('npi_irr', 'Irritability/lability')]),
 ('Affective',           [('npi_depd', 'Depression/dysphoria'), ('npi_anx', 'Anxiety')]),
 ('Somatic/Behavioural', [('npi_apa', 'Apathy/indifference'), ('npi_mot', 'Aberrant motor'), ('npi_nite', 'Night-time behaviour'), ('npi_app', 'Appetite/eating')]),
 ('Other',               [('npi_elat', 'Euphoria/elation')]),
]
C_DEM = '#a60000'; C_CN = '#1d5f22'; C_SIG = '#E69F00'; C_NORM = '#0072B2'

df = pd.read_csv(DATA_PATH)
df = df[df['clinical_diagnosis'].isin(['AD', 'FTD', 'CN'])].copy()
df = df[df['demo_resid'].isin([1, 3, 4, 5, 6])].copy()
df['DEM'] = df['clinical_diagnosis'].isin(['AD', 'FTD'])
df = df[df['npi_total'].notna()].copy()
df['SIG'] = df['npi_total'] >= 3

def prevalence(mask, item):
    sub = pd.to_numeric(df.loc[mask, item], errors='coerce').dropna()
    return 100 * (sub > 0).mean() if len(sub) > 0 else np.nan

def sev_contrib(mask, item):
    """Contribution to the NPI-Q total: severity (npi_*sev, 1-3) with absent = 0,
       averaged over the subgroup. Summed across the 12 symptoms = group's mean
       NPI-Q total score."""
    d = df.loc[mask]
    pres = pd.to_numeric(d[item], errors='coerce')                 # presence 0/1
    sev = pd.to_numeric(d.get(item + 'sev'), errors='coerce') if (item + 'sev') in d.columns else pd.Series(np.nan, index=d.index)
    valid = pres.notna()
    contrib = pd.Series(np.where(pres.eq(1), sev.fillna(0.0), 0.0), index=d.index)[valid]
    return contrib.mean() if len(contrib) > 0 else np.nan

# ordered symptom positions with small gaps between domains
labels = []; ypos = []; items = []; dom_centers = []
y = 0.0; gap = 0.8
for dom, syms in DOMAINS:
    ys = []
    for it, lab in syms:
        labels.append(lab); items.append(it); ypos.append(y); ys.append(y); y += 1
    dom_centers.append((dom, sum(ys) / len(ys))); y += gap
ypos = np.array(ypos); ymax = y

# masks
m_sig = df['SIG']; m_norm = ~df['SIG']
m_dem_sig = df['DEM'] & df['SIG']; m_cn_sig = (~df['DEM']) & df['SIG']
n_dem_sig = int(m_dem_sig.sum()); n_cn_sig = int(m_cn_sig.sum())

prev_sig = [prevalence(m_sig, it) for it in items]
prev_norm = [prevalence(m_norm, it) for it in items]
sev_dsig = [sev_contrib(m_dem_sig, it) for it in items]
sev_csig = [sev_contrib(m_cn_sig, it) for it in items]

# ---- panel configs: (letter, title, valsA, valsB, colA, colB, legA, legB, metric) ----
CFG = [
 ('a', 'Severity contribution to significant NPS', sev_dsig, sev_csig, C_DEM, C_CN,
      f'Dementia (n={n_dem_sig})', f'Controls (n={n_cn_sig})', 'sev'),
 ('b', 'Significant vs Normal NPS', prev_sig, prev_norm, C_SIG, C_NORM, 'Significant NPS', 'Normal NPS', 'prev'),
]
fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 120/25.4), sharey=False)
bh = 0.38
for ax, (letter, title, pA, pB, cA, cB, gA, gB, metric) in zip(axes, CFG):
    ax.barh(ymax-1-ypos+bh/2, pA, height=bh, color=cA, zorder=3)
    ax.barh(ymax-1-ypos-bh/2, pB, height=bh, color=cB, zorder=3)
    ax.set_ylim(-1, ymax); ax.set_yticks(ymax-1-ypos)
    if ax is axes[0]:
        ax.set_yticklabels(labels, fontsize=6, va='center')
    else:
        ax.set_yticklabels([]); ax.tick_params(axis='y', length=0)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    if metric == 'prev':
        ax.set_xlim(0, 70); ax.set_xticks([0, 20, 40, 60]); ax.set_xlabel('Prevalence (%)', fontsize=6)
    else:
        ax.set_xlim(0, 1.3); ax.set_xticks([0, 0.4, 0.8, 1.2])
        ax.set_xlabel('Mean NPI-Q severity contribution', fontsize=6)
    ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontsize=8, fontweight='bold', ha='left', va='bottom')
    ax.text(0.075, 1.02, title, transform=ax.transAxes, fontsize=6.5, ha='left', va='bottom')
    handles = [Patch(facecolor=cA, edgecolor='black', linewidth=0.6, label=gA),
               Patch(facecolor=cB, edgecolor='black', linewidth=0.6, label=gB)]
    ax.legend(handles=handles, loc='lower right', fontsize=5.5, frameon=False)

# domain labels on the far left of panel a (offset scaled to panel a's x-range)
axL = axes[0]
dom_x = -0.72 * axL.get_xlim()[1]     # robust to the axis range (severity 0-1.3 or prevalence 0-70)
for dom, yc in dom_centers:
    dom_disp = dom.replace('/', '/\n') if len(dom) > 12 else dom
    axL.text(dom_x, ymax-1-yc, dom_disp, fontsize=6, ha='center', va='center', rotation=90, clip_on=False)

plt.subplots_adjust(left=0.26, right=0.99, top=0.90, bottom=0.10, wspace=0.08)
for ext in ['pdf', 'png', 'svg']:
    fig.savefig(f'Fig_symptom_prevalence.{ext}', dpi=450, facecolor='white')
print("Saved: Fig_symptom_prevalence.{pdf,png,svg}")
print(f"n significant: dementia={n_dem_sig}, controls={n_cn_sig}")
