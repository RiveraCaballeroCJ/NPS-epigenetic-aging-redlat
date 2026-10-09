# ============================================================================
# Mediation (causal, Imai) - FULL SAMPLE only (WITH country)
#   Mediator = Brain EAA (the only FDR-significant clock in the full sample).
#   Standardized ACME, 5000 bootstrap resamples. BH-FDR within the block.
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: Mediation_NPS_ACME_FULL_country.xlsx
# ============================================================================
import pandas as pd, numpy as np
import statsmodels.api as sm
from statsmodels.stats.mediation import Mediation
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings('ignore')

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'
N_BOOT = 5000; SEED = 2025; MIN_N = 20

df = pd.read_csv(DATA_PATH)
df = df[df['clinical_diagnosis'].isin(['AD','FTD','CN'])].copy()
df = df[df['demo_resid'].isin([1,3,4,5,6])].copy()
df['NPS'] = (df['npi_total'] >= 3).astype(float); df.loc[df['npi_total'].isna(), 'NPS'] = np.nan
df['DX'] = df['clinical_diagnosis'].isin(['AD','FTD']).astype(float)
df['sexm'] = df['demo_sex'] - 1

MED_FULL = [('Brain_EAA', 'Brain EAA')]
OUTCOMES = [('mmse_total', 'MMSE (cognition)'), ('cdr_cdrtot', 'CDR (dementia severity)'),
            ('t_adlq_tot', 'T-ADLQ (function)')]
COV = 'demo_age + sexm + demo_smoke + cog_ed + meds_nps + DX + C(demo_resid)'
cov_vars = ['demo_age', 'sexm', 'demo_smoke', 'cog_ed', 'meds_nps', 'DX', 'demo_resid']

def mediate_one(data, mcol, ycol, COV, cov_vars):
    d = data.dropna(subset=['NPS', mcol, ycol] + cov_vars).copy()
    if len(d) < MIN_N:
        return None
    d['M'] = (d[mcol] - d[mcol].mean()) / d[mcol].std(ddof=0)
    d['Y'] = (d[ycol] - d[ycol].mean()) / d[ycol].std(ddof=0)
    med_model = sm.OLS.from_formula(f'M ~ NPS + {COV}', data=d)
    out_model = sm.OLS.from_formula(f'Y ~ NPS + M + {COV}', data=d)
    np.random.seed(SEED)
    med = Mediation(out_model, med_model, 'NPS', 'M').fit(method='bootstrap', n_rep=N_BOOT)
    s = med.summary()
    def g(row, col):
        try:
            return float(s.loc[row, col])
        except Exception:
            return np.nan
    return dict(N=len(d),
        ACME=g('ACME (average)', 'Estimate'), ACME_lo=g('ACME (average)', 'Lower CI bound'),
        ACME_hi=g('ACME (average)', 'Upper CI bound'), ACME_p=g('ACME (average)', 'P-value'),
        ADE=g('ADE (average)', 'Estimate'), ADE_p=g('ADE (average)', 'P-value'),
        TOT=g('Total effect', 'Estimate'), PROP=g('Prop. mediated (average)', 'Estimate'))

rows = []
print(f"\n################ Full sample (with country)  |  COV: {COV} ################")
for mcol, mname in MED_FULL:
    if mcol not in df.columns:
        print(f"  [missing] {mcol}"); continue
    for ycol, yname in OUTCOMES:
        r = mediate_one(df, mcol, ycol, COV, cov_vars)
        if r is None:
            print(f"  [skipped] {mname} -> {yname}: N<{MIN_N}"); continue
        rows.append({'Block':'Full sample', 'Mediator':mname, 'Outcome':yname, 'N':r['N'],
            'ACME (std)':round(r['ACME'],4), '95% CI':f"[{r['ACME_lo']:.4f}, {r['ACME_hi']:.4f}]",
            'ACME p':round(r['ACME_p'],4), 'ADE (std)':round(r['ADE'],4), 'ADE p':round(r['ADE_p'],4),
            'Total (std)':round(r['TOT'],4),
            'Prop. mediated':round(r['PROP'],3) if pd.notna(r['PROP']) else 'NA', '_p':r['ACME_p']})
        print(f"  {mname} -> {yname}: N={r['N']} ACMEstd={r['ACME']:.4f} "
              f"CI[{r['ACME_lo']:.4f},{r['ACME_hi']:.4f}] p={r['ACME_p']:.4f}")

res = pd.DataFrame(rows)
res['ACME P-FDR'] = multipletests(res['_p'].values, method='fdr_bh')[1].round(4)
res = res.drop(columns=['_p'])
print("\n================= FULL SAMPLE, WITH COUNTRY (BH-FDR within block) =================")
print(res.to_string(index=False))
res.to_excel('Mediation_NPS_ACME_FULL_country.xlsx', index=False)
print("\nSaved: Mediation_NPS_ACME_FULL_country.xlsx")
