# ============================================================================
# Mediation (causal, Imai) - FTD (WITH country, no Dx)
#   Mediators = 10 EAA measures that survived FDR in the FTD regression (with country).
#   Standardized ACME, 5000 bootstrap resamples.
#   BH-FDR WITHIN EACH CLOCK FAMILY (consistent with the regression):
#     First/Second (Horvath), Next Generation, Organs/Systems.
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: Mediation_NPS_ACME_FTD_country.xlsx
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
df['sexm'] = df['demo_sex'] - 1
df_ftd = df[df['clinical_diagnosis'] == 'FTD'].copy()

# (mediator_col, label, clock family)
MED_FTD = [('Horvath_EAA', 'Horvath', 'First/Second'),
           ('Retroclock_EAA', 'Retroclock', 'Next'), ('DNAmFitAge_EAA', 'DNAmFitAge', 'Next'),
           ('AdaptAge_EAA', 'AdaptAge', 'Next'), ('IntrinsicCapacity', 'IntrinsicCapacity', 'Next'),
           ('DamAge_EAA', 'DamAge', 'Next'),
           ('Blood_EAA', 'Blood', 'Organs'), ('Brain_EAA', 'Brain', 'Organs'),
           ('Liver_EAA', 'Liver', 'Organs'), ('Metabolic_EAA', 'Metabolic', 'Organs')]
OUTCOMES = [('mmse_total', 'MMSE (cognition)'), ('cdr_cdrtot', 'CDR (dementia severity)'),
            ('t_adlq_tot', 'T-ADLQ (function)')]
COV = 'demo_age + sexm + demo_smoke + cog_ed + meds_nps + C(demo_resid)'
cov_vars = ['demo_age', 'sexm', 'demo_smoke', 'cog_ed', 'meds_nps', 'demo_resid']

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
print(f"\n################ FTD (with country, no Dx)  |  COV: {COV} ################")
for mcol, mname, fam in MED_FTD:
    if mcol not in df_ftd.columns:
        print(f"  [missing] {mcol}"); continue
    for ycol, yname in OUTCOMES:
        r = mediate_one(df_ftd, mcol, ycol, COV, cov_vars)
        if r is None:
            print(f"  [skipped] {mname} -> {yname}: N<{MIN_N}"); continue
        rows.append({'Block':'FTD', 'Family':fam, 'Mediator':mname, 'Outcome':yname, 'N':r['N'],
            'ACME (std)':round(r['ACME'],4), '95% CI':f"[{r['ACME_lo']:.4f}, {r['ACME_hi']:.4f}]",
            'ACME p':round(r['ACME_p'],4), 'ADE (std)':round(r['ADE'],4), 'ADE p':round(r['ADE_p'],4),
            'Total (std)':round(r['TOT'],4),
            'Prop. mediated':round(r['PROP'],3) if pd.notna(r['PROP']) else 'NA', '_p':r['ACME_p']})
        print(f"  [{fam}] {mname} -> {yname}: N={r['N']} ACMEstd={r['ACME']:.4f} "
              f"CI[{r['ACME_lo']:.4f},{r['ACME_hi']:.4f}] p={r['ACME_p']:.4f}")

res = pd.DataFrame(rows)
# BH-FDR WITHIN EACH CLOCK FAMILY
res['ACME P-FDR'] = np.nan
for fam in res['Family'].unique():
    m = res['Family'] == fam
    res.loc[m, 'ACME P-FDR'] = multipletests(res.loc[m, '_p'].values, method='fdr_bh')[1]
res['ACME P-FDR'] = res['ACME P-FDR'].round(4)
res = res.drop(columns=['_p'])
print("\n========== FTD, WITH COUNTRY (BH-FDR within clock family) ==========")
print(res.to_string(index=False))
res.to_excel('Mediation_NPS_ACME_FTD_country.xlsx', index=False)
print("\nSaved: Mediation_NPS_ACME_FTD_country.xlsx")
