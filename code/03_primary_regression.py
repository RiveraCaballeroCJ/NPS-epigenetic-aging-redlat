# ============================================================================
# EAA ~ NPS STATUS (binary) - PRIMARY MODEL: OLS + country fixed effects
# Covariates: Age, Sex, Smoking, Education, Medication_NPS, Diagnosis (dementia vs control)
# Country = fixed effect (dummies, Colombia reference). NO mixed model.
#   _complete   = Table 2 (beta std | 95% CI std | P | P-FDR | Cohen's d | Partial R2)
#   _covariates = SUPPLEMENTARY: every predictor with beta std | 95% CI std | P | P-FDR
#                 | Cohen's f2 | Partial R2 + model R2  (Country = block f2/Partial R2)
#   PART 1 VIF | PART 2 PRIMARY full | PART 3 STRATIFIED AD/FTD/CN | FDR within family
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: EAA_primary_OLS_Dx.xlsx
# ============================================================================
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant
from scipy.stats import t as t_dist
import warnings; warnings.filterwarnings('ignore')

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'
VALID_COUNTRIES = [1, 3, 4, 5, 6]

df = pd.read_csv(DATA_PATH)
df = df[df['clinical_diagnosis'].isin(['AD','FTD','CN'])].copy()
df = df[df['demo_resid'].isin(VALID_COUNTRIES)].copy()
country_map = {1:'Argentina', 3:'Chile', 4:'Colombia', 5:'Mexico', 6:'Peru'}
df['country_name'] = df['demo_resid'].map(country_map)
df['dx_binary'] = df['clinical_diagnosis'].isin(['AD','FTD']).astype(int)
df['nps_binary'] = df['npi_total'].apply(lambda x: (1 if x >= 3 else 0) if pd.notna(x) else np.nan)
df_ad  = df[df['clinical_diagnosis']=='AD'].copy()
df_ftd = df[df['clinical_diagnosis']=='FTD'].copy()
df_cn  = df[df['clinical_diagnosis']=='CN'].copy()
clock_families = {
    'First/Second Generation Clocks': [
        ('Horvath_EAA','Horvath'), ('Hannum_EAA','Hannum'), ('cAge_EAA','cAge'),
        ('PhenoAge_EAA','PhenoAge'), ('PCGrimAge_EAA','PCGrimAge')],
    'Next Generation Clocks': [
        ('IntrinClock_EAA','IntrinClock'), ('Retroclock_EAA','Retroclock'),
        ('Retroclockv2_EAA','Retroclock v2'), ('OMICmAge_EAA','OMICmAge'), ('DNAmFitAge_EAA','DNAmFitAge'),
        ('AdaptAge_EAA','AdaptAge'), ('DamAge_EAA','DamAge'),
        ('IntrinsicCapacity','IntrinsicCapacity'), ('DunedinPACE','DunedinPACE')],
    'Organs/Systems Clocks': [
        ('Blood_EAA','Blood'), ('Brain_EAA','Brain'), ('Heart_EAA','Heart'), ('Lung_EAA','Lung'),
        ('Liver_EAA','Liver'), ('Kidney_EAA','Kidney'), ('Immune_EAA','Immune'),
        ('Inflammation_EAA','Inflammation'), ('Hormone_EAA','Hormone'), ('Metabolic_EAA','Metabolic'),
        ('MusculoSkeletal_EAA','MusculoSkeletal'), ('SystemsAge_EAA','SystemsAge')],
}
clock_names = [c for fam in clock_families.values() for c, _ in fam]
clock_to_family = {c: fam for fam, lst in clock_families.items() for c, _ in lst}
predictor    = 'nps_binary'
COVARS_FULL  = ['demo_age','demo_sex','demo_smoke','cog_ed','meds_nps','dx_binary']
COVARS_STRAT = ['demo_age','demo_sex','demo_smoke','cog_ed','meds_nps']
CTRY = "C(country_name, Treatment(reference='Colombia'))"
pretty = {'demo_age':'Age','demo_sex':'Sex','demo_smoke':'Smoking','cog_ed':'Education',
          'meds_nps':'Medication_NPS','dx_binary':'Diagnosis','nps_binary':'NPS_status'}
def pretty_term(name):
    if name.startswith('C(country_name'):
        for code, nm in country_map.items():
            if f'T.{nm}' in name: return f'Country: {nm} (vs Colombia)'
        return 'Country'
    return pretty.get(name, name)
def p_star(p):
    if pd.isna(p): return 'NA'
    s = '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else ''
    return f"{p:.3f}{s}"
def pooled_sd(a, b):
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2: return np.nan
    return np.sqrt(((n1-1)*a.var(ddof=1) + (n2-1)*b.var(ddof=1)) / (n1+n2-2))
def fdr_within_family(dfr):
    dfr = dfr.copy(); dfr['p_FDR'] = np.nan
    dfr['Family'] = dfr['Clock'].map(clock_to_family)
    for fam in dfr['Family'].dropna().unique():
        idx = dfr.index[dfr['Family']==fam]
        pv = dfr.loc[idx,'p_nominal'].values
        if len(pv): dfr.loc[idx,'p_FDR'] = multipletests(pv, method='fdr_bh')[1]
    return dfr
def _formula(clock, fixed):
    return f"{clock} ~ " + " + ".join(fixed) + f" + {CTRY}"
def fit_ols(data, clock, fixed):
    cols = list(fixed) + [clock, 'country_name']
    sub = data[cols].copy(); sub['demo_sex'] = sub['demo_sex'] - 1
    for c in list(fixed) + [clock]:
        sub[c] = pd.to_numeric(sub[c], errors='coerce')
    sub = sub.dropna()
    if len(sub) < 20 or sub[predictor].nunique() < 2:
        return None, None
    return smf.ols(_formula(clock, fixed), data=sub).fit(), sub
def compute_es(sub, clock, fixed, m):
    """Per-term Cohen's f2 and partial R2; Country evaluated as a block."""
    r2f = m.rsquared; es = {}
    for v in fixed:
        red = smf.ols(f"{clock} ~ " + " + ".join([f for f in fixed if f != v]) + f" + {CTRY}", data=sub).fit()
        partial = r2f - red.rsquared
        es[pretty.get(v, v)] = {'f2': partial/(1-r2f) if (1-r2f) > 0 else np.nan, 'pr2': partial}
    red = smf.ols(f"{clock} ~ " + " + ".join(fixed), data=sub).fit()   # drop country block
    partial = r2f - red.rsquared
    es['Country'] = {'f2': partial/(1-r2f) if (1-r2f) > 0 else np.nan, 'pr2': partial}
    return es
def run_all(data, covars):
    fixed = [predictor] + covars
    rows, coefs, es_by_clock = [], [], {}
    for clock in clock_names:
        if clock not in data.columns: continue
        m, sub = fit_ols(data, clock, fixed)
        if m is None: continue
        b, se, p = m.params[predictor], m.bse[predictor], m.pvalues[predictor]
        ci = m.conf_int().loc[predictor]; ci_l, ci_u = float(ci[0]), float(ci[1])
        sd_y = sub[clock].std(ddof=1)
        scale = (sub[predictor].std(ddof=1)/sd_y) if sd_y > 0 else np.nan
        beta_std = b*scale; beta_ci_l, beta_ci_u = ci_l*scale, ci_u*scale
        sp = pooled_sd(sub.loc[sub[predictor]==1, clock], sub.loc[sub[predictor]==0, clock])
        cohen_d = b/sp if (sp and sp > 0) else np.nan
        es = compute_es(sub, clock, fixed, m); es_by_clock[clock] = es
        rec = es['NPS_status']
        rows.append({'Clock':clock, 'N':len(sub),
                     'beta_std':beta_std, 'beta_CI_L':beta_ci_l, 'beta_CI_U':beta_ci_u,
                     'Cohen_d':cohen_d, 'p_nominal':p, 'R2':m.rsquared, 'R2_adj':m.rsquared_adj,
                     'f2_NPS':rec['f2'], 'partialR2_NPS':rec['pr2']})
        # ---- standardized beta + standardized CI for EVERY term (incl. country dummies) ----
        exog_names = list(m.model.exog_names)
        exog = np.asarray(m.model.exog)
        tcrit = t_dist.ppf(0.975, int(m.df_resid))
        for name in m.params.index:
            if name == 'Intercept': continue
            j = exog_names.index(name)
            sd_x = np.std(exog[:, j], ddof=1)
            sc = (sd_x/sd_y) if (sd_y > 0 and np.isfinite(sd_x)) else np.nan
            bb, ss = m.params[name], m.bse[name]
            coefs.append({'Clock':clock, 'Predictor':pretty_term(name),
                          'beta_std':bb*sc,
                          'beta_CI_L_std':(bb - tcrit*ss)*sc,
                          'beta_CI_U_std':(bb + tcrit*ss)*sc,
                          'p_num':m.pvalues[name], 'R2':round(m.rsquared,4)})
    dfr = pd.DataFrame(rows)
    if len(dfr): dfr = fdr_within_family(dfr)
    coefs_df = pd.DataFrame(coefs)
    if len(coefs_df):
        # FDR WITHIN EACH CLOCK FAMILY, per predictor (consistent with the main table)
        coefs_df['Family'] = coefs_df['Clock'].map(clock_to_family)
        coefs_df['p_FDR'] = np.nan
        for nm in coefs_df['Predictor'].unique():
            for fam in coefs_df['Family'].dropna().unique():
                idx = coefs_df.index[(coefs_df['Predictor']==nm) & (coefs_df['Family']==fam)]
                if len(idx):
                    coefs_df.loc[idx,'p_FDR'] = multipletests(
                        coefs_df.loc[idx,'p_num'].values, method='fdr_bh')[1]
    return dfr, coefs_df, es_by_clock
# ---- Table 2 (trimmed): beta std | 95% CI std | P | P-FDR | Cohen's d | Partial R2 ----
def build_complete(dfr):
    d = dfr.set_index('Clock'); out = []
    for fam, clocks in clock_families.items():
        out.append({'Epigenetic Age Acceleration': fam})
        for clock, lab in clocks:
            if clock not in d.index: continue
            r = d.loc[clock]
            out.append({'Epigenetic Age Acceleration':lab,
                        'β':round(r['beta_std'],3),
                        '95% CI (std)':f"[{r['beta_CI_L']:.3f}, {r['beta_CI_U']:.3f}]",
                        'P-value':p_star(r['p_nominal']), 'P-FDR':p_star(r['p_FDR']),
                        "Cohen's d":round(r['Cohen_d'],3),
                        'Partial R2':round(r['partialR2_NPS'],4) if pd.notna(r['partialR2_NPS']) else np.nan})
    return pd.DataFrame(out)
# ---- SUPPLEMENTARY = covariates: beta std | 95% CI std | P | P-FDR | f2 | Partial R2 + model R2 ----
def build_covariates(coefs_df, es):
    if not len(coefs_df): return coefs_df
    name_map = {c: lab for fam in clock_families.values() for c, lab in fam}
    out = coefs_df.copy()
    out['Epigenetic Age Acceleration'] = out['Clock'].map(name_map).fillna(out['Clock'])
    out['β (std)'] = out['beta_std'].round(3)
    out['95% CI (std)'] = out.apply(lambda r: f"[{r['beta_CI_L_std']:.3f}, {r['beta_CI_U_std']:.3f}]", axis=1)
    out['P-value'] = out['p_num'].apply(p_star); out['P-FDR'] = out['p_FDR'].apply(p_star)
    def es_get(r, key):
        rec = es.get(r['Clock'], {})
        rec = rec.get('Country') if str(r['Predictor']).startswith('Country:') else rec.get(r['Predictor'])
        if not rec: return ''
        v = rec.get(key)
        return '' if (v is None or (isinstance(v, float) and pd.isna(v))) else round(v, 4)
    out["Cohen's f2"] = out.apply(lambda r: es_get(r, 'f2'), axis=1)
    out['Partial R2'] = out.apply(lambda r: es_get(r, 'pr2'), axis=1)
    out = out[['Epigenetic Age Acceleration','Predictor','β (std)','95% CI (std)','P-value','P-FDR',
               "Cohen's f2",'Partial R2','R2']].copy()
    out = out.rename(columns={'R2':'R2 model'})
    dup = out['Epigenetic Age Acceleration'].duplicated()
    out['Epigenetic Age Acceleration'] = out['Epigenetic Age Acceleration'].astype(object)
    out['R2 model'] = out['R2 model'].astype(object)
    out.loc[dup, 'Epigenetic Age Acceleration'] = ''; out.loc[dup, 'R2 model'] = ''
    return out
def vif_table(data, covars):
    fixed = [predictor] + covars
    d = data[fixed + ['country_name']].dropna().copy(); d['demo_sex'] = d['demo_sex'] - 1
    for c in fixed: d[c] = pd.to_numeric(d[c], errors='coerce')
    d = d.dropna(); design = d[fixed].astype(float).copy(); lab = dict(pretty)
    for code, nm in country_map.items():
        if nm != 'Colombia' and (d['country_name']==nm).any():
            design[f'C_{nm}'] = (d['country_name']==nm).astype(float); lab[f'C_{nm}'] = f'Country: {nm} (vs Colombia)'
    X = add_constant(design); out = []
    for i, c in enumerate(X.columns):
        if c == 'const': continue
        try:
            v = variance_inflation_factor(X.values, i); v = round(float(v),3) if np.isfinite(v) else np.nan
        except Exception:
            v = np.nan
        out.append({'Predictor':lab.get(c,c), 'VIF':v})
    return pd.DataFrame(out)
# PART 1 - VIF
print("="*95 + "\nPART 1 - VIF\n" + "="*95)
vif_all = []
for g, dd, cv in [('Full',df,COVARS_FULL),('AD',df_ad,COVARS_STRAT),
                  ('FTD',df_ftd,COVARS_STRAT),('CN',df_cn,COVARS_STRAT)]:
    v = vif_table(dd, cv); v.insert(0,'Group',g); vif_all.append(v)
vif_all = pd.concat(vif_all, ignore_index=True); print(vif_all.to_string(index=False), "\n")
# PART 2 - PRIMARY (full sample)
print("="*95 + "\nPART 2 - PRIMARY (full sample; OLS + country FE + Dx; FDR within family)\n" + "="*95)
res_full, coef_full, es_full = run_all(df, COVARS_FULL)
t = res_full.copy(); t['P-value']=t['p_nominal'].apply(p_star); t['P-FDR']=t['p_FDR'].apply(p_star)
print(t.sort_values('p_nominal')[['Clock','Family','N','beta_std','beta_CI_L','beta_CI_U','f2_NPS','Cohen_d','P-value','P-FDR']].round(4).to_string(index=False))
print(f"\n>> FDR-sig: {(res_full['p_FDR']<.05).sum()} -> " + ", ".join(res_full[res_full['p_FDR']<.05]['Clock'].tolist()) + "\n")
# PART 3 - STRATIFIED (Dx dropped: constant within stratum)
print("="*95 + "\nPART 3 - STRATIFIED (OLS + country FE, no Dx; FDR within family)\n" + "="*95)
strat = {'AD':df_ad, 'FTD':df_ftd, 'CN':df_cn}
res_strat, coef_strat, es_strat = {}, {}, {}
for g, dd in strat.items():
    print(f"\n--- {g} (n={len(dd)}) ---")
    r, c, e = run_all(dd, COVARS_STRAT)
    res_strat[g], coef_strat[g], es_strat[g] = r, c, e
    if len(r):
        tt = r.copy(); tt['P-value']=tt['p_nominal'].apply(p_star); tt['P-FDR']=tt['p_FDR'].apply(p_star)
        print(tt.sort_values('p_nominal')[['Clock','N','beta_std','beta_CI_L','beta_CI_U','Cohen_d','P-value','P-FDR']].round(4).to_string(index=False))
        print(f">> FDR-sig: {(r['p_FDR']<.05).sum()} -> " + ", ".join(r[r['p_FDR']<.05]['Clock'].tolist()))
with pd.ExcelWriter('EAA_primary_OLS_Dx.xlsx', engine='openpyxl') as w:
    vif_all.to_excel(w, sheet_name='VIF_all', index=False)
    build_complete(res_full).to_excel(w, sheet_name='PRIMARY_complete', index=False)
    build_covariates(coef_full, es_full).to_excel(w, sheet_name='PRIMARY_covariates', index=False)
    for g in strat:
        if len(res_strat[g]):
            build_complete(res_strat[g]).to_excel(w, sheet_name=f'S_{g}_complete', index=False)
            build_covariates(coef_strat[g], es_strat[g]).to_excel(w, sheet_name=f'S_{g}_covariates', index=False)
print("\nSaved: EAA_primary_OLS_Dx.xlsx")
