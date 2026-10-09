# ============================================================================
# Sociodemographic and clinical characterization by NPS status
# NPS = clinically significant neuropsychiatric symptoms (NPI-Q total >= 3)
#
# - ReDLat countries: [1, 3, 4, 5, 6] (Argentina, Chile, Colombia, Mexico, Peru)
# - Global cohort summary (total N) + breakdown by country
# - Characterization tables, Significant vs Normal NPS:
#       Full sample and diagnosis-stratified (AD, FTD, CN)
# - Shapiro-Wilk + Levene reported for reproducibility; Welch t = primary test
# - NPI-Q: descriptive only, no test (grouping derived from NPI = circular)
# - Categorical: chi-square; Fisher (2x2) or Fisher-Freeman-Halton (rxc, Monte
#   Carlo) when any expected cell < 5
# - Available (non-missing) N reported per variable
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: NPS_characterization_tables.xlsx
# ============================================================================
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2_contingency, ttest_ind, levene, shapiro, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 1. LOAD AND PREPARE
# ============================================================================
print("="*80)
print("CHARACTERIZATION BY NPS STATUS — SIGNIFICANT vs NORMAL")
print("="*80 + "\n")

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'
df = pd.read_csv(DATA_PATH)
df = df[df['clinical_diagnosis'].isin(['AD', 'FTD', 'CN'])].copy()

# ReDLat countries (cohort definition used across all analyses)
VALID_COUNTRIES = [1, 3, 4, 5, 6]
df = df[df['demo_resid'].isin(VALID_COUNTRIES)].copy()

# Groups and binary variables
df['dementia'] = df['clinical_diagnosis'].isin(['AD', 'FTD']).astype(int)  # 1 = dementia, 0 = control

def categorize_nps(npi_score):
    if pd.isna(npi_score):
        return np.nan
    return 'Significant' if npi_score >= 3 else 'Normal'
df['nps_group'] = df['npi_total'].apply(categorize_nps)

# Label maps
countries = {1: 'Argentina', 3: 'Chile', 4: 'Colombia', 5: 'Mexico', 6: 'Peru'}
sex_map = {1: 'Female', 2: 'Male'}
smoking_map = {0: 'No', 1: 'Yes'}
df['country'] = df['demo_resid'].map(countries)
df['sex'] = df['demo_sex'].map(sex_map)
df['cognitive_status'] = df['clinical_diagnosis'].map({'AD': 'AD', 'FTD': 'FTD', 'CN': 'CN'})
df['smoking_status'] = df['demo_smoke'].map(smoking_map)

# Diagnostic subgroups (for counts)
df_ad = df[df['clinical_diagnosis'] == 'AD'].copy()
df_ftd = df[df['clinical_diagnosis'] == 'FTD'].copy()
df_cn = df[df['clinical_diagnosis'] == 'CN'].copy()
print(f"Total cohort (N): {len(df)}")
print(f"  AD: {len(df_ad)} | FTD: {len(df_ftd)} | CN: {len(df_cn)}\n")

# ============================================================================
# 2. GLOBAL COHORT SUMMARY (TOTAL N) + BREAKDOWN BY COUNTRY
# ============================================================================
print("="*80)
print("GLOBAL COHORT SUMMARY")
print("="*80 + "\n")
N_total = len(df)
age_mean = df['demo_age'].mean()
age_sd = df['demo_age'].std()
n_female = (df['sex'] == 'Female').sum()
pct_female = n_female / N_total * 100
n_dementia = (df['dementia'] == 1).sum()
pct_dementia = n_dementia / N_total * 100
print(f"Total N: {N_total}")
print(f"Age: {age_mean:.1f} +/- {age_sd:.1f} years (mean +/- SD)")
print(f"Women: {n_female} ({pct_female:.1f}%)")
print(f"Dementia: {n_dementia} ({pct_dementia:.1f}%)\n")

print("-"*80)
print("BREAKDOWN BY COUNTRY")
print("-"*80 + "\n")
global_rows = []
for code in VALID_COUNTRIES:
    cname = countries[code]
    sub = df[df['demo_resid'] == code]
    n_c = len(sub)
    if n_c == 0:
        continue
    pct_c = n_c / N_total * 100
    n_dem_c = (sub['dementia'] == 1).sum()
    n_fem_c = (sub['sex'] == 'Female').sum()
    age_m_c = sub['demo_age'].mean()
    age_sd_c = sub['demo_age'].std()
    edu_m_c = sub['cog_ed'].mean()
    edu_sd_c = sub['cog_ed'].std()
    print(f"{cname}:")
    print(f"  N = {n_c} ({pct_c:.1f}% of total)")
    print(f"  Dementia (n): {n_dem_c}")
    print(f"  Women (n): {n_fem_c}")
    print(f"  Age: {age_m_c:.1f} +/- {age_sd_c:.1f} years")
    print(f"  Education: {edu_m_c:.1f} +/- {edu_sd_c:.1f} years\n")
    global_rows.append({
        'Country': cname, 'N': n_c, '%_of_total': round(pct_c, 1),
        'Dementia_n': int(n_dem_c), 'Female_n': int(n_fem_c),
        'Age_mean': round(age_m_c, 1), 'Age_SD': round(age_sd_c, 1),
        'Education_mean': round(edu_m_c, 1), 'Education_SD': round(edu_sd_c, 1),
    })
df_global = pd.DataFrame(global_rows)
df_global_total = pd.DataFrame([{
    'Country': 'TOTAL', 'N': N_total, '%_of_total': 100.0,
    'Dementia_n': int(n_dementia), 'Female_n': int(n_female),
    'Age_mean': round(age_mean, 1), 'Age_SD': round(age_sd, 1),
    'Education_mean': round(df['cog_ed'].mean(), 1),
    'Education_SD': round(df['cog_ed'].std(), 1),
}])
df_global = pd.concat([df_global, df_global_total], ignore_index=True)
print("Global summary table:")
print(df_global.to_string(index=False))
print()

# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================
def cohens_d(g1, g2):
    n1, n2 = len(g1), len(g2)
    if n1 < 2 or n2 < 2:
        return np.nan
    v1, v2 = np.var(g1, ddof=1), np.var(g2, ddof=1)
    pooled = np.sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2))
    return np.nan if pooled == 0 else (np.mean(g1) - np.mean(g2)) / pooled

def cramers_v(chi2, n, min_dim):
    if n == 0 or min_dim <= 1:
        return np.nan
    return np.sqrt(chi2 / (n * (min_dim - 1)))

def pval_stars(p):
    if pd.isna(p):
        return "N/A"
    s = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
    return f"{p:.3f}{s}"

def categorical_test(table):
    """
    Return (stat_str, p, effect_str, test_name).
    Chi-square or Fisher depending on expected cells < 5.
    - 2x2: Fisher exact.
    - rxc: Fisher-Freeman-Halton via chi2 with Monte Carlo simulation.
    """
    if table.size == 0 or table.shape[0] < 2 or table.shape[1] < 2:
        return "N/A", np.nan, "N/A", "none"
    chi2, p_chi, dof, expected = chi2_contingency(table)
    n_total = table.values.sum()
    min_dim = min(table.shape)
    v = cramers_v(chi2, n_total, min_dim)
    min_expected = expected.min()
    if min_expected >= 5:
        return f"chi2 = {chi2:.3f}, df = {dof}", p_chi, f"V = {v:.3f}", "chi2"
    # Some expected cell < 5 -> exact test
    if table.shape == (2, 2):
        _, p_fisher = fisher_exact(table.values)
        return f"Fisher exact (min E={min_expected:.2f})", p_fisher, f"V = {v:.3f}", "fisher"
    else:
        # rxc: Fisher-Freeman-Halton approximated by Monte Carlo
        try:
            res = stats.chi2_contingency(table, method="monte_carlo", n_resamples=10000)
            p_mc = res.pvalue
            return (f"chi2={chi2:.3f} (MC, min E={min_expected:.2f})",
                    p_mc, f"V = {v:.3f}", "montecarlo")
        except Exception:
            return (f"chi2={chi2:.3f}, df={dof} (min E={min_expected:.2f}, approx.)",
                    p_chi, f"V = {v:.3f}", "chi2_warn")

def continuous_test(g1, g2):
    """
    Run Shapiro (per group) and Levene (reproducibility); report Welch t as the
    primary test. Returns a dict with all quantities.
    """
    g1 = pd.Series(g1).dropna()
    g2 = pd.Series(g2).dropna()
    out = {'n1': len(g1), 'n2': len(g2)}

    def _shapiro(x):
        if 3 <= len(x) <= 5000:
            try:
                return shapiro(x).pvalue
            except Exception:
                return np.nan
        return np.nan
    out['shapiro_p_g1'] = _shapiro(g1)
    out['shapiro_p_g2'] = _shapiro(g2)
    try:
        out['levene_p'] = levene(g1, g2).pvalue
    except Exception:
        out['levene_p'] = np.nan
    try:
        t, p = ttest_ind(g1, g2, equal_var=False)  # Welch (primary)
        out['t'], out['p'] = t, p
    except Exception:
        out['t'], out['p'] = np.nan, np.nan
    out['d'] = cohens_d(g1.values, g2.values)
    return out

# ============================================================================
# 4. CHARACTERIZATION TABLE — Significant vs Normal NPS
# ============================================================================
def build_row(char, g1v='', g2v='', ts='', p='', es='', ndisp=''):
    return {'Characteristic': char, 'Significant': g1v, 'Normal': g2v,
            'Test': ts, 'p-value': p, 'Effect Size': es, 'N_available': ndisp}

def continuous_block(rows, label, g1, g2, do_test=True):
    """Add Mean (SD) for a continuous variable. do_test=False for NPI (circular)."""
    s1, s2 = g1.dropna(), g2.dropna()
    rows.append(build_row(label))
    header_idx = len(rows) - 1
    rows.append(build_row('Mean (SD)',
                          f'{s1.mean():.1f} ({s1.std():.1f})' if len(s1) else 'NA',
                          f'{s2.mean():.1f} ({s2.std():.1f})' if len(s2) else 'NA',
                          ndisp=f'{len(s1)}/{len(s2)}'))
    if do_test:
        res = continuous_test(g1, g2)
        rows[-1]['Test'] = f"t = {res['t']:.3f}"
        rows[-1]['p-value'] = pval_stars(res['p'])
        rows[-1]['Effect Size'] = f"d = {res['d']:.3f}"
        rows[header_idx]['Test'] = (f"Shapiro p: {res['shapiro_p_g1']:.3f}/"
                                    f"{res['shapiro_p_g2']:.3f}; "
                                    f"Levene p: {res['levene_p']:.3f}")
    else:
        rows[-1]['Test'] = 'by design (no test)'
    return rows

def categorical_block(rows, label, df_sub, group_col, colname, categories):
    rows.append(build_row(label))
    header_idx = len(rows) - 1
    dfg = df_sub[df_sub[group_col].isin(['Significant', 'Normal'])]
    n1 = (dfg[group_col] == 'Significant').sum()
    n2 = (dfg[group_col] == 'Normal').sum()
    for cat in categories:
        c1 = ((dfg[group_col] == 'Significant') & (dfg[colname] == cat)).sum()
        c2 = ((dfg[group_col] == 'Normal') & (dfg[colname] == cat)).sum()
        rows.append(build_row(f'{cat} (N, %)',
                              f'{c1} ({c1/n1*100:.1f}%)' if n1 else '0',
                              f'{c2} ({c2/n2*100:.1f}%)' if n2 else '0'))
    table = pd.crosstab(dfg[group_col], dfg[colname])
    n_avail = int(dfg[colname].notna().sum())
    ts, p, es, _ = categorical_test(table)
    rows[header_idx]['Test'] = ts
    rows[header_idx]['p-value'] = pval_stars(p)
    rows[header_idx]['Effect Size'] = es
    rows[header_idx]['N_available'] = f'{n_avail}'
    return rows

def create_characterization_table(df_subset, title):
    print(f"\n{'='*80}\n{title}\n{'='*80}")
    dfg = df_subset[df_subset['nps_group'].isin(['Significant', 'Normal'])].copy()
    a = dfg[dfg['nps_group'] == 'Significant']
    n = dfg[dfg['nps_group'] == 'Normal']
    na, nn = len(a), len(n)
    print(f"Significant: N={na} | Normal: N={nn}\n")
    rows = [build_row(f'Total (N={na+nn})',
                      f'{na} ({na/(na+nn)*100:.1f}%)' if (na+nn) else '0',
                      f'{nn} ({nn/(na+nn)*100:.1f}%)' if (na+nn) else '0')]
    continuous_block(rows, 'NPI-Q Total Score', a['npi_total'], n['npi_total'], do_test=False)
    categorical_block(rows, 'Sex', dfg, 'nps_group', 'sex', ['Female', 'Male'])
    continuous_block(rows, 'Age (years)', a['demo_age'], n['demo_age'], do_test=True)
    categorical_block(rows, 'Cognitive Status', dfg, 'nps_group', 'cognitive_status', ['AD', 'FTD', 'CN'])
    continuous_block(rows, 'Education (years)', a['cog_ed'], n['cog_ed'], do_test=True)
    continuous_block(rows, 'MMSE Total Score', a['mmse_total'], n['mmse_total'], do_test=True)
    continuous_block(rows, 'T-ADLQ Functionality', a['t_adlq_tot'], n['t_adlq_tot'], do_test=True)
    categorical_block(rows, 'Smoking History', dfg, 'nps_group', 'smoking_status', ['Yes', 'No'])
    present_countries = [countries[c] for c in VALID_COUNTRIES if c in dfg['demo_resid'].values]
    categorical_block(rows, 'Country', dfg, 'nps_group', 'country', present_countries)
    dfg['med_binary'] = dfg['meds_nps'].apply(lambda x: 'Yes' if pd.notna(x) and x >= 1
                                              else ('No' if x == 0 else np.nan))
    categorical_block(rows, 'NPS Medication', dfg, 'nps_group', 'med_binary', ['Yes', 'No'])
    table = pd.DataFrame(rows)
    print(table.to_string(index=False))
    return table

# ============================================================================
# 5. BUILD TABLES AND SAVE
#    Full sample + diagnosis-stratified (AD, FTD, CN). No dementia/controls split.
# ============================================================================
excel_file = 'NPS_characterization_tables.xlsx'
writer = pd.ExcelWriter(excel_file, engine='openpyxl')

df_global.to_excel(writer, sheet_name='Global_by_Country', index=False)

t_full = create_characterization_table(
    df, 'FULL SAMPLE - Clinically Significant vs Normal NPS')
t_full.to_excel(writer, sheet_name='Full_sample', index=False)

t_ad = create_characterization_table(
    df[df['clinical_diagnosis'] == 'AD'], 'AD - Clinically Significant vs Normal NPS')
t_ad.to_excel(writer, sheet_name='AD', index=False)

t_ftd = create_characterization_table(
    df[df['clinical_diagnosis'] == 'FTD'], 'FTD - Clinically Significant vs Normal NPS')
t_ftd.to_excel(writer, sheet_name='FTD', index=False)

t_cn = create_characterization_table(
    df[df['clinical_diagnosis'] == 'CN'], 'CN - Clinically Significant vs Normal NPS')
t_cn.to_excel(writer, sheet_name='CN', index=False)

for sh in writer.sheets:
    ws = writer.sheets[sh]
    ws.column_dimensions['A'].width = 26
    ws.column_dimensions['B'].width = 18
    ws.column_dimensions['C'].width = 18
    ws.column_dimensions['D'].width = 40
    ws.column_dimensions['E'].width = 14
    ws.column_dimensions['F'].width = 16
    ws.column_dimensions['G'].width = 14
writer.close()

print("\n" + "="*80)
print(f"Saved: {excel_file}")
print("   Sheets: Global_by_Country, Full_sample, AD, FTD, CN")
print("="*80)
print("\nNOTES:")
print("  - NPI-Q: Mean (SD) only, no test (grouping derived from NPI = circular).")
print("  - Shapiro/Levene reported; Welch t = primary test.")
print("  - Categorical: Fisher (2x2) or Monte Carlo (rxc) if any expected cell < 5.")
print("  - N_available = non-missing n per variable (missing-data handling, STROBE 12c).")
