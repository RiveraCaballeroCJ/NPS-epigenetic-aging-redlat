# ============================================================================
# NPI-Q symptom characterization by NPS status (Clinically Significant vs Normal)
# ----------------------------------------------------------------------------
# For each of the 12 NPI-Q symptoms, organized by the four Gonzalez domains,
# reports SEPARATELY for the clinically significant and the normal NPS groups:
#   - Presence: n (%) of participants with the symptom (npi_* == 1)
#   - Mean severity (s.d.): among those WHO HAVE the symptom (npi_*sev, 1-3)
#   - Severity breakdown: n (%) mild (1) / moderate (2) / severe (3)
# Two separate tables (one per group). Complete-case per symptom.
#
# Input : analysis dataset (CSV). Set DATA_PATH below.
# Output: NPI_symptom_characterization_NPS.xlsx
# ============================================================================
import pandas as pd, numpy as np
import warnings; warnings.filterwarnings('ignore')

DATA_PATH = '/content/drive/MyDrive/Doctorado/Epi_cr3x.csv'

# Domains (Gonzalez et al.) and display names
DOMAINS = {
 'Psychosis':           ['npi_del', 'npi_hall'],
 'Activation':          ['npi_agit', 'npi_disn', 'npi_irr'],
 'Affective':           ['npi_depd', 'npi_anx'],
 'Somatic/Behavioural': ['npi_apa', 'npi_mot', 'npi_nite', 'npi_app'],
 'Other':               ['npi_elat'],   # euphoria/elation: not part of the four Gonzalez domains
}
NAME = {'npi_del':'Delusions', 'npi_hall':'Hallucinations', 'npi_agit':'Agitation/aggression',
 'npi_disn':'Disinhibition', 'npi_irr':'Irritability/lability', 'npi_depd':'Depression/dysphoria',
 'npi_anx':'Anxiety', 'npi_apa':'Apathy/indifference', 'npi_mot':'Aberrant motor behaviour',
 'npi_nite':'Night-time behaviour', 'npi_app':'Appetite/eating', 'npi_elat':'Euphoria/elation'}

df = pd.read_csv(DATA_PATH)
df = df[df['clinical_diagnosis'].isin(['AD', 'FTD', 'CN'])].copy()
df = df[df['demo_resid'].isin([1, 3, 4, 5, 6])].copy()

# NPS status from the NPI-Q total score (>=3 = clinically significant; 0-2 = normal)
df = df[df['npi_total'].notna()].copy()
df['NPS_status'] = np.where(df['npi_total'] >= 3, 'Significant NPS', 'Normal NPS')

SUBSETS = {
    'Significant NPS': df['npi_total'] >= 3,
    'Normal NPS':      df['npi_total'] < 3,
}

def sev_col(pres_col):
    """Severity column name for a presence column (npi_del -> npi_delsev)."""
    return pres_col + 'sev'

def characterize(sub):
    """Build the characterization rows for one group (a dataframe subset)."""
    rows = []
    for dom, items in DOMAINS.items():
        rows.append({'Domain / Symptom': dom, '_is_header': True})
        for it in items:
            if it not in sub.columns:
                continue
            pres = pd.to_numeric(sub[it], errors='coerce')
            n_valid = pres.notna().sum()
            n_pres = int((pres == 1).sum())
            pct_pres = 100 * n_pres / n_valid if n_valid > 0 else np.nan
            # severity among those who HAVE the symptom
            sc = sev_col(it)
            if sc in sub.columns:
                sev = pd.to_numeric(sub[sc], errors='coerce')
                sev_have = sev[(pres == 1) & sev.notna()]
                if len(sev_have) > 0:
                    msev = f"{sev_have.mean():.2f} ({sev_have.std():.2f})"
                    n1 = int((sev_have == 1).sum()); n2 = int((sev_have == 2).sum()); n3 = int((sev_have == 3).sum())
                    tot = len(sev_have)
                    mild = f"{n1} ({100*n1/tot:.1f}%)"; mod = f"{n2} ({100*n2/tot:.1f}%)"; sev3 = f"{n3} ({100*n3/tot:.1f}%)"
                else:
                    msev = '-'; mild = mod = sev3 = '-'
            else:
                msev = 'NA'; mild = mod = sev3 = 'NA'
            rows.append({'Domain / Symptom': NAME.get(it, it),
                'Present, n (%)': f"{n_pres} ({pct_pres:.1f}%)" if n_valid > 0 else 'NA',
                'Mean severity (s.d.)': msev,
                'Mild, n (%)': mild, 'Moderate, n (%)': mod, 'Severe, n (%)': sev3,
                '_is_header': False})
    return pd.DataFrame(rows)

tables = {}
for g, mask in SUBSETS.items():
    sub = df[mask]
    t = characterize(sub)
    tables[g] = t
    print("="*80)
    print(f"{g}  (N = {len(sub)})")
    print("="*80)
    show = t.drop(columns=['_is_header']).fillna('')
    print(show.to_string(index=False))
    print()

# ---- save: two sheets, domain headers kept as separator rows ----
with pd.ExcelWriter('NPI_symptom_characterization_NPS.xlsx', engine='openpyxl') as w:
    for g in ['Significant NPS', 'Normal NPS']:
        cols = ['Domain / Symptom', 'Present, n (%)', 'Mean severity (s.d.)',
                'Mild, n (%)', 'Moderate, n (%)', 'Severe, n (%)']
        out = tables[g].copy()
        for c in cols:
            if c not in out.columns:
                out[c] = ''
        out = out[cols].fillna('')
        out.to_excel(w, sheet_name=g, index=False)

print("Saved: NPI_symptom_characterization_NPS.xlsx (sheets: Significant NPS, Normal NPS)")
print("\nNote: severity is summarized only among participants who have each symptom.")
print("Euphoria/elation is shown under 'Other' (not part of the four Gonzalez domains).")
