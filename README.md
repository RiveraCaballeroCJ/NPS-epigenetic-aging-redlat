# NPS-epigenetic-aging-redlat
Code for: Associations between Neuropsychiatric symptoms and epigenetic age acceleration in ReDLat

# Neuropsychiatric symptoms and epigenetic age acceleration (ReDLat)

Code for the manuscript *"Neuropsychiatric Symptoms Linked to Accelerated Brain and Epigenetic Aging Across Healthy Aging and Dementia in Latin America"* (<autores>, <año>).
We test whether clinically significant neuropsychiatric symptoms (NPS, NPI-Q ≥ 3)
are associated with epigenetic age acceleration (EAA) across 26 DNA-methylation
clocks in a multi-country Latin American cohort (ReDLat), in the full sample and
by diagnosis (Alzheimer's disease, behavioural frontotemporal dementia, and
cognitively normal controls).

## Repository structure
- `code/` — analysis scripts (Python).
- `figures/` — code to reproduce the figures (`.py` for Python, `.R` for R).
- `requirements.txt` — Python dependencies.
- `install_packages.R` — R dependencies (Figure 2).

## Software and dependencies
**Python** (analyses and Figure 3):
```bash
pip install -r requirements.txt
```
pandas, numpy, statsmodels, scipy, factor_analyzer, matplotlib, openpyxl.

**R / RStudio** (Figure 2 — chord diagrams and forest plots):
```r
install.packages(c("readxl","circlize","dplyr","ggplot2","stringr","svglite"))
```
`grid` ships with base R.

## Analyses
- Primary regression: EAA ~ NPS, OLS with country fixed effects (Colombia reference), full sample and diagnosis-stratified (AD / FTD / CN). FDR within each clock family.
- Exploratory causal mediation (ACME): NPS → EAA → clinical outcomes (MMSE, CDR, T-ADLQ).
- Sensitivity analyses: composite aging factors (EFA), symptom domains (data-driven EFA and theory-based), sex-stratified, leave-one-country-out, country-level random-effects meta-analysis, repeated split-sample, bootstrap, and continuous NPS.

## How to run
Each script sets the input data path at the top and writes result tables or figures.
Install the dependencies above, set the path, and run the script.

## Data availability
Individual-level ReDLat data are **not** included in this repository. They are
available under controlled access from the ReDLat consortium upon reasonable
request and with the appropriate ethical approvals.

## Citation
<cita completa del paper cuando esté publicado>

## License
MIT
