# Replication Materials — Blinded Review Version

**Manuscript:** Synthetic Data Can Mislead Forecasting Benchmarks: Regime-Dependent Performance in Real vs. Simulated U.S. Retail Sales, 2009–2023

**Submitted to:** Journal of Retailing and Consumer Services (Elsevier, SSCI Q1, IF=13.1) — double-anonymized review

This repository contains all code needed to reproduce the results in the paper.
Full author details and a permanent DOI will be added upon acceptance.

---

## Repository Contents

| File | Description |
|---|---|
| `analysis_v2.R` | Main analysis — downloads FRED data, fits all models, runs DM tests, produces all figures |
| `analysis_synthetic_appendix.R` | Synthetic data generation and appendix figures (figA1, figA2) |
| `README.md` | This file |

> Note for reviewers: `retail_data_processed.csv` and `walkforward_predictions.csv` are generated automatically when you run `analysis_v2.R`. They are not included in the repo to avoid file size issues.

---

## Key Results

| Regime | Best model | DM stat | p-value |
|---|---|---|---|
| Post-GFC Recovery (2009–2012) | SARIMA | −1.62 | 0.133 |
| Stable Growth (2013–2019) | **SARIMA** | **−2.26** | **0.027** |
| COVID-19 Disruption (2020–2021) | SARIMA | +0.59 | 0.563 |
| Post-COVID Rebound (2022–2023) | SARIMA | −1.77 | 0.090 |

Switch strategy: −0.5% RMSE vs SARIMA alone — honest negative result reported.

**Variable importance reversal (real vs synthetic):**

| Predictor | Synthetic data | Real FRED data |
|---|---|---|
| Seasonal index | 48.5% (dominant) | 1.4% (last) |
| AR(1) | 5.6% | 17.4% (dominant) |
| Short AR lags combined | 6.8% | 43.7% |

---

## Data Sources (all FRED, public domain)

All data downloaded automatically by `analysis_v2.R` via the `fredr` package.

| Series ID | Description |
|---|---|
| `MRTSSM44000USS` | U.S. Retail Sales: Retail Trade, NSA |
| `CPIAUCSL` | Consumer Price Index, All Urban Consumers |
| `UNRATE` | Civilian Unemployment Rate |
| `DCOILWTICO` | WTI Crude Oil Price |

Get a free FRED API key at: https://fred.stlouisfed.org/docs/api/api_key.html

---

## How to Reproduce

```r
# 1. Set your free FRED API key
Sys.setenv(FRED_API_KEY = "your_key_here")

# 2. Run main analysis (~15 min on standard laptop)
source("analysis_v2.R")

# 3. Run synthetic appendix figures
source("analysis_synthetic_appendix.R")
```

**Requirements:** R 4.x · `set.seed(42)` for full numerical reproducibility

**Key packages** (installed automatically): `fredr`, `tidyverse`, `randomForest`, `forecast`, `strucchange`, `sandwich`, `lmtest`, `Metrics`, `scales`

---

*Note: This repository is intentionally anonymized for double-anonymized peer review. All analytical decisions, code, and results are the author's own. Full attribution will be added upon acceptance.*
