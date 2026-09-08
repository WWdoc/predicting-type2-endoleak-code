# predicting-type2-endoleak-code
Source code for manuscript on predicting Type 2 endoleak


# Reproducible analysis code for predicting type II endoleak

This repository contains generic R scripts corresponding to the tables and figures in the manuscript. Patient-level data are intentionally not included. The scripts use a public, standardized data contract and do not contain manuscript results, personal file paths, or manually entered patient values.

## Repository layout

```text
config.R                         Analysis settings and standardized column names
run_all.R                        Optional batch runner
R/_common.R                     Shared import, validation, formatting, and export helpers
R/01_table1_baseline.R          Table 1: baseline characteristics
R/02_table2_procedural.R        Table 2: procedural and hospitalization details
R/03_table3_outcomes.R          Table 3: clinical outcomes
R/04_table4_anatomy.R           Table 4: aneurysm anatomy measurements
R/05_table5_cox_models.R        Table 5: sequential Cox models
R/06_figures3_4_survival.R      Figures 3 and 4: five Kaplan-Meier curves and sac-status panel
R/07_figure5_roc.R              Figure 5: ROC curves for four anatomical markers
R/08_figure6_rcs.R              Figure 6: restricted cubic spline analyses
R/09_supp_figure1_agreement.R   Supplementary Figure 1: Bland-Altman plots
R/10_supp_figure2_pairplot.R    Supplementary Figure 2: Spearman pair plot
R/11_supp_figure3_cox_lasso.R   Supplementary Figure 3: Cox-LASSO sensitivity analysis
R/12_supp_figure4_bootstrap.R   Supplementary Figure 4: bootstrap AUC and threshold stability
R/13_supp_figure5_distributions.R Supplementary Figure 5: four threshold comparisons
R/14_supp_table1_agreement.R    Supplementary Table 1: reproducibility statistics
R/15_supp_table2_ph_tests.R     Supplementary Table 2: proportional-hazards tests
R/16_supp_table3_diagnostics.R  Supplementary Table 3: 17 diagnostic methods
```

## Data privacy

Place a de-identified analysis file at `data/analysis_data.csv` and, for the observer study, a de-identified file at `data/operator_variability.csv`. The `.gitignore` file excludes everything in `data/` except its documentation. Before a public upload, run `git status --ignored` and confirm that no patient-level file is staged.

The expected standardized columns are listed in `DATA_DICTIONARY.md`. If the local dataset uses different names, either rename a de-identified copy or edit only the mappings in `config.R`.

## Software

R 4.3 or later is recommended. Install the required packages once:

```r
install.packages(c(
  "broom", "dplyr", "flextable", "GGally", "ggplot2", "ggpubr",
  "glmnet", "gtsummary", "openxlsx", "pROC", "patchwork", "readr",
  "readxl", "rlang", "scales", "survival", "survminer", "tidyr"
))
```

## Running the analyses

Run scripts from the repository root. Each script accepts an optional input path and output directory:

```bash
Rscript R/01_table1_baseline.R data/analysis_data.csv outputs
Rscript R/06_figures3_4_survival.R data/analysis_data.csv outputs
Rscript R/09_supp_figure1_agreement.R data/operator_variability.csv outputs
```

To run all analyses after the two de-identified inputs have been prepared:

```bash
Rscript run_all.R
```

Tables are written as CSV and editable Word documents. Figures are written as vector PDF and 300-dpi PNG. Stochastic procedures use a fixed seed of 42. Bootstrap repetitions default to 1000 and can be changed in `config.R`.

## Statistical notes

- Continuous summaries and tests are paired: mean (SD) with Welch t tests, or median (IQR) with Mann-Whitney tests.
- Categorical comparisons use chi-square tests unless sparse expected counts require Fisher or simulated-exact inference.
- Cox models report hazard ratios with 95% confidence intervals, fitted sample size, event count, concordance, and Schoenfeld-residual tests.
- ROC analyses use DeLong confidence intervals and paired DeLong comparisons; sensitivity and specificity use Wilson confidence intervals.
- The Cox-LASSO analysis uses stratified 10-fold cross-validation and bootstrap selection frequencies.
- All analyses are observational and should be interpreted as associations, not causal effects.

## Citation and license

Please cite the associated article when reusing this workflow. The repository owner may apply the license selected when the GitHub repository is created.
